package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;
import edu.pitt.dbmi.data.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDataFileReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import static java.lang.Math.*;

/**
 * @author Madelyn Glymour 7/21/18
 */
public class DemixerNongaussian {
    private final DataSet data;
    private final TetradMatrix X;
    private final int numComponents;
    private final int numVars;
    private final TetradMatrix weights;
    private final TetradMatrix[] S;
    private final TetradMatrix[] K;
    private final TetradMatrix posteriorProbs;
    private final double learningRate;

    private int T;

    private TetradMatrix[] bias;
    private TetradMatrix[] W;


    private DemixerNongaussian(DataSet data, int numComponents, double learningRate) {
        this.data = data;
        this.X = data.getDoubleData();
        this.T = X.rows();
        this.numVars = X.columns();
        this.numComponents = numComponents;
        this.bias = new TetradMatrix[this.numComponents];
        this.weights = new TetradMatrix(1, this.numComponents);
        this.posteriorProbs = new TetradMatrix(T, this.numComponents);
        this.W = new TetradMatrix[this.numComponents];
        this.K = new TetradMatrix[this.numComponents];
        this.S = new TetradMatrix[this.numComponents];
        this.learningRate = learningRate;

        for (int t = 0; t < T; t++) {
            for (int k = 0; k < numComponents; k++) {
                posteriorProbs.set(t, k, 1.0 / numComponents);
            }
        }
    }

    private MixtureModelNongaussian demix() {
        for (int k = 0; k < numComponents; k++) {
            weights.set(0, k, 1);//RandomUtil.getInstance().nextUniform(0, 1));
        }

        normalize(weights, 0);

        for (int k = 0; k < numComponents; k++) {
            FastIca ica = new FastIca(X, numVars);
            FastIca.IcaResult result = ica.findComponents();

            TetradMatrix _W = result.getW();

//            System.out.println("W = " + _W);
//
            W[k] = _W;
//            W[k] = whiten(_W, 0.001);

            TetradMatrix _bias = new TetradMatrix(1, numVars);

            for (int i = 0; i < numVars; i++) {
                _bias.set(0, i, 0);// RandomUtil.getInstance().nextNormal(0, 0.01));
            }

            bias[k] = _bias;
            K[k] = new TetradMatrix(numVars, numVars);
            S[k] = X.minus(repmat(_bias.getRow(0), T)).times(W[k].transpose());
        }

        for (int k = 0; k < numComponents; k++) {
            for (int v = 0; v < numVars; v++) {
                double sechSum = 0;
                double _sum = 0;
                double tanhSum = 0;

                for (int t = 0; t < T; t++) {
                    sechSum += 1.0 / pow(cosh(S[k].get(t, v)), 2.0);
                    _sum += pow(S[k].get(t, v), 2.0);
                    tanhSum += tanh(S[k].get(t, v)) * S[k].get(t, v);
                }

                sechSum /= T;
                _sum /= T;
                tanhSum /= T;

                double kurtosis = signum(sechSum * _sum - tanhSum);

                K[k].set(v, v, kurtosis);
            }
        }

        TetradMatrix tempWeights = new TetradMatrix(weights);

        boolean weightsUnequal = true;
        ArrayList<Double> diffsList;
        int i = 0;

        System.out.println("Weights: " + weights);
        maximization();

        while (weightsUnequal) {
            expectation();
            maximization();

            if (i++ > 1) break;

//            System.out.println("Weights: " + weights);

//            diffsList = new ArrayList<>();
//
//            for (int k = 0; k < numComponents; k++) {
//                diffsList.add(abs(weights.get(0, k) - tempWeights.get(0, k)));
//            }
//
//            Collections.sort(diffsList);
//
//            if (diffsList.get(numComponents - 1) < 0.001 || i > 100) {
//                weightsUnequal = false;
//            }
//
//            for (int k = 0; k < numComponents; k++) {
//                tempWeights.set(0, k, weights.get(0, k));
//            }
//
//            i++;
        }

        return new MixtureModelNongaussian(data, posteriorProbs, invertWMatrices(W), S, bias, weights);
    }

    private TetradMatrix whiten(TetradMatrix w, double v) {
        w = new TetradMatrix(w);

        for (int r = 0; r < w.rows(); r++) {
            for (int c = 0; c < w.columns(); c++) {
                w.set(r, c, w.get(r, c) + RandomUtil.getInstance().nextNormal(0, v));
            }
        }

        return w;
    }

    private TetradMatrix repmat(TetradVector _bias, int rows) {
        TetradMatrix rep = X.like();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < _bias.size(); c++) {
                rep.set(r, c, _bias.get(c));
            }
        }

        return rep;
    }

    /*
     * Inline updating of the mixing matrix A, updated every time a new posterior is determined
     */
    private TetradMatrix[] adaptMixingMatrices() {

        for (int k = 0; k < numComponents; k++) {
            TetradMatrix _tanhS = X.like(); // NxM

            for (int t = 0; t < T; t++) {
                for (int i = 0; i < numVars; i++) {
                    _tanhS.set(t, i, Math.tanh(S[k].get(t, i)));
                }
            }

            TetradMatrix skskt = S[k].transpose().times((S[k])); // Mx1 * 1xM = MxM
            TetradMatrix tanhSk = _tanhS.transpose().times(S[k]); // Mx1 * 1xM = MxM
            TetradMatrix kurtosisSquare = K[k].times(tanhSk);
            TetradMatrix smketc = TetradMatrix.identity(numVars).minus(kurtosisSquare);
            TetradMatrix minusSquare = smketc.minus(skskt);

            TetradMatrix probs = posteriorProbs.column(k);
            TetradMatrix newW = new TetradMatrix(numVars, numVars);

            for (int t = 0; t < T; t++) {
                double p = probs.get(t, 0);
                TetradMatrix tempWMatrix = minusSquare.times(W[k]).scalarMult(p * learningRate / (T * numComponents));
                newW = newW.plus(tempWMatrix);
            }

            final double det = newW.det();

            if (abs(det) > 0) {
                W[k] = W[k].plus(newW);
            }
        }

        return W;
    }

    private TetradMatrix[] invertWMatrices(TetradMatrix[] m) {
        TetradMatrix[] n = new TetradMatrix[m.length];

        for (int k = 0; k < m.length; k++) {
            n[k] = m[k].inverse();
        }

        return n;
    }

    private TetradMatrix[] adaptBiasVectors(int t) {
        for (int k = 0; k < numComponents; k++) {
            TetradMatrix _W = W[k];
            TetradMatrix sourceVector = S[k];
            double gamma = posteriorProbs.get(t, k);

            TetradMatrix _bias = new TetradMatrix(1, numVars);
            double sum = 0;
            double det = Math.log(Math.abs(_W.det()));

            if (Double.isNaN(det)) {
                return bias;
            }

            for (int i = 0; i < numVars; i++) {
                TetradMatrix bMatrix = _W.row(i);
                _bias = _bias.plus((bMatrix.scalarMult(Math.tanh(sourceVector.get(t, i)))).plus(bMatrix.scalarMult(sourceVector.get(t, i))));

                double tempSum = -.5 * Math.log(2.0 * Math.PI) - (K[k].get(i, i) * log(cosh(S[k].get(t, i))))
                        - (pow(S[k].get(t, i), 2.0) / 2.0);

                if (K[k].get(i, i) > 0) {
                    tempSum = tempSum - Math.log(0.7413);
                }

                sum += tempSum;
            }

            double L = -sum - det;

            _bias = _bias.scalarMult(gamma * L);
            bias[k] = _bias;
        }

        return bias;
    }

    /*
     * Find posteriors of observations based on maximally likely values of mixing matrices, bias vectors, source vectors, and weights
     */
    private boolean expectation() {

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(T, numComponents); // Nxk
        double sum;
        double det;
        double _sum = 0.0;

        for (int k = 0; k < numComponents; k++) {
            det = Math.log(Math.abs(W[k].det()));

            if (Double.isNaN(det)) {
                return true;
            }

            for (int t = 0; t < X.rows(); t++) {
                sum = 0;

                for (int i = 0; i < X.columns(); i++) {
                    double _l = -.5 * Math.log(2.0 * Math.PI) - (K[k].get(i, i) * log(cosh(S[k].get(t, i))))
                            - (pow(S[k].get(t, i), 2.0) / 2.0);

                    if (K[k].get(i, i) > 0) {
                        _l = _l - Math.log(0.7413);
                    }

                    sum += _l;
                }

                double L = sum - det;

                likelihoods.set(t, k, L);
                sum += L;
            }
        }

        System.out.println("Likelihoods = " + likelihoods);
        System.out.println("L = " + _sum);
        System.out.println(" Done");

        for (int t = 0; t < T; t++) {
            double[] probs = new double[numComponents];

            for (int k = 0; k < numComponents; k++) {
                probs[k] = exp(likelihoods.get(t, k)) * weights.get(0, k);
            }

            for (int k = 0; k < numComponents; k++) {
                posteriorProbs.set(t, k, probs[k]);
            }
            
            normalize(posteriorProbs, t);

//            System.out.println("t = " + t + " " + Arrays.toString(_weightedProbs)
//                    + " " + likelihoods.row(t) + " weights = " + weights
//                    + " " + posteriorProbs.row(t));

            System.out.println("PROBS = " + posteriorProbs.row(t));

//        System.out.println("## posterior probs = " + posteriorProbs);

            bias = adaptBiasVectors(t);
        }

        W = adaptMixingMatrices();
        return true;
    }

    private void maximization() {

        // find values of weights
        double sum;

        for (int k = 0; k < numComponents; k++) {
            sum = 0;
            int count = 0;

            for (int r = 0; r < T; r++) {
                final double prob = posteriorProbs.get(r, k);

                if (!Double.isNaN(prob)) {
                    sum += prob;
                    count++;
                }
            }

            weights.set(0, k, sum / count);
        }

        normalize(weights, 0);

        System.out.println("weights = " + weights);

        for (int k = 0; k < numComponents; k++) {
            S[k] = X.minus(repmat(bias[k].getRow(0), T)).times(W[k].transpose());
        }
    }

    private void normalize(TetradMatrix weights, int row) {
        double _sum = 0.0;

        for (int k = 0; k < numComponents; k++) {
            _sum += weights.get(row, k);
        }

        for (int k = 0; k < numComponents; k++) {
            weights.set(0, k, weights.get(row, k) / _sum);
        }
    }


    private static DataSet loadData(String path) {
        try {
            ContinuousTabularDataFileReader dataReader = new ContinuousTabularDataFileReader(
                    new File(path), Delimiter.COMMA);
            dataReader.setHasHeader(true);
            return (DataSet) DataConvertUtils.toDataModel(dataReader.readInData());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public double getLearningRate() {
        return learningRate;
    }

    public static void main(String... args) {

        DataSet dataSet = loadData("/Users/user/Downloads/mixfile1.csv");

        long startTime = System.currentTimeMillis();

        DemixerNongaussian pedro = new DemixerNongaussian(dataSet, 2, 0.001);
        MixtureModelNongaussian model = pedro.demix();

        long elapsed = System.currentTimeMillis() - startTime;

        DataSet[] datasets = model.getDemixedData();

        System.out.println(model.getMixingMatrices()[0]);
        System.out.println();
        System.out.println(model.getBias()[0]);
        System.out.println();
        System.out.println("\n\nDatasets:");

        int count = 1;

        for (DataSet _dataSet : datasets) {
            System.out.println("#" + count++ + ": rows = " + _dataSet.getNumRows() + " cols = " + _dataSet.getNumColumns());
        }

        System.out.println();

        System.out.println("Elapsed: " + (elapsed / 1000));
    }

}