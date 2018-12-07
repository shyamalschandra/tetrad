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
    private TetradMatrix X;
    private int C;
    private TetradMatrix weights;
    //    private TetradMatrix[] A;
    private TetradMatrix[] W;
    private TetradVector[] bias;
    private TetradMatrix[] S;
    private TetradMatrix[] K;
    private TetradMatrix posteriorProbs;
    private double learningRate = .001;

    public DemixerNongaussian(DataSet data) {
        X = data.getDoubleData();
    }

    public MixtureModelNongaussian demix(DataSet data, int numComponents) {
        this.C = numComponents;


        weights = new TetradMatrix(1, C);
        W = new TetradMatrix[C];
        K = new TetradMatrix[C];
        S = new TetradMatrix[C];
        posteriorProbs = new TetradMatrix(X.rows(), C);
        bias = new TetradVector[C];

        double sum = 0.0;

        for (int k = 0; k < C; k++) {
            weights.set(0, k, RandomUtil.getInstance().nextUniform(0, 1));
            sum += weights.get(0, k);
        }

        for (int k = 0; k < C; k++) {
            weights.set(0, k, weights.get(0, k) / sum);
        }

        for (int k = 0; k < C; k++) {
            weights.set(0, k, 1.0 / C);

            FastIca ica = new FastIca(data.getDoubleData(), X.columns());
            FastIca.IcaResult result = ica.findComponents();

            X = result.getX();

//            TetradMatrix _A = result.getA();
            TetradMatrix _W = result.getW();

            System.out.println("W = " + _W);

            W[k] = _W;

            TetradVector _bias = new TetradVector(X.columns());

            for (int i = 0; i < X.columns(); i++) {
                _bias.set(i, .1);
            }

            bias[k] = _bias;
            K[k] = new TetradMatrix(X.columns(), X.columns());
            S[k] = X.minus(repmat(_bias, X.rows())).times(W[k].transpose());
        }

        for (int k = 0; k < C; k++) {
            for (int i = 0; i < X.columns(); i++) {
                double sechSum = 0;
                double _sum = 0;
                double tanhSum = 0;

                for (int r = 0; r < X.rows(); r++) {
                    sechSum += 1.0 / pow(cosh(S[k].get(r, i)), 2.0);
                    _sum += pow(S[k].get(r, i), 2.0);
                    tanhSum += tanh(S[k].get(r, i)) * S[k].get(r, i);
                }

                sechSum /= X.rows();
                _sum /= X.rows();
                tanhSum /= X.rows();

                double kurtosis = signum(sechSum * _sum - tanhSum);

                K[k].set(i, i, kurtosis);
            }
        }

        TetradMatrix tempWeights = new TetradMatrix(weights);

        boolean weightsUnequal = true;
        ArrayList<Double> diffsList;
        int iterCounter = 0;

        System.out.println("Weights: " + weights);
//        maximization();


        while (weightsUnequal) {
            if (!expectation()) {
//                break;
            }
            maximization();

            System.out.println("Weights: " + weights);

            diffsList = new ArrayList<>();

            for (int k = 0; k < C; k++) {
                diffsList.add(abs(weights.get(0, k) - tempWeights.get(0, k)));
            }

            Collections.sort(diffsList);

            if (diffsList.get(C - 1) < 0.001 || iterCounter > 100) {
                weightsUnequal = false;
            }

            for (int k = 0; k < C; k++) {
                tempWeights.set(0, k, weights.get(0, k));
            }

            iterCounter++;
        }

        return new MixtureModelNongaussian(data, posteriorProbs, invertWMatrices(W), S, bias, weights);
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

        for (int c = 0; c < C; c++) {
            TetradMatrix _tanhS = X.like(); // NxM

            for (int t = 0; t < X.rows(); t++) {
                for (int i = 0; i < X.columns(); i++) {
                    _tanhS.set(t, i, Math.tanh(S[c].get(t, i)));
                }
            }

            TetradMatrix skskt = S[c].transpose().times((S[c])); // Mx1 * 1xM = MxM
            TetradMatrix tanhSk = _tanhS.transpose().times(S[c]); // Mx1 * 1xM = MxM
            TetradMatrix kurtosisSquare = K[c].times(tanhSk);
            TetradMatrix smketc = TetradMatrix.identity(X.columns()).minus(kurtosisSquare);
            TetradMatrix minusSquare = smketc.minus(skskt);

            TetradMatrix probs = posteriorProbs.column(c);
            TetradMatrix newW = new TetradMatrix(X.columns(), X.columns());

            learningRate = 0.001;

            for (int t = 0; t < X.rows(); t++) {
                double p = probs.get(t, 0);
                TetradMatrix tempWMatrix = minusSquare.times(W[c]).scalarMult(p * getLearningRate() / X.rows());
                newW = newW.plus(tempWMatrix);
            }

            final double det = newW.det();

            if (Double.isNaN(det)) {
                return W;
            }

            System.out.println("newW = " + newW + " det = " + det);

            if (abs(det) > 0.01) {
                W[c] = W[c].plus(newW);
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

    private TetradVector[] adaptBiasVectors(int t) {
        for (int k = 0; k < C; k++) {
            TetradMatrix _W = W[k];
            TetradMatrix sourceVector = S[k];
            double gamma = posteriorProbs.get(t, k);

            TetradVector _bias = new TetradVector(X.columns());
            double sum = 0;
            double det = Math.log(Math.abs(_W.det()));

            if (Double.isNaN(det)) {
                return bias;
            }

            for (int i = 0; i < X.columns(); i++) {
                TetradVector bMatrix = _W.getRow(i);
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
     * Helper method for matrix multiplication
     */
    private TetradMatrix toColumnMatrix(TetradVector vector) {
        TetradMatrix m = new TetradMatrix(vector.size(), 1);

        for (int i = 0; i < vector.size(); i++) {
            m.set(i, 0, vector.get(i));
        }

        return m;
    }

    /*
     * Helper method for matrix multiplication
     */
    private TetradMatrix toRowMatrix(TetradVector vector) {
        TetradMatrix m = new TetradMatrix(1, vector.size());

        for (int i = 0; i < vector.size(); i++) {
            m.set(0, i, vector.get(i));
        }

        return m;
    }


    /*
     * Find posteriors of observations based on maximally likely values of mixing matrices, bias vectors, source vectors, and weights
     */
    private boolean expectation() {

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(X.rows(), C); // Nxk
        double sum;
        double det;

        double L = 0;

        for (int c = 0; c < C; c++) {
            det = Math.log(Math.abs(W[c].det()));

            if (Double.isNaN(det)) {
                return true;
            }

            for (int t = 0; t < X.rows(); t++) {
                sum = 0;

                for (int i = 0; i < X.columns(); i++) {
                    double _l = -.5 * Math.log(2.0 * Math.PI) - (K[c].get(i, i) * log(cosh(S[c].get(t, i))))
                            - (pow(S[c].get(t, i), 2.0) / 2.0);

                    if (K[c].get(i, i) > 0) {
                        _l = _l - Math.log(0.7413);
                    }

                    sum += _l;
                }

                double l = sum - det;

                likelihoods.set(t, c, l);

                L += l;
            }
        }

        System.out.println("L = " + L);

        for (int t = 0; t < X.rows(); t++) {
            double[] _weightedProbs = new double[C];
            double _sum = 0.0;

            for (int c = 0; c < C; c++) {
                _weightedProbs[c] = exp(likelihoods.get(t, c)) * weights.get(0, c);
                _sum += _weightedProbs[c];
            }

            for (int c = 0; c < C; c++) {
                double prob = _weightedProbs[c] / _sum;
                posteriorProbs.set(t, c, prob);
            }

//            System.out.println("t = " + t + " " + Arrays.toString(_weightedProbs)
//                    + " " + likelihoods.row(t) + " weights = " + weights
//                    + " " + posteriorProbs.row(t));

//            System.out.println("PROBS = " + posteriorProbs.row(t));
        }

//        System.out.println("## posterior probs = " + posteriorProbs);

        for (int t = 0; t < X.rows(); t++) {
            W = adaptMixingMatrices();
            bias = adaptBiasVectors(t);
        }

        return true;
    }

    private void maximization() {

        TetradMatrix weights = new TetradMatrix(1, C);

        // find values of weights
        double sum;

        for (int c = 0; c < C; c++) {
            sum = 0;

            for (int r = 0; r < X.rows(); r++) {
                sum += posteriorProbs.get(r, c);
            }

            weights.set(0, c, sum / X.rows());
        }

        this.weights = weights;

        System.out.println("weights = " + weights);

        for (int k = 0; k < C; k++) {
            S[k] = X.minus(repmat(bias[k], X.rows())).times(W[k].transpose());
        }
    }

    public static void main(String... args) {
        DataSet dataSet = loadData("/Users/user/Downloads/mixfile1.csv");

        DemixerNongaussian pedro = new DemixerNongaussian(dataSet);
        long startTime = System.currentTimeMillis();
        MixtureModelNongaussian model = pedro.demix(dataSet, 2);
        long elapsed = System.currentTimeMillis() - startTime;

        DataSet[] datasets = model.getDemixedData();

        System.out.println(model.getMixingMatrices()[0]);
        System.out.println();
        System.out.println(model.getBiasVectors()[0]);
        System.out.println();
        System.out.println("\n\nDatasets:");

        int count = 1;

        for (DataSet _dataSet : datasets) {
            System.out.println("#" + count++ + ": rows = " + _dataSet.getNumRows() + " cols = " + _dataSet.getNumColumns());
        }

        System.out.println();

        System.out.println("Elapsed: " + (elapsed / 1000));
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

    public void setLearningRate(double learningRate) {
        this.learningRate = learningRate;
    }
}