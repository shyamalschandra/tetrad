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
    private double[] weights;
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


        weights = new double[C];
//        A = new TetradMatrix[C];
        W = new TetradMatrix[C];
        K = new TetradMatrix[C];
        S = new TetradMatrix[C];
        posteriorProbs = new TetradMatrix(X.rows(), C);
        bias = new TetradVector[C];

        double _sum = 0.0;

        for (int k = 0; k < C; k++) {
            weights[k] = RandomUtil.getInstance().nextUniform(0, 1);
            _sum += weights[k];
        }

        for (int k = 0; k < C; k++) {
            weights[k] /= _sum;
        }

        for (int k = 0; k < C; k++) {
//            weights[k] = 1.0 / C;

            FastIca ica = new FastIca(data.getDoubleData(), X.columns());
            FastIca.IcaResult result = ica.findComponents();

            X = result.getX();

            TetradMatrix _A = result.getA();
            TetradMatrix _W = result.getW();

            System.out.println("A = " + _A);


//            A[k] = _A;
            W[k] = _W;

            TetradVector _bias = new TetradVector(X.columns());

            for (int i = 0; i < X.columns(); i++) {
                _bias.set(i, 0);
            }

            bias[k] = _bias;
            K[k] = new TetradMatrix(X.columns(), X.columns());
            S[k] = X.minus(repmat(_bias, X.rows())).times(W[k].transpose());
        }

        for (int k = 0; k < C; k++) {
            for (int i = 0; i < X.columns(); i++) {
                double sechSum = 0;
                double sum = 0;
                double tanhSum = 0;

                for (int r = 0; r < X.rows(); r++) {
                    sechSum += 1.0 / pow(cosh(S[k].get(r, i)), 2.0);
                    sum += pow(S[k].get(r, i), 2.0);
                    tanhSum += tanh(S[k].get(r, i)) * S[k].get(r, i);
                }

                sechSum /= X.rows();
                sum /= X.rows();
                tanhSum /= X.rows();

                double kurtosis = signum(sechSum * sum - tanhSum);

                K[k].set(i, i, kurtosis);
            }
        }

        double[] tempWeights = Arrays.copyOf(weights, weights.length);

        boolean weightsUnequal = true;
        ArrayList<Double> diffsList;
        int iterCounter = 0;

        System.out.println("Weights: " + Arrays.toString(weights));
        maximization();



        while (weightsUnequal) {
            expectation();
            maximization();

            System.out.println("Weights: " + Arrays.toString(weights));

            diffsList = new ArrayList<>();

            for (int k = 0; k < C; k++) {
                diffsList.add(abs(weights[k] - tempWeights[k]));
            }

            Collections.sort(diffsList);

            if (diffsList.get(C - 1) < 0.001 || iterCounter > 100) {
                weightsUnequal = false;
            }

            for (int k = 0; k < C; k++) {
                tempWeights[k] = weights[k];
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
    private TetradMatrix[] adaptMixingMatrices(int t) {

        for (int c = 0; c < C; c++) {
            double _prob = posteriorProbs.get(t, c); // scalar
            TetradMatrix _tanhS = X.like(); // NxM

            for (int i = 0; i < X.columns(); i++) {
                _tanhS.set(t, i, Math.tanh(S[c].get(t, i)));
            }

            TetradMatrix skskt = S[c].row(t).transpose().times((S[c].row(t))); // Mx1 * 1xM = MxM
            TetradMatrix tanhSk = _tanhS.row(t).transpose().times(S[c].row(t)); // Mx1 * 1xM = MxM
            TetradMatrix kurtosisSquare = K[c].times(tanhSk);
            TetradMatrix smketc = TetradMatrix.identity(X.columns()).minus(kurtosisSquare);
            TetradMatrix minusSquare = smketc.minus(skskt);

//            TetradMatrix tempWMatrix = minusSquare.times(W[c]).scalarMult(-_prob * learningRate); // MxM * (MxM - (MxM * MxM - MxM)) = MxM
            TetradMatrix tempWMatrix = minusSquare.times(W[c]).scalarMult(_prob * getLearningRate()); // MxM * (MxM - (MxM * MxM - MxM)) = MxM

            TetradMatrix newW = W[c].plus(tempWMatrix);

            if (newW.det() == 0) {
                newW = W[c];
            }

            W[c] = newW;
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
            // bias[k] = bias[k].plus(bias);
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
    private void expectation() {

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(X.rows(), C); // Nxk
        double sum;
        double det;

        double L = 0;

        for (int c = 0; c < C; c++) {
            det = Math.log(Math.abs(W[c].det()));

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
                _weightedProbs[c] = exp(likelihoods.get(t, c)) * weights[c];
                _sum += _weightedProbs[c];
            }

            for (int c = 0; c < C; c++) {
                double prob = _weightedProbs[c] / _sum;
                posteriorProbs.set(t, c, prob);
            }

            System.out.println("t = " + t + " " + Arrays.toString(_weightedProbs)
                    + " " + likelihoods.row(t) + " weights = " + Arrays.toString(weights)
                    + " " + posteriorProbs.row(t));

            System.out.println("PROBS = " + posteriorProbs.row(t));
        }

        System.out.println(posteriorProbs);

        for (int t = 0; t < X.rows(); t++) {
            W = adaptMixingMatrices(t);
//            A = invertWMatrices(W);
            bias = adaptBiasVectors(t);
        }
    }

    /*
     * Find maximally likely values for source vectors, bias vectors, and weights of components
     */
    private void maximization1() {

        // find values of weights
        double sum;

        for (int c = 0; c < C; c++) {
            sum = 0;

            for (int r = 0; r < X.rows(); r++) {
                sum += posteriorProbs.get(r, c);
            }

            weights[c] = sum / X.rows();
        }

        for (int k = 0; k < C; k++) {
            S[k] = X.minus(repmat(bias[k], X.rows())).times(W[k].transpose());
        }
    }

    private void maximization() {

        double[] weights = new double[C];

        // find values of weights
        double sum;

        for (int c = 0; c < C; c++) {
            sum = 0;

            for (int r = 0; r < X.rows(); r++) {
                sum += posteriorProbs.get(r, c);
            }

            weights[c] = 1. / C;// sum / X.rows();
        }

        this.weights = weights;

//        sourceVector = X.like();
//
//        for (int k = 0; k < C; k++) {
//            // find values of source vectors
//            for (int r = 0; r < X.rows(); r++) {
//                TetradVector vector = X.getRow(r).minus(bias[k]);
//                TetradVector tempVector = A[k].transpose().inverse().times(vector);
//
//                sourceVector.assignRow(r, tempVector);
//            }
//        }

        System.out.println("W[0] = " + W[0]);

        System.out.println("W[1] = " + W[1]);

        for (int k = 0; k < C; k++) {
            S[k] = X.minus(repmat(bias[k], X.rows())).times(W[k].transpose());
        }

        // find values of bias vectors
        //  TetradVector biasVector;
        //  double biasVectorDenom;
        // sourceVector = new TetradMatrix(numCases, numVars);

        // I:
        // for (int I = 0; I < K; I++) {
        //    biasVector = new TetradVector(numVars);
        //    for (int c = 0; c < numVars; c++) {
        //       biasVector.set(c, 0);
        //    }

        //    biasVectorDenom = 0;

        //   for (int r = 0; r < numCases; r++) {
        //       biasVector = biasVector.plus(dataMatrix.getRow(r)).scalarMult(gammas.get(r, I));
        //       biasVectorDenom += gammas.get(r, I);
        //  }

        //   biasVector = biasVector.scalarMult(1.0 / biasVectorDenom);
        //   biasVector = stats.getBiasVectors()[I];

        //debugging
            /* for (int b = 0; b < biasVector.size(); b++) {
                double val = biasVector.get(b);

                if (Double.isNaN(val)) {
                    throw new IllegalArgumentException("NaN in Bias Vector " + I);
//                    continue I;
                }
            } */

    }

    public static void main(String... args) {
        DataSet dataSet = loadData("/Users/user/Downloads/16.txt");

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
                    new File(path), Delimiter.WHITESPACE);
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