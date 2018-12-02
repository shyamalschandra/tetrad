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
    private TetradMatrix[] A;
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
        this.X = data.getDoubleData();

        FastIca ica = new FastIca(X, X.columns());
        FastIca.IcaResult result = ica.findComponents();

        TetradMatrix _A = result.getA();
        TetradMatrix _W = result.getW();

        weights = new double[C];
        A = new TetradMatrix[C];
        W = new TetradMatrix[C];
        K = new TetradMatrix[C];
        S = new TetradMatrix[C];
        posteriorProbs = new TetradMatrix(X.rows(), C);
        bias = new TetradVector[C];

        for (int k = 0; k < C; k++) {
            weights[k] = 1.0 / C;

            A[k] = _A;
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

        return new MixtureModelNongaussian(data, posteriorProbs, A, S, bias, weights);
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

        for (int k = 0; k < C; k++) {
            double _gamma = posteriorProbs.get(t, k); // scalar
            TetradMatrix _tanhS = X.like(); // NxM

            for (int i = 0; i < X.columns(); i++) {
                _tanhS.set(t, i, Math.tanh(S[k].get(t, i)));
            }

            TetradMatrix skskt = toColumnMatrix(S[k].getRow(t)).times(toRowMatrix(S[k].getRow(t))); // Mx1 * 1xM = MxM
            TetradMatrix tanhSk = toColumnMatrix(_tanhS.getRow(t)).times(toRowMatrix(S[k].getRow(t))); // Mx1 * 1xM = MxM
            TetradMatrix kurtosisSquare = K[k].times(tanhSk);
            TetradMatrix smketc = TetradMatrix.identity(X.columns()).minus(kurtosisSquare);
            TetradMatrix minusSquare = smketc.minus(skskt);

            // tempMixingMatrix = minusSquare.times(wMatrix).scalarMult(-_gamma * learningRate); // MxM * (MxM - (MxM * MxM - MxM)) = MxM
            TetradMatrix tempWMatrix = minusSquare.times(W[k]).scalarMult(_gamma * getLearningRate()); // MxM * (MxM - (MxM * MxM - MxM)) = MxM

            TetradMatrix newW = W[k].plus(tempWMatrix);

            if (newW.det() == 0) {
                newW = W[k];
            }

            W[k] = newW;
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

        for (int k = 0; k < C; k++) {
            det = Math.log(Math.abs(A[k].det()));

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
            }
        }

        for (int t = 0; t < X.rows(); t++) {
            double[] _weightedLIkelihoods = new double[C];
            double _sum = 0.0;

            for (int c = 0; c < C; c++) {
                _weightedLIkelihoods[c] = exp(likelihoods.get(t, c)) * weights[c];
                _sum += _weightedLIkelihoods[c];
            }

            for (int c = 0; c < C; c++) {
                double prob = _weightedLIkelihoods[c] / _sum;
                posteriorProbs.set(t, c, prob);
            }
        }

        for (int t = 0; t < X.rows(); t++) {
            W = adaptMixingMatrices(t);
            A = invertWMatrices(W);
            bias = adaptBiasVectors(t);
        }
    }

    /*
     * Find maximally likely values for source vectors, bias vectors, and weights of components
     */
    private void maximization() {

        // find values of weights
        double sum;

        for (int i = 0; i < C; i++) {
            sum = 0;

            for (int r = 0; r < X.rows(); r++) {
                sum += posteriorProbs.get(r, i);
            }

            weights[i] = sum / X.rows();
        }

        for (int k = 0; k < C; k++) {
            S[k] = X.minus(repmat(bias[k], X.rows())).times(W[k].transpose());
        }
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
