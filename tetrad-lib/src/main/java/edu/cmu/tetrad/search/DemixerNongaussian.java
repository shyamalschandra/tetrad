package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;
import edu.pitt.dbmi.data.Delimiter;
import edu.pitt.dbmi.data.reader.DataReader;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDataFileReader;
import org.apache.commons.math3.linear.SingularMatrixException;

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
    private int numCases;
    private int numVars;
    private TetradMatrix dataMatrix;

    public DemixerNongaussian(DataSet data) {
        numCases = data.getNumRows();
        numVars = data.getNumColumns();
        dataMatrix = data.getDoubleData();
    }

    public MixtureModelNongaussian demix(DataSet data, int numComponents) {

        double[] weights = new double[numComponents];
        TetradMatrix[] A = new TetradMatrix[numComponents];
        TetradMatrix[] W = new TetradMatrix[numComponents];
        TetradVector[] bias = new TetradVector[numComponents];
        TetradMatrix[] S = new TetradMatrix[numComponents];
        TetradMatrix[] K = new TetradMatrix[numComponents];
        TetradMatrix _S;
        TetradMatrix _K;
        TetradVector _bias;

        TetradMatrix dataCopy = dataMatrix.copy();

        FastIca ica = new FastIca(dataCopy, numVars);
        FastIca.IcaResult result = ica.findComponents();

        TetradMatrix _A = result.getA();
        TetradMatrix _W = result.getW();

        for (int k = 0; k < numComponents; k++) {
            weights[k] = 1.0 / numComponents;
            _K = new TetradMatrix(numVars, numVars); // MxM

            _bias = new TetradVector(numVars); // 1xM

            for (int i = 0; i < numVars; i++) {
                _bias.set(i, 0);
            }

            A[k] = _A;
            W[k] = _W;

            System.out.println(A[k]);

            bias[k] = _bias;
            K[k] = _K;

            _S = dataMatrix.minus(repmat(_bias)).times(W[k].transpose());

            S[k] = _S;
//            System.out.println(_S);
        }

        double sechSum;
        double sum;
        double tanhSum;
        double kurtosis;

        for (int k = 0; k < numComponents; k++) {
            _K = K[k]; // MxM
            _S = S[k]; // NxM

            for (int i = 0; i < numVars; i++) {
                sechSum = 0;
                sum = 0;
                tanhSum = 0;

                for (int r = 0; r < numCases; r++) {
                    sechSum += 1.0 / pow(cosh(_S.get(r, i)), 2.0);
                    sum += pow(_S.get(r, i), 2.0);
                    tanhSum += tanh(_S.get(r, i)) * _S.get(r, i);
                }

                sechSum /= numCases;
                sum /= numCases;
                tanhSum /= numCases;
                kurtosis = signum(sechSum * sum - tanhSum);

                _K.set(i, i, kurtosis);
            }

            System.out.println("K = " + _K);
        }

        double[] tempWeights = Arrays.copyOf(weights, weights.length);

        boolean weightsUnequal = true;
        ArrayList<Double> diffsList;
        int iterCounter = 0;

        System.out.println("Weights: " + Arrays.toString(weights));

        DeterminingStats stats = new DeterminingStats();
        stats.setSourceVectors(S);
        stats.setWeights(weights);
        stats.setMixingMatrices(A);
        stats.setBiasVectors(bias);
        stats.setK(numComponents);
        stats.setKurtosisMatrices(K);
        stats.setWMatrices(W);

        while (weightsUnequal) {
            expectation(stats);
            maximization(stats);

            System.out.println("Weights: " + Arrays.toString(stats.getWeights()));

          /*  double likely;

            for(int k = 0; k < numComponents; k++){
                for(int t = 0; t < numCases; t++){
                    likely = stats.getGammas().get(t, k) * stats.getWeights()[k];
                    System.out.println(likely);
                }
            }*/

            diffsList = new ArrayList<>();

            for (int k = 0; k < numComponents; k++) {
                diffsList.add(abs(stats.getWeights()[k] - tempWeights[k]));
            }

            Collections.sort(diffsList);

            if (diffsList.get(numComponents - 1) < 0.001 || iterCounter > 100) {
                weightsUnequal = false;
            }

            for (int k = 0; k < numComponents; k++) {
                tempWeights[k] = stats.getWeights()[k];
            }

            iterCounter++;
        }

        return new MixtureModelNongaussian(data, stats.getPosteriorProbs(), stats.getMixingMatrices(),
                stats.getSourceVectors(), stats.getBiasVectors(), stats.getWeights());
    }

    private TetradMatrix repmat(TetradVector _bias) {
        TetradMatrix rep = new TetradMatrix(numCases, numVars);

        for (int r = 0; r < numCases; r++) {
            for (int c = 0; c < _bias.size(); c++) {
                rep.set(r, c, _bias.get(c));
            }
        }
        return rep;
    }

    private TetradMatrix whiten(TetradMatrix A) {
        TetradMatrix S = A.like();

        for (int r = 0; r < A.rows(); r++) {
            for (int c = 0; c < A.columns(); c++) {
                S.set(r, c, A.get(r, c) + RandomUtil.getInstance().nextNormal(0, 1));
            }
        }

        return S;
    }

    /*
     * Inline updating of the mixing matrix A, updated every time a new posterior is determined
     */
    private TetradMatrix[] adaptMixingMatrices(DeterminingStats stats, int t) {

        TetradMatrix[] A = stats.getMixingMatrices();
        TetradMatrix[] W = stats.getWMatrices();
        TetradMatrix[] K = stats.getKurtosisMatrices();
        TetradMatrix[] S = stats.getSourceVectors();
        TetradMatrix gammas = stats.getPosteriorProbs();

        for (int k = 0; k < stats.getNumComponents(); k++) {
            TetradMatrix _W = W[k];
            TetradMatrix _K = K[k]; // MxM
            TetradMatrix _S = S[k]; // NxM
            double _gamma = gammas.get(t, k); // scalar
            TetradMatrix _tanhS = new TetradMatrix(numCases, numVars); // NxM

            for (int i = 0; i < numVars; i++) {
                _tanhS.set(t, i, Math.tanh(_S.get(t, i)));
            }

            TetradMatrix skskt = toColumnMatrix(_S.getRow(t)).times(toRowMatrix(_S.getRow(t))); // Mx1 * 1xM = MxM
            TetradMatrix tanhSk = toColumnMatrix(_tanhS.getRow(t)).times(toRowMatrix(_S.getRow(t))); // Mx1 * 1xM = MxM

            TetradMatrix kurtosisSquare = _K.times(tanhSk);

            TetradMatrix smketc = TetradMatrix.identity(numVars).minus(kurtosisSquare);

            TetradMatrix minusSquare = smketc.minus(skskt);

            final double learningRate = .01;
            // tempMixingMatrix = minusSquare.times(wMatrix).scalarMult(-_gamma * learningRate); // MxM * (MxM - (MxM * MxM - MxM)) = MxM
            TetradMatrix tempWMatrix = minusSquare.times(_W).scalarMult(_gamma * learningRate); // MxM * (MxM - (MxM * MxM - MxM)) = MxM

            TetradMatrix newW = _W.plus(tempWMatrix);

            if (newW.det() == 0) {
                newW = _W;
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

    private TetradVector[] adaptBiasVectors(DeterminingStats stats, int t) {

        TetradMatrix[] W = stats.getWMatrices();
        TetradMatrix[] S = stats.getSourceVectors();
        TetradVector[] bias = stats.getBiasVectors();
        TetradMatrix gammas = stats.getPosteriorProbs();
        TetradMatrix[] K = stats.getKurtosisMatrices();

        for (int k = 0; k < stats.getNumComponents(); k++) {
            TetradMatrix _W = W[k];
            TetradMatrix sourceVector = S[k];
            double gamma = gammas.get(t, k);

            TetradVector _bias = new TetradVector(numVars);
            double sum = 0;
            double det = Math.log(Math.abs(_W.det()));

            for (int i = 0; i < numVars; i++) {
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
     * Find maximally likely values for source vectors, bias vectors, and weights of components
     */
    private void maximization(DeterminingStats stats) {

        int K = stats.getNumComponents();
//        TetradMatrix[] A = stats.getMixingMatrices();
        TetradMatrix[] W = stats.getWMatrices();
        TetradMatrix[] S = stats.getSourceVectors();
        TetradVector[] bias = stats.getBiasVectors();
        TetradMatrix gammas = stats.getPosteriorProbs();

        double[] weights = new double[K];

        TetradMatrix sourceVector;

        // find values of weights
        double sum;

        for (int i = 0; i < K; i++) {
            sum = 0;

            for (int r = 0; r < numCases; r++) {
                sum += gammas.get(r, i);
            }

            weights[i] = sum / numCases;
        }

        stats.setWeights(weights);

        sourceVector = new TetradMatrix(numCases, numVars);
        for (int k = 0; k < stats.getNumComponents(); k++) {
            final TetradMatrix inverse = W[k].transpose();

            // find values of source vectors
            for (int r = 0; r < numCases; r++) {
                TetradVector vector = dataMatrix.getRow(r).minus(bias[k]);
                TetradVector tempVector = inverse.times(vector);

                sourceVector.assignRow(r, tempVector);
            }

            S[k] = sourceVector;
        }

        stats.setSourceVectors(S);

//        System.out.println(sourceVector);
    }

    /*
     * Find posteriors of observations based on maximally likely values of mixing matrices, bias vectors, source vectors, and weights
     */
    private void expectation(DeterminingStats stats) {

        TetradMatrix[] A = stats.getMixingMatrices();
        TetradMatrix[] W = stats.getWMatrices();
        TetradMatrix[] K = stats.getKurtosisMatrices();
        TetradMatrix[] S = stats.getSourceVectors();
        TetradVector[] bias = stats.getBiasVectors();
        double[] weights = stats.getWeights();

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(numCases, stats.getNumComponents()); // Nxk
        double sum;
        double det;

        for (int k = 0; k < stats.getNumComponents(); k++) {
            det = Math.log(Math.abs(A[k].det()));

            for (int t = 0; t < numCases; t++) {
                sum = 0;

                for (int i = 0; i < numVars; i++) {
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

        // determine posteriors
        TetradMatrix gammas = new TetradMatrix(numCases, stats.getNumComponents()); // N x K

        for (int t = 0; t < numCases; t++) {
            double[] _weightedLIkelihoods = new double[stats.getNumComponents()];
            double _sum = 0.0;

            for (int c = 0; c < stats.getNumComponents(); c++) {
                _weightedLIkelihoods[c] = exp(likelihoods.get(t, c)) * weights[c];
                _sum += _weightedLIkelihoods[c];
            }

            for (int c = 0; c < stats.getNumComponents(); c++) {
                double gamma = _weightedLIkelihoods[c] / _sum;
                gammas.set(t, c, gamma);
            }
        }

        stats.setGammas(gammas);

        for (int t = 0; t < numCases; t++) {
            W = adaptMixingMatrices(stats, t);
            A = invertWMatrices(W);
            bias = adaptBiasVectors(stats, t);
        }

        stats.setMixingMatrices(A);
        stats.setWMatrices(W);
        stats.setBiasVectors(bias);
    }

    /*
     * This class holds all of the information that we're going to update as we go
     */
    private class DeterminingStats {
        private TetradMatrix[] mixingMatrices;
        private TetradMatrix[] wMatrices;
        private TetradVector[] biasVectors;
        private TetradMatrix[] sourceVectors;
        private TetradMatrix[] kurtosisMatrices;
        private TetradMatrix gammas;
        private double[] weights;
        private int K;

        void setMixingMatrices(TetradMatrix[] mixingMatrices) {
            this.mixingMatrices = mixingMatrices;
        }
        void setWMatrices(TetradMatrix[] wMatrices) {
            this.wMatrices = wMatrices;
        }
        void setBiasVectors(TetradVector[] biasVectors) {
            this.biasVectors = biasVectors;
        }
        void setSourceVectors(TetradMatrix[] sourceVectors) {
            this.sourceVectors = sourceVectors;
        }
        void setKurtosisMatrices(TetradMatrix[] kurtosisMatrices) {
            this.kurtosisMatrices = kurtosisMatrices;
        }
        void setGammas(TetradMatrix gammas) {
            this.gammas = gammas;
        }
        void setWeights(double[] weights) {
            this.weights = weights;
        }
        void setK(int K) {
            this.K = K;
        }
        TetradMatrix[] getMixingMatrices() {
            return mixingMatrices;
        }
        TetradMatrix[] getWMatrices() {
            return wMatrices;
        }
        TetradVector[] getBiasVectors() {
            return biasVectors;
        }
        TetradMatrix[] getSourceVectors() {
            return sourceVectors;
        }
        TetradMatrix[] getKurtosisMatrices() {
            return kurtosisMatrices;
        }
        TetradMatrix getPosteriorProbs() {
            return gammas;
        }
        double[] getWeights() {
            return weights;
        }
        int getNumComponents() {
            return K;
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

}
