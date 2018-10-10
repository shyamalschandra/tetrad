package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;

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

    public MixtureModelNongaussian demix(DataSet data, int K) {

        double[] weights = new double[K];
        TetradMatrix[] mixingMatrices = new TetradMatrix[K];
        TetradVector[] biasVectors = new TetradVector[K];
        TetradMatrix[] sourceVectors = new TetradMatrix[K];
        TetradMatrix[] kurtosisMatrices = new TetradMatrix[K];
        TetradMatrix matrix;
        TetradMatrix kurtosisMatrix;
        TetradVector vector;

        FastIca ica = new FastIca(dataMatrix, numVars);
        FastIca.IcaResult result = ica.findComponents();
        TetradMatrix A = result.getA();

        for (int I = 0; I < K; I++) {
            weights[I] = 1.0 / K;

            matrix = new TetradMatrix(numVars, numVars); // MxM
            vector = new TetradVector(numVars); // 1xM
            kurtosisMatrix = new TetradMatrix(numVars, numVars); // MxM

            for (int r = 0; r < A.rows(); r++) {
                for (int c = 0; c < A.columns(); c++) {
                    matrix.set(r, c, A.get(r, c) + RandomUtil.getInstance().nextNormal(0, 1));
                }
            }

            for (int c = 0; c < A.columns(); c++) {
                vector.set(c, 0);
            }

            for (int c = 0; c < numVars; c++) {
                kurtosisMatrix.set(I, c, 0);
            }

            vector.set(I, 0);

            mixingMatrices[I] = matrix;
            biasVectors[I] = vector;
            kurtosisMatrices[I] = kurtosisMatrix;

            matrix = new TetradMatrix(numCases, numVars); // NxM

            for (int r = 0; r < numCases; r++) {
                TetradVector v = dataMatrix.getRow(r).minus(vector); // Mx1
                v = mixingMatrices[I].inverse().times(v); // MxM * Mx1 = Mx1
                matrix.assignRow(r, v);
            }

            sourceVectors[I] = matrix;
        }

        double sechSum;
        double sum;
        double tanhSum;
        double kurtosis;

        for (int I = 0; I < K; I++) {
            kurtosisMatrix = kurtosisMatrices[I]; // MxM
            matrix = sourceVectors[I]; // NxM

            for (int c = 0; c < numVars; c++) {
                sechSum = 0;
                sum = 0;
                tanhSum = 0;

                for (int r = 0; r < numCases; r++) {
                    sechSum += 1.0 / pow(cosh(matrix.get(r, c)), 2.0);
                    sum += pow(matrix.get(r, c), 2.0);
                    tanhSum += tanh(matrix.get(r, c)) * matrix.get(r, c);
                }

                sechSum /= numCases;
                sum /= numCases;
                tanhSum /= numCases;
                kurtosis = signum(sechSum * sum - tanhSum);

                kurtosisMatrix.set(I, c, kurtosis);
            }

            kurtosisMatrices[I] = kurtosisMatrix; // MxM
        }

        double[] tempWeights = Arrays.copyOf(weights, weights.length);

        boolean weightsUnequal = true;
        ArrayList<Double> diffsList;
        int iterCounter = 0;

        System.out.println("Weights: " + Arrays.toString(weights));

        DeterminingStats stats = new DeterminingStats();
        stats.setSourceVectors(sourceVectors);
        stats.setWeights(weights);
        stats.setMixingMatrices(mixingMatrices);
        stats.setBiasVectors(biasVectors);
        stats.setK(K);
        stats.setKurtosisMatrices(kurtosisMatrices);

        while (weightsUnequal) {
            expectation(stats);
            maximization(stats);

            System.out.println("Weights: " + Arrays.toString(stats.getWeights()));

            diffsList = new ArrayList<>();

            for (int I = 0; I < K; I++) {
                diffsList.add(abs(stats.getWeights()[I] - tempWeights[I]));
            }

            Collections.sort(diffsList);

            if (diffsList.get(K - 1) < 0.001 || iterCounter > 100) {
                weightsUnequal = false;
            }

            for (int I = 0; I < K; I++) {
                tempWeights[I] = stats.getWeights()[I];
            }

            iterCounter++;
        }

        return new MixtureModelNongaussian(data, stats.getGammas(), stats.getMixingMatrices(),
                stats.getSourceVectors(), stats.getBiasVectors(), stats.getWeights());
    }

    /*
     * Inline updating of the mixing matrix A, updated every time a new posterior is determined
     */
    private TetradMatrix[] adaptMixingMatrices(DeterminingStats stats, int r) {

        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] kurtosisMatrices = stats.getKurtosisMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        TetradMatrix gammas = stats.getGammas();

        TetradMatrix mixingMatrix;
        TetradMatrix kurtosisMatrix;
        TetradMatrix sourceVector;
        TetradMatrix tanhSourceVector;
        double gamma;

        I:
        for (int I = 0; I < mixingMatrices.length; I++) {
            mixingMatrix = mixingMatrices[I]; // MxM
            kurtosisMatrix = kurtosisMatrices[I]; // MxM
            sourceVector = sourceVectors[I]; // NxM
            gamma = gammas.get(r, I); // scalar
            tanhSourceVector = new TetradMatrix(numCases, numVars); // NxM

            for (int c = 0; c < numVars; c++) {
                tanhSourceVector.set(r, c, Math.tanh(sourceVector.get(r, c)));
            }

            TetradMatrix sourceSquare = toColumnMatrix(sourceVector.getRow(r)).times(toRowMatrix(sourceVector.getRow(r))); // Mx1 * 1xM = MxM
            TetradMatrix tanhSquare = toColumnMatrix(tanhSourceVector.getRow(r)).times(toRowMatrix(sourceVector.getRow(r))); // Mx1 * 1xM = MxM

            TetradMatrix kurtosisSquare = kurtosisMatrix.times(tanhSquare);

            TetradMatrix identitySquare = TetradMatrix.identity(numVars).minus(kurtosisSquare);

            TetradMatrix minusSquare = identitySquare.minus(sourceSquare);

            final double learningRate = 1e-8;
            TetradMatrix tempMixingMatrix = mixingMatrix.times((minusSquare)).scalarMult(-learningRate * gamma); // MxM * (MxM - (MxM * MxM - MxM)) = MxM

            System.out.println(tempMixingMatrix);

            //debugging
            for (int c1 = 0; c1 < mixingMatrices[I].columns(); c1++) {
                for (int c2 = 0; c2 < mixingMatrices[I].columns(); c2++) {
                    double val = mixingMatrices[I].get(c1, c2);

                    if (Double.isNaN(val)) {
                        System.out.println("NaN in Mixing Matrix " + I);
                        continue I;
                    } else if (Double.isInfinite(val)) {
                        System.out.println("Inf in Mixing Matrix " + I);
                        continue I;
                    }
                }
            }

            mixingMatrices[I] = mixingMatrix.plus(tempMixingMatrix);
        }

        return mixingMatrices;
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

        int K = stats.getK();
        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        TetradVector[] biasVectors = new TetradVector[K];
        TetradMatrix gammas = stats.getGammas();

        double[] weights = new double[K];

        TetradMatrix sourceVector;

        // find values of weights
        double sum;

        for (int I = 0; I < K; I++) {
            sum = 0;

            for (int r = 0; r < numCases; r++) {
                for (int c = 0; c < numVars; c++) {
                    sum += gammas.get(r, I);
                }
            }

            weights[I] = sum / numCases;

            //debugging
            if (Double.isNaN(weights[I])) {
                throw new NullPointerException("NaN in Weight " + Integer.toString(I));
            }
        }

        stats.setWeights(weights);

        // find values of bias vectors
        TetradVector biasVector;
        double biasVectorDenom;
        sourceVector = new TetradMatrix(numCases, numVars);

        I:
        for (int I = 0; I < K; I++) {
            biasVector = new TetradVector(numVars);
            for (int c = 0; c < numVars; c++) {
                biasVector.set(c, 0);
            }

            biasVectorDenom = 0;

            for (int r = 0; r < numCases; r++) {
                biasVector = biasVector.plus(dataMatrix.getRow(r)).scalarMult(gammas.get(r, I));
                biasVectorDenom += gammas.get(r, I);
            }

            biasVector = biasVector.scalarMult(1.0 / biasVectorDenom);

            //debugging
            for (int b = 0; b < biasVector.size(); b++) {
                double val = biasVector.get(b);

                if (Double.isNaN(val)) {
                    System.out.println("NaN in Bias Vector " + I);
                    continue I;
                }
            }

            // find values of source vectors
            for (int r = 0; r < numCases; r++) {
                TetradVector tempVector = mixingMatrices[I].transpose().times(dataMatrix.getRow(r).minus(biasVector));
                sourceVector.assignRow(r, tempVector);
            }

            //debugging
            for (int r = 0; r < sourceVector.rows(); r++) {
                for (int c = 0; c < sourceVector.columns(); c++) {
                    double val = sourceVector.get(r, c);
                    if (Double.isNaN(val)) {
                        System.out.println("NaN in Source Vector " + I);
                        continue I;
                    }
                }
            }

            biasVectors[I] = biasVector;
            sourceVectors[I] = sourceVector;
        }

        stats.setBiasVectors(biasVectors);
        stats.setSourceVectors(sourceVectors);
    }

    /*
     * Find posteriors of observations based on maximally likely values of mixing matrices, bias vectors, source vectors, and weights
     */
    private void expectation(DeterminingStats stats) {

        int K = stats.getK();
        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] kurtosisMatrices = stats.getKurtosisMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        double[] weights = stats.getWeights();

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(numCases, K); // Nxk
        double sum;
        double det;

        for (int I = 0; I < K; I++) {
            det = Math.log(Math.abs(mixingMatrices[I].det()));

            for (int r = 0; r < numCases; r++) {
                sum = 0;

                for (int c = 0; c < numVars; c++) {
                    double tempSum = (kurtosisMatrices[I].get(I, c) * log(cosh(sourceVectors[I].get(r, c)))
                            - pow(sourceVectors[I].get(I, c), 2.0) / 2.0);
                    sum += tempSum;
                }

                double likelihood = -sum - det;

                if (Double.isInfinite(likelihood)) {
                    System.out.println("HEY!!! likelihood = " + likelihood);
                    continue;
                } else if (Double.isNaN(likelihood)) {
                    System.out.println("HEY!!! likelihood = " + likelihood);
                    continue;
                }

                likelihoods.set(r, I, likelihood);
            }
        }

        // determine posteriors
        double gamma;
        double divisor;
        TetradMatrix gammas = new TetradMatrix(numCases, K); // N x K

        for (int r = 0; r < numCases; r++) {
            for (int c1 = 0; c1 < K; c1++) {
                gamma = likelihoods.get(r, c1) * weights[c1];
                divisor = 0;

                for (int c2 = 0; c2 < K; c2++) {
                    divisor += likelihoods.get(r, c2) * weights[c2];
                }

                //System.out.println(Double.toString(gamma) + "\t" + Double.toString(divisor));
                gamma = gamma / divisor;

                //debugging
                if (Double.isNaN(gamma)) {
                    System.out.println("NaN in Gamma " + Integer.toString(c1) + ", " + Integer.toString(r));
                    continue;
                }

                gammas.set(r, c1, gamma);
            }
        }

        stats.setGammas(gammas);

        for (int r = 0; r < numCases; r++) {
            mixingMatrices = adaptMixingMatrices(stats, r);
        }

        stats.setMixingMatrices(mixingMatrices);
    }

    /*
     * This class holds all of the information that we're going to update as we go
     */
    private class DeterminingStats {
        private TetradMatrix[] mixingMatrices;
        private TetradVector[] biasVectors;
        private TetradMatrix[] sourceVectors;
        private TetradMatrix[] kurtosisMatrices;
        private TetradMatrix gammas;
        private double[] weights;
        private int K;

        public DeterminingStats() {

        }

        public void setMixingMatrices(TetradMatrix[] mixingMatrices) {
            this.mixingMatrices = mixingMatrices;
        }

        public void setBiasVectors(TetradVector[] biasVectors) {
            this.biasVectors = biasVectors;
        }

        public void setSourceVectors(TetradMatrix[] sourceVectors) {
            this.sourceVectors = sourceVectors;
        }

        public void setKurtosisMatrices(TetradMatrix[] kurtosisMatrices) {
            this.kurtosisMatrices = kurtosisMatrices;
        }

        public void setGammas(TetradMatrix gammas) {
            this.gammas = gammas;
        }

        public void setWeights(double[] weights) {
            this.weights = weights;
        }

        public void setK(int K) {
            this.K = K;
        }

        public TetradMatrix[] getMixingMatrices() {
            return mixingMatrices;
        }

        public TetradVector[] getBiasVectors() {
            return biasVectors;
        }

        public TetradMatrix[] getSourceVectors() {
            return sourceVectors;
        }

        public TetradMatrix[] getKurtosisMatrices() {
            return kurtosisMatrices;
        }

        public TetradMatrix getGammas() {
            return gammas;
        }

        public double[] getWeights() {
            return weights;
        }

        public int getK() {
            return K;
        }
    }

    public static void main(String... args) {
        edu.cmu.tetrad.data.DataReader dataReader = new edu.cmu.tetrad.data.DataReader();
        dataReader.setVariablesSupplied(false);
        dataReader.setDelimiter(DelimiterType.WHITESPACE);

        DataSet dataSet = null;
        try {
            dataSet = dataReader.parseTabular(new File("/Users/user/Downloads/gamma_3000_4var.txt"));
        } catch (IOException e) {
            e.printStackTrace();
        }

        DemixerNongaussian pedro = new DemixerNongaussian(dataSet);
        long startTime = System.currentTimeMillis();
        MixtureModelNongaussian model = pedro.demix(dataSet, 3);
        long elapsed = System.currentTimeMillis() - startTime;

        double[] weights = model.getWeights();
        for (int i = 0; i < weights.length; i++) {
            System.out.print(Double.toString(weights[i]) + "\t");
        }

        DataSet[] datasets = model.getDemixedData();

        System.out.println("\n\nDatasets:");

        int count = 1;

        for (DataSet _dataSet : datasets) {
            System.out.println("#" + count++ + ": rows = " + _dataSet.getNumRows() + " cols = " + _dataSet.getNumColumns());
        }

        System.out.println();

        System.out.println("Elapsed: " + Long.toString(elapsed / 1000));
    }

}
