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

    public MixtureModelNongaussian demix(DataSet data, int numComponents) {

        double[] weights = new double[numComponents];
        TetradMatrix[] mixingMatrices = new TetradMatrix[numComponents];
        TetradVector[] biasVectors = new TetradVector[numComponents];
        TetradMatrix[] sourceVectors = new TetradMatrix[numComponents];
        TetradMatrix[] kurtosisMatrices = new TetradMatrix[numComponents];
        TetradMatrix sourceVectorMatrix;
        TetradMatrix kurtosisMatrix;
        TetradVector vector;

        FastIca ica = new FastIca(dataMatrix, numVars);
        FastIca.IcaResult result = ica.findComponents();
        TetradMatrix A = result.getA();

        for (int k = 0; k < numComponents; k++) {
            weights[k] = 1.0 / numComponents;

            sourceVectorMatrix = new TetradMatrix(numVars, numVars); // MxM
            vector = new TetradVector(numVars); // 1xM
            kurtosisMatrix = new TetradMatrix(numVars, numVars); // MxM

            for (int r = 0; r < A.rows(); r++) {
                for (int i = 0; i < A.columns(); i++) {
                    sourceVectorMatrix.set(r, i, A.get(r, i) + RandomUtil.getInstance().nextNormal(0, 1));
                }
            }

            for (int i = 0; i < numVars; i++) {
                vector.set(i, 0);
            }

            mixingMatrices[k] = sourceVectorMatrix;
            biasVectors[k] = vector;
            kurtosisMatrices[k] = kurtosisMatrix;

            sourceVectorMatrix = new TetradMatrix(numCases, numVars); // NxM

            for (int r = 0; r < numCases; r++) {
                TetradVector v = dataMatrix.getRow(r).minus(vector); // Mx1
                v = mixingMatrices[k].inverse().times(v); // MxM * Mx1 = Mx1
                sourceVectorMatrix.assignRow(r, v);
            }

            sourceVectors[k] = sourceVectorMatrix;
        }

        double sechSum;
        double sum;
        double tanhSum;
        double kurtosis;

        for (int k = 0; k < numComponents; k++) {
            kurtosisMatrix = kurtosisMatrices[k]; // MxM
            sourceVectorMatrix = sourceVectors[k]; // NxM

            for (int i = 0; i < numVars; i++) {
                sechSum = 0;
                sum = 0;
                tanhSum = 0;

                for (int r = 0; r < numCases; r++) {
                    sechSum += 1.0 / pow(cosh(sourceVectorMatrix.get(r, i)), 2.0);
                    sum += pow(sourceVectorMatrix.get(r, i), 2.0);
                    tanhSum += tanh(sourceVectorMatrix.get(r, i)) * sourceVectorMatrix.get(r, i);
                }

                sechSum /= numCases;
                sum /= numCases;
                tanhSum /= numCases;
                kurtosis = signum(sechSum * sum - tanhSum);

                kurtosisMatrix.set(i, i, kurtosis);
            }

            kurtosisMatrices[k] = kurtosisMatrix; // MxM
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
        stats.setK(numComponents);
        stats.setKurtosisMatrices(kurtosisMatrices);

        while (weightsUnequal) {
            expectation(stats);
            maximization(stats);

            System.out.println("Weights: " + Arrays.toString(stats.getWeights()));

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

        return new MixtureModelNongaussian(data, stats.getGammas(), stats.getMixingMatrices(),
                stats.getSourceVectors(), stats.getBiasVectors(), stats.getWeights());
    }

    /*
     * Inline updating of the mixing matrix A, updated every time a new posterior is determined
     */
    private TetradMatrix[] adaptMixingMatrices(DeterminingStats stats, int t) {

        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] kurtosisMatrices = stats.getKurtosisMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        TetradMatrix gammas = stats.getGammas();

        TetradMatrix mixingMatrix;
        TetradMatrix K;
        TetradMatrix sourceVector;
        TetradMatrix tanhSourceVector;
        double gamma;


        I:
        for (int k = 0; k < stats.getNumComponents(); k++) {
            mixingMatrix = mixingMatrices[k]; // MxM
            K = kurtosisMatrices[k]; // MxM
            sourceVector = sourceVectors[k]; // NxM
            gamma = gammas.get(t, k); // scalar
            tanhSourceVector = new TetradMatrix(numCases, numVars); // NxM

            for (int i = 0; i < numVars; i++) {
                tanhSourceVector.set(t, i, Math.tanh(sourceVector.get(t, i)));
            }

            TetradMatrix skskt = toColumnMatrix(sourceVector.getRow(t)).times(toRowMatrix(sourceVector.getRow(t))); // Mx1 * 1xM = MxM
            TetradMatrix tanhSk = toColumnMatrix(tanhSourceVector.getRow(t)).times(toRowMatrix(sourceVector.getRow(t))); // Mx1 * 1xM = MxM

            TetradMatrix kurtosisSquare = K.times(tanhSk);

            TetradMatrix smketc = TetradMatrix.identity(numVars).minus(kurtosisSquare);

            TetradMatrix minusSquare = smketc.minus(skskt);

            final double learningRate = 1e-3;
            TetradMatrix tempMixingMatrix = mixingMatrix.times((minusSquare)).scalarMult(-gamma * learningRate); // MxM * (MxM - (MxM * MxM - MxM)) = MxM

            System.out.println(tempMixingMatrix);

            //debugging
            for (int c1 = 0; c1 < mixingMatrices[k].columns(); c1++) {
                for (int c2 = 0; c2 < mixingMatrices[k].columns(); c2++) {
                    double val = mixingMatrices[k].get(c1, c2);

                    if (Double.isNaN(val)) {
                        throw new IllegalArgumentException("NaN in Mixing Matrix " + k);
//                        continue I;
                    } else if (Double.isInfinite(val)) {
                        throw new IllegalArgumentException("Inf in Mixing Matrix " + k);
//                        continue I;
                    }
                }
            }

            mixingMatrices[k] = mixingMatrix.plus(tempMixingMatrix);
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

        int K = stats.getNumComponents();
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
                    throw new IllegalArgumentException("NaN in Bias Vector " + I);
//                    continue I;
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

        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] kurtosisMatrices = stats.getKurtosisMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        double[] weights = stats.getWeights();

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(numCases, stats.getNumComponents()); // Nxk
        double sum;
        double det;

        for (int k = 0; k < stats.getNumComponents(); k++) {
            det = Math.log(Math.abs(mixingMatrices[k].det()));

            for (int t = 0; t < numCases; t++) {
                sum = 0;

                for (int i = 0; i < numVars; i++) {
                    double tempSum = (kurtosisMatrices[k].get(i, i) * log(cosh(sourceVectors[k].get(t, i)))
                            - pow(sourceVectors[k].get(t, i), 2.0) / 2.0);

//                    double d = log(cosh(sourceVectors[k].get(t, i)));
                    double d =  tempSum;

                    if (Double.isInfinite(d)) d = 0.0;
                    if (Double.isNaN(d)) d = 0.0;

                    if (Double.isNaN(d) || Double.isInfinite(d)) {
                        continue;
//                        throw new IllegalArgumentException("this thing = " + d + " sum = " + sum + " det = " + det);
                    }

                    sum += tempSum;
                }

                double likelihood = -sum - det;

                if (Double.isNaN(likelihood) || Double.isInfinite(likelihood)) {
                    throw new IllegalArgumentException("HEY!!! likelihood = " + likelihood + " sum = " + sum + " det = " + det);
                }

                likelihoods.set(t, k, likelihood);
            }
        }

        // determine posteriors
        double gamma;
        double divisor;
        TetradMatrix gammas = new TetradMatrix(numCases, stats.getNumComponents()); // N x K

        for (int r = 0; r < numCases; r++) {
            for (int c1 = 0; c1 < stats.getNumComponents(); c1++) {
                gamma = likelihoods.get(r, c1) * weights[c1];
                divisor = 0;

                for (int c2 = 0; c2 < stats.getNumComponents(); c2++) {
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

        public int getNumComponents() {
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
