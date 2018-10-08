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
 * Created by user on 7/21/18.
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

    public MixtureModelNongaussian demix(DataSet data, int k) {

        double[] weights = new double[k];
        TetradMatrix[] mixingMatrices = new TetradMatrix[k];
        TetradVector[] biasVectors = new TetradVector[k];
        TetradMatrix[] sourceVectors = new TetradMatrix[k];
        TetradMatrix[] kurtosisMatrices = new TetradMatrix[k];
        TetradMatrix matrix;
        TetradMatrix kurtosisMatrix;
        TetradVector vector;

        FastIca ica = new FastIca(dataMatrix, numVars);
        FastIca.IcaResult result = ica.findComponents();
        TetradMatrix baseMatrix = result.getA();

        for (int i = 0; i < k; i++) {
            weights[i] = 1.0 / k;

            matrix = new TetradMatrix(numVars, numVars); // MxM
            vector = new TetradVector(numVars); // 1xM
            kurtosisMatrix = new TetradMatrix(numVars, numVars); // MxM
            for (int j = 0; j < numVars; j++) {
                for (int w = 0; w < numVars; w++) {
                    matrix.set(j, w, baseMatrix.get(j, w) + RandomUtil.getInstance().nextNormal(0, 1));
                    kurtosisMatrix.set(j, w, 0);
                }
                vector.set(j, 0);
            }

            // System.out.println(Double.toString(matrix.det()));

            mixingMatrices[i] = matrix;
            biasVectors[i] = vector;
            kurtosisMatrices[i] = kurtosisMatrix;

            matrix = new TetradMatrix(numCases, numVars); // NxM
            for (int j = 0; j < numCases; j++) {

                TetradVector v = dataMatrix.getRow(j).minus(vector); // Mx1

                v = mixingMatrices[i].inverse().times(v); // MxM * Mx1 = Mx1
                matrix.assignRow(j, v);
            }
            sourceVectors[i] = matrix;
        }

        double sechSum;
        double sum;
        double tanhSum;
        double kurtosis;
        for (int i = 0; i < k; i++) {
            kurtosisMatrix = kurtosisMatrices[i]; // MxM
            matrix = sourceVectors[i]; // NxM
            for (int j = 0; j < numVars; j++) {
                sechSum = 0;
                sum = 0;
                tanhSum = 0;

                for (int z = 0; z < numCases; z++) {
                    sechSum += pow(cosh(matrix.get(z, j)) * 4, 2.0);
                    sum += pow(matrix.get(z, j), 2.0);
                    tanhSum += tanh(matrix.get(z, j)) * matrix.get(z, j);
                }

                sechSum /= numCases;
                sum /= numCases;
                tanhSum /= numCases;
                kurtosis = signum(sechSum * sum - tanhSum);

                kurtosisMatrix.set(i, j, kurtosis);
            }

            kurtosisMatrices[i] = kurtosisMatrix; // MxM
        }

        double[] tempWeights = new double[k];

        for (int i = 0; i < k; i++) {
            tempWeights[i] = weights[i];
        }

        boolean weightsUnequal = true;
        ArrayList<Double> diffsList;
        int iterCounter = 0;

        System.out.println("Weights: " + Arrays.toString(weights));

        DeterminingStats stats = new DeterminingStats();
        stats.setSourceVectors(sourceVectors);
        stats.setWeights(weights);
        stats.setMixingMatrices(mixingMatrices);
        stats.setBiasVectors(biasVectors);
        stats.setK(k);
        stats.setKurtosisMatrices(kurtosisMatrices);

        while (weightsUnequal) {
            expectation(stats);
            maximization(stats);

            System.out.println("Weights: " + Arrays.toString(stats.getWeights()));


            diffsList = new ArrayList<>();
            for (int i = 0; i < k; i++) {
                diffsList.add(Math.abs(stats.getWeights()[i] - tempWeights[i]));
            }
            Collections.sort(diffsList);
            if (diffsList.get(k - 1) < 0.001 || iterCounter > 100) {
                weightsUnequal = false;
            }

            for (int i = 0; i < k; i++) {
                tempWeights[i] = stats.getWeights()[i];
            }

            iterCounter++;
        }

        MixtureModelNongaussian model = new MixtureModelNongaussian(data, stats.getGammas(), stats.getMixingMatrices(), stats.getSourceVectors(), stats.getBiasVectors(), stats.getWeights());
        //   demixed = true;

        return model;
    }

    /*
     * Inline updating of the mixing matrix A, updated every time a new posterior is determined
     */
    private TetradMatrix[] adaptMixingMatrices(DeterminingStats stats, int observation) {

        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] kurtosisMatrices = stats.getKurtosisMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        TetradMatrix gammas = stats.getGammas();

        TetradMatrix mixingMatrix;
        TetradMatrix kurtosisMatrix;
        TetradMatrix sourceVector;
        TetradMatrix tanhSourceVector;
        double gamma;

        for (int i = 0; i < mixingMatrices.length; i++) {
            mixingMatrix = mixingMatrices[i]; // MxM
            kurtosisMatrix = kurtosisMatrices[i]; // MxM
            sourceVector = sourceVectors[i]; // NxM
            gamma = gammas.get(observation, i); // scalar
            tanhSourceVector = new TetradMatrix(numCases, numVars); // NxM
            for (int j = 0; j < numCases; j++) {
                for (int w = 0; w < numVars; w++) {
                    tanhSourceVector.set(j, w, Math.tanh(sourceVector.get(j, w)));
                }
            }

            TetradMatrix sourceSquare = toColumnMatrix(sourceVector.getRow(observation)).times(toRowMatrix(sourceVector.getRow(observation))); // Mx1 * 1xM = MxM
            TetradMatrix tanhSquare = toColumnMatrix(tanhSourceVector.getRow(observation)).times(toRowMatrix(sourceVector.getRow(observation))); // Mx1 * 1xM = MxM

            TetradMatrix kurtosisSquare = kurtosisMatrix.times(tanhSquare);

            TetradMatrix identitySquare = TetradMatrix.identity(numVars).minus(kurtosisSquare);

            TetradMatrix minusSquare = identitySquare.minus(sourceSquare);

            TetradMatrix tempMixingMatrix = mixingMatrix.times((minusSquare)).scalarMult(-1 * gamma * .002); // MxM * (MxM - (MxM * MxM - MxM)) = MxM

            System.out.println(tempMixingMatrix);

            //debugging
            for (int r = 0; r < mixingMatrices[i].rows(); r++) {
                for (int c = 0; c < mixingMatrices[i].columns(); c++) {
                    double val = mixingMatrices[i].get(r, c);
                    if (Double.isNaN(val)) {
                        throw new IllegalArgumentException("NaN in Mixing Matrix " + i);
                    } else if (Double.isInfinite(val)) {
                        throw new IllegalArgumentException("Inf in Mixing Matrix " + i);
                    }
                }
            }

            mixingMatrices[i] = mixingMatrix.plus(tempMixingMatrix);
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
    private DeterminingStats maximization(DeterminingStats stats) {

        int k = stats.getK();
        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        TetradVector[] biasVectors = new TetradVector[k];
        TetradMatrix gammas = stats.getGammas();

        double[] weights = new double[k];

        TetradMatrix sourceVector;

        // find values of weights
        double sum;
        for (int i = 0; i < k; i++) {
            sum = 0;
            for (int j = 0; j < numCases; j++) {
                for (int w = 0; w < numVars; w++) {
                    sum += gammas.get(j, i);
                }
            }
            weights[i] = sum / numCases;

            //debugging
            if (Double.isNaN(weights[i])) {
                throw new NullPointerException("NaN in Weight " + Integer.toString(i));
            }
        }

        stats.setWeights(weights);

        // find values of bias vectors
        TetradVector biasVector;
        double biasVectorDenom;
        sourceVector = new TetradMatrix(numCases, numVars);

        for (int i = 0; i < k; i++) {
            biasVector = new TetradVector(numVars);
            for (int j = 0; j < numVars; j++) {
                biasVector.set(j, 0);
            }

            biasVectorDenom = 0;

            for (int j = 0; j < numCases; j++) {
                biasVector = biasVector.plus(dataMatrix.getRow(j)).scalarMult(gammas.get(j, i));
                biasVectorDenom += gammas.get(j, i);

            }

            biasVector = biasVector.scalarMult(1.0 / biasVectorDenom);

            //debugging
            for (int b = 0; b < biasVector.size(); b++) {
                double val = biasVector.get(b);
                if (Double.isNaN(val)) {
                    continue;
//                    throw new NullPointerException("NaN in Bias Vector " + Integer.toString(i));
                }
            }

            // find values of source vectors
            for (int j = 0; j < numCases; j++) {
                TetradVector tempVector = mixingMatrices[i].transpose().times(dataMatrix.getRow(j).minus(biasVector));
                sourceVector.assignRow(j, tempVector);
            }

            //debugging
            for (int r = 0; r < sourceVector.rows(); r++) {
                for (int c = 0; c < sourceVector.columns(); c++) {
                    double val = sourceVector.get(r, c);
                    if (Double.isNaN(val)) {
                        continue;
//                        throw new NullPointerException("NaN in Source Vector " + Integer.toString(i));
                    }
                }
            }

            biasVectors[i] = biasVector;
            sourceVectors[i] = sourceVector;
        }

        stats.setBiasVectors(biasVectors);
        stats.setSourceVectors(sourceVectors);

        return stats;
    }

    /*
     * Find posteriors of observations based on maximally likely values of mixing matrices, bias vectors, source vectors, and weights
     */
    private void expectation(DeterminingStats stats) {

        int k = stats.getK();
        TetradMatrix[] mixingMatrices = stats.getMixingMatrices();
        TetradMatrix[] kurtosisMatrices = stats.getKurtosisMatrices();
        TetradMatrix[] sourceVectors = stats.getSourceVectors();
        double[] weights = stats.getWeights();

        // determine log-likelihoods
        TetradMatrix likelihoods = new TetradMatrix(numCases, k); // Nxk
        double sum;
        double det;
        for (int w = 0; w < k; w++) {
            det = Math.log(Math.abs(mixingMatrices[w].det()));
            for (int i = 0; i < numCases; i++) {
                sum = 0;
                for (int j = 0; j < numVars; j++) {

                    double tempSum = (kurtosisMatrices[w].get(w, j) * log(cosh(sourceVectors[w].get(i, j))) - pow(sourceVectors[w].get(i, j), 2.0) / 2.0);

                    sum += tempSum;
                }

                double likelihood = -sum - det;

                if (Double.isInfinite(likelihood)) {
                    throw new IllegalArgumentException("HEY!!! likelihood = " + likelihood);
                } else if (Double.isNaN(likelihood)) {
                    throw new IllegalArgumentException("HEY!!! likelihood = " + likelihood);
                } else {
                    System.out.println("likelihood = " + likelihood + " det = " + det);
                }

                likelihoods.set(i, w, likelihood);
            }
        }

        // determine posteriors
        double gamma;
        double divisor;
        TetradMatrix gammas = new TetradMatrix(numCases, k); // Nxk
        for (int i = 0; i < numCases; i++) {
            for (int w = 0; w < k; w++) {
                gamma = likelihoods.get(i, w) * weights[w];
                divisor = 0;
                for (int j = 0; j < k; j++) {
                    divisor += likelihoods.get(i, j) * weights[j];
                }

                //System.out.println(Double.toString(gamma) + "\t" + Double.toString(divisor));
                gamma = gamma / divisor;
                //debugging
                if (Double.isNaN(gamma)) {
                    throw new NullPointerException("NaN in Gamma " + Integer.toString(w) + ", " + Integer.toString(i));
                }

                gammas.set(i, w, gamma);
            }
            stats.setGammas(gammas);

            mixingMatrices = adaptMixingMatrices(stats, i);
            stats.setMixingMatrices(mixingMatrices);
        }

        stats.setGammas(gammas);
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
        private int k;

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

        public void setK(int k) {
            this.k = k;
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
            return k;
        }
    }

    public static void main(String... args) {
        edu.cmu.tetrad.data.DataReader dataReader = new edu.cmu.tetrad.data.DataReader();
        dataReader.setVariablesSupplied(true);
        dataReader.setDelimiter(DelimiterType.TAB);

        DataSet dataSet = null;
        try {
            dataSet = dataReader.parseTabular(new File("/Users/user/Downloads/gamma_3000_4var.txt"));
        } catch (IOException e) {
            e.printStackTrace();
        }

        DemixerNongaussian pedro = new DemixerNongaussian(dataSet);
        long startTime = System.currentTimeMillis();
        MixtureModelNongaussian model = pedro.demix(dataSet, 2);
        long elapsed = System.currentTimeMillis() - startTime;

        double[] weights = model.getWeights();
        for (int i = 0; i < weights.length; i++) {
            System.out.print(Double.toString(weights[i]) + "\t");
        }

      /*  try {
            FileWriter writer = new FileWriter(new File("/Users/user/Documents/Demix_Testing/NonGaussian/bruce3_cases.txt"));
            BufferedWriter bufferedWriter = new BufferedWriter(writer);

            for (int i = 0; i < dataSet.getNumRows(); i++) {
                bufferedWriter.write(Integer.toString(model.getDistribution(i)) + "\n");
            }
            bufferedWriter.flush();
            bufferedWriter.close();

            DataSet[] dataSets = model.getDemixedData();

            for(int i = 0; i < dataSets.length; i++){
                writer = new FileWriter(new File("/Users/user/Documents/Demix_Testing/bruce3_demixed_" + Integer.toString(i+1) + ".txt"));
                bufferedWriter = new BufferedWriter(writer);
                bufferedWriter.write(dataSets[i].toString());
                bufferedWriter.flush();
                bufferedWriter.close();
            }
        }catch(Exception e){
            e.printStackTrace();
        } */

        System.out.println("Elapsed: " + Long.toString(elapsed / 1000));

    }

}
