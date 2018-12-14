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

import static java.lang.Math.*;

/**
 * @author Madelyn Glymour 7/21/18
 */
public class DemixerNongaussian {

    // The following fields are fixed.

    // The given data, N rows, numVars vars.
    private final DataSet data;

    // The data as an N x numVars matrix.
    private final TetradMatrix X;

    // Sample size
    private final int N;

    // The number of variables.
    private final int numVars;

    // The number of components. This must be specified up front.
    private final int numComponents;

    // Learning rate, by default 0.001.
    private final double learningRate;

    // The following fields are all updated.

    // S matrices - N x numVars, numComponents of them. This is updated.
    private final TetradMatrix[] S;

    // Diagonal kurtosis matrices, numComponents of them. This is updated.
    private final TetradMatrix[] K;

    // The weight of each component. This is updated.
    private final TetradMatrix weights;

    // Likelihood of each record, per component. This is updated but is rewritten every time and is only
    // made a field to save memory allocation.
    private final TetradMatrix likelihoods;

    // Posterior probability of each record, per component. This is updated.
    private final TetradMatrix posteriorProbs;

    // Biases. This is updated.
    private TetradMatrix[] bias;

    // Unmixing matrices. These are updated.
    private TetradMatrix[] W;

    private DemixerNongaussian(DataSet data, int numComponents, double learningRate) {

        // Fixed.
        this.data = data;
        this.X = data.getDoubleData();
        this.N = X.rows();
        this.numVars = X.columns();
        this.numComponents = numComponents;
        this.learningRate = learningRate;

        // Updated.
        this.bias = new TetradMatrix[this.numComponents];
        this.weights = new TetradMatrix(1, this.numComponents);
        this.likelihoods = new TetradMatrix(N, this.numComponents);
        this.posteriorProbs = new TetradMatrix(N, this.numComponents);
        this.W = new TetradMatrix[this.numComponents];
        this.K = new TetradMatrix[this.numComponents];
        this.S = new TetradMatrix[this.numComponents];
    }

    private MixtureModelNongaussian demix() {
        for (int k = 0; k < numComponents; k++) {
            weights.set(0, k,  RandomUtil.getInstance().nextUniform(0, 1));
        }

        normalize(weights, 0);

        // Initializing posterior probability of each record to be uniform.
        for (int n = 0; n < N; n++) {
            for (int k = 0; k < numComponents; k++) {
                posteriorProbs.set(n, k,  weights.get(0, k));
            }
        }

        FastIca ica = new FastIca(X, numVars);
        FastIca.IcaResult result = ica.findComponents();

        TetradMatrix _W = result.getW();

        for (int k = 0; k < numComponents; k++) {
            W[k] = _W;
            W[k] = whiten(_W, 0.2);
        }

        for (int k = 0; k < numComponents; k++) {
            TetradMatrix _bias = new TetradMatrix(1, numVars);

            for (int i = 0; i < numVars; i++) {
                _bias.set(0, i, RandomUtil.getInstance().nextNormal(0, 0.01));
            }

            bias[k] = _bias;
        }

        for (int k = 0; k < numComponents; k++) {
            adaptSourceMatrices(k);
        }

        initializeKurtosisMatrices();

        double[] _likelihoods = new double[numComponents];
        for (int k = 0; k < numComponents; k++) _likelihoods[k] = Double.NEGATIVE_INFINITY;

        boolean[] converged = new boolean[numComponents]; // All false.
        for (int k = 0; k < numComponents; k++) converged[k] = false;

        while (true) {
            boolean allConverged = tryConvergingUntilAllCoverged(_likelihoods, converged);
            maximization(); // each time.

            System.out.println("Weights: " + weights);

            if (allConverged) {
                break;
            }
        }

        return new MixtureModelNongaussian(data, posteriorProbs, invertWMatrices(W), bias, weights);
    }

    private boolean tryConvergingUntilAllCoverged(double[] _likelihoods, boolean[] converged) {
        for (int k = 0; k < numComponents; k++) {
            if (!converged[k]) {
                double l = expectation(k);

                // likelihood just for component k..
                if (l > _likelihoods[k]) {
                    _likelihoods[k] = l;
                } else {
                    converged[k] = true;
                }
            }
        }

        printLikelihoodfromLikelihoodArray();
        printLikelihoodFromPosteriorsArray();

        boolean allConverged = true;

        for (int k = 0; k < numComponents; k++) {
            if (!converged[k]) {
                allConverged = false;
            }
        }

        return allConverged;
    }

    private void printLikelihoodfromLikelihoodArray() {
        double l = 0;

        for (int n = 0; n < N; n++) {
            for (int k = 0; k < numComponents; k++) {
                l += likelihoods.get(n, k);
            }
        }

        System.out.println("L from likelihoods array = " + l);
    }

    private void printLikelihoodFromPosteriorsArray() {
        double l = 0;

        for (int n = 0; n < N; n++) {
            for (int k = 0; k < numComponents; k++) {
                l += log(posteriorProbs.get(n, k));
            }
        }

        System.out.println("L from posteriors array = " + l);
    }

    private void initializeKurtosisMatrices() {
        for (int k = 0; k < numComponents; k++) {
            if (K[k] == null) K[k] = new TetradMatrix(numVars, numVars);

            for (int v = 0; v < numVars; v++) {
                double sechSum = 0;
                double _sum = 0;
                double tanhSum = 0;

                for (int n = 0; n < N; n++) {
                    sechSum += 1.0 / pow(cosh(S[k].get(n, v)), 2.0);
                    _sum += pow(S[k].get(n, v), 2.0);
                    tanhSum += tanh(S[k].get(n, v)) * S[k].get(n, v);
                }

                sechSum /= N;
                _sum /= N;
                tanhSum /= N;

                double kurtosis = signum(sechSum * _sum - tanhSum);

                K[k].set(v, v, kurtosis);
            }
        }
    }


    /*
     * Inline updating of the mixing matrix A, updated every time a new posterior is determined
     */
    private void adaptMixingMatrices(int k) {

        TetradMatrix _tanhS = X.like(); // NxM

        for (int n = 0; n < N; n++) {
            for (int i = 0; i < numVars; i++) {
                _tanhS.set(n, i, Math.tanh(S[k].get(n, i)));
            }
        }

        TetradMatrix skskt = S[k].transpose().times((S[k])); // Mx1 * 1xM = MxM
        TetradMatrix tanhSk = _tanhS.transpose().times(S[k]); // Mx1 * 1xM = MxM
        TetradMatrix kurtosisSquare = K[k].times(tanhSk);
        TetradMatrix smketc = TetradMatrix.identity(numVars).minus(kurtosisSquare);
        TetradMatrix minusSquare = smketc.minus(skskt);

        TetradMatrix delta = new TetradMatrix(numVars, numVars);

        for (int n = 0; n < N; n++) {
            double p = posteriorProbs.get(n, k);
            delta = delta.plus(minusSquare.times(W[k]).scalarMult(p * learningRate / (N * numComponents)));
        }

        final double det = delta.det();

        if (abs(det) > 0) {
            W[k] = W[k].plus(delta);
            }
    }

    private void adaptBiasVectors(int k) {
        TetradMatrix _W = W[k];
        TetradMatrix sourceVector = S[k];
        TetradMatrix _bias = new TetradMatrix(1, numVars);

        for (int n = 0; n < N; n++) {
            double prob = posteriorProbs.get(n, k);

            double sum = 0;
            double det = Math.log(Math.abs(_W.det()));

            for (int i = 0; i < numVars; i++) {
                TetradMatrix bMatrix = _W.row(i);
                _bias = _bias.plus((bMatrix.scalarMult(Math.tanh(sourceVector.get(n, i)))).plus(bMatrix.scalarMult(sourceVector.get(n, i))));

                double tempSum = -.5 * Math.log(2.0 * Math.PI) - (K[k].get(i, i) * log(cosh(S[k].get(n, i))))
                        - (pow(S[k].get(n, i), 2.0) / 2.0);

                if (K[k].get(i, i) > 0) {
                    tempSum = tempSum - Math.log(0.7413);
                }

                sum += tempSum;
            }

            double L = (sum - det);

            _bias = _bias.scalarMult(prob * L);
        }

        bias[k] = bias[k].plus(_bias);
    }

    /*
     * Find posteriors of observations based on maximally likely values of mixing matrices, bias vectors,
     * source vectors, and weights
     */
    private double expectation(int k) {
        calculateLikelihoods(k);
        return calculatePosteriors(k);
    }

    private void maximization() {
        for (int k = 0; k < numComponents; k++) {
            adaptMixingMatrices(k);
        }

        for (int k = 0; k < numComponents; k++) {
            adaptBiasVectors(k);
        }

        adaptWeights();

        for (int k = 0; k < numComponents; k++) {
            adaptSourceMatrices(k);
        }
    }

    private void calculateLikelihoods(int k) {
        double det = Math.log(Math.abs(W[k].det()));

        if (det == 0 || Double.isNaN(det) || Double.isInfinite(det)) return;

        for (int n = 0; n < N; n++) {
            double sum = 0;

            for (int c = 0; c < numVars; c++) {
                double l = -.5 * Math.log(2.0 * Math.PI) - (K[k].get(c, c) * log(cosh(S[k].get(n, c))))
                        - (pow(S[k].get(n, c), 2.0) / 2.0);

                if (K[k].get(c, c) > 0) {
                    l = l - Math.log(0.7413);
                }

                sum += l;
            }

            double L = (sum - det);

            if (Double.isNaN(L) || Double.isInfinite(L)) {
                continue;
            }

            likelihoods.set(n, k, L);
        }
    }

//    // Returns the likelihood just for component k.
//    private double calculatePosteriors(int k) {
//        double likelihood = 0.0;
//
//        for (int n = 0; n < N; n++) {
//            double prob = exp(likelihoods.get(n, k));// * weights.get(0, k);
//            posteriorProbs.set(n, k, prob);
//            likelihood += likelihoods.get(n, k);
//        }
//
//        return likelihood;
//    }

    // Returns the likelihood just for component k.
    private double calculatePosteriors(int k) {
        double likelihood = 0.0;

        for (int n = 0; n < N; n++) {
            double prob = exp(likelihoods.get(n, k)) * weights.get(0, k);

            double divisor = 0;
            for (int i = 0; i < numComponents; i++){
                divisor += exp(likelihoods.get(n, i)) * weights.get(0, i);
            }

            prob = prob/divisor;

            posteriorProbs.set(n, k, prob);
            likelihood += likelihoods.get(n, k);
        }

        return likelihood;
    }


    private void adaptWeights() {
        double sum;

        for (int k = 0; k < numComponents; k++) {
            sum = 0;

            for (int n = 0; n < N; n++) {
                sum += posteriorProbs.get(n, k);
            }

            weights.set(0, k, sum);
        }

        normalize(weights, 0);
    }

    // Needs bias and W.
    private void adaptSourceMatrices(int k) {
        S[k] = X.minus(repmat(bias[k].getRow(0), N)).times(W[k].transpose());
    }

    public static void main(String... args) {

        DataSet dataSet = loadData("/Users/user/Downloads/mixfile1.csv");

        long startTime = System.currentTimeMillis();

        DemixerNongaussian pedro = new DemixerNongaussian(dataSet, 2, 0.001);
        MixtureModelNongaussian model = pedro.demix();

        long elapsed = System.currentTimeMillis() - startTime;

        System.out.println("\n================ WHAT I LEARNED ===================");

        for (int k = 0; k < model.getMixingMatrices().length; k++) {
            System.out.println("W = " + model.getMixingMatrices()[0]);
            System.out.println();
            System.out.println("Weights = " + model.getWeights());
            System.out.println();
        }

        DataSet[] datasets = model.getDemixedData();

        System.out.println("\n\nDatasets:");

        int count = 1;

        for (DataSet _dataSet : datasets) {
            System.out.println("#" + count++ + ": rows = " + _dataSet.getNumRows() + " cols = " + _dataSet.getNumColumns());
        }

        System.out.println();

        System.out.println("Elapsed: " + (elapsed / 1000));
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


    private TetradMatrix[] invertWMatrices(TetradMatrix[] m) {
        TetradMatrix[] n = new TetradMatrix[m.length];

        for (int k = 0; k < m.length; k++) {
            n[k] = m[k].inverse();
        }

        return n;
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
}