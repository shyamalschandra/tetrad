package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.util.TetradMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by user on 7/21/18.
 */
public class MixtureModelNongaussian {

    private DataSet data;
    private TetradMatrix[] mixingMatrices;
    private TetradMatrix weights;
    private TetradMatrix[] bias;
    private int[] cases;
    private int[] caseCounts;
    private double[][] dataArray;
    private TetradMatrix posteriorProbabilities;

    public MixtureModelNongaussian(DataSet data, TetradMatrix posteriorProbabilities, TetradMatrix[] mixingMatrices,
                                   TetradMatrix[] biasVectors, TetradMatrix weights) {

        this.data = data;
        this.dataArray = data.getDoubleData().toArray();
        this.posteriorProbabilities = posteriorProbabilities;
        this.mixingMatrices = mixingMatrices;
        this.bias = biasVectors;
        this.weights = weights;

        this.cases = new int[data.getNumRows()];

        for (int t = 0; t < data.getNumRows(); t++) {
            cases[t] = argMax(t);
        }

        this.caseCounts = new int[weights.columns()];

        for (int i = 0; i < weights.columns(); i++) {
            caseCounts[i] = 0;
        }

        for (int i = 0; i < cases.length; i++) {
            for (int j = 0; j < weights.columns(); j++) {
                if (cases[i] == j) {
                    caseCounts[j]++;
                    break;
                }
            }
        }
    }

    /*
     * Classifies a given case into a model, based on which model has the highest gamma value for that case.
     */
    public int argMax(int t) {

        // hard classification
        int component = -1;
        double highest = -1;

        for (int k = 0; k < weights.columns(); k++) {
            final double prob = posteriorProbabilities.get(t, k) / weights.get(0, k);
            if (prob > highest) {
                highest = prob;
                component = k;
            }
        }

        return component;
    }

    public DataSet[] getDemixedData() {
        int k = weights.columns();

        List<List<Integer>> indices = new ArrayList<>();

        for (int _k = 0; _k < getMixingMatrices().length; _k++) {
            indices.add(new ArrayList<>());
        }

        System.out.println();

        for (int i = 0; i < cases.length; i++) {
            int _case = cases[i];
            System.out.print(_case);

            if (_case == -1) continue;

            indices.get(_case).add(i);
        }

        System.out.println();
        System.out.println();

        DoubleDataBox[] dataBoxes = new DoubleDataBox[mixingMatrices.length];

        for (int _k = 0; _k < mixingMatrices.length; _k++) {
            dataBoxes[_k] = new DoubleDataBox(indices.get(_k).size(), data.getNumColumns());

            for (int i = 0; i < indices.get(_k).size(); i++) {
                for (int j = 0; j < data.getNumColumns(); j++) {
                    dataBoxes[_k].set(i, j, data.getDouble(i, j));
                }
            }
        }

        DataSet[] dataSets = new DataSet[k];

        for (int i = 0; i < k; i++) {
            dataSets[i] = new BoxDataSet(dataBoxes[i], data.getVariables());
        }

        return dataSets;
    }

    public double[][] getData() {
        return dataArray;
    }

    public TetradMatrix getWeights() {
        return weights;
    }

    public int[] getCases() {
        return cases;
    }

    public TetradMatrix[] getMixingMatrices() { return mixingMatrices; }

    public TetradMatrix[] getBias() {return bias; }
}