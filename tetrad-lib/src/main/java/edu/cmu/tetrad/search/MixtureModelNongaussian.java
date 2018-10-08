package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.BoxDataSet;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DoubleDataBox;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;

/**
 * Created by user on 7/21/18.
 */
public class MixtureModelNongaussian {

    private DataSet data;
    private TetradMatrix gammas;
    private TetradMatrix[] mixingMatrices;
    private TetradMatrix[] sourceVectors;
    private TetradVector[] biasVectors;
    private double[] weights;
    private int[] cases;
    private int[] caseCounts;
    private double[][] dataArray;
    private double[][] gammaArray;

    public MixtureModelNongaussian(DataSet data, TetradMatrix gammas, TetradMatrix[] mixingMatrices, TetradMatrix[] sourceVectors, TetradVector[] biasVectors, double[] weights){

        this.data = data;
        this.dataArray = data.getDoubleData().toArray();
        this.gammas = gammas;
        this.gammaArray = gammas.toArray();
        this.mixingMatrices = mixingMatrices;
        this.sourceVectors = sourceVectors;
        this.biasVectors = biasVectors;
        this.weights = weights;

        this.cases = new int[data.getNumRows()];

        for(int i = 0; i < cases.length; i++){
            cases[i] = getDistribution(i);
        }

        this.caseCounts = new int[weights.length];

        for(int i = 0; i < weights.length; i++){
            caseCounts[i] = 0;
        }

        for(int i = 0; i < cases.length; i++){
            for(int j = 0; j < weights.length; j++){
                if(cases[i] == j){
                    caseCounts[j]++;
                    break;
                }
            }
        }
    }

    /*
     * Classifies a given case into a model, based on which model has the highest gamma value for that case.
     */
    public int getDistribution(int caseNum) {

        // hard classification
        int dist = 0;
        double highest = 0;
        for (int i = 0; i < weights.length; i++) {
            if (gammas.get(caseNum, i) > highest) {
                highest = gammas.get(caseNum, i);
                dist = i;
            }

        }

        return dist;
    }

    public DataSet[] getDemixedData(){
        int k = weights.length;
        DoubleDataBox[] dataBoxes = new DoubleDataBox[k];
        int[] caseIndices = new int[k];

        for(int i = 0; i < k; i++){
            dataBoxes[i] = new DoubleDataBox(caseCounts[i], data.getNumColumns());
            caseIndices[i] = 0;
        }

        int index;
        DoubleDataBox box;
        int count;
        for(int i = 0; i < cases.length; i++){
            index = cases[i];
            box = dataBoxes[index];
            count = caseIndices[index];
            for(int j = 0; j < data.getNumColumns(); j++){
                box.set(count, j, data.getDouble(i, j));
            }
            dataBoxes[index] = box;
            caseIndices[index] = count+1;
        }

        DataSet[] dataSets = new DataSet[k];
        for(int i = 0; i < k; i++){
            dataSets[i] = new BoxDataSet(dataBoxes[i],data.getVariables());
        }

        return dataSets;
    }

    public double[][] getData(){
        return dataArray;
    }

    public double[][] getGammas(){
        return gammaArray;
    }

    public double[] getWeights(){
        return weights;
    }

    public int[] getCaseCounts() { return caseCounts; }

    public int[] getCases() { return cases; }
}
