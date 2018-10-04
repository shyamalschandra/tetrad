///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.TetradMatrix;

import java.util.*;

import static java.lang.Math.log;

/**
 * Implements a discrete Gaussian BIC score for FGS.
 *
 * @author Bryan Andrews
 */
public class DiscreteGaussianScore implements Score {

    private int N;

    private double structurePrior = 1.0;

    private DataSet mixedDataSet;

    private List<Node> mixedVariables;

    private DataSet dataSet;

    private TetradMatrix cov;

    private Map<Integer, List<Integer>> map;

    private static double LOG2PI = log(2.0 * Math.PI);

    /**
     * Constructs the score using a covariance matrix.
     */
    public DiscreteGaussianScore(DataSet dataSet, double structurePrior) {
        if (dataSet == null) {
            throw new NullPointerException();
        }

        this.N = dataSet.getNumRows();
        this.structurePrior = structurePrior;
        this.mixedDataSet = dataSet;
        this.mixedVariables = dataSet.getVariables();
        this.map = new HashMap<>();

        List<Node> A = new ArrayList<>();
        List<double[]> B = new ArrayList<>();

        int i = 0;
        int i_ = 0;
        while (i_ < dataSet.getNumColumns()) {

            Node v = dataSet.getVariable(i_);

            if (v instanceof DiscreteVariable) {

                Map<String, Integer> keys = new HashMap<>();
                for (int j = 0; j < N; j++) {
                    String key = v.getName().concat("_");
                    key = key.concat(Integer.toString(dataSet.getInt(j, i_)));
                    if (!keys.containsKey(key)) {
                        keys.put(key, i);
                        Node v_ = new ContinuousVariable(key);
                        A.add(v_);
                        B.add(new double[N]);
                        i++;
                    }
                    B.get(keys.get(key))[j] = 1;
                }

                i--;
                keys.remove(A.get(i).getName());
                A.remove(i);
                B.remove(i);

                List<Integer> index = new ArrayList<>();
                index.addAll(keys.values());
                map.put(i_, index);
                i_++;

            } else {

                A.add(v);
                double[] b = new double[N];
                for (int j = 0; j < N; j++) {
                    b[j] = dataSet.getDouble(j,i_);
                }
                B.add(b);

                List<Integer> index = new ArrayList<>();
                index.add(i);
                map.put(i_, index);
                i++;
                i_++;
            }
        }

        double[][] B_ = new double[N][B.size()];
        for (int j = 0; j < B.size(); j++) {
            for (int k = 0; k < N; k++) {
                B_[k][j] = B.get(j)[k];
            }
        }


        TetradMatrix D = new TetradMatrix(B_);
        TetradMatrix D_ = DataUtils.standardizeData(D);
        this.dataSet = new BoxDataSet(new DoubleDataBox(D_.toArray()),A);
        this.cov = this.dataSet.getCovarianceMatrix();

    }


    /**
     * Calculates the sample likelihood and BIC score.
     */
    public double localScore(int i, int... parents) {

        List<Integer> num = new ArrayList<>();
        List<Integer> den = new ArrayList<>();
        List current = map.get(i);
        for (int j = 0; j < current.size(); j++) {
            num.add((Integer) current.get(j));
        }
        for (int j : parents) {
            current = map.get(j);
            for (int k = 0; k < current.size(); k++) {
                num.add((Integer) current.get(k));
                den.add((Integer) current.get(k));
            }
        }
        int[] num_ = new int[num.size()];
        int[] den_ = new int[den.size()];
        for (int j = 0; j < num.size(); j++) {
            num_[j] = num.get(j);
        }
        for (int j = 0; j < den.size(); j++) {
            den_[j] = den.get(j);
        }

        double lik = getLik(num_, den_);
        int dof = getDof(num_.length, den_.length);

        return lik - dof / 2.0 * Math.log(N) + getStructurePrior(parents.length);
    }

    private double getLik(int[] num, int[] den) {

        double num_ = -N / 2.0 * log(cov.getSelection(num, num).det()) + num.length * (1 + LOG2PI);
        double den_ = -N / 2.0 * log(cov.getSelection(den, den).det()) + den.length * (1 + LOG2PI);

        return num_ - den_;
    }

    private int getDof(int num, int den) {

        int num_= (num * (num + 1)) / 2;
        int den_ = (den * (den + 1)) / 2;

        return num_ - den_;
    }

    private double getStructurePrior(int k) {
        if (structurePrior <= 0) {
            return 0;
        } else {
            int n = mixedVariables.size() - 1;
            double p = structurePrior / (double) n;
            return k * Math.log(p) + (n - k) * Math.log(1.0 - p);
        }
    }

//    private TetradMatrix cov(TetradMatrix x) {
//        return new TetradMatrix(new Covariance(x.getRealMatrix(), true).getCovarianceMatrix());
//    }

    private int[] append(int[] parents, int extra) {
        int[] all = new int[parents.length + 1];
        System.arraycopy(parents, 0, all, 0, parents.length);
        all[parents.length] = extra;
        return all;
    }

    public double localScoreDiff(int x, int y, int[] z) {
        return localScore(y, append(z, x)) - localScore(y, z);
    }

    @Override
    public double localScoreDiff(int x, int y) {
        return localScore(y, x) - localScore(y);
    }

    /**
     * Specialized scoring method for a single parent. Used to speed up the effect edges search.
     */
    public double localScore(int i, int parent) {
        return localScore(i, new int[]{parent});
    }

    /**
     * Specialized scoring method for no parents. Used to speed up the effect edges search.
     */
    public double localScore(int i) {
        return localScore(i, new int[0]);
    }

    public int getSampleSize() {
        return N;
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return bump > 0;
    }

    @Override
    public List<Node> getVariables() {
        return mixedVariables;
    }

    @Override
    public Node getVariable(String targetName) {
        for (Node node : mixedVariables) {
            if (node.getName().equals(targetName)) {
                return node;
            }
        }
        return null;
    }

    @Override
    public int getMaxDegree() {
        return (int) Math.ceil(Math.log(dataSet.getNumRows()));
    }

    @Override
    public boolean determines(List<Node> z, Node y) {
        return false;
    }

}



