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

package edu.cmu.tetrad.test;

import edu.cmu.tetrad.data.CovarianceMatrixOnTheFly;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.Fask;
import edu.cmu.tetrad.search.Fask_C;
import edu.cmu.tetrad.search.SemBicScore;
import edu.cmu.tetrad.util.*;
import edu.pitt.dbmi.data.Delimiter;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDataFileReader;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.*;
import static edu.cmu.tetrad.util.StatUtils.correlation;
import static java.lang.Math.*;

/**
 * Pulling this test out for Madelyn.
 *
 * @author jdramsey
 */
public class TestFask {

    private NormalDistribution normal = new NormalDistribution(0, 1);
    private double cutoff;

    @Test
    public void test1() {
        try {
            DataSet dataSet = loadData("/Users/user/Downloads/data2.txt", Delimiter.TAB, true);
            List<Node> variables = dataSet.getVariables();
            dataSet = DataUtils.standardizeData(dataSet);
            int N = dataSet.getNumRows();

            double[][] data = dataSet.getDoubleData().transpose().toArray();

            Graph graph = new EdgeListGraph(variables);

            double alpha = 0.00000001;

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    final double leftright1 = leftright(data[i], data[j]);

                    if (leftright1 > 0) {
                        System.out.println("Adding " + Edges.directedEdge(variables.get(i), variables.get(j)));
                        graph.addDirectedEdge(variables.get(i), variables.get(j));
                    }

//                    if (Double.isNaN(leftright1)) continue;;

                    System.out.println((i + 1) + " " + (j + 1)
                            + " lr " + leftright1);

                    for (int k = 0; k < 3; k++) {
                        if (i == k) continue;
                        final double leftright = leftright(data[i], data[j], new double[][]{data[i]});

                        if (leftright > 0) {
                            System.out.println("Adding " + Edges.directedEdge(variables.get(i), variables.get(j)));
                            graph.addDirectedEdge(variables.get(i), variables.get(j));
                        }

//                        if (Double.isNaN(leftright)) continue;;

                        System.out.println((i + 1) + " " + (j + 1) + " " + (k + 1)
                                + " lr " + leftright);
                    }
                }
            }

            System.out.println("\nGraph = " + graph);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private DataSet loadData(String path, Delimiter delimiter, boolean header) throws IOException {
        ContinuousTabularDataFileReader dataReader = new ContinuousTabularDataFileReader(new File(path), delimiter);
        dataReader.setHasHeader(header);

        return (DataSet) DataConvertUtils.toDataModel(dataReader.readInData());
    }


    @Test
    public void test2() {
        try {
            DataSet dataSet = loadData("/Users/user/Box/data/Sachs/data.logged.txt", Delimiter.TAB, true);
//            DataSet dataSet = loadData("/Users/user/Downloads/data1.txt", Delimiter.TAB, true);
            List<Node> variables = dataSet.getVariables();
            dataSet = DataUtils.center(dataSet);
            int N = dataSet.getNumRows();

            double[][] data = dataSet.getDoubleData().transpose().toArray();

            Graph graph = new EdgeListGraph(variables);

            for (int i = 0; i < variables.size(); i++) {
                for (int j = 0; j < variables.size(); j++) {
                    if (i == j) continue;
                    ;

                    Node X = variables.get(i);
                    Node Y = variables.get(j);

                    double[] x = data[i];
                    double[] y = data[j];

                    System.out.println("X = " + X + " Y = " + Y);

                    if (leftright(x, y) > 0 && leftright(y, x) < 0) {
                        graph.addDirectedEdge(X, Y);
                    }
                }
            }

            for (int i = 0; i < variables.size(); i++) {
                for (int j = 0; j < variables.size(); j++) {
                    if (i == j) continue;
                    ;

                    Node X = variables.get(i);
                    Node Y = variables.get(j);

                    if (graph.isAdjacentTo(X, Y)) continue;

                    double[] x = data[i];
                    double[] y = data[j];

                    List<Node> ZZ = graph.getParents(Y);
                    ZZ.remove(X);
                    ZZ.remove(Y);

                    for (Node Z : ZZ) {
                        double[] z = data[variables.indexOf(Z)];

                        System.out.println("X = " + X + " Y = " + Y + " | Z = " + Z);

                        if (leftright(x, y, z) > 0 && leftright(y, x, z) < 0) {
                            graph.addDirectedEdge(X, Y);
                        }
                    }
                }
            }

            System.out.println("\nGraph = " + graph);

            GraphUtils.saveGraph(graph, new File("/Users/user/Downloads/mygraph.txt"), false);


        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void test2a() {
        try {
            setCutoff(0.00001);
            double margin = .0;
            int depth = 5;

//            DataSet dataSet = loadData("/Users/user/Box/data/4cellLineData/4celllines.noindicators.txt", Delimiter.TAB, true);
//            DataSet dataSet = loadData("/Users/user/Box/data/Sachs/data.logged.txt", Delimiter.TAB, true);
            DataSet dataSet = loadData("/Users/user/Downloads/data2.txt", Delimiter.TAB, true);
            List<Node> variables = dataSet.getVariables();
            dataSet = DataUtils.center(dataSet);
            int N = dataSet.getNumRows();

            double[][] data = dataSet.getDoubleData().transpose().toArray();

            Graph graph = new EdgeListGraph(variables);

            {
                boolean added = true;

                for (int d = 0; d <= 0; d++) {
                    if (!added) break;

                    added = false;

                    for (int i = 0; i < variables.size(); i++) {
                        for (int j = 0; j < variables.size(); j++) {
                            if (i == j) continue;

                            Node X = variables.get(i);
                            Node Y = variables.get(j);

                            double[] x = data[i];
                            double[] y = data[j];

                            List<Node> ZZ = graph.getParents(Y);
                            ZZ.remove(X);

                            if (ZZ.size() < d) continue;

                            ChoiceGenerator gen = new ChoiceGenerator(ZZ.size(), d);
                            int[] choice;

                            while ((choice = gen.next()) != null) {
                                List<Node> Z = GraphUtils.asList(choice, ZZ);

                                double[][] z = new double[Z.size()][];

                                for (int t = 0; t < Z.size(); t++) {
                                    z[t] = data[variables.indexOf(Z.get(t))];
                                }

                                System.out.println("X = " + X + " Y = " + Y + " | Z = " + Z);

                                if (consistent(x, y, z, margin)) {
                                    graph.removeEdges(X, Y);
                                    graph.addDirectedEdge(X, Y);
                                    added = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            System.out.println("\n===Removing===");

//            {
//                boolean removed = true;
//
//                for (int d = 1; d <= depth; d++) {
//                    if (!removed) break;
//
//                    removed = false;
//
//                    for (Edge edge : graph.getEdges()) {
//                        if (graph.getEdges(edge.getNode1(), edge.getNode2()).size() == 2) continue;;
//
//                        Node X = Edges.getDirectedEdgeTail(edge);
//                        Node Y = Edges.getDirectedEdgeHead(edge);
//
//                        int i = variables.indexOf(X);
//                        int j = variables.indexOf(Y);
//
//                        double[] x = data[i];
//                        double[] y = data[j];
//
//                        List<Node> ZZ = graph.getParents(Y);
//                        ZZ.remove(X);
//
//                        if (ZZ.size() < d) continue;
//
//                        ChoiceGenerator gen = new ChoiceGenerator(ZZ.size(), d);
//                        int[] choice;
//
//                        while ((choice = gen.next()) != null) {
//                            List<Node> Z = GraphUtils.asList(choice, ZZ);
//
//                            double[][] z = new double[Z.size()][];
//
//                            for (int t = 0; t < Z.size(); t++) {
//                                z[t] = data[variables.indexOf(Z.get(t))];
//                            }
//
//                            System.out.println("X = " + X + " Y = " + Y + " | Z = " + Z);
//
//                            if (leftright(x, y, z) < margin || leftright(y, x, z) > -margin) {
//                                graph.removeEdges(X, Y);
//                                removed = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//            }

//            for (int i = 0; i < variables.size(); i++) {
//                for (int j = i + 1; j < variables.size(); j++) {
//                    Node X = variables.get(i);
//                    Node Y = variables.get(j);
//
//                    double[] x = data[i];
//                    double[] y = data[j];
//
//                    if (bidirected(x, y, graph, X, Y, dataSet, data)) {
//                        graph.addDirectedEdge(X, Y);
//                        graph.addDirectedEdge(Y, X);
//                    }
//                }
//            }

            System.out.println("\nGraph = " + graph);
            GraphUtils.saveGraph(graph, new File("/Users/user/Downloads/mygraph.txt"), false);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private boolean consistent(double[] x, double[] y, double[][] z, double margin) {
        return leftright(x, y, z) > margin && leftright(y, x, z) < -margin;
    }

    @Test
    public void test3() {
        try {
            DataSet dataSet = loadData("/Users/user/Box/data/Sachs/data.logged.txt", Delimiter.TAB, true);
            List<Node> variables = dataSet.getVariables();
            dataSet = DataUtils.standardizeData(dataSet);
            int N = 2000;//dataSet.getNumRows();

            double[][] data = dataSet.getDoubleData().transpose().toArray();

            Graph graph = new EdgeListGraph(variables);

            for (int i = 0; i < variables.size(); i++) {
                for (int j = 0; j < variables.size(); j++) {
                    if (i == j) continue;

                    System.out.println("i  = " + i + " j = " + j);

                    Node X = variables.get(i);
                    Node Y = variables.get(j);

                    System.out.println("\nX = " + X + " Y = " + Y);

                    // Centered
                    final double[] x = data[i];
                    final double[] y = data[j];

                    double alpha = .0000000000000000000001;

                    if (adj(x, y)) {
//                        double lrxy ;
//                        try {
//                            lrxy = leftright(x, y, N, alpha);
//                        } catch (Exception e) {
//                            lrxy = Double.NaN;
//                        }

                        if (leftright(x, y) > 0 && leftright(y, x) < 0) {
                            graph.addDirectedEdge(X, Y);
                        }
                    }

                    for (int k = 0; k < variables.size(); k++) {
                        if (k == i) continue;

                        if (!adj(y, data[k])) {
                            continue;
                        }

                        double lrxy;

                        try {
                            lrxy = leftright(x, y, data[k]);
                        } catch (Exception e) {
                            lrxy = Double.NaN;
                        }

                        if (lrxy > 0) {
                            graph.addDirectedEdge(X, Y);
                        }
                    }
                }
            }

            System.out.println(graph);

            GraphUtils.saveGraph(graph, new File("/Users/user/Downloads/mygraph.txt"), false);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test4() {
        try {
//            DataSet dataSet = loadData("/Users/user/Box/data/4cellLineData/4celllines.noindicators.txt", Delimiter.TAB, true);
//            DataSet dataSet = loadData("/Users/user/Box/data/Sachs/data.logged.txt", Delimiter.TAB, true);
            DataSet dataSet = loadData("/Users/user/Downloads/data2.txt", Delimiter.TAB, true);


//            Graph initial = GraphUtils.loadGraphTxt(new File("/Users/user/Box/data/4cellLineData/goldNet2_35N_no_circles.graph.txt"));

//            System.out.println(initial);

            Fask_C fask = new Fask_C(dataSet, new SemBicScore(new CovarianceMatrixOnTheFly(dataSet)));
//            fask.setInitialGraph(initial);
//            fask.setDelta(-.2);

            Graph graph = fask.search();

            System.out.println("\nGraph = " + graph);
            GraphUtils.saveGraph(graph, new File("/Users/user/Downloads/mygraph.txt"), false);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private boolean adj(double[] x, double[] y) {
        double c1 = StatUtils.cov(x, y, x, 0, +1)[0];
        double c2 = StatUtils.cov(x, y, y, 0, +1)[0];

        double d1 = (covariance(x, y) / variance(x)) * (StatUtils.cov(x, x, x, 0, +1)[0]
                - StatUtils.cov(x, x, y, 0, +1)[0]);
        double d2 = (covariance(x, y) / variance(y)) * (StatUtils.cov(y, y, y, 0, +1)[0]
                - StatUtils.cov(y, y, x, 0, +1)[0]);

        System.out.println("d1 = " + d1 + " d2 = " + d2);

        double d3 = Math.max(abs(d1), abs(d2));

        return (Math.abs(c1 - c2)) > d3;
    }

    private double leftright(double[] x, double[] y, double[]... z) {
        double[][] cond = new double[z.length + 1][];
        cond[0] = x;
        System.arraycopy(z, 0, cond, 1, z.length);
        double[] ry = residuals(y, cond);
        double a = covariance(x, y) / variance(x);
        final double lr = E(a, x, ry, y, +1) - E(a, x, ry, y, -1);
        System.out.println(lr);
        return lr;
    }

    private double getP(double N, double alpha, double lr, double r) {
        double _z = 0.5 * sqrt(N) * (log(1 + r) - log(1 - r));
        final double v = normal.cumulativeProbability(abs(_z));
        double p1 = (1.0 - v);
        System.out.println("  corr(x, ry^2) = " + r + " p = " + p1 + " lr = " + lr);
        System.out.println("  significant = " + (p1 < alpha));
        return p1;
    }


//    private double leftright(double[] x, double[] y, double[] w, double N, double alpha) {
//        x = Arrays.copyOf(x, x.length);
//        y = Arrays.copyOf(y, y.length);
//
//        double[] ry = residuals(y, new double[][]{x, w});
//
//        double a = covariance(x, y) / variance(x);
//        double lr = E(a, x, ry, y, -1) - E(a, x, ry, y, +1);
//
//        if (!Double.isNaN(lr)) {
//            System.out.println();
////            System.out.println("corr(x, ry) = " + correlation(x, ry));
////            System.out.println("  corr(x^2, ry) = " + correlation(square(x), ry));
//            final double r = correlation(x, square(ry));
//            double z = 0.5 * sqrt(N) * (log(1 + r) - log(1 - r));
//            final double v = normal.cumulativeProbability(abs(z));
//            double p = (1.0 - v);
//            System.out.println("  corr(x, ry^2) = " + r + " p = " + p + " lr = " + lr);
//            System.out.println("  significant = " + (p < alpha));
//            if (p < 0.05) lr = Double.NaN;
//        }
//
//        return lr;
//    }

    private double[] residuals(double[] _y, double[][] _x) {
        TetradMatrix y = new TetradMatrix(new double[][]{_y}).transpose();
        TetradMatrix x = new TetradMatrix(_x).transpose();

        TetradMatrix xT = x.transpose();
        TetradMatrix xTx = xT.times(x);
        TetradMatrix xTxInv = xTx.inverse();
        TetradMatrix xTy = xT.times(y);
        TetradMatrix b = xTxInv.times(xTy);

        TetradMatrix yHat = x.times(b);
        if (yHat.columns() == 0) yHat = y.copy();

        return y.minus(yHat).getColumn(0).toArray();
    }

    private static double E(double a, double[] x, double[] ry, double[] y, double dir) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            final double _x = x[k];
            final double _ry = ry[k];
            final double _y = abs(a) * _x + _ry;

            if (_x * dir > 0 && _y * dir < 0) {
                exy += _x * _ry;
                n++;
            }
        }

        return exy / n;
    }

    private double[] square(double[] x) {
        double[] square = Arrays.copyOf(x, x.length);
        for (int i = 0; i < square.length; i++) square[i] = square[i] * square[i];
        return square;
    }

    private double leftRightMinnesota(double[] x, double[] y) {
//        x = correctSkewness(x);
//        y = correctSkewness(y);

        final double cxyx = cov(x, y, x, 0, 1)[0];
        final double cxyy = cov(x, y, y, 0, 1)[0];
        final double cxxx = cov(x, x, x, 0, 1)[0];
        final double cyyx = cov(y, y, x, 0, 1)[0];
        final double cxxy = cov(x, x, y, 0, 1)[0];
        final double cyyy = cov(y, y, y, 0, 1)[0];

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;
        double b1 = cxyy / cyyy;
        double b2 = cxyx / cyyx;

        double Q = (a2 > 0) ? a1 / a2 : a2 / a1;
        double R = (b2 > 0) ? b1 / b2 : b2 / b1;

        double lr = Q - R;

//        if (StatUtils.correlation(x, y) < 0) lr += delta;

//        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));
//
//        if (sk_ey < 0) {
//            lr *= -1;
//        }

//        final double a = correlation(x, y);
//
//        if (a < 0 && sk_ey > -.2) {
//            lr *= -1;
//        }

        return lr;
    }

    public static double[] cov(double[] x, double[] y, double[] condition, double threshold, double direction) {
        double exy = 0.0;
        double exx = 0.0;
        double eyy = 0.0;

        double ex = 0.0;
        double ey = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (direction > threshold) {
                if (condition[k] > threshold) {
                    exy += x[k] * y[k];
                    exx += x[k] * x[k];
                    eyy += y[k] * y[k];
                    ex += x[k];
                    ey += y[k];
                    n++;
                }
            } else if (direction < threshold) {
                if (condition[k] > threshold) {
                    exy += x[k] * y[k];
                    exx += x[k] * x[k];
                    eyy += y[k] * y[k];
                    ex += x[k];
                    ey += y[k];
                    n++;
                }
            }
        }

        exy /= n;
        exx /= n;
        eyy /= n;
        ex /= n;
        ey /= n;

        double sxy = exy - ex * ey;
        double sx = exx - ex * ex;
        double sy = eyy - ey * ey;

        return new double[]{sxy, sxy / sqrt(sx * sy), sx, sy, (double) n, ex, ey, sxy / sx};
    }

    int depth = 5;

    private boolean bidirected(double[] x, double[] y, Graph G0, Node X, Node Y, DataSet dataSet, double[][] data) {

        Set<Node> adjSet = new HashSet<>(G0.getAdjacentNodes(X));
        adjSet.addAll(G0.getAdjacentNodes(Y));
        List<Node> adj = new ArrayList<>(adjSet);
        adj.remove(X);
        adj.remove(Y);

        DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), Math.min(depth, adj.size()));
        int[] choice;

        while ((choice = gen.next()) != null) {
            List<Node> _adj = GraphUtils.asList(choice, adj);
            double[][] _Z = new double[_adj.size()][];

            for (int f = 0; f < _adj.size(); f++) {
                Node _z = _adj.get(f);
                int column = dataSet.getColumn(_z);
                _Z[f] = data[column];
            }

            double pc = 0;
            double pc1 = 0;
            double pc2 = 0;

            try {
                pc = partialCorrelation(x, y, _Z, x, Double.NEGATIVE_INFINITY, +1);
                pc1 = partialCorrelation(x, y, _Z, x, 0, +1);
                pc2 = partialCorrelation(x, y, _Z, y, 0, +1);
            } catch (SingularMatrixException e) {
                System.out.println("Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                TetradLogger.getInstance().log("info", "Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                continue;
            }

            int nc = StatUtils.getRows(x, Double.NEGATIVE_INFINITY, +1).size();
            int nc1 = StatUtils.getRows(x, 0, +1).size();
            int nc2 = StatUtils.getRows(y, 0, +1).size();

            double z = 0.5 * (log(1.0 + pc) - log(1.0 - pc));
            double z1 = 0.5 * (log(1.0 + pc1) - log(1.0 - pc1));
            double z2 = 0.5 * (log(1.0 + pc2) - log(1.0 - pc2));

            double zv1 = (z - z1) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc1 - 3)));
            double zv2 = (z - z2) / sqrt((1.0 / ((double) nc - 3) + 1.0 / ((double) nc2 - 3)));

            boolean rejected1 = abs(zv1) > cutoff;
            boolean rejected2 = abs(zv2) > cutoff;

            boolean possibleTwoCycle = false;

            if (zv1 < 0 && zv2 > 0 && rejected1) {
                possibleTwoCycle = true;
            } else if (zv1 > 0 && zv2 < 0 && rejected2) {
                possibleTwoCycle = true;
            } else if (rejected1 && rejected2) {
                possibleTwoCycle = true;
            }

            if (!possibleTwoCycle) {
                return false;
            }
        }

        return true;
    }

    private double partialCorrelation(double[] x, double[] y, double[][] z, double[] condition, double threshold, double direction) throws SingularMatrixException {
        double[][] cv = StatUtils.covMatrix(x, y, z, condition, threshold, direction);
        TetradMatrix m = new TetradMatrix(cv).transpose();
        return StatUtils.partialCorrelation(m);
    }

    /**
     * Sets the significance level at which independence judgments should be made.  Affects the cutoff for partial
     * correlations to be considered statistically equal to zero.
     */
    private void setCutoff(double alpha) {
        if (alpha < 0.0 || alpha > 1.0) {
            throw new IllegalArgumentException("Significance out of range: " + alpha);
        }

        this.cutoff = StatUtils.getZForAlpha(alpha);
    }

    public static void main(String... args) {
        new TestFask().test1();
    }
}




