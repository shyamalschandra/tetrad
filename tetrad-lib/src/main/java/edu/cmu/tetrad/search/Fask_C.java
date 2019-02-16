///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (c) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
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
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.util.*;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.io.File;
import java.util.List;
import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.max;
import static edu.cmu.tetrad.util.StatUtils.min;
import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Math.max;
import static java.lang.Math.*;

/**
 * Fast adjacency search followed by robust skew orientation. Checks are done for adding
 * two-cycles. The two-cycle checks do not require non-Gaussianity. The robust skew
 * orientation of edges left or right does.
 *
 * @author Joseph Ramsey
 */
public final class Fask_C implements GraphSearch {

    // The score to be used for the FAS adjacency search.
    private IndependenceTest test;

    // An initial graph to orient, skipping the adjacency step.
    private Graph initialGraph = null;

    // The colData sets being analyzed. They must all have the same variables and the same
    // number of records.
    private DataSet dataSet;

    // For the Fast Adjacency Search.
    private int depth = 1000;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // Data as a double[][].
    private final double[][] colData;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

    // True if skew adjacencies should be included in the output.
    private boolean useSkewAdjacencies = true;

    // A threshold for including extra adjacencies due to skewness.
    private double skewEdgeAlpha = 0.05;

    // Alpha for orienting 2-cycles.
    private double twoCycleAlpha = 0.05;

    // The list of variables.
    private final List<Node> variables;

    // True iff verbose output should be printed.
    private boolean verbose = false;
    private double cutoff;

    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public Fask_C(DataSet dataSet, Score score) {
//        final SemBicScore score = new SemBicScore(new CovarianceMatrixOnTheFly(dataSet));
//        score.setPenaltyDiscount(5);
        test = new IndTestScore(score);

        this.dataSet = dataSet;

        dataSet = DataUtils.center(dataSet);
        colData = dataSet.getDoubleData().transpose().toArray();
        this.variables = dataSet.getVariables();

        if (isVerbose()) {
            for (Node variable : variables) {
                TetradLogger.getInstance().forceLogMessage(variable + " skewness = " + skewness(variable));
            }
        }
    }

    //======================================== PUBLIC METHODS ====================================//

    /**
     * Runs the search on the concatenated data, returning a graph, possibly cyclic, possibly with
     * two-cycles. Runs the fast adjacency search (FAS, Spirtes et al., 2000) follows by a modification
     * of the robust skew rule (Pairwise Likelihood Ratios for Estimation of Non-Gaussian Structural
     * Equation Models, Smith and Hyvarinen), together with some heuristics for orienting two-cycles.
     *
     * @return the graph. Some of the edges may be undirected (though it shouldn't be many in most cases)
     * and some of the adjacencies may be two-cycles.
     */
    public Graph search() {
        setCutoff(twoCycleAlpha);

        TetradLogger.getInstance().forceLogMessage("\nStarting FASK-B Algorithm");

        TetradLogger.getInstance().forceLogMessage("\nSmoothly skewed:");

        for (Node node : dataSet.getVariables()) {
            if (smoothlySkewed(node)) {
                TetradLogger.getInstance().forceLogMessage(node.getName());
            }
        }

        TetradLogger.getInstance().forceLogMessage("\nNot smoothly skewed:");

        for (Node node : dataSet.getVariables()) {
            if (!smoothlySkewed(node)) {
                TetradLogger.getInstance().forceLogMessage(node.getName());
            }
        }

        TetradLogger.getInstance().forceLogMessage("");

        List<Node> variables = dataSet.getVariables();
        Graph graph = new EdgeListGraph(variables);
        dataSet = DataUtils.center(dataSet);
        double[][] data = dataSet.getDoubleData().transpose().toArray();

        {
            if (initialGraph == null) {

                FasStable fas = new FasStable(test);
                fas.setDepth(5);
                fas.setVerbose(false);
                fas.setKnowledge(knowledge);
                Graph graph2 = fas.search();

                List<Node> nodes = dataSet.getVariables();

                for (int i = 0; i < nodes.size(); i++) {
                    for (int j = i + 1; j < nodes.size(); j++) {
                        double[] x = data[i];
                        double[] y = data[j];

                        if (isAdj(x, y)) {
                            graph2.addUndirectedEdge(nodes.get(i), nodes.get(j));
                        }
                    }
                }

                initialGraph = graph2;
            } else {
                initialGraph = GraphUtils.undirectedGraph(initialGraph);
                initialGraph = GraphUtils.replaceNodes(initialGraph, dataSet.getVariables());

                if (initialGraph == null) throw new NullPointerException("Initial graph is null.");
            }
        }

//        for (Edge edge : initialGraph.getEdges()) {
//            if (RandomUtil.getInstance().nextDouble() > 0.5) {
//                initialGraph.removeEdge(edge);
//                initialGraph.addEdge(Edges.directedEdge(edge.getNode1(), edge.getNode2()));
//            } else {
//                initialGraph.removeEdge(edge);
//                initialGraph.addEdge(Edges.directedEdge(edge.getNode2(), edge.getNode1()));
//            }
//        }

        {

            Set<Node> changed1 = new HashSet<>(variables);
            Set<Node> changed2 = new HashSet<>(variables);
            for (Edge edge : initialGraph.getEdges()) graph.addEdge(edge);

            for (int d = 0; d < 20; d++) {
//            while (!changed1.isEmpty()) {

                if (changed1.isEmpty()) break;
//
                changed1 = changed2;
                changed2 = new HashSet<>();

                for (Edge edge : graph.getEdges()) {
                    Node X, Y;

                    if (Edges.isUndirectedEdge(edge)) {
                        X = edge.getNode1();
                        Y = edge.getNode2();
                    } else {
                        X = Edges.getDirectedEdgeTail(edge);
                        Y = Edges.getDirectedEdgeHead(edge);
                    }

                    if (!changed1.contains(Y)) continue;

                    if (edgeForbiddenByKnowledge(X, Y)) {
                        // Don't add an edge.
                    } else if (knowledgeOrients(X, Y)) {
                    } else if (knowledgeOrients(Y, X)) {
                    } else {

                        int i = variables.indexOf(X);
                        int j = variables.indexOf(Y);

                        double[] x = data[i];
                        double[] y = data[j];

                        final List<Node> Z = graph.getParents(Y);
                        Z.remove(X);

                        double[][] z = new double[Z.size()][];

                        for (int t = 0; t < Z.size(); t++) {
                            final Node V = Z.get(t);
                            z[t] = data[variables.indexOf(V)];
                        }

                        System.out.println("X = " + X + " Y = " + Y + " | Z = " + Z);

                        final boolean cxy = leftright(x, y, z) > 0;
                        final boolean cyx = leftright(y, x, z) > 0;

                        if (cxy && !cyx && !(graph.getEdges(X, Y).size() == 1 && graph.getEdge(X, Y).pointsTowards(Y))) {
                            graph.removeEdges(X, Y);
                            graph.addDirectedEdge(X, Y);
                            changed2.add(Y);
                        } else if (cyx && !cxy && !(graph.getEdges(X, Y).size() == 1 && graph.getEdge(Y, X).pointsTowards(X))) {
                            graph.removeEdges(Y, X);
                            graph.addDirectedEdge(Y, X);
                            changed2.add(X);
                        }
                    }
                }
            }
        }

        return graph;
    }

    private boolean isAdj(double[] x, double[] y) {
        double c1 = StatUtils.cov(x, y, x, 0, +1)[1];
        double c2 = StatUtils.cov(x, y, y, 0, +1)[1];

        double d1 = (covariance(x, y) / variance(x)) * (StatUtils.cov(x, x, x, 0, +1)[0]
                - StatUtils.cov(x, x, y, 0, +1)[0]);
        double d2 = (covariance(x, y) / variance(y)) * (StatUtils.cov(y, y, y, 0, +1)[0]
                - StatUtils.cov(y, y, x, 0, +1)[0]);

        System.out.println("d1 = " + d1 + " d2 = " + d2);

        double d3 = Math.max(abs(d1), abs(d2));

        return abs(c1 - c2) > .3;//  d3 + 0;
    }

    private double leftright(double[] x, double[] y, double[]... z) {
        double[][] cond = new double[z.length + 1][];
        cond[0] = x;
        System.arraycopy(z, 0, cond, 1, z.length);
        double[] ry = residuals(y, cond);
        double a = covariance(x, y) / variance(x);

        return E(a, x, ry, +1) - E(a, x, ry, -1);
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

    public Graph getInitialGraph() {
        return initialGraph;
    }

    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    public void setSkewEdgeAlpha(double skewEdgeAlpha) {
        this.skewEdgeAlpha = skewEdgeAlpha;
    }

    public boolean isUseFasAdjacencies() {
        return useFasAdjacencies;
    }

    public void setUseSkewAdjacencies(boolean useSkewAdjacencies) {
        this.useSkewAdjacencies = useSkewAdjacencies;
    }

    public void setUseFasAdjacencies(boolean useFasAdjacencies) {
        this.useFasAdjacencies = useFasAdjacencies;
    }

    public void setTwoCycleAlpha(double twoCycleAlpha) {
        this.twoCycleAlpha = twoCycleAlpha;
    }

    /**
     * @return The depth of search for the Fast Adjacency Search (FAS).
     */

    public int getDepth() {
        return depth;
    }

    /**
     * @param depth The depth of search for the Fast Adjacency Search (S). The default is -1.
     *              unlimited. Making this too high may results in statistical errors.
     */
    public void setDepth(int depth) {
        this.depth = depth;
    }

    /**
     * @return The elapsed time in milliseconds.
     */
    public long getElapsedTime() {
        return (long) 0;
    }

    /**
     * @return the current knowledge.
     */
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    /**
     * @param knowledge Knowledge of forbidden and required edges.
     */
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /////////////////////////////////////// PRIVATE METHODS ////////////////////////////////////

    private boolean smoothlySkewed(Node X) {
        double[] x = colData[variables.indexOf(X)];
        return smoothlySkewed(x);
    }

    private boolean smoothlySkewed(double[] x) {
        x = Arrays.copyOf(x, x.length);
        Arrays.sort(x);

        double min = min(x);
        double max = max(x);

        double _max = max(abs(min), abs(max));

        boolean smoothPositive = true;
        boolean smoothNegative = true;

        int numIntervals = 5;

        for (int i = 0; i < numIntervals; i++) {
            double t = (i * _max) / numIntervals;

            double a1 = Integrator.getArea(x1 -> 0, -t, 0, 5);
            double a2 = Integrator.getArea(x2 -> 0, 0, t, 5);

            if (a1 < a2) {
                smoothPositive = false;
            }

            if (a1 > a2) {
                smoothNegative = false;
            }
        }

        return smoothPositive || smoothNegative;
    }

    private double skewness(Node X) {
        double[] x = colData[variables.indexOf(X)];
        return StatUtils.skewness(x);
    }

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

    private static double E(double a, double[] x, double[] ry, double dir) {
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

    private boolean knowledgeOrients(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) || knowledge.isRequired(left.getName(), right.getName());
    }

    private boolean edgeForbiddenByKnowledge(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) && knowledge.isForbidden(left.getName(), right.getName());
    }
}







