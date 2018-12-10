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

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.util.*;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.awt.*;
import java.util.*;
import java.util.List;

import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Math.*;

/**
 * Fast adjacency search followed by robust skew orientation. Checks are done for adding
 * two-cycles. The two-cycle checks do not require non-Gaussianity. The robust skew
 * orientation of edges left or right does.
 *
 * @author Joseph Ramsey
 */
public final class Fask_B implements GraphSearch {

    // The score to be used for the FAS adjacency search.
    private IndependenceTest test;

    // An initial graph to orient, skipping the adjacency step.
    private Graph initialGraph = null;

    // Elapsed time of the search, in milliseconds.
    private long elapsed = 0;

    // The colData sets being analyzed. They must all have the same variables and the same
    // number of records.
    private DataSet dataSet;

    // For the Fast Adjacency Search.
    private int depth = -1;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // Data as a double[][].
    private final double[][] colData;

    // A threshold for including extra adjacencies due to skewness.
    private double skewEdgeAlpha = 0.05;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

    // True if skew adjacencies should be included in the output.
    private boolean useSkewAdjacencies = true;

    // Alpha for orienting 2-cycles.
    private double twoCycleAlpha = 0.05;

    // Cutoff for orienting 2-cycles, calculated from twoCycleAlpha,
    private double twoCycleCutoff;

    // Regression procedure (linear) for conditioning.
    private RegressionDataset regressionDataset;

    // The list of variables.
    private final List<Node> variables;

    // The mask threshold, if used. If this is set to zero, it will not be used and ordinary conditional reasoning
    // will ensue. If set to non-zero, the mask will be applied, and then ordinary conditional reasoning will enssure.
    private double maskThreshold = 0.3;

    // For identifying small (negative) correlations beyond which judgment of direction should be reversed. In the
    // range (-1, 0).
    private double delta = -0.2;

    // True iff errors ey are assumed to be positively skewed. If false, a negative skew check will be performed.
    private boolean assumeErrorsPositivelySkewed = true;

    // True iff verbose output should be printed.
    private boolean verbose = false;

    // True if a mask should be use to estimate adjacencies.
    private boolean useMask = false;

    // True if the data is multiplicative; in this case, the judgment is reversed.
    private boolean multiplicative = false;

    // True iff skewnesses of each variable should be made positive by multiply by the signum of its skewness.
    private boolean correctSkews = true;

    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public Fask_B(DataSet dataSet, IndependenceTest test) {
        this.dataSet = dataSet;
        this.test = test;

        dataSet = DataUtils.center(dataSet);
        colData = dataSet.getDoubleData().transpose().toArray();
        this.variables = dataSet.getVariables();

        if (isVerbose()) {
            for (Node variable : variables) {
                TetradLogger.getInstance().forceLogMessage(variable + " skewness = " + skewness(variable));
            }
        }

        regressionDataset = new RegressionDataset(dataSet);
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
        long start = System.currentTimeMillis();

        TetradLogger.getInstance().forceLogMessage("\nStarting FASK-B Algorithm");

        setCutoff();

        Graph fasGraph = new EdgeListGraph(dataSet.getVariables());

        if (isUseFasAdjacencies()) {
            if (getInitialGraph() != null) {
                TetradLogger.getInstance().forceLogMessage("\nUsing initial graph.");
                System.out.println("Initial graph: " + initialGraph);

                Graph g1 = new EdgeListGraph(getInitialGraph().getNodes());

                for (Edge edge : getInitialGraph().getEdges()) {
                    Node X = edge.getNode1();
                    Node Y = edge.getNode2();

                    if (!g1.isAdjacentTo(X, Y)) g1.addUndirectedEdge(X, Y);
                }

                g1 = GraphUtils.replaceNodes(g1, dataSet.getVariables());
                Graph graph = new EdgeListGraph(g1.getNodes());

                for (Edge edge : g1.getEdges()) {
                    Node X = edge.getNode1();
                    Node Y = edge.getNode2();

                    if (getLrxy(X, Y) > leftRight(Y, X)) {
                        graph.addDirectedEdge(X, Y);
                    } else {
                        graph.addDirectedEdge(Y, X);
                    }
                }

                for (Edge edge : g1.getEdges()) {
                    Node X = edge.getNode1();
                    Node Y = edge.getNode2();

                    if (graph.isAdjacentTo(X, Y) && !edgeForbiddenByKnowledge(X, Y)
                            && !knowledgeOrients(X, Y) && !knowledgeOrients(Y, X)
                            && !Edges.isBidirectedEdge(graph.getEdge(X, Y))
                            && twocycle(X, Y, graph)) {
                        final double lrxy = leftRight(X, Y);
                        final double lryx = leftRight(Y, X);

                        if (!(lrxy < 0 && lryx < 0)) {
                            graph.removeEdges(X, Y);
                            graph.addDirectedEdge(X, Y);
                            graph.addDirectedEdge(Y, X);
                        }
                    }

                }

                return graph;
            } else {
                TetradLogger.getInstance().forceLogMessage("\nFAS");

                FasStable fas = new FasStable(test);
                fas.setDepth(getDepth());
                fas.setVerbose(false);
                fas.setKnowledge(knowledge);
                fasGraph = fas.search();
            }
        }

        TetradLogger.getInstance().forceLogMessage("\nOrienting required edges");
        SearchGraphUtils.pcOrientbk(knowledge, fasGraph, fasGraph.getNodes());

        Graph graph = new EdgeListGraph(variables);

        TetradLogger.getInstance().forceLogMessage("\nOrienting skeleton");
        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                if (edgeForbiddenByKnowledge(X, Y)) continue;

                final double[] x = colData[i];
                final double[] y = colData[j];

                double c1 = StatUtils.cov(x, y, x, 0, +1)[1];
                double c2 = StatUtils.cov(x, y, y, 0, +1)[1];

                if ((isUseFasAdjacencies() && fasGraph.isAdjacentTo(X, Y))
                        || ((isUseSkewAdjacencies() && (isUseMask() ?
                        Math.abs(c1 - c2) > getMaskThreshold() : skewAdjacent(X, Y, Collections.emptyList()))))) {
                    orientEdge(graph, X, Y);
                }
            }
        }

        System.out.println(graph);

        Graph graph2 = new EdgeListGraph(graph);

        TetradLogger.getInstance().forceLogMessage("\nRemoving extra edges by conditioning");
        removeExtraEdges(graph);

        TetradLogger.getInstance().forceLogMessage("\nOrienting 2-cycles");
        for (int i = 0; i < variables.size(); i++) {
            for (int j = 0; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                if (graph.isAdjacentTo(X, Y) && !edgeForbiddenByKnowledge(X, Y)
                        && !knowledgeOrients(X, Y) && !knowledgeOrients(Y, X)
                        && !Edges.isBidirectedEdge(graph.getEdge(X, Y))
                        && twocycle(X, Y, graph)) {
                    final double lrxy = leftRight(X, Y);
                    final double lryx = leftRight(Y, X);

                    if (!(lrxy < 0 && lryx < 0)) {
                        graph.removeEdges(X, Y);
                        graph.addDirectedEdge(X, Y);
                        graph.addDirectedEdge(Y, X);
                    }
                }
            }
        }

        TetradLogger.getInstance().forceLogMessage("\nAdding edges back in and re-running removal of extra edges (given 2-cycles)");

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                final double[] x = colData[i];
                final double[] y = colData[j];

                double c1 = StatUtils.cov(x, y, x, 0, +1)[1];
                double c2 = StatUtils.cov(x, y, y, 0, +1)[1];

                if (!graph.isAdjacentTo(X, Y) && graph2.isAdjacentTo(X, Y)) {
                    if (Edges.isBidirectedEdge(graph2.getEdge(X, Y))) continue;

                    if ((isUseFasAdjacencies() && fasGraph.isAdjacentTo(X, Y))
                            || (isUseSkewAdjacencies() && getMaskThreshold() != 0 ?
                            Math.abs(c1 - c2) > getMaskThreshold() : ((skewAdjacent(X, Y, Collections.emptyList()))))) {
                        orientEdge(graph, X, Y);
                    }
                }
            }
        }

        for (Edge edge : graph.getEdges()) {
            Node X = edge.getNode1();
            Node Y = edge.getNode2();
            final double[] x = colData[variables.indexOf(X)];
            final double[] y = colData[variables.indexOf(Y)];
            double r = correlation(x, y);
            if (abs(r) < 0.03) {
                edge.setLineColor(Color.ORANGE);
            }
        }

        removeExtraEdges(graph);

        TetradLogger.getInstance().forceLogMessage("\n\nFinal graph: \n\n" + graph);

        TetradLogger.getInstance().forceLogMessage("\nDone");

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

    private double getLrxy(Node x, Node y) {
        return leftRight(x, y);
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

    private boolean isUseSkewAdjacencies() {
        return useSkewAdjacencies;
    }

    public void setUseSkewAdjacencies(boolean useSkewAdjacencies) {
        this.useSkewAdjacencies = useSkewAdjacencies;
    }

    public void setUseFasAdjacencies(boolean useFasAdjacencies) {
        this.useFasAdjacencies = useFasAdjacencies;
    }

    public void setTwoCycleAlpha(double twoCycleAlpha) {
        this.twoCycleAlpha = twoCycleAlpha;
        setCutoff();
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
        return elapsed;
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

    public void setMaskThreshold(double maskThreshold) {
        this.maskThreshold = maskThreshold;
    }

    public double getDelta() {
        return delta;
    }

    public void setDelta(double delta) {
        this.delta = delta;
    }

//    public boolean isAssumeErrorsPositivelySkewed() {
//        return assumeErrorsPositivelySkewed;
//    }
//
//    public void setAssumeErrorsPositivelySkewed(boolean assumeErrorsPositivelySkewed) {
//        this.assumeErrorsPositivelySkewed = assumeErrorsPositivelySkewed;
//    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public boolean isUseMask() {
        return useMask;
    }

    public void setUseMask(boolean useMask) {
        this.useMask = useMask;
    }

    public void setMultiplicative(boolean multiplicative) {
        this.multiplicative = multiplicative;
    }

    public boolean isCorrectSkews() {
        return correctSkews;
    }

    public void setCorrectSkews(boolean correctSkews) {
        this.correctSkews = correctSkews;
    }

    /////////////////////////////////////// PRIVATE METHODS ////////////////////////////////////

    /**
     * Sets the significance level at which independence judgments should be made.  Affects the cutoff for partial
     * correlations to be considered statistically equal to zero.
     */
    private void setCutoff() {
        this.twoCycleCutoff = StatUtils.getZForAlpha(twoCycleAlpha);
    }

    private double getMaskThreshold() {
        return maskThreshold;
    }

    private void orientEdge(Graph graph, Node X, Node Y) {
        if (graph.getEdges(X, Y).size() == 2) return;
        graph.removeEdges(X, Y);

        if (knowledgeOrients(X, Y)) {
            graph.addDirectedEdge(X, Y);
        } else if (knowledgeOrients(Y, X)) {
            graph.addDirectedEdge(Y, X);
        } else {
            if (!edgeForbiddenByKnowledge(X, Y)) {
                if (knowledgeOrients(X, Y)) {
                    graph.addDirectedEdge(X, Y);
                } else if (knowledgeOrients(Y, X)) {
                    graph.addDirectedEdge(Y, X);
                } else {
                    if (leftRight(X, Y) > leftRight(Y, X)) {
                        graph.addDirectedEdge(X, Y);
                    } else {
                        graph.addDirectedEdge(Y, X);
                    }
                }
            }
        }
    }

    private void removeExtraEdges(Graph graph) {
        final int depth2 = this.depth == -1 ? 1000 : this.depth;
        boolean existsAnother = true;

        for (int d = 0; d < depth2; d++) {
            if (!existsAnother) break;
            existsAnother = false;

            for (Node head : variables) {
                List<Node> allParents1 = graph.getParents(head);
                List<Node> remove = new ArrayList<>();

                TAIL:
                for (Node parent : allParents1) {
                    Edge edge3 = graph.getDirectedEdge(parent, head);
                    Node tail = Edges.getDirectedEdgeTail(edge3);

                    if (!graph.isAdjacentTo(tail, head)) continue;

                    List<Node> c = getRelevantParents(graph, tail, head);

                    if (d > c.size()) continue;
                    existsAnother = true;

                    ChoiceGenerator gen = new ChoiceGenerator(c.size(), d);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        if (choice.length == 0) continue;
                        List<Node> Z = GraphUtils.asList(choice, c);

                        if (!skewAdjacent(tail, head, Z)) {
                            remove.add(tail);
                            TetradLogger.getInstance().forceLogMessage("Removing " + Edges.directedEdge(tail, head) + " | " + Z);
                            continue TAIL;
                        }
                    }
                }

                if (!new HashSet<>(remove).equals(new HashSet<>(allParents1))) {
                    for (Node tail1 : remove) {
                        graph.removeEdges(tail1, head);
                    }
                }
            }
        }
    }

    private List<Node> getRelevantParents(Graph graph, Node tail, Node head) {
        List<Node> parents1 = graph.getParents(head);
        parents1.addAll(graph.getParents(tail));
        parents1.remove(head);
        parents1.remove(tail);
        List<Node> c = new ArrayList<>();

        for (Node p1 : new ArrayList<>(parents1)) {

            if (graph.getEdges(p1, head).size() == 2) {
                continue;
            }

            if (graph.getEdges(p1, tail).size() == 2) {
                continue;
            }

            Edge edge = graph.getEdge(p1, head);

            if (edge != null && Edges.isUndirectedEdge(edge)) {
                continue;
            }

            if (edge != null && Edges.isBidirectedEdge(edge)) {
                continue;
            }

            Edge edge2 = graph.getEdge(p1, tail);

            if (edge2 != null && Edges.isUndirectedEdge(edge2)) {
                continue;
            }

            if (edge2 != null && Edges.isBidirectedEdge(edge2)) {
                continue;
            }

            if (!directedPathsFromTo(graph, tail, p1, 3).isEmpty() && !knowledge.isForbidden(p1.getName(), head.getName())) {
                if (!c.contains(p1)) {
                    c.add(p1);
                }
            }

            if (!directedPathsFromTo(graph, head, p1, 3).isEmpty() && !knowledge.isForbidden(p1.getName(), tail.getName())) {
                if (!c.contains(p1)) {
                    c.add(p1);
                }
            }
        }

        return c;
    }

    private double skewness(Node X) {
        double[] x = colData[variables.indexOf(X)];
        return StatUtils.skewness(x);
    }

    private boolean skewAdjacent(Node X, Node Y, List<Node> Z) {
        boolean b1 = false, b2 = false;

        try {
            Fask_B.E hx = new E(X, Y, Z, X).invoke();
            Fask_B.E hy = new E(X, Y, Z, Y).invoke();

            double[] dx = hx.getR();
            double[] dy = hy.getR();

            int nx = hx.getRows().size();
            int ny = hy.getRows().size();

            // Unequal variances, unequal sample sizes, T test, 2-sided
            double exyy = variance(dy) / ((double) ny);
            double exyx = variance(dx) / ((double) nx);
            double t = (mean(dy) - mean(dx)) / sqrt(exyy + exyx);
            double df = ((exyy + exyx) * (exyy + exyx)) / ((exyy * exyy) / (ny - 1)) + ((exyx * exyx) / (nx - 1));

            double p = new TDistribution(df).cumulativeProbability(t);
            b1 = p < skewEdgeAlpha;
        } catch (Exception e) {
            //
        }

        try {
            Fask_B.E hy = new E(Y, X, Z, Y).invoke();
            Fask_B.E hx = new E(Y, X, Z, X).invoke();

            double[] dx = hx.getR();
            double[] dy = hy.getR();

            int nx = hx.getRows().size();
            int ny = hy.getRows().size();

            // Unequal variances, unequal sample sizes, T test, 2-sided
            double exyy = variance(dy) / ((double) ny);
            double exyx = variance(dx) / ((double) nx);
            double t = (mean(dx) - mean(dy)) / sqrt(exyy + exyx);
            double df = ((exyy + exyx) * (exyy + exyx)) / ((exyy * exyy) / (ny - 1)) + ((exyx * exyx) / (nx - 1));

            double p = new TDistribution(df).cumulativeProbability(t);
            b2 = p < skewEdgeAlpha;
        } catch (Exception e) {
            //
        }

        return b1 || b2;
    }

    private double leftRight(Node X, Node Y) {
        double[] x = colData[variables.indexOf(X)];
        double[] y = colData[variables.indexOf(Y)];

        final double sx = StatUtils.skewness(x);
        final double sy = StatUtils.skewness(y);

        if (isCorrectSkews()) {
            x = times(x, signum(sx));
            y = times(y, signum(sy));
        }

        double lr = E(x, y, y, -1) / E(x, x, y, -1) - E(x, y, y, +1) / E(x, x, y, +1);

        if (correlation(x, y) < getDelta()) {
            lr *= -1;
        }

        if (isVerbose()) {
            TetradLogger.getInstance().forceLogMessage(
                    Edges.directedEdge(X, Y)
                            + " X = " + X.getName()
                            + " Y = " + Y.getName()
                            + " LR = " + lr
                            + " sx = " + sx
                            + " sy = " + sy
                            + " corr = " + correlation(x, y)
            );
        }

        return lr;
    }

    private double leftRight2(Node X, Node Y) {
        double[] x = colData[variables.indexOf(X)];
        double[] y = colData[variables.indexOf(Y)];

        final double sx = StatUtils.skewness(x);
        final double sy = StatUtils.skewness(y);

        if (isCorrectSkews()) {
            x = times(x, signum(sx));
            y = times(y, signum(sy));
        }

        final double left = E(x, y, x) / sqrt(E(x, x, x) * E(y, y, x));
        final double right = E(x, y, y) / sqrt(E(x, x, y) * E(y, y, y));

        double lr = (left - right);

        if (correlation(x, y) < getDelta()) {
            lr *= -1;
        }

        if (isVerbose()) {
            TetradLogger.getInstance().forceLogMessage(
                    Edges.directedEdge(X, Y)
                            + " X = " + X.getName()
                            + " Y = " + Y.getName()
                            + " LR = " + lr
                            + " sx = " + sx
                            + " sy = " + sy
                            + " corr = " + correlation(x, y)
            );
        }

        return lr;
    }

    private double[] times(double[] data, double r) {
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * r;
        return data2;
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

    private static double E(double[] x, double[] y, double[] condition) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (condition[k] > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
    }

    private double E(double[] x, double[] y, double[] condition, double dir) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (dir * condition[k] > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
    }

    private static int N(double[] condition) {
        int n = 0;

        for (int k = 0; k < condition.length; k++) {
            if (condition[k] > 0) {
                n++;
            }
        }

        return n;
    }

    private boolean twocycle(Node X, Node Y, Graph graph) {
        if (multiplicative) return false;

        double[] x = colData[variables.indexOf(X)];
        double[] y = colData[variables.indexOf(Y)];

        final List<Node> adjX = graph.getAdjacentNodes(X);
        final List<Node> adjY = graph.getAdjacentNodes(Y);

        for (Node a : new ArrayList<>(adjX)) {
            if (graph.getEdges(a, X).size() == 2 || Edges.isUndirectedEdge(graph.getEdge(a, X))) {
                adjX.remove(a);
            }
        }

        for (Node a : new ArrayList<>(adjY)) {
            if (graph.getEdges(a, Y).size() == 2 || Edges.isUndirectedEdge(graph.getEdge(a, Y))) {
                adjY.remove(a);
            }
        }

        Set<Node> adjSet = new HashSet<>(adjX);
        adjSet.addAll(adjY);
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
                _Z[f] = colData[column];
            }

            double r;
            double rx;
            double ry;

            try {
                r = partialCorrelation(x, y, _Z, x, Double.NEGATIVE_INFINITY, +1);
                rx = partialCorrelation(x, y, _Z, x, 0, +1);
                ry = partialCorrelation(x, y, _Z, y, 0, +1);
            } catch (SingularMatrixException e) {
                System.out.println("Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                TetradLogger.getInstance().log("info", "Singularity X = " + X + " Y = " + Y + " adj = " + adj);
                continue;
            }

            int n = StatUtils.getRows(x, Double.NEGATIVE_INFINITY, +1).size();
            int nx = StatUtils.getRows(x, 0, +1).size();
            int ny = StatUtils.getRows(y, 0, +1).size();

            double z = 0.5 * (log(1.0 + r) - log(1.0 - r));
            double zx = 0.5 * (log(1.0 + rx) - log(1.0 - rx));
            double zy = 0.5 * (log(1.0 + ry) - log(1.0 - ry));

            double zvx = (z - zx) / sqrt((1.0 / ((double) n - 3) + 1.0 / ((double) nx - 3)));
            double zvy = (z - zy) / sqrt((1.0 / ((double) n - 3) + 1.0 / ((double) ny - 3)));

            boolean rejectedx = abs(zvx) > twoCycleCutoff;
            boolean rejectedy = abs(zvy) > twoCycleCutoff;

            boolean possibleCycle = false;

            if (zvx < 0 && zvy > 0 && rejectedx) {
                possibleCycle = true;
            } else if (zvx > 0 && zvy < 0 && rejectedy) {
                possibleCycle = true;
            } else if (rejectedx && rejectedy) {
                possibleCycle = true;
            }

            if (!possibleCycle) {
                return false;
            }
        }

        return true;
    }

    private double partialCorrelation(double[] x, double[] y, double[][] z, double[] condition,
                                      double threshold, double direction) throws SingularMatrixException {
        double[][] cv = covMatrix(x, y, z, condition, threshold, direction);
        TetradMatrix m = new TetradMatrix(cv).transpose();
        return StatUtils.partialCorrelation(m);
    }

    private boolean knowledgeOrients(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) || knowledge.isRequired(left.getName(), right.getName());
    }

    private boolean edgeForbiddenByKnowledge(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) && knowledge.isForbidden(left.getName(), right.getName());
    }

    private static double[][] covMatrix(double[] x, double[] y, double[][] z, double[] condition,
                                        double threshold, double direction) {
        List<Integer> rows = getRows(x, condition, threshold, direction);

        double[][] allData = new double[z.length + 2][];

        allData[0] = x;
        allData[1] = y;

        System.arraycopy(z, 0, allData, 2, z.length);

        double[][] subdata = new double[allData.length][rows.size()];

        for (int c = 0; c < allData.length; c++) {
            for (int i = 0; i < rows.size(); i++) {
                try {
                    subdata[c][i] = allData[c][rows.get(i)];
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        double[][] cov = new double[z.length + 2][z.length + 2];

        for (int i = 0; i < z.length + 2; i++) {
            for (int j = i; j < z.length + 2; j++) {
                double c = StatUtils.covariance(subdata[i], subdata[j]);
                cov[i][j] = c;
                cov[j][i] = c;
            }
        }

        return cov;
    }

    private static List<Integer> getRows(double[] x, double[] condition, double threshold, double direction) {
        List<Integer> rows = new ArrayList<>();

        for (int k = 0; k < x.length; k++) {
            if (direction > threshold) {
                if (condition[k] > threshold) {
                    rows.add(k);
                }
            } else if (direction < threshold) {
                if (condition[k] > threshold) {
                    rows.add(k);
                }
            } else {
                if (condition[k] > threshold) {
                    rows.add(k);
                }
            }
        }

        return rows;
    }

    private class E {
        private Node x;
        private Node y;
        private List<Node> z;
        private Node condition;
        private List<Integer> rows;
        private double[] rxy_over_erxx;

        public E(Node X, Node Y, List<Node> Z, Node condition) {
            x = X;
            y = Y;
            z = Z;
            this.condition = condition;
        }

        public List<Integer> getRows() {
            return rows;
        }

        public double[] getR() {
            return rxy_over_erxx;
        }

        public Fask_B.E invoke() {
            final double[] _w = colData[variables.indexOf(condition)];

            rows = StatUtils.getRows(_w, 0, +1);

            int[] _rows = new int[rows.size()];
            for (int i = 0; i < rows.size(); i++) _rows[i] = rows.get(i);

            regressionDataset.setRows(_rows);
            double[] rx = regressionDataset.regress(x, z).getResiduals().toArray();
            double[] ry = regressionDataset.regress(y, z).getResiduals().toArray();

            double[] rxy = new double[rows.size()];

            for (int i = 0; i < rows.size(); i++) {
                rxy[i] = rx[i] * ry[i];
            }

            double[] rxx = new double[rows.size()];

            for (int i = 0; i < rows.size(); i++) {
                rxx[i] = rx[i] * rx[i];
            }

            rxy_over_erxx = Arrays.copyOf(rxy, rxy.length);
            double erxx = mean(rxx);

            for (int i = 0; i < rxy_over_erxx.length; i++) rxy_over_erxx[i] /= erxx;

            return this;
        }

    }

    private static List<List<Node>> directedPathsFromTo(Graph graph, Node node1, Node node2, int maxLength) {
        List<List<Node>> paths = new LinkedList<>();
        directedPathsFromToVisit(graph, node1, node2, new LinkedList<Node>(), paths, maxLength);
        return paths;
    }

    private static void directedPathsFromToVisit(Graph graph, Node node1, Node node2,
                                                 LinkedList<Node> path, List<List<Node>> paths, int maxLength) {
        if (maxLength != -1 && path.size() > maxLength - 2) {
            return;
        }

        int witnessed = 0;

        for (Node node : path) {
            if (node == node1) {
                witnessed++;
            }
        }

        if (witnessed > 1) {
            return;
        }

        path.addLast(node1);

        for (Edge edge : graph.getEdges(node1)) {
            Node child = traverseDirected(node1, edge);

            if (child == null) {
                continue;
            }

            if (child == node2) {
                LinkedList<Node> _path = new LinkedList<>(path);
                _path.add(child);
                paths.add(_path);
                continue;
            }

            if (path.contains(child)) {
                continue;
            }

            directedPathsFromToVisit(graph, child, node2, path, paths, maxLength);
        }

        path.removeLast();
    }

    /**
     * For A -> B, given A, returns B; otherwise returns null.
     */
    private static Node traverseDirected(Node node, Edge edge) {
        if (node == edge.getNode1()) {
            if ((edge.getEndpoint1() == Endpoint.TAIL || edge.getEndpoint1() == Endpoint.CIRCLE) &&
                    (edge.getEndpoint2() == Endpoint.ARROW)) {
                return edge.getNode2();
            }
        } else if (node == edge.getNode2()) {
            if ((edge.getEndpoint2() == Endpoint.TAIL || edge.getEndpoint2() == Endpoint.CIRCLE) &&
                    (edge.getEndpoint1() == Endpoint.ARROW)) {
                return edge.getNode1();
            }
        }

        return null;
    }
}







