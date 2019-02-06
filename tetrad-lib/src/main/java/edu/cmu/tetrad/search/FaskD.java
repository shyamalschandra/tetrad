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
import org.apache.commons.math3.linear.SingularMatrixException;

import java.util.*;

import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Math.*;

/**
 * Fast adjacency search followed by robust skew orientation. Checks are done for adding
 * two-cycles. The two-cycle checks do not require non-Gaussianity. The robust skew
 * orientation of edges left or right does.
 *
 * @author Joseph Ramsey
 */
public final class FaskD implements GraphSearch {

    // The score to be used for the FAS adjacency search.
    private final Score score;
    private final List<Node> variables;

    // An initial graph to orient, skipping the adjacency step.
    private Graph initialGraph = null;

    // Elapsed time of the search, in milliseconds.
    private long elapsed = 0;

    // The data sets being analyzed. They must all have the same variables and the same
    // number of records.
    private DataSet dataSet = null;

    // For the Fast Adjacency Search.
    private int depth = -1;

    // For the SEM BIC score, for the Fast Adjacency Search.
    private double penaltyDiscount = 1;

    // Alpha for orienting 2-cycles. Usually needs to be low.
    private double alpha = 1e-6;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // Data as a double[][].
    private final double[][] data;

    // Cutoff for T tests for 2-cycle tests.
    private double cutoff;

    // A threshold for including extra adjacencies due to skewness.
    private double skewEdgeThreshold = 0.3;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

    // True if skew adjacencies should be included in the output.
    private boolean useSkewAdjacencies = true;

//    // Threshold for reversing casual judgments for negative coefficients.
//    private double delta = -0.2;

    private RegressionDataset regressionDataset;


    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public FaskD(DataSet dataSet, Score score) {
        this.dataSet = dataSet;
        this.score = score;

        data = dataSet.getDoubleData().transpose().toArray();
        this.variables = dataSet.getVariables();

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

        setCutoff(alpha);

        dataSet = DataUtils.getNonparanormalTransformed(dataSet);
//        DataSet dataSet = DataUtils.standardizeData(this.dataSet);

        List<Node> variables = dataSet.getVariables();
        double[][] colData = dataSet.getDoubleData().transpose().toArray();
        Graph G0;

        if (getInitialGraph() != null) {
            Graph g1 = new EdgeListGraph(getInitialGraph().getNodes());

            for (Edge edge : getInitialGraph().getEdges()) {
                Node x = edge.getNode1();
                Node y = edge.getNode2();

                if (!g1.isAdjacentTo(x, y)) g1.addUndirectedEdge(x, y);
            }

            g1 = GraphUtils.replaceNodes(g1, dataSet.getVariables());

            G0 = g1;
        } else {
            IndependenceTest test = new IndTestScore(score, dataSet);
            System.out.println("FAS");

            FasStable fas = new FasStable(test);
            fas.setDepth(getDepth());
            fas.setVerbose(false);
            fas.setKnowledge(knowledge);
            G0 = fas.search();
        }

        SearchGraphUtils.pcOrientbk(knowledge, G0, G0.getNodes());

        System.out.println("Orientation");

        Graph graph = new EdgeListGraph(variables);

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                System.out.println("\nX = " + X + " Y = " + Y);

                // Centered
                final double[] x = colData[i];
                final double[] y = colData[j];

                boolean adj = isAdj(x, y);

                final boolean initialadj = initialGraph != null && initialGraph.isAdjacentTo(X, Y);

//                double sd1 = sd(residuals(x, new double[][]{y}));
//                double sd2 = sd(residuals(y, new double[][]{x}));

                if (initialadj || (isUseFasAdjacencies() && G0.isAdjacentTo(X, Y))
                        || (isUseSkewAdjacencies() && adj)) {//(Math.abs(c1 - c2)) > d3 + getSkewEdgeThreshold())) {
                    if (edgeForbiddenByKnowledge(X, Y)) {
                        // Don't add an edge.
                    } else if (knowledgeOrients(X, Y)) {
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        graph.addDirectedEdge(Y, X);
                    }
//                    else if (bidirected(x, y, G0, X, Y)) {
//                        Edge edge1 = Edges.directedEdge(X, Y);
//                        Edge edge2 = Edges.directedEdge(Y, X);
//                        graph.addEdge(edge1);
//                        graph.addEdge(edge2);
//                    }
                    else {
                        final double lrxy = leftright(x, y);
                        final double lryx = leftright(y, x);

                        if (lrxy> 0) {
                            graph.addDirectedEdge(X, Y);
                        } else if (lryx > 0) {
                            graph.addDirectedEdge(Y, X);
                        }

//                        if (lrxy < 0 && lryx < 0) {
////                            graph.addDirectedEdge(X, Y);
////                            graph.addDirectedEdge(Y, X);
//                            graph.addBidirectedEdge(X, Y);
//                        } else if (lrxy > 0 && lryx > 0) {
////                            graph.addDirectedEdge(X, Y);
////                            graph.addDirectedEdge(Y, X);
//                            graph.addBidirectedEdge(X, Y);
//                        } else if (lrxy > 0) {
//                            graph.addDirectedEdge(X, Y);
//                        } else if (lryx > 0) {
//                            graph.addDirectedEdge(Y, X);
//                        }
                    }
                }
            }
        }

//        removeExtraEdges(graph);
//        orientTwoCycles(variables, colData, G0, graph);
//        removeExtr    aEdges(graph);

        System.out.println();
        System.out.println("Done");

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

    private boolean isAdj(double[] x, double[] y) {
        double c1 = StatUtils.cov(x, y, x, 0, +1)[0];
        double c2 = StatUtils.cov(x, y, y, 0, +1)[0];

        double d1 = (covariance(x, y) / variance(x)) * (StatUtils.cov(x, x, x, 0, +1)[0]
                - StatUtils.cov(x, x, y, 0, +1)[0]);
        double d2 = (covariance(x, y) / variance(y)) * (StatUtils.cov(y, y, y, 0, +1)[0]
                - StatUtils.cov(y, y, x, 0, +1)[0]);

        System.out.println("d1 = " + d1 + " d2 = " + d2);

        double d3 = Math.max(abs(d1), abs(d2));

        return abs(c1 - c2) > d3 + getSkewEdgeThreshold();
    }

    private void orientTwoCycles(List<Node> variables, double[][] colData, Graph g0, Graph graph) {
        for (Edge edge : graph.getEdges()) {
            Node X = edge.getNode1();
            Node Y = edge.getNode2();
            int i = variables.indexOf(X);
            int j = variables.indexOf(Y);

            if (edgeForbiddenByKnowledge(X, Y)) continue;

            final double[] x = colData[i];
            final double[] y = colData[j];

            if (graph.isAdjacentTo(X, Y)
                    && !knowledgeOrients(X, Y) && !knowledgeOrients(Y, X)
                    && bidirected(x, y, g0, X, Y)) {
                graph.removeEdges(X, Y);
                graph.addDirectedEdge(X, Y);
                graph.addDirectedEdge(Y, X);
            }
        }
    }

    private void removeExtraEdges(Graph graph) {
        final int depth2 = this.depth == -1 ? 1000 : this.depth;

        boolean existsAnother = true;

        for (int d = 1; d <= depth2; d++) {
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
                    if (graph.getEdges(tail, head).size() == 2) continue;

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
//        parents1.addAll(graph.getParents(tail));
        parents1.remove(head);
//        parents1.remove(tail);
        Set<Node> c = new HashSet<>();

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

            // indirect path
            if (existsPath(graph, tail, p1, head) && existsPath(graph, p1, head, tail)) {
                c.add(p1);
            }

            // cycle back to tail
//            if (existsPath(graph, head, p1, tail) && existsPath(graph, p1, tail, head)) {
//                c.add(p1);
//            }

            // head
            if (existsPath(graph, p1, tail, head) && existsPath(graph, p1, head, tail)) {
                c.add(p1);
            }
        }

        return new ArrayList<>(c);
    }

    private boolean skewAdjacent(Node X, Node Y, List<Node> Z) {
        boolean b1 = false, b2 = false;

        try {
            FaskD.E hx = new FaskD.E(X, Y, Z, X).invoke();

            double[] rxx = hx.getRx();
            double[] ryx = hx.getRy();

            FaskD.E hy = new FaskD.E(X, Y, Z, Y).invoke();

            double[] rxy = hx.getRx();
            double[] ryy = hx.getRy();

            double[] x = data[variables.indexOf(X)];
            double[] y = data[variables.indexOf(Y)];

            double c1 = StatUtils.cov(rxx, ryx, rxx, Double.NEGATIVE_INFINITY, +1)[0];
            double c2 = StatUtils.cov(rxy, ryy, rxy, Double.NEGATIVE_INFINITY, +1)[0];

            double d1 = (covariance(rxx, rxx) / variance(rxx)) * (StatUtils.cov(rxx, rxx, rxx, Double.NEGATIVE_INFINITY, +1)[0]
                    - StatUtils.cov(rxx, rxx, rxx, Double.NEGATIVE_INFINITY, +1)[0]);
            double d2 = (covariance(rxy, ryy) / variance(ryy)) * (StatUtils.cov(ryy, ryy, ryy, Double.NEGATIVE_INFINITY, +1)[0]
                    - StatUtils.cov(ryy, ryy, rxy, Double.NEGATIVE_INFINITY, +1)[0]);

            System.out.println("c1 = " + c1 + " c2 = " + c2 + " d1 = " + d1 + " d2 = " + d2);

            double d3 = Math.max(abs(d1), abs(d2));

            b1 = abs(c1 - c2) > d3 - getSkewEdgeThreshold();
        } catch (Exception e) {
            e.printStackTrace();
        }

//        try {
//            Fask.E hx = new Fask.E(X, Y, Z, Y).invoke();
//            Fask.E hy = new Fask.E(X, Y, Z, X).invoke();
//
//            double[] x = hx.getRx();
//            double[] y = hy.getRy();
//
//            double c1 = StatUtils.cov(x, y, x, 0, +1)[0];
//            double c2 = StatUtils.cov(x, y, y, 0, +1)[0];
//
//            double d1 = (covariance(x, y) / variance(x)) * (StatUtils.cov(x, x, x, 0, +1)[0]
//                    - StatUtils.cov(x, x, y, 0, +1)[0]);
//            double d2 = (covariance(x, y) / variance(y)) * (StatUtils.cov(y, y, y, 0, +1)[0]
//                    - StatUtils.cov(y, y, x, 0, +1)[0]);
//
//            System.out.println("d1 = " + d1 + " d2 = " + d2);
//
//            double d3 = Math.max(abs(d1), abs(d2));
//
//            b2 = abs(c1 - c2) > d3 + getSkewEdgeThreshold();
//
//
//        } catch (Exception e) {
//            //
//        }

        return b1;//|| b2;

    }

    private class E {
        private Node x;
        private Node y;
        private List<Node> z;
        private Node condition;
        private List<Integer> rows;
        private double[] rx;
        private double[] ry;

        E(Node X, Node Y, List<Node> Z, Node condition) {
            x = X;
            y = Y;
            z = Z;

            if (z.contains(x) || z.contains(y)) {
                throw new IllegalArgumentException("Z should not contain X or Y.");
            }

            this.condition = condition;
        }

        public List<Integer> getRows() {
            return rows;
        }

        public double[] getRx() {
            return rx;
        }
        public double[] getRy() {
            return ry;
        }

        public FaskD.E invoke() {

            if (condition != null) {
                final double[] _w = data[variables.indexOf(condition)];
                rows = StatUtils.getRows(_w, 0, +1);
            } else {
                rows = new ArrayList<>();
                for (int i = 0; i < dataSet.getNumRows(); i++) rows.add(i);
            }


            int[] _rows = new int[rows.size()];
            for (int i = 0; i < rows.size(); i++) _rows[i] = rows.get(i);

            regressionDataset.setRows(_rows);

            double[] rx = regressionDataset.regress(x, z).getResiduals().toArray();
            double[] ry = regressionDataset.regress(y, z).getResiduals().toArray();

            this.rx = rx;
            this.ry = ry;

            return this;
        }

    }

    private boolean existsPath(Graph graph, Node X, Node Y, Node Z) {
        final List<List<Node>> paths = directedPathsFromTo(graph, X, Y);
        boolean found = false;

        for (List<Node> path : paths) {
            if (path.contains(Z)) continue;
            if (path.isEmpty()) continue;
            found = true;
            break;

        }

        return found && !knowledge.isForbidden(X.getName(), Y.getName());
    }

    private static List<List<Node>> directedPathsFromTo(Graph graph, Node node1, Node node2) {
        List<List<Node>> paths = new LinkedList<>();
        directedPathsFromToVisit(graph, node1, node2, new LinkedList<>(), paths, 3);
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


    private boolean bidirected(double[] x, double[] y, Graph G0, Node X, Node Y) {

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

//    private boolean leftright(double[] x, double[] y) {
//        double left = cu(x, y, x) / (sqrt(cu(x, x, x) * cu(y, y, x)));
//        double right = cu(x, y, y) / (sqrt(cu(x, x, y) * cu(y, y, y)));
//        double lr = left - right;
//
//        double r = StatUtils.correlation(x, y);
//        double sx = StatUtils.skewness(x);
//        double sy = StatUtils.skewness(y);
//
//        r *= signum(sx) * signum(sy);2
//        lr *= signum(r);
//        if (r < getDelta()) lr *= -1;
//
//        return lr > 0;
//    }

    // If X->Y then this number should be > 0.
    private double leftright(double[] x, double[] y) {
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);

        double[] ry = residuals(y, new double[][]{x});

        double a = covariance(x, y) / variance(x);
        double lr = E(a, x, ry, y, -1) - E(a, x, ry, y, +1);

        System.out.println("corr(x, ry) = " + correlation(x, ry));
        System.out.println("  corr(x^2, ry) = " + correlation(square(x), ry));
        System.out.println("  corr(x, ry^2) = " + correlation(x, square(ry)));

        return lr;
    }

    private double[] square(double[] x) {
        double[] square = Arrays.copyOf(x, x.length);
        for (int i = 0; i < square.length; i++) square[i] = square[i] * square[i];
        return square;
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

    private static double E(double a, double[] x, double[] ry, double[] y, double dir) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            final double _x = x[k];
            final double _ry = ry[k];
            final double _y = abs(a) * _x + _ry;

            if (_x * dir < 0 && _y * dir > 0) {
                exy += _x * _ry;
                n++;
            }
        }

        return exy / n;
    }

    private static double cu(double[] x, double[] y, double[] condition) {
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

    private static double cu(double[] x, double[] y, double[] condition, double dir) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (condition[k] * dir > 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
    }

    private static double cu2(double[] x, double[] _x, double[] _y) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (_x[k] > 0 && _y[k] > 0) {
                exy += x[k] * x[k];
                n++;
            }
        }

        return exy / n;
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
     * @return Returns the penalty discount used for the adjacency search. The default is 1,
     * though a higher value is recommended, say, 2, 3, or 4.
     */
    public double getPenaltyDiscount() {
        return penaltyDiscount;
    }

    /**
     * @param penaltyDiscount Sets the penalty discount used for the adjacency search.
     *                        The default is 1, though a higher value is recommended, say,
     *                        2, 3, or 4.
     */
    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }

    /**
     * @param alpha Alpha for orienting 2-cycles. Needs to be on the low side usually. Default 1e-6.
     */
    public void setAlpha(double alpha) {
        this.alpha = alpha;
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

    //======================================== PRIVATE METHODS ====================================//

    private boolean knowledgeOrients(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) || knowledge.isRequired(left.getName(), right.getName());
    }

    private boolean edgeForbiddenByKnowledge(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) && knowledge.isForbidden(left.getName(), right.getName());
    }

    public Graph getInitialGraph() {
        return initialGraph;
    }

    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    public double getSkewEdgeThreshold() {
        return skewEdgeThreshold;
    }

    public void setSkewEdgeThreshold(double skewEdgeThreshold) {
        this.skewEdgeThreshold = skewEdgeThreshold;
    }

    public boolean isUseFasAdjacencies() {
        return useFasAdjacencies;
    }

    public void setUseFasAdjacencies(boolean useFasAdjacencies) {
        this.useFasAdjacencies = useFasAdjacencies;
    }

    public boolean isUseSkewAdjacencies() {
        return useSkewAdjacencies;
    }

    public void setUseSkewAdjacencies(boolean useSkewAdjacencies) {
        this.useSkewAdjacencies = useSkewAdjacencies;
    }
//
//    public double getDelta() {
//        return delta;
//    }
//
//    public void setDelta(double delta) {
//        this.delta = delta;
//    }
}





