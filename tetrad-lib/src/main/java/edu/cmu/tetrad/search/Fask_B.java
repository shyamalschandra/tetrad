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

import java.util.*;
import java.util.List;

import static edu.cmu.tetrad.util.StatUtils.*;
import static java.lang.Math.*;
import static java.lang.Math.min;

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

    // Data with skew corrected.
    private final double[][] skewCorrected;

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

    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public Fask_B(DataSet dataSet, IndependenceTest test) {
        this.dataSet = dataSet;
        this.test = test;

        dataSet = DataUtils.standardizeData(dataSet);
        colData = dataSet.getDoubleData().transpose().toArray();
        this.variables = dataSet.getVariables();

        for (Node variable : variables) {
            System.out.println(variable + " skewness = " + skewness(variable));
        }

        skewCorrected = new double[colData.length][];

        for (int i = 0; i < skewCorrected.length; i++) {
            skewCorrected[i] = correctSkewness(colData[i]);
        }

        regressionDataset = new RegressionDataset(dataSet);
    }

    //======================================== PUBLIC METHODS ====================================//

    public Graph getInitialGraph() {
        return initialGraph;
    }

    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    public double getSkewEdgeAlpha() {
        return skewEdgeAlpha;
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
     * Sets the significance level at which independence judgments should be made.  Affects the cutoff for partial
     * correlations to be considered statistically equal to zero.
     */
    private void setCutoff() {
        this.twoCycleCutoff = StatUtils.getZForAlpha(twoCycleAlpha);
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

        setCutoff();

        Graph fasGraph;

        if (getInitialGraph() != null) {
            Graph g1 = new EdgeListGraph(getInitialGraph().getNodes());

            for (Edge edge : getInitialGraph().getEdges()) {
                Node X = edge.getNode1();
                Node Y = edge.getNode2();

                if (!g1.isAdjacentTo(X, Y)) g1.addUndirectedEdge(X, Y);
            }

            g1 = GraphUtils.replaceNodes(g1, dataSet.getVariables());

            fasGraph = g1;
        } else {
            System.out.println("FAS");

            FasStable fas = new FasStable(test);
            fas.setDepth(getDepth());
            fas.setVerbose(false);
            fas.setKnowledge(knowledge);
            fasGraph = fas.search();

//            SemBicScore score = new SemBicScore(new CovarianceMatrixOnTheFly((DataSet) test.getData()));
//            score.setPenaltyDiscount(2);
//            Fges fges = new Fges(score);
//            fges.setMaxDegree(getDepth());
//            fges.setKnowledge(knowledge);
//            fasGraph = fges.search();
        }

        SearchGraphUtils.pcOrientbk(knowledge, fasGraph, fasGraph.getNodes());

        System.out.println("Orientation");

        Graph graph = new EdgeListGraph(variables);

        for (int i = 0; i < variables.size(); i++) {
            for (int j = i + 1; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                final double[] x = colData[i];
                final double[] y = colData[j];

                double c1 = StatUtils.cov(x, y, x, 0, +1)[1];
                double c2 = StatUtils.cov(x, y, y, 0, +1)[1];

                if ((isUseFasAdjacencies() && fasGraph.isAdjacentTo(X, Y))
                        || (isUseSkewAdjacencies() && getMaskThreshold() != 0 ?
                        Math.abs(c1 - c2) > getMaskThreshold() : ((skewAdjacent(X, Y, Collections.emptyList()))))) {
                    orientEdge(getDelta(), graph, X, Y, false);
                }
            }
        }

        removeExtraEdges(graph);

        for (int i = 0; i < variables.size(); i++) {
            for (int j = 0; j < variables.size(); j++) {
                Node X = variables.get(i);
                Node Y = variables.get(j);

                if (graph.isAdjacentTo(X, Y)) {
                    graph.removeEdges(X, Y);
                    orientEdge(getDelta(), graph, X, Y, true);
                }
            }
        }

        System.out.println();
        System.out.println("Done");

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

    private void orientEdge(double delta, Graph graph, Node X, Node Y, boolean useTwoCycleRule) {
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
                } else if (useTwoCycleRule && twocycle(X, Y, graph)) {
                    graph.addDirectedEdge(X, Y);
                    graph.addDirectedEdge(Y, X);
                } else {
                    final boolean lrxy = leftRightMinnesota(X, Y, delta, false);
                    final boolean lryx = leftRightMinnesota(Y, X, delta, false);

                    if (!lrxy && !lryx) {
                        graph.addBidirectedEdge(X, Y);
                    } else if (lrxy & lryx) {
                        if (useTwoCycleRule) {
                            graph.addNondirectedEdge(X, Y);
                        } else {
                            graph.addDirectedEdge(X, Y);
                            graph.addDirectedEdge(Y, X);
                        }
                    } else if (lrxy) {
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
                            System.out.println("Removing " + Edges.directedEdge(tail, head) + " | " + Z);
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

            if (!directedPathsFromTo(graph, tail, p1, 3).isEmpty()
                    || !directedPathsFromTo(graph, head, p1, 3).isEmpty()) {
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

    //======================================== PRIVATE METHODS ====================================//

    private boolean skewAdjacent(Node X, Node Y, List<Node> Z) {
        boolean b1, b2;

        {
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
        }

        {
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
        }

        return b1 || b2;
    }

    private boolean leftRight0(Node X, Node Y, double delta, boolean printStuff) {
        double[] x = skewCorrected[variables.indexOf(X)];
        double[] y = skewCorrected[variables.indexOf(Y)];

        final double cxyx = cu(x, y, x);
        final double cxyy = cu(x, y, y);

        double left = cxyx / sqrt(cu(x, x, x) * cu(y, y, x));
        double right = cxyy / sqrt(cu(x, x, y) * cu(y, y, y));

        double lr = left - right;

        double r = StatUtils.correlation(x, y);
        double sx = StatUtils.skewness(x);
        double sy = StatUtils.skewness(y);

        r *= signum(sx) * signum(sy);
        lr *= signum(r);
        if (r < -.2) lr *= -1;

        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));

        if (sk_ey < 0) {
            lr *= -1;
        }

        return lr > 0;
    }

    // If x->y, returns true
    private boolean leftRightMinnesota(Node X, Node Y, double delta, boolean printStuff) {
        double[] x = skewCorrected[variables.indexOf(X)];
        double[] y = skewCorrected[variables.indexOf(Y)];

        final double cxyx = cov(x, y, x);
        final double cxyy = cov(x, y, y);
        final double cxxx = cov(x, x, x);
        final double cyyx = cov(y, y, x);
        final double cxxy = cov(x, x, y);
        final double cyyy = cov(y, y, y);

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;
        double b1 = cxyy / cyyy;
        double b2 = cxyx / cyyx;

        double Q = (a2 > 0) ? a1 / a2 : a2 / a1;
        double R = (b2 > 0) ? b1 / b2 : b2 / b1;

        double lr = Q - R;

        if (StatUtils.correlation(x, y) < 0) lr += delta;

        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));

        if (sk_ey < 0) {
            lr *= -1;
        }

        return lr > 0;
    }

    private static double cov(double[] x, double[] y, double[] condition) {
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

    // If x->y, returns true
    private boolean leftRight2(Node X, Node Y, double delta, boolean printStuff) {
        double[] x = skewCorrected[variables.indexOf(X)];
        double[] y = skewCorrected[variables.indexOf(Y)];

        final double cxyx = e(x, y, x, true);
        final double cxyy = e(x, y, y, true);
        final double cxxx = e(x, x, x, true);
        final double cxxy = e(x, x, y, true);

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;

        double lr = (a1 - a2);

        final double sk_ey = StatUtils.skewness(residuals(y, new double[][]{x}));

        final double a = correlation(x, y);

        if (a < 0 && sk_ey > delta) {// &&*/ abs(sk_ey) > abs(skewness(X))) {
            lr *= -1;
        }

        if (sk_ey < 0) {
            lr *= -1;
        }

        if (printStuff) {
            System.out.println(Edges.directedEdge(X, Y) + " sk(X) = "
                    + StatUtils.skewness(x) + " sk(Y) = " + StatUtils.skewness(y)
                    + " corr(X, Y) = " + correlation(x, y) + " sk_ey = " + sk_ey + " LR = " + lr);
        }

        return lr > 0;
    }

    private boolean twocycle(Node X, Node Y, Graph graph) {
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

    private static double e(double[] x, double[] y, double[] condition, boolean positive) {
        double exy = 0.0;

        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (positive && condition[k] > 0) {
                exy += x[k] * y[k];
                n++;
            } else if (!positive && condition[k] < 0) {
                exy += x[k] * y[k];
                n++;
            }
        }

        return exy / n;
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

    private double[] correctSkewness(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
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

    private double coef(double[] _x, double[][] _y) {
        TetradMatrix x = new TetradMatrix(new double[][]{_x}).transpose();
        TetradMatrix y = new TetradMatrix(_y).transpose();

        TetradMatrix b = x.transpose().times(x).inverse().times(x.transpose().times(y));

        return b.get(0, 0);
    }

    public void setMaskThreshold(double maskThreshold) {
        this.maskThreshold = maskThreshold;
    }

    private double getMaskThreshold() {
        return maskThreshold;
    }

    public double getDelta() {
        return delta;
    }

    public void setDelta(double delta) {
        this.delta = delta;
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







