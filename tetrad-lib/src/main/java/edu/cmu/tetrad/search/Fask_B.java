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

import java.util.*;

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
    private final IndependenceTest test;

    // An initial graph to orient, skipping the adjacency step.
    private Graph initialGraph = null;

    // Elapsed time of the search, in milliseconds.
    private long elapsed = 0;

    // The colData sets being analyzed. They must all have the same variables and the same
    // number of records.
    private DataSet dataSet;

    // For the Fast Adjacency Search.
    private int depth = -1;

    // For the SEM BIC score, for the Fast Adjacency Search.
    private double penaltyDiscount = 1;

    // Alpha for orienting 2-cycles. Usually needs to be low.
    private double alpha = 1e-6;

    // Knowledge the the search will obey, of forbidden and required edges.
    private IKnowledge knowledge = new Knowledge2();

    // Data as a double[][].
    private final double[][] colData;

    // A threshold for including extra adjacencies due to skewness.
    private double skewEdgeAlpha = 0.3;

    // True if FAS adjacencies should be included in the output.
    private boolean useFasAdjacencies = true;

    // True if skew adjacencies should be included in the output.
    private boolean useSkewAdjacencies = true;

    // Threshold for reversing casual judgments for negative coefficients.
    private double delta = -0.2;

    private double twoCycleAlpha = 1e-6;
    private double twoCycleCutoff;

    private RegressionDataset regressionDataset;
    private final List<Node> variables;

    /**
     * @param dataSet These datasets must all have the same variables, in the same order.
     */
    public Fask_B(DataSet dataSet, IndependenceTest test) {
        this.dataSet = dataSet;
        this.test = test;

        dataSet = DataUtils.standardizeData(dataSet);
        colData = dataSet.getDoubleData().transpose().toArray();
        this.variables = dataSet.getVariables();
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

    public double getDelta() {
        return delta;
    }

    public void setDelta(double delta) {
        this.delta = delta;
    }

    public void setUseFasAdjacencies(boolean useFasAdjacencies) {
        this.useFasAdjacencies = useFasAdjacencies;
    }

    public void setTwoCycleAlpha(double twoCycleAlpha) {
        this.twoCycleAlpha = twoCycleAlpha;
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

                if ((isUseFasAdjacencies() && G0.isAdjacentTo(X, Y))
                        || (isUseSkewAdjacencies() && (skewAdjacent(X, Y, Collections.emptyList())))) {
                    if (!edgeForbiddenByKnowledge(X, Y)) {
                        if (knowledgeOrients(X, Y)) {
                            graph.addDirectedEdge(X, Y);
                        } else if (knowledgeOrients(Y, X)) {
                            graph.addDirectedEdge(Y, X);
                        } else {
                            final boolean leftrightxy = leftright(X, Y);
                            final boolean leftrightyx = leftright(Y, X);

                            if (leftrightxy && leftrightyx) {
                                graph.addDirectedEdge(X, Y);
                                graph.addDirectedEdge(Y, X);
                            } else if (!leftrightxy && !leftrightyx) {
                                graph.addDirectedEdge(X, Y);
                                graph.addDirectedEdge(Y, X);
                            } else if (leftrightxy) {
                                graph.addDirectedEdge(X, Y);
                            } else {
                                graph.addDirectedEdge(Y, X);
                            }
                        }
                    }
                }
            }
        }

        Set<Edge> edges = graph.getEdges();

        EDGES:
        for (Edge edge : new ArrayList<>(edges)) {
            if (edge.isDirected()) {
                Node tail = Edges.getDirectedEdgeTail(edge);
                Node head = Edges.getDirectedEdgeHead(edge);

                List<Node> parents = graph.getParents(head);
                parents.remove(tail);

                for (Node p : new ArrayList<>(parents)) {
                    if (GraphUtils.directedPathsFromTo(graph, tail, p, 3).isEmpty()) {
                        parents.remove(p);
                    }
                }

                final int depth = this.depth == -1 ? 1000 : this.depth;
                final int min = min(depth, parents.size());

                ChoiceGenerator gen = new ChoiceGenerator(parents.size(), min);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    if (choice.length == 0) continue;
                    List<Node> Z = GraphUtils.asList(choice, parents);

                    if (!skewAdjacent(tail, head, Z)) {
                        graph.removeEdge(edge);
                        continue EDGES;
                    }
                }
            }
        }

        for (Edge edge : edges) {
            Node X = edge.getNode1();
            Node Y = edge.getNode2();

            if (!(knowledgeOrients(X, Y) || knowledgeOrients(Y, X))) {
                if (bidirected(X, Y, G0)) {
                    Edge edge1 = Edges.directedEdge(X, Y);
                    Edge edge2 = Edges.directedEdge(Y, X);
                    graph.addEdge(edge1);
                    graph.addEdge(edge2);
                }
            }
        }

        System.out.println();
        System.out.println("Done");

        long stop = System.currentTimeMillis();
        this.elapsed = stop - start;

        return graph;
    }

    //======================================== PRIVATE METHODS ====================================//

    private boolean skewAdjacent(Node X, Node Y, List<Node> Z) {
        final double[] _x = colData[variables.indexOf(X)];
        final double[] _y = colData[variables.indexOf(Y)];

        List<Integer> rowsx = StatUtils.getRows(_x, 0, +1);
        List<Integer> rowsy = StatUtils.getRows(_y, 0, +1);

        int n1 = rowsx.size();
        int n2 = rowsy.size();

        int[] _rowsx = new int[rowsx.size()];
        for (int i = 0; i < rowsx.size(); i++) _rowsx[i] = rowsx.get(i);

        int[] _rowsy = new int[rowsy.size()];
        for (int i = 0; i < rowsy.size(); i++) _rowsy[i] = rowsy.get(i);

        regressionDataset.setRows(_rowsx);
        double[] rxzx = regressionDataset.regress(X, Z).getResiduals().toArray();
        double[] ryzx = regressionDataset.regress(Y, Z).getResiduals().toArray();

        regressionDataset.setRows(_rowsy);
        double[] rxzy = regressionDataset.regress(X, Z).getResiduals().toArray();
        double[] ryzy = regressionDataset.regress(Y, Z).getResiduals().toArray();

        double[] d1 = new double[rowsx.size()];

        for (int i = 0; i < rowsx.size(); i++) {
            d1[i] = rxzx[i] * ryzx[i];
        }

        double[] d2 = new double[rowsy.size()];

        for (int i = 0; i < rowsy.size(); i++) {
            d2[i] = rxzy[i] * ryzy[i];
        }

        // Unequal variances, uequal sample sizes, T test, 2-sided
        double e2 = variance(d2) / ((double) n2);
        double e1 = variance(d1) / ((double) n1);
        double t = abs(mean(d2) - mean(d1)) / sqrt(e2 + e1);
        double df = ((e2 + e1) * (e2 + e2)) / ((e2 * e2) / (n2 - 1)) + ((e1 * e1) / (n1 - 1));
        double cutoff = new TDistribution(df).inverseCumulativeProbability(1.0 - alpha / 2.0);
        return t > cutoff;
    }


    // If x->y, returns true
    private boolean leftright1(Node X, Node Y) {
        final double[] x = colData[variables.indexOf(X)];
        final double[] y = colData[variables.indexOf(Y)];

        final double cxyx = e(x, y, x);
        final double cxyy = e(x, y, y);
        final double cxxx = e(x, x, x);
        final double cyyx = e(y, y, x);
        final double cxxy = e(x, x, y);
        final double cyyy = e(y, y, y);

        double a1 = cxyx / cxxx;
        double a2 = cxyy / cxxy;
        double b1 = cxyy / cyyy;
        double b2 = cxyx / cyyx;

        double Q = (a2 > 0) ? a1 / a2 : a2 / a1;
        double R = (b2 > 0) ? b1 / b2 : b2 / b1;

        double lr = Q - R;

        if (correlation(x, y) < 0) lr += getDelta();
        return lr > 0;
    }

    // If x->y, returns true
    private boolean leftright(Node X, Node Y) {
        final double[] x = colData[variables.indexOf(X)];
        final double[] y = colData[variables.indexOf(Y)];

        final double cxyx = e(x, y, x);
        final double cxyy = e(x, y, y);
        final double cxxx = e(x, x, x);
        final double cxxy = e(x, x, y);

        double exeyy = (cxyy / cxxy - cxyx / cxxx) * cxxy;

        if (correlation(x, y) < 0 && correlation(x, y) > getDelta()) exeyy *= -1;

        return exeyy < 0;
    }

    private boolean bidirected(Node X, Node Y, Graph G0) {
        final double[] x = colData[variables.indexOf(X)];
        final double[] y = colData[variables.indexOf(Y)];

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

    private static double e(double[] x, double[] y, double[] condition) {
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

    private double[] correctSkewnesses(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
        return data2;
    }
}






