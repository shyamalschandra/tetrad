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

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradLogger;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static edu.cmu.tetrad.util.StatUtils.correlation;
import static edu.cmu.tetrad.util.StatUtils.mean;
import static edu.cmu.tetrad.util.StatUtils.variance;
import static java.lang.Math.*;

/**
 * Implements the PC ("Peter/Clark") algorithm, as specified in Chapter 6 of Spirtes, Glymour, and Scheines, "Causation,
 * Prediction, and Search," 2nd edition, with a modified rule set in step D due to Chris Meek. For the modified rule
 * set, see Chris Meek (1995), "Causal inference and causal explanation with background knowledge."
 *
 * @author Joseph Ramsey.
 */
public class SkewSearch implements GraphSearch {

    private final DataSet dataSet;
    private final double[][] data;
    private final List<Node> nodes;

    /**
     * Forbidden and required edges for the search.
     */
    private IKnowledge knowledge = new Knowledge2();

    /**
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    private int depth = 1000;

    /**
     * True if cycles are to be aggressively prevented. May be expensive for large graphs (but also useful for large
     * graphs).
     */
    private boolean aggressivelyPreventCycles = false;

    private boolean verbose = false;

    private boolean fdr = false;
    private double delta = -0.01;
    private RegressionDataset regressionDataset;
    private double twoCycleAlpha = 0.05;
    private double twoCycleAlphaCutoff;

    //=============================CONSTRUCTORS==========================//

    /**
     * Constructs a new PC search using the given independence test as oracle.
     */
    public SkewSearch(DataSet dataSet) {
        dataSet = DataUtils.standardizeData(dataSet);
        this.dataSet = dataSet;
        this.data = dataSet.getDoubleData().transpose().toArray();
        this.nodes = dataSet.getVariables();
        regressionDataset = new RegressionDataset(dataSet);
    }

    private double[] correctSkewnesses(double[] data) {
        double skewness = StatUtils.skewness(data);
        double[] data2 = new double[data.length];
        for (int i = 0; i < data.length; i++) data2[i] = data[i] * Math.signum(skewness);
        return data2;
    }

    //==============================PUBLIC METHODS========================//

    /**
     * @return true iff edges will not be added if they would create cycles.
     */
    public boolean isAggressivelyPreventCycles() {
        return this.aggressivelyPreventCycles;
    }

    /**
     * @param aggressivelyPreventCycles Set to true just in case edges will not be addeds if they would create cycles.
     */
    public void setAggressivelyPreventCycles(boolean aggressivelyPreventCycles) {
        this.aggressivelyPreventCycles = aggressivelyPreventCycles;
    }

    /**
     * @return the knowledge specification used in the search. Non-null.
     */
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    /**
     * Sets the knowledge specification to be used in the search. May not be null.
     */
    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) {
            throw new NullPointerException();
        }

        this.knowledge = knowledge;
    }

    /**
     * @return the current depth of search--that is, the maximum number of conditioning nodes for any conditional
     * independence checked.
     */
    public int getDepth() {
        return depth;
    }

    /**
     * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
     * checked.
     *
     * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
     *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
     *              machines.
     */
    public void setDepth(int depth) {
        if (depth < -1) {
            throw new IllegalArgumentException("Depth must be -1 or >= 0: " + depth);
        }

        this.depth = depth;
    }

    /**
     * Runs PC starting with a complete graph over all nodes of the given conditional independence test, using the given
     * independence test and knowledge and returns the resultant graph. The returned graph will be a pattern if the
     * independence information is consistent with the hypothesis that there are no latent common causes. It may,
     * however, contain cycles or bidirected edges if this assumption is not born out, either due to the actual presence
     * of latent common causes, or due to statistical errors in conditional independence judgments.
     */
    @Override
    public Graph search() {
        setCutoffs();

        Graph graph = new EdgeListGraph(nodes);
        List<Node> empty = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                if (i == j) continue;

                Node X = nodes.get(i);
                Node Y = nodes.get(j);

                if (!edgeForbiddenByKnowledge(X, Y)) {
                    if (knowledgeOrients(X, Y)) {
                        graph.addDirectedEdge(X, Y);
                    } else if (knowledgeOrients(Y, X)) {
                        graph.addDirectedEdge(Y, X);
                    } else {
                        if (existsEdge(X, Y, empty)) {
                            final double[] x = data[nodes.indexOf(X)];
                            final double[] y = data[nodes.indexOf(Y)];

                            if (!leftright(y, x)) {
                                graph.addDirectedEdge(X, Y);
                            } else if (!leftright(x, y)) {
                                graph.addDirectedEdge(Y, X);
                            }
                        }
                    }
                }
            }
        }

        for (Edge edge : new ArrayList<>(graph.getEdges())) {
            final Node X = Edges.getDirectedEdgeTail(edge);
            final Node Y = Edges.getDirectedEdgeHead(edge);

            final double[] x = data[nodes.indexOf(X)];
            final double[] y = data[nodes.indexOf(Y)];


            if (!knowledgeOrients(X, Y) && !knowledgeOrients(Y, X)
                    && graph.getEdges(edge.getNode1(), edge.getNode2()).size() != 2
                    && bidirected(x, y, graph, X, Y)) {
                graph.removeEdges(X, Y);

                Edge edge1 = Edges.directedEdge(X, Y);
                Edge edge2 = Edges.directedEdge(Y, X);
                edge1.setLineColor(Color.GREEN);
                edge2.setLineColor(Color.GREEN);
                edge1.setBold(true);
                edge2.setBold(true);
                graph.addEdge(edge1);
                graph.addEdge(edge2);
            }
        }


        EDGES:
        for (Edge edge : new ArrayList<>(graph.getEdges())) {
            if (edge.isDirected()) {
                Node tail = Edges.getDirectedEdgeTail(edge);
                Node head = Edges.getDirectedEdgeHead(edge);

                List<Node> parents = graph.getParents(head);
                parents.remove(tail);

                final int depth = this.depth == -1 ? 1000 : this.depth;
                final int min = min(depth, parents.size());

                DepthChoiceGenerator gen = new DepthChoiceGenerator(parents.size(), min);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    if (choice.length == 0) continue;

                    List<Node> Z = GraphUtils.asList(choice, parents);

                    if (!existsEdge(tail, head, Z)) {
                        graph.removeEdge(edge);
                        continue EDGES;
                    }
                }
            }
        }

        return graph;
    }

    @Override
    public long getElapsedTime() {
        return 0;
    }

    double alpha = 0.01;

    private boolean existsEdge(Node x, Node y, List<Node> Z) {
        final double[] _x = data[nodes.indexOf(x)];
        final double[] _y = data[nodes.indexOf(y)];

        List<Integer> rowsx = StatUtils.getRows(_x, 0, +1);
        int[] _rowsx = new int[rowsx.size()];
        for (int i = 0; i < rowsx.size(); i++) _rowsx[i] = rowsx.get(i);

        List<Integer> rowsy = StatUtils.getRows(_y, 0, +1);
        int[] _rowsy = new int[rowsy.size()];
        for (int i = 0; i < rowsy.size(); i++) _rowsy[i] = rowsy.get(i);

        regressionDataset.setRows(_rowsx);
        double[] rxzx = regressionDataset.regress(x, Z).getResiduals().toArray();
        double[] ryzx = regressionDataset.regress(y, Z).getResiduals().toArray();

        regressionDataset.setRows(_rowsy);
        double[] rxzy = regressionDataset.regress(x, Z).getResiduals().toArray();
        double[] ryzy = regressionDataset.regress(y, Z).getResiduals().toArray();

        double[] sxyx = new double[rxzx.length];

        for (int i = 0; i < rxzx.length; i++) {
            sxyx[i] = rxzx[i] * ryzx[i];
        }

        double[] sxyy = new double[rxzy.length];

        for (int i = 0; i < rxzy.length; i++) {
            sxyy[i] = rxzy[i] * ryzy[i];
        }

        double zv3 = (mean(sxyx)) / sqrt(variance(sxyx) / sxyx.length);
        double zv4 = (mean(sxyy)) / sqrt(variance(sxyy) / sxyy.length);

        double c3 = new TDistribution(sxyx.length).cumulativeProbability(abs(zv3));

        double c4 = new TDistribution(sxyy.length).cumulativeProbability(abs(zv4));

        boolean b3 = 2 * (1.0 - c3) < alpha;
        boolean b4 = 2 * (1.0 - c4) < alpha;

        return b3 && b4;
    }


    // If x->y, returns true
    // If x->y, returns true
    private boolean leftright(double[] x, double[] y) {
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

        if (StatUtils.correlation(x, y) < 0) lr += getDelta();
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

    private void setCutoffs() {
        if (twoCycleAlpha < 0.0 || twoCycleAlpha > 1.0) {
            throw new IllegalArgumentException("Significance out of range for two-cycles: " + twoCycleAlpha);
        }

        this.twoCycleAlphaCutoff = StatUtils.getZForAlpha(twoCycleAlpha);
    }

    private double partialCorrelation(double[] x, double[] y, double[][] z, double[] condition, double threshold,
                                      double direction) throws SingularMatrixException {
        double[][] cv = covMatrix(x, y, z, condition, threshold, direction);
        TetradMatrix m = new TetradMatrix(cv).transpose();
        return StatUtils.partialCorrelation(m);
    }

    public static double[][] covMatrix(double[] x, double[] y, double[][] z, double[] condition, double threshold,
                                       double direction) {
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

    public static List<Integer> getRows(double[] x, double[] condition, double threshold, double direction) {
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

    public double getDelta() {
        return delta;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * True iff the algorithm should be run with False Discovery Rate tests.
     */
    public boolean isFdr() {
        return fdr;
    }

    public void setFdr(boolean fdr) {
        this.fdr = fdr;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setDelta(double delta) {
        this.delta = delta;
    }

    private boolean edgeForbiddenByKnowledge(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) && knowledge.isForbidden(left.getName(), right.getName());
    }

    private boolean knowledgeOrients(Node left, Node right) {
        return knowledge.isForbidden(right.getName(), left.getName()) || knowledge.isRequired(left.getName(), right.getName());
    }

    private boolean bidirected(double[] x, double[] y, Graph G0, Node X, Node Y) {
        Set<Node> adjSet = new HashSet<>(G0.getAdjacentNodes(X));
        adjSet.addAll(G0.getAdjacentNodes(Y));
        List<Node> adj = new ArrayList<>(adjSet);
        adj.remove(X);
        adj.remove(Y);

        DepthChoiceGenerator gen = new DepthChoiceGenerator(adj.size(), min(depth, adj.size()));
        int[] choice;

        while ((choice = gen.next()) != null) {
            List<Node> _adj = GraphUtils.asList(choice, adj);
            double[][] _Z = new double[_adj.size()][];

            for (int f = 0; f < _adj.size(); f++) {
                Node _z = _adj.get(f);
                int column = dataSet.getColumn(_z);
                _Z[f] = data[column];
            }

            double pc;
            double pc1;
            double pc2;

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

            boolean rejected1 = abs(zv1) > twoCycleAlphaCutoff;
            boolean rejected2 = abs(zv2) > twoCycleAlphaCutoff;

            boolean possibleCycle = false;

            if (zv1 < 0 && zv2 > 0 && rejected1) {
                possibleCycle = true;
            } else if (zv1 > 0 && zv2 < 0 && rejected2) {
                possibleCycle = true;
            } else if (rejected1 && rejected2) {
                possibleCycle = true;
            }

            if (!possibleCycle) {
                return false;
            }
        }

        return true;
    }

    public void setTwoCycleAlpha(double twoCycleAlpha) {
        this.twoCycleAlpha = twoCycleAlpha;
    }
}




