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
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.PermutationGenerator;
import edu.cmu.tetrad.util.TetradMatrix;

import java.util.*;

import static java.lang.StrictMath.abs;

/**
 * Implements the LiNGAM algorithm in Shimizu, Hoyer, Hyvarinen, and Kerminen, A linear nongaussian acyclic model for
 * causal discovery, JMLR 7 (2006). Largely follows the Matlab code.
 * <p>
 * We use FGES with knowledge of causal order for the pruning step.
 *
 * @author Joseph Ramsey
 */
public class Lingam {
    private double penaltyDiscount = 2;

    //================================CONSTRUCTORS==========================//

    /**
     * Constructs a new LiNGAM algorithm with the given alpha level (used for pruning).
     */
    public Lingam() {
    }

    public Graph search(DataSet data) {
        data = DataUtils.center(data);

        DataSet dataSet = data;
//        dataSet = DataUtils.standardizeData(dataSet);
        TetradMatrix X = dataSet.getDoubleData();
        FastIca fastIca = new FastIca(X, 30);

//        TetradMatrix winit = new TetradMatrix(X.columns(), X.columns());
//        for (int i1 = 0; i1 < X.columns(); i1++) {
//            for (int j = 0; j < X.columns(); j++) {
//                winit.set(i1, j, );
//            }
//        }

//        fastIca.setWInit(winit);
        fastIca.setVerbose(false);
        FastIca.IcaResult result11 = fastIca.findComponents();
        TetradMatrix W = result11.getW().transpose();

        System.out.println("W = " + W);

        PermutationGenerator gen1 = new PermutationGenerator(W.rows());
        int[] perm1 = new int[0];
        double sum1 = Double.NEGATIVE_INFINITY;
        int[] choice1;

        while ((choice1 = gen1.next()) != null) {
            double sum = 0.0;

            for (int i1 = 0; i1 < W.rows(); i1++) {
                final double c = W.get(i1, choice1[i1]);
                sum += 1.0 / abs(c);
            }

            if (sum > sum1) {
                sum1 = sum;
                perm1 = Arrays.copyOf(choice1, choice1.length);
            }
        }

        TetradMatrix WTilde = W.getSelection(perm1, perm1);

        System.out.println("WTilde before normalization = " + WTilde);

        for (int j = 0; j < WTilde.columns(); j++) {
            for (int i = j; i > 0; i--) {
                WTilde.set(i, j, WTilde.get(i, j) * WTilde.get(j, j));
            }
        }

        System.out.println("WTilde after normalization = " + WTilde);

        final int m = dataSet.getNumColumns();
        TetradMatrix b = TetradMatrix.identity(m).minus(WTilde.transpose());
        System.out.println("b = " + b);

        final List<Node> variables = data.getVariables();

        if (true) {

            PermutationGenerator gen2 = new PermutationGenerator(b.rows());
            int[] perm2 = new int[0];
            double sum2 = Double.POSITIVE_INFINITY;
            int[] choice2;

            while ((choice2 = gen2.next()) != null) {
                double sum = 0.0;

                for (int i1 = 0; i1 < W.rows(); i1++) {
                    for (int j = i1; j < W.rows(); j++) {
                        final double c = b.get(choice2[i1], choice2[j]);
                        sum += c * c;
                    }
                }

                if (sum < sum2) {
                    sum2 = sum;
                    perm2 = Arrays.copyOf(choice2, choice2.length);
                }
            }

            TetradMatrix BTilde = b.getSelection(perm2, perm2);

            System.out.println("BTilde = " + BTilde);

            final SemBicScore score = new SemBicScore(new CovarianceMatrixOnTheFly(data));
            score.setPenaltyDiscount(penaltyDiscount);
            Fges fges = new Fges(score);

            IKnowledge knowledge = new Knowledge2();

            for (int i = 0; i < variables.size(); i++) {
                knowledge.addToTier(i + 1, variables.get(perm2[i]).getName());
            }

            fges.setKnowledge(knowledge);

            final Graph graph = fges.search();
            System.out.println("graph Returning this graph: " + graph);
            return graph;
        } else {
            List<Double> elements = new ArrayList<>();
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    elements.add(abs(b.get(i, j)));
                }
            }

            Collections.sort(elements);
            double cutoff = elements.get(m * (m + 1) / 2);

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    if (abs(b.get(i, j)) <= cutoff) {
                        b.set(i, j, 0);
                    }
                }
            }

            int[] perm3 = null;
            int extra = 0;
            int[] choice2;

            LOOP:
            while (true) {
                PermutationGenerator gen2 = new PermutationGenerator(b.rows());

                PERM:
                while (((choice2 = gen2.next()) != null)) {
                    for (int i = 0; i < m; i++) {
                        for (int j = i; j < m; j++) {
                            if (b.get(choice2[i], choice2[j]) != 0) {
                                System.out.println(b.getSelection(choice2, choice2));

                                System.out.println("Rejecting " + Arrays.toString(choice2));
                                continue PERM;
                            }
                        }
                    }

                    System.out.println("PP = " + b);

                    perm3 = Arrays.copyOf(choice2, choice2.length);
                    break LOOP;
                }

                extra++;

                cutoff = elements.get(m * (m + 1) / 2 + extra);

                for (int q = 0; q < m; q++) {
                    for (int r = 0; r < m; r++) {
                        if (abs(b.get(q, r)) <= cutoff) {
                            b.set(q, r, 0);

                            System.out.println("BB = " + b);
                        }
                    }
                }
            }


            System.out.println("B after zeros: " + b);

            Graph graph = new EdgeListGraph(variables);

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    if (b.get(j, i) != 0) {
                        graph.addDirectedEdge(variables.get(perm3[j]), variables.get(perm3[i]));
                    }
                }
            }

            return graph;
        }
    }

    //================================PUBLIC METHODS========================//

    private int[] algorithmC(TetradMatrix b) {
        int m = b.rows();

        List<Double> elements = new ArrayList<>();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                elements.add(abs(b.get(i, j)));
            }
        }

        Collections.sort(elements);
        double cutoff = elements.get(m * (m + 1) / 2);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (abs(b.get(i, j)) <= cutoff) {
                    b.set(i, j, 0);
                }
            }
        }

        int[] perm3 = null;
        int extra = 0;
        int[] choice2;

        PermutationGenerator gen2 = new PermutationGenerator(b.rows());

        PERM:
        while (((choice2 = gen2.next()) != null)) {
            for (int i = 0; i < m; i++) {
                for (int j = i; j < m; j++) {
                    final double c = b.get(choice2[i], choice2[j]);

                    if (c != 0) {
                        extra++;

                        cutoff = elements.get(m * (m + 1) / 2 + extra);

                        for (int q = 0; q < m; q++) {
                            for (int r = 0; r < m; r++) {
                                if (abs(b.get(q, r)) <= cutoff) {
                                    b.set(q, r, 0);
                                }
                            }
                        }

                        continue PERM;
                    }
                }
            }

            perm3 = Arrays.copyOf(choice2, choice2.length);
            break;
        }

        return perm3;
    }

    public void setPenaltyDiscount(double penaltyDiscount) {
        this.penaltyDiscount = penaltyDiscount;
    }

    public static class CausalOrder {
        private int[] perm;

        CausalOrder(int[] perm) {
            this.perm = perm;
        }

        int[] getPerm() {
            return perm;
        }
    }

}

