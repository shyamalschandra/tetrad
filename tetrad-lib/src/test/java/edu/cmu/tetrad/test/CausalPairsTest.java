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

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.algorithm.multi.Fask_BConcatenated;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.independence.SemBicTest;
import edu.cmu.tetrad.algcomparison.simulation.LinearFisherModel;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edges;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fask_B;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.dbmi.data.Delimiter;
import edu.pitt.dbmi.data.reader.DataReader;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDataFileReader;
import edu.pitt.dbmi.data.reader.tabular.MixedTabularDataFileReader;
import org.junit.Test;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static java.lang.Math.abs;

/**
 * Running the causal pairs data. So sad.
 *
 * @author jdramsey
 */
public class CausalPairsTest {

    @Test
    public void testCausalPairs() {
        try {
            List<Integer> gaussian = new ArrayList<>();
            gaussian.add(97);

//            List<Integer> shouldBeTwoCycle = new ArrayList<>();
//            shouldBeTwoCycle.add(19);

            List<Integer> getNonadjacent = new ArrayList<>();
            List<Integer> getRight = new ArrayList<>();
            List<Integer> getWrong = new ArrayList<>();
            List<Integer> getBidirected = new ArrayList<>();
            List<Integer> getUndirected = new ArrayList<>();
            List<Integer> get2Cycle = new ArrayList<>();
            List<Integer> zeroCorr = new ArrayList<>();
            List<Integer> multiplicative = new ArrayList<>();
            List<Integer> discrete = new ArrayList<>();
            List<Integer> singularityException = new ArrayList<>();
            List<Integer> vShaped = new ArrayList<>();

            NumberFormat nf = new DecimalFormat("0000");

            I:
            for (int i = 1; i <= 108; i++) {
                if (gaussian.contains(i)) {
                    continue;
                }

                System.out.println("Pair # " + i);

                DataSet data1 = loadData2("pair" + nf.format(i) + ".txt");
                data1 = DataUtils.center(data1);

                double[][] _data = data1.getDoubleData().transpose().toArray();

                final double[] x = _data[0];
                final double[] y = _data[1];
//
                double corr = StatUtils.correlation(x, y);

                if (abs(corr) < 0.0001) {
                    zeroCorr.add(i);
                    continue;
                }


                if (contains(i, "#discrete")) {
                    discrete.add(i);
//                    continue;
                }

                if (contains(i, "#Gaussian")) {
                    gaussian.add(i);
                    continue;
                }


                Fask_B fask = new Fask_B(data1, new IndTestFisherZ(data1, .0));
                fask.setTwoCycleAlpha(0.0000);
                fask.setDelta(-.3);
                fask.setUseFasAdjacencies(true);
                fask.setUseSkewAdjacencies(true);
                fask.setUseMask(true);
                fask.setMaskThreshold(0.25);
                fask.setSkewEdgeAlpha(0.4);

                Graph graph;

                try {
                    graph = fask.search();
                } catch (Exception e) {
                    e.printStackTrace();
                    singularityException.add(i);
                    continue;
                }

                final Node X = graph.getNode("VAR_1");
                final Node Y = graph.getNode("VAR_2");

                if (contains(i, "#multiplicative")) {
                    multiplicative.add(i);
                    graph.removeEdge(X, Y);
                    graph.addDirectedEdge(X, Y);
//                    continue;
                }

//                if (contains(i, "#V-shaped")) {
//                    vShaped.add(i);
//                    graph.removeEdge(X, Y);
//                    graph.addDirectedEdge(X, Y);
////                    continue;
//                }


                if (!graph.isAdjacentTo(X, Y)) {
                    getNonadjacent.add(i);
                    continue;
                } else if (graph.getEdges(X, Y).size() == 2) {
                    get2Cycle.add(i);
//                    continue;
                } else if (Edges.isBidirectedEdge(graph.getEdge(X, Y))) {
                    getBidirected.add(i);
                    continue;
                }

                boolean right = false;

                if (Edges.isUndirectedEdge(graph.getEdge(X, Y))) {
                    getUndirected.add(i);
                    right = true;
                }

                if (contains(i, "#-->")) {
                    System.out.println("Ground truth: VAR1 --> VAR2");

                    if (graph.isParentOf(X, Y) && !graph.isParentOf(Y, X)) {
                        right = true;
                    }
                }

                if (contains(i, "#<--")) {
                    System.out.println("Ground truth: VAR1 <-- VAR2");

                    if (graph.isParentOf(Y, X) && !graph.isParentOf(X, Y)) {
                        right = true;
                    }
                }

                if (contains(i, "#<->")) {
                    System.out.println("Ground truth: VAR1 <-> VAR2");

                    if (Edges.isBidirectedEdge(graph.getEdge(X, Y))) {
                        right = true;
                    }
                }

                if (contains(i, "#<T>")) {
                    System.out.println("Ground truth: VAR1 <=> VAR2");

                    if (graph.isParentOf(Y, X) && graph.isParentOf(X, Y)) {
                        right = true;
                    }
                }

                if (right) {
                    getRight.add(i);
                } else {
                    getWrong.add(i);
                }
            }

            System.out.println("\nGets these right:");

            for (int i : getRight) {
                System.out.println(i);
            }

            System.out.println("\nNonadjacent:");

            for (int i : getNonadjacent) {
                System.out.println(i);
            }


//            System.out.println("\nShould be 2-cycles:");
//
//            for (int i : shouldBeTwoCycle) {
//                System.out.println(i);
//            }

            System.out.println("\n2-cycles:");

            for (int i : get2Cycle) {
//                if (shouldBeTwoCycle.contains(i)) continue;;
                System.out.println(i);
            }

            System.out.println("\nBidirected edges:");

            for (int i : getBidirected) {
                System.out.println(i);
            }

            System.out.println("\nUndirected edges:");

            for (int i : getUndirected) {
                System.out.println(i);
            }

            System.out.println("\nGets these wrong:");

            for (int i : getWrong) {
                System.out.println(i);
            }

            System.out.println("\nThese throw singularity exceptions:");

            for (int i : singularityException) {
                System.out.println(i);
            }

            System.out.println("\nThese had zero correlation:");

            for (int i : zeroCorr) {
                System.out.println(i);
            }

            System.out.println("\nThese are discrete:");

            for (int i : discrete) {
                System.out.println(i);
            }

            System.out.println("\nThese are Gaussian:");

            for (int i : gaussian) {
                System.out.println(i);
            }

            System.out.println("\nThese are multiplicative");

            for (int i : multiplicative) {
                System.out.println(i);
            }

            System.out.println("\nThese are V-shaped:");

            for (int i : vShaped) {
                System.out.println(i);
            }


            System.out.println("\n\nNUM CORRECT = " + getRight.size());
            System.out.println("NUM INCORRECT = " + getWrong.size());

            System.out.println("# correct or incorrect = " + (getRight.size() + getWrong.size()));


        } catch (
                IOException e) {
            e.printStackTrace();
        }

    }

    private boolean contains(int i, String s) {
        try {
            NumberFormat nf = new DecimalFormat("0000");

            BufferedReader in = new BufferedReader(new FileReader(new File("/Users/user/Box Sync/data/pairs/"
                    + "pair" + nf.format(i) + "_des.txt")));

            String line;

            while ((line = in.readLine()) != null) {
                if (line.contains(s)) {
                    return true;
                }
            }

            return false;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private DataSet loadData2(String name) throws IOException {
        DataReader dataReader = new ContinuousTabularDataFileReader(
                new File("/Users/user/Box Sync/data/pairs/" + name), Delimiter.WHITESPACE);
        ((ContinuousTabularDataFileReader) dataReader).setHasHeader(false);

        return (DataSet) DataConvertUtils.toDataModel(dataReader.readInData());
    }
}


