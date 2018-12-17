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

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edges;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fask_B;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.util.DataConvertUtils;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.dbmi.data.Delimiter;
import edu.pitt.dbmi.data.reader.DataReader;
import edu.pitt.dbmi.data.reader.tabular.ContinuousTabularDataFileReader;
import org.junit.Test;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

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

                final boolean xNormal = new AndersonDarlingTest(x).getP() > 0.05;
                final boolean yNormal = new AndersonDarlingTest(y).getP() > 0.05;
                System.out.println("x normal = " + xNormal + " y normal = " + yNormal);
//
                double corr = StatUtils.correlation(x, y);

                if (abs(corr) < 0.0001) {
                    zeroCorr.add(i);
//                    continue;
                }


                if (contains(i, "#discrete")) {
                    discrete.add(i);
//                    continue;
                }

//                if (contains(i, "#Gaussian")) {
//                    gaussian.add(i);
////                    continue;
//                }

                if (xNormal && yNormal) {
                    gaussian.add(i);
//                    continue;
                }


                Fask_B fask = new Fask_B(data1, new IndTestFisherZ(data1, .001));
                fask.setTwoCycleAlpha(0.000);
//                fask.setDelta(-.05);
                fask.setUseFasAdjacencies(true);
                fask.setUseSkewAdjacencies(true);
                fask.setUseMask(false);
                fask.setMaskThreshold(0.001);
                fask.setSkewEdgeAlpha(.05);
                fask.setVerbose(true);

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
                    graph.removeEdges(X, Y);
                    graph.addDirectedEdge(X, Y);
//                    continue;
                }

                if (contains(i, "#V-shaped")) {
                    vShaped.add(i);
                    graph.removeEdges(X, Y);
                    graph.addDirectedEdge(X, Y);
//                    continue;
                }


                if (!graph.isAdjacentTo(X, Y)) {
                    getNonadjacent.add(i);
                    continue;
                } else if (graph.getEdges(X, Y).size() == 2) {
                    get2Cycle.add(i);
//                    continue;
                } else if (Edges.isBidirectedEdge(graph.getEdge(X, Y))) {
                    getBidirected.add(i);
//                    continue;
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

            printList("Gets these right:", getRight);
            printList("Nonadjacent:", getNonadjacent);
            printList("2-cycles:", get2Cycle);
            printList("Bidirected edges:", getBidirected);
            printList("Undirected edges:", getUndirected);
            printList("Gets these wrong:", getWrong);
            printList("These throw singularity exceptions:", singularityException);
            printList("These had zero correlation:", zeroCorr);
            printList("These are discrete:", discrete);
            printList("These are Gaussian:", gaussian);
            printList("These are multiplicative", multiplicative);
            printList("These are V-shaped:", vShaped);

            System.out.println("\n\nNUM CORRECT = " + getRight.size());
            System.out.println("NUM INCORRECT = " + getWrong.size());

            System.out.println("# correct or incorrect = " + (getRight.size() + getWrong.size()));


        } catch (
                IOException e) {
            e.printStackTrace();
        }

    }

    private void printList(String s, List<Integer> vShaped) {
        System.out.println("\n" + s);

        for (int i : vShaped) {
            System.out.println(i);
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


