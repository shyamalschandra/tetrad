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
import edu.cmu.tetrad.algcomparison.independence.SemBicTest;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
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
import java.util.*;

import static java.lang.Math.abs;
import static java.lang.Math.log;

/**
 * Pulling this test out for Madelyn.
 *
 * @author jdramsey
 */
public class TestSimulatedFmri {

    private void task(boolean testing) {
        Parameters parameters = new Parameters();
        parameters.set("penaltyDiscount", 1);
        parameters.set("twoCycleAlpha", .0000001);
        parameters.set("faskDelta", -.2);
        parameters.set("depth", 5);
        parameters.set("extraEdgeThreshold", 10);
        parameters.set("maskThreshold", 10);
        parameters.set("skewEdgeAlpha", 0.01);
        parameters.set("errorsPositivelySkewed", true);


        testing = false;

        parameters.set("numRuns", 5);
        parameters.set("randomSelectionSize", 3);

        parameters.set("useFasAdjacencies", true);
        parameters.set("useCorrDiffAdjacencies", true);

        parameters.set("Structure", "Placeholder");

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn("Structure"));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
//        statistics.add(new MathewsCorrAdj());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadPrecisionCommonAdjacencies());
        statistics.add(new ArrowheadRecall());
        statistics.add(new TwoCyclePrecision());
        statistics.add(new TwoCycleRecall());
        statistics.add(new TwoCycleFalsePositive());
        statistics.add(new TwoCycleFalseNegative());
        statistics.add(new TwoCycleTruePositive());
        statistics.add(new ElapsedTime());
        statistics.setWeight("AHR", 1.0);
        statistics.setWeight("2CP", 1.0);
        statistics.setWeight("2CR", 1.0);
        statistics.setWeight("2CFP", 1.0);

        Simulations simulations = new Simulations();

        if (!testing) {
            String dir = "/Users/user/Downloads/Cycles_Data_fMRI/";
            String subdir = "data_fslfilter";

            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network1_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network2_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network3_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network4_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_contr", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_contr_p2n6", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_contr_p6n2", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network6_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network6_contr", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network7_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network7_contr", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network8_amp_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network8_amp_contr", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network8_contr_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network9_amp_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network9_amp_contr", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network9_contr_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Diamond", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Markov_Complex_1", subdir));
        } else {

            String dir = "/Users/user/Downloads/CyclesTestingData/";
            String subdir = "data_fslfilter";

            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network1_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network2_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network3_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network4_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_cont", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_cont_p3n7", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network5_cont_p7n3", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network6_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network6_cont", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network7_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network7_cont", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network8_amp_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network8_amp_cont", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network8_cont_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network9_amp_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network9_amp_cont", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Network9_cont_amp", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Diamond", subdir));
            simulations.add(new LoadContinuousDataAndSingleGraph(
                    dir + "Markov_Complex_1", subdir));
        }

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new FaskConcatenated(new SemBicScore()));
        algorithms.add(new Fask_BConcatenated(new SemBicTest()));
//        algorithms.add(new SkewSearchConcatenated(new FisherZSkew()));
//        algorithms.add(new FaskConcatenated(new SemBicScore()));
//        algorithms.add(new Fask_BConcatenated(new SemBicTest()));
//        algorithms.add(new FaskConcatenated(new SemBicScore()));
//        algorithms.add(new FaskGfciConcatenated(new SemBicTest()));

//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.RSkew));

        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(true);
        comparison.setSortByUtility(false);
        comparison.setShowUtilities(false);
        comparison.setParallelized(false);
        comparison.setSaveGraphs(false);
        comparison.setTabDelimitedTables(false);
        comparison.setSaveGraphs(true);

        String directory;

        if (!testing) {
            directory = "comparison_training";
        } else {
            directory = "comparison_testing";
        }

        comparison.compareFromSimulations(directory, simulations, algorithms, statistics, parameters);
    }

    @Test
    public void task2() {
        Parameters parameters = new Parameters();
        parameters.set("penaltyDiscount", 1);
        parameters.set("twoCycleAlpha", .0000001);
        parameters.set("faskDelta", -.2);
        parameters.set("depth", 5);
        parameters.set("extraEdgeThreshold", 10);
        parameters.set("maskThreshold", 10);
        parameters.set("skewEdgeAlpha", 0.01);
        parameters.set("errorsPositivelySkewed", true);

        parameters.set("numRuns", 5);
        parameters.set("randomSelectionSize", 5);

        parameters.set("Structure", "Placeholder");

        Statistics statistics = new Statistics();

        statistics.add(new ParameterColumn("Structure"));
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
//        statistics.add(new MathewsCorrAdj());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadPrecisionCommonAdjacencies());
        statistics.add(new ArrowheadRecall());
        statistics.add(new TwoCyclePrecision());
        statistics.add(new TwoCycleRecall());
        statistics.add(new TwoCycleFalsePositive());
        statistics.add(new TwoCycleFalseNegative());
        statistics.add(new TwoCycleTruePositive());
        statistics.add(new ElapsedTime());
        statistics.setWeight("AHR", 1.0);
        statistics.setWeight("2CP", 1.0);
        statistics.setWeight("2CR", 1.0);
        statistics.setWeight("2CFP", 1.0);

        Simulations simulations = new Simulations();

        Algorithms algorithms = new Algorithms();

        for (int i = 1; i <= 28; i++) {
            if (i == 21) continue;
//            simulations.add(new LoadContinuousDataSmithSim("/Users/user/Downloads/smithsim/", i));
//            simulations.add(new LoadContinuousDataPwdd7("/Users/user/Downloads/pwdd7/", i, "50_BOLDdemefilt1"));
            simulations.add(new LoadContinuousDataPwdd7("/Users/user/Downloads/pwdd7/", i, "50_BOLDnoise"));
        }

//        algorithms.add(new LofsConcatenated(Lofs2.Rule.FASKLR));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.R1));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.R3));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.RSkew));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.RSkewE));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.Skew));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.SkewE));
//        algorithms.add(new LofsConcatenated(Lofs2.Rule.Patel));

//        algorithms.add(new FaskConcatenated(new SemBicScore()));
        algorithms.add(new Fask_BConcatenated(new SemBicTest()));
//        algorithms.add(new SkewSearchConcatenated());

//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.R1));
//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.Patel));
//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.R3));
//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.Skew));
//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.RSkew));
//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.RSkewE));
//        algorithms.add(new FasLofsConcatenated(Lofs2.Rule.SkewE));

        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(true);
        comparison.setSortByUtility(false);
        comparison.setShowUtilities(false);
        comparison.setParallelized(false);
        comparison.setSaveGraphs(false);
        comparison.setTabDelimitedTables(false);
        comparison.setSaveGraphs(true);

        String directory = "smithsim";

        comparison.compareFromSimulations(directory, simulations, algorithms, statistics, parameters);
    }

    @Test
    public void test4celllinesdata() {
        try {
            DataSet data1 = loadData("BT20_excluded1.csv");
            DataSet data2 = loadData("BT549_excluded1.csv");
            DataSet data3 = loadData("MCF7_excluded1.csv");
            DataSet data4 = loadData("UACC812_excluded1.csv");

            List<Integer> continuousColumns = new ArrayList<>();
            for (Node node : data1.getVariables()) {
                if (node instanceof ContinuousVariable) {
                    continuousColumns.add(data1.getColumn(node));
                }
            }
            int[] _continuouscolumns = new int[continuousColumns.size()];
            for (int i = 0; i < continuousColumns.size(); i++) _continuouscolumns[i] = continuousColumns.get(i);

            DataSet data1a = data1.subsetColumns(_continuouscolumns);
            DataSet data2a = data2.subsetColumns(_continuouscolumns);
            DataSet data3a = data3.subsetColumns(_continuouscolumns);
            DataSet data4a = data4.subsetColumns(_continuouscolumns);

            DataSet concatenated = DataUtils.concatenate(data1a, data2a, data3a, data4a);

            List<String> allCategories = new ArrayList<>();

            allCategories.addAll(addCategories(data1));
            allCategories.addAll(addCategories(data2));
            allCategories.addAll(addCategories(data3));
            allCategories.addAll(addCategories(data4));

            System.out.println(allCategories);

            List<Node> indicators = new ArrayList<>();
            Map<String, ContinuousVariable> hash = new HashMap<>();

            for (String category : allCategories) {
                final ContinuousVariable e = new ContinuousVariable(category);
                indicators.add(e);
                if (concatenated.getVariable(e.getName()) == null) {
                    concatenated.addVariable(e);
                    hash.put(category, e);
                }
            }

            for (String s : hash.keySet()) {
                Node var = concatenated.getVariable(s);
                int c = concatenated.getColumn(var);

                for (int r = 0; r < concatenated.getNumRows(); r++) {
                    concatenated.setDouble(r, c, 0.0 + RandomUtil.getInstance().nextNormal(0, 0.001));
                }
            }

            int r = -1;

            r = addIndicatorData(data1, concatenated, 4, hash, r);

            r = addIndicatorData(data2, concatenated, 4, hash, r);

            r = addIndicatorData(data3, concatenated, 4, hash, r);

            r = addIndicatorData(data4, concatenated, 4, hash, r);

            System.out.println(concatenated);

            PrintStream out = new PrintStream("/Users/user/Box Sync/data/4cellLineData/4celllines.with.indicators.jittered.txt");
            out.println(concatenated);
            out.flush();

            String[] delays = {"5min", "15min", "30min", "60min", "2hr", "4hr"};
            List<DataSet> allDelays = new ArrayList<>();

            for (String delay : delays) {
                Node _delay = hash.get(delay);
                int hr = concatenated.getColumn(_delay);
                List<Integer> rows = new ArrayList<>();

                for (int r2 = 0; r2 < concatenated.getNumRows(); r2++) {
                    if (abs(concatenated.getDouble(r2, hr) - 1.0) < 0.1) {
                        rows.add(r2);
                    }
                }

                int[] _rows = new int[rows.size()];
                for (int i = 0; i < rows.size(); i++) _rows[i] = rows.get(i);

                DataSet subset = concatenated.subsetRows(_rows);

                allDelays.add(subset);

                System.out.println(subset.getNumRows());

                PrintStream out2 = new PrintStream("/Users/user/Box Sync/data/4cellLineData/4celllines.with.indicators." + delay + ".jittered.txt");
                out2.println(subset);
                out2.flush();


            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testCausalPairs() {
        try {
            List<Integer> multiplicativeError = new ArrayList<>();
//            multiplicativeError.add(1); // multiplicative error
            multiplicativeError.add(2); // multiplicative error
//            multiplicativeError.add(4); // multiplicative error
            multiplicativeError.add(5); // multiplicative error
            multiplicativeError.add(6); // multiplicative error
            multiplicativeError.add(7); // multiplicative error
            multiplicativeError.add(9); // multiplicative error
            multiplicativeError.add(10); // multiplicative error
            multiplicativeError.add(11); // multiplicative error
            multiplicativeError.add(12); // multiplicative error
            multiplicativeError.add(17); // multiplicative error
            multiplicativeError.add(25); // multiplicative error ?
            multiplicativeError.add(80); // multiplicative error ?

            List<Integer> gaussian = new ArrayList<>();
            gaussian.add(97);

//            List<Integer> shouldBeTwoCycle = new ArrayList<>();
//            shouldBeTwoCycle.add(19);

            List<Integer> getNonadjacent = new ArrayList<>();
            List<Integer> getRight = new ArrayList<>();
            List<Integer> getWrong = new ArrayList<>();
            List<Integer> getBidirected = new ArrayList<>();
            List<Integer> get2Cycle = new ArrayList<>();
            List<Integer> deterministic = new ArrayList<>();
            List<Integer> zeroCorr = new ArrayList<>();
            List<Integer> imbalanced = new ArrayList<>();

            I:
            for (int i = 1; i <= 108; i++) {
//                if (multiplicativeError.contains(i)) {
//                    continue;
//                }

                if (gaussian.contains(i)) {
                    continue;
                }

                System.out.println("Pair # " + i);

                NumberFormat nf = new DecimalFormat("0000");

                DataSet data1 = loadData2("pair" + nf.format(i) + ".txt");
                data1 = DataUtils.center(data1);

                double[][] _data = data1.getDoubleData().transpose().toArray();

//                for (int r = 0; r < data1.getNumRows(); r++) {
//                    for (int c = 0; c < data1.getNumColumns(); c++) {
//                        data1.setDouble(r, c, log(data1.getDouble(r, c)));
//                    }
//                }

                final double[] x = _data[0];
                final double[] y = _data[1];
//
//                if (n(x) < 0.4 * x.length) {
//                    imbalanced.add(i);
//                    continue;
//                }
//                if (n(y) < 0.4 * x.length) {
//                    imbalanced.add(i);
//                    continue;
//                }

                double corr = StatUtils.correlation(x, y);

                if (abs(corr) < 0.05) {
                    zeroCorr.add(i);
                    continue;
                }


                Fask_B fask = new Fask_B(data1, new IndTestFisherZ(data1, 0.05));
                fask.setTwoCycleAlpha(0.000);
                fask.setAssumeErrorsPositivelySkewed(false);
                fask.setDelta(-.2);

                Graph graph;
                try {
                    graph = fask.search();
                } catch (Exception e) {
                    e.printStackTrace();
                    deterministic.add(i);
                    continue;
                }

                final Node X = graph.getNode("VAR_1");
                final Node Y = graph.getNode("VAR_2");

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

                if (multiplicativeError.contains(i)) {
                    Edge edge = graph.getEdge(X, Y);
                    graph.removeEdge(edge);
                    graph.addEdge(edge.reverse());
                }

                BufferedReader in = new BufferedReader(new FileReader(new File("/Users/user/Box Sync/data/pairs/"
                        + "pair" + nf.format(i) + "_des.txt")));

                String line;

                while ((line = in.readLine()) != null) {

                    if (line.contains("#-->")) {
                        System.out.println("Ground truth: VAR1 --> VAR2");

                        if (graph.isParentOf(X, Y)) {
                            getRight.add(i);
                            continue I;
                        }
                    } else if (line.contains("#<--")) {
                        System.out.println("Ground truth: VAR1 <-- VAR2");

                        if (graph.isParentOf(Y, X)) {
                            getRight.add(i);
                            continue I;
                        }
                    } else if (line.contains("#<=>")) {
                        System.out.println("Ground truth: VAR1 <=> VAR2");

                        if (graph.isParentOf(Y, X) && graph.isParentOf(X, Y)) {
                            getRight.add(i);
                            continue I;
                        }
                    }
                }

                getWrong.add(i);
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

            System.out.println("\nGets these wrong:");

            for (int i : getWrong) {
                System.out.println(i);
            }

            System.out.println("\nThese were deterministic:");

            for (int i : deterministic) {
                System.out.println(i);
            }

            System.out.println("\nThese had zero correlation:");

            for (int i : zeroCorr) {
                System.out.println(i);
            }

            System.out.println("\nImbalanced centering:");

            for (int i : imbalanced) {
                System.out.println(i);
            }

            System.out.println("\nMultiplicative error:");

            for (int i : multiplicativeError) {
                System.out.println(i);
            }

            System.out.println("\nGaussian:");

            for (int i : gaussian) {
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

    private static double n(double[] x) {
        int n = 0;

        for (int k = 0; k < x.length; k++) {
            if (x[k] > 0) {
                n++;
            }
        }

        return n;
    }

    private int addIndicatorData(DataSet data1, DataSet concatenated, int discreteColumns, Map<String, ContinuousVariable> hash, int r) {
        for (int i = 0; i < data1.getNumRows(); i++) {
            ++r;

            for (int j = 0; j < discreteColumns; j++) {
                int v = data1.getInt(i, j);

                if (v != -99) {
                    DiscreteVariable var = (DiscreteVariable) data1.getVariable(j);
                    String cat = var.getCategory(v);
                    Node n = hash.get(cat);
                    int cc = concatenated.getColumn(n);
                    concatenated.setDouble(r, cc, 1.0 + RandomUtil.getInstance().nextNormal(0, 0.001));
                }
            }

        }

        return r;
    }

    private List<String> addCategories(DataSet data1) {
        List<String> allCategories = new ArrayList<>();
        for (Node node : data1.getVariables()) {
            if (node instanceof DiscreteVariable) {
                allCategories.addAll(((DiscreteVariable) node).getCategories());
            }
        }
        return allCategories;
    }

    private DataSet loadData(String name) throws IOException {
        DataReader dataReader = new MixedTabularDataFileReader(20
                , new File("/Users/user/Box Sync/data/4cellLineData/" + name), Delimiter.COMMA);
        ((MixedTabularDataFileReader) dataReader).setHasHeader(true);

        return (DataSet) DataConvertUtils.toDataModel(dataReader.readInData());
    }

    private DataSet loadData2(String name) throws IOException {
        DataReader dataReader = new ContinuousTabularDataFileReader(
                new File("/Users/user/Box Sync/data/pairs/" + name), Delimiter.WHITESPACE);
        ((ContinuousTabularDataFileReader) dataReader).setHasHeader(false);

        return (DataSet) DataConvertUtils.toDataModel(dataReader.readInData());
    }

    public static void main(String... args) {
        new TestSimulatedFmri().task(false);
    }
}



