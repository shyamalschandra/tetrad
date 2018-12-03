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
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Fci;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Fci_Sachs;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Gfci;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Gfci_Sachs;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.Fges_Sachs;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll_Sachs;
import edu.cmu.tetrad.algcomparison.independence.*;
import edu.cmu.tetrad.algcomparison.score.*;
import edu.cmu.tetrad.algcomparison.simulation.LoadDataAndGraphs;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.search.PcAll;
import edu.cmu.tetrad.util.Parameters;
import org.junit.Test;

/**
 * Pulling this test out for Madelyn.
 *
 * @author jdramsey
 */
public class TestSachs {

    private void task() {
        Parameters parameters = new Parameters();
        parameters.set("structurePrior", 1);
        parameters.set("discretize", false);
        parameters.set("penaltyDiscount", 1, 2, 4, 8);
        parameters.set("alpha", 0.1, 0.01, 0.001, 0.0001);
        parameters.set("depth", -1);
        parameters.set("maxDistinctValuesDiscrete", 2);

        // PC-Stable
        parameters.set("stableFAS", true);
        parameters.set("concurrentFAS", false);
        parameters.set("colliderDiscoveryRule", 1);
        parameters.set("conflictRule", 3);

        // CCI
        parameters.set("cciAlpha", 0.001);
        parameters.set("cciScoreAlpha", 0.001);

        // KCI
        parameters.set("kciAlpha", 0.001);
        parameters.set("kciScoreAlpha", 0.001);

        // Resampling
        parameters.set("numberResampling", 100);
        parameters.set("resampleSize", 7466);
        parameters.set("resamplingWithReplacement", true);
        parameters.set("resamplingEnsemble", 0);

        Statistics statistics = new Statistics();

        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadPrecisionCommonAdjacencies());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ElapsedTime());

        Simulations simulations = new Simulations();

        String dir = "/home/bandrews/Desktop/fask/Sachs/";
        String subdir = "data";

        simulations.add(new LoadDataAndGraphs(dir + "mixed"));
//        simulations.add(new LoadContinuousDataAndSingleGraph(dir + "jittered", subdir));
//        simulations.add(new LoadContinuousDataAndSingleGraph(dir + "jittered_subsampled", subdir));

        Algorithms algorithms = new Algorithms();

        algorithms.add(new Fges_Sachs(new ConditionalGaussianBicScore()));
//        algorithms.add(new Fges_Sachs(new SemBicScore()));
//        algorithms.add(new Fges_Sachs(new CciScore()));
//        algorithms.add(new Fges_Sachs(new KciMatlabScore()));
        algorithms.add(new PcAll_Sachs(new ConditionalGaussianLRT()));
//        algorithms.add(new PcAll_Sachs(new SemBicTest()));
//        algorithms.add(new PcAll_Sachs(new FisherZ()));
//        algorithms.add(new PcAll_Sachs(new CciTest()));
//        algorithms.add(new PcAll_Sachs(new KciMatlab()));

        algorithms.add(new Gfci_Sachs(new ConditionalGaussianLRT(), new ConditionalGaussianBicScore()));
//        algorithms.add(new Gfci_Sachs(new SemBicTest(), new SemBicScore()));
//        algorithms.add(new Gfci_Sachs(new FisherZ(), new SemBicScore()));
//        algorithms.add(new Gfci_Sachs(new CciTest(), new CciScore()));
//        algorithms.add(new Gfci_Sachs(new KciMatlab(), new KciMatlabScore()));
        algorithms.add(new Fci_Sachs(new ConditionalGaussianLRT()));
//        algorithms.add(new Fci_Sachs(new SemBicTest()));
//        algorithms.add(new Fci_Sachs(new FisherZ()));
//        algorithms.add(new Fci_Sachs(new CciTest()));
//        algorithms.add(new Fci_Sachs(new KciMatlab()));



        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(false);
        comparison.setSortByUtility(false);
        comparison.setShowUtilities(false);
        comparison.setParallelized(false);
        comparison.setSaveGraphs(false);
        comparison.setTabDelimitedTables(false);
        comparison.setSaveGraphs(true);

        String directory = "comparison_sachs";

        comparison.compareFromSimulations(directory, simulations, algorithms, statistics, parameters);
    }


    public static void main(String... args) {
        new TestSachs().task();
    }
}



