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

package edu.cmu.tetrad.algcomparison.examples;

import edu.cmu.tetrad.algcomparison.Comparison;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithms;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.PcAll;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.*;
import edu.cmu.tetrad.algcomparison.simulation.GeneralSemSimulationExample4;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.util.Parameters;

/**
 * An example script to simulate data and run a comparison analysis on it.
 *
 * @author jdramsey
 */
public class Example4 {

    public static void main(String... args) {
        Parameters parameters = new Parameters();

        parameters.set("numRuns", 1);
        parameters.set("numMeasures", 20);
        parameters.set("avgDegree", 2);
        parameters.set("numLatents", 0);
        parameters.set("differentGraphs", true);
        parameters.set("sampleSize", 1000);
        parameters.set("fastFDR", true);
        parameters.set("verbose", true);
        parameters.set("alpha", 0.001);
        parameters.set("kciAlpha", 0.05);
        parameters.set("penaltyDiscount", 1);
        parameters.set("maxDegree", 4);
        parameters.set("colliderDiscoveryRule", 1);
        parameters.set("conflictRule", 1);
        parameters.set("useMaxPOrientationHeuristic", false);
        parameters.set("maxPOrientationMaxPathLength", 3);
        parameters.set("depth", -1);
        parameters.set("stableFAS", false);
        parameters.set("concurrentFAS", false);
        parameters.set("cciScoreAlpha", .01);
        parameters.set("numBasisFunctions", 30);
        parameters.set("kernelType", 2);
        parameters.set("kernelMultiplier", 1.0);
        parameters.set("basisType", 2);
        parameters.set("kernelRegressionSampleSize", 100);
        parameters.set("numDependenceSpotChecks", 0);

        String function = "TSUM($)";
        parameters.set("generalSemFunctionTemplateMeasured", function);
        parameters.set("generalSemFunctionTemplateLatent", function);
        parameters.set("generalSemErrorTemplate", "U(-.2, .2)");

//        parameters.set("percentDiscrete", 50);

        Statistics statistics = new Statistics();

//        statistics.add(new ParameterColumn("colliderDiscoveryRule"));
//        statistics.add(new ParameterColumn("conflictRule"));
//
//        statistics.add(new ParameterColumn("alpha"));
//        statistics.add(new ParameterColumn("cciScoreAlpha"));
//        statistics.add(new ParameterColumn("conflictRule"));
//        statistics.add(new ParameterColumn("numBasisFunctions"));
//        statistics.add(new ParameterColumn("kernelType"));
//        statistics.add(new ParameterColumn("kernelMultiplier"));
//        statistics.add(new ParameterColumn("basisType"));

//        statistics.add(new AdjacencyPrecision());
//        statistics.add(new AdjacencyRecall());
//        statistics.add(new ArrowheadPrecision());
//        statistics.add(new ArrowheadRecall());

        statistics.add(new AdjacencyTP());
        statistics.add(new AdjacencyFP());
        statistics.add(new AdjacencyFN());
        statistics.add(new ArrowheadTP());
        statistics.add(new ArrowheadFP());
        statistics.add(new ArrowheadFN());

//        statistics.add(new MathewsCorrAdj());
//        statistics.add(new MathewsCorrArrow());
//        statistics.add(new F1Adj());
//        statistics.add(new F1Arrow());
        statistics.add(new SHD());
        statistics.add(new ElapsedTime());

//        statistics.setWeight("AP", 1.0);
//        statistics.setWeight("AR", 0.5);
//        statistics.setWeight("AHP", 1.0);
//        statistics.setWeight("AHR", 0.5);
//        statistics.setWeight("AHFP", 1.0);
//        statistics.setWeight("AHFN", 0.5);
        statistics.setWeight("SHD", 1.0);
//        statistics.setWeight("AHP", .5);


        //==============ALGORIHTHMS============//
        Algorithms algorithms = new Algorithms();

//        algorithms.add(new PcAll(new Kci()));
//        algorithms.add(new PcAll(new KciMatlab()));

//        algorithms.add(new PcAll(new ResidualCITMatlab()));
//        algorithms.add(new/Library/Frameworks/R.Framework/Libraries PcAll(new FcitJRI()));
////
//        algorithms.add(new PcAll(new RcotJRI()));
//        algorithms.add(new PcAll(new RcitJRI()));
        algorithms.add(new PcAll(new CciTest()));
//        algorithms.add(new Fges(new CciScore()));
//        algorithms.add(new PcAll(new FcitJRI()));

////        algorithms.add(new PcAll(new ConditionalGaussianLRT()));
////        algorithms.add(new PcAll(new FisherZ()));
////        algorithms.add(new PcAll(new SemBicTest()));
//
////        algorithms.add(new Fci(new CciTest()));
//
        Simulations simulations = new Simulations();
        simulations.add(new GeneralSemSimulationExample4(new RandomForward()));
        Comparison comparison = new Comparison();

        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(false);
        comparison.setSortByUtility(true);
        comparison.setShowUtilities(false);
        comparison.setParallelized(false);
        comparison.setComparisonGraph(Comparison.ComparisonGraph.true_DAG);

        comparison.setSaveGraphs(true);
//        comparison.setSavePatterns(true);
        comparison.setSavePags(true);

        final String dir = "example4";
//        comparison.saveToFiles(dir, simulations.getSimulations().get(0), parameters);
//        comparison.compareFromFiles(dir, dir, algorithms, statistics, parameters);
//
        comparison.compareFromSimulations(dir, simulations, algorithms, statistics, parameters);
    }
}




