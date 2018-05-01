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
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.*;
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.independence.ConditionalGaussianLRT;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.score.ConditionalGaussianBicScore;
import edu.cmu.tetrad.algcomparison.score.SemBicScore;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulation;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.intervention.CGISimulation;
import edu.cmu.tetrad.intervention.algorithms.*;
import edu.cmu.tetrad.util.Parameters;

/**
 * An example script to simulate data and run a comparison analysis on it.
 *
 * @author jdramsey
 */
public class ExampleCompareInterventionSimulation {
    public static void main(String... args) {
        Parameters parameters = new Parameters();
        parameters.set("numRuns", 10);
        parameters.set("numMeasures", 100);
        parameters.set("numLatents", 0);
        parameters.set("sampleSize", 200);
        parameters.set("percentDiscrete", 0);
        parameters.set("minCategories", 2);
        parameters.set("maxCategories", 5);
        parameters.set("differentGraphs", true);

        parameters.set("avgDegree", 2);
        parameters.set("maxDegree", 100);
        parameters.set("maxIndegree", 100);
        parameters.set("maxOutdegree", 100);

        parameters.set("meanLow", 0);
        parameters.set("meanHigh", 1);

        parameters.set("interventionSize", 200);
        parameters.set("numInterventions", 5);
        parameters.set("percentIDiscrete", 100);
        parameters.set("minICategories", 1);
        parameters.set("maxICategories", 3);
        parameters.set("minEffected", 1);
        parameters.set("maxEffected", 1);

        parameters.set("minPotency", 0.8);
        parameters.set("maxPotency", 1.0);

//        parameters.set("minPotency", 0);
//        parameters.set("maxPotency", 0);

        parameters.set("percentCInfluence", 0);
        parameters.set("discreteCInfluence", 0);
        parameters.set("continuousCInfluence", 0);

//        parameters.set("percentCInfluence", 0.1);
//        parameters.set("discreteCInfluence", 0.5);
//        parameters.set("continuousCInfluence", 1.0);

        parameters.set("observationCondition", 2);
        parameters.set("hand", 0);
        parameters.set("removeContext", false);

        parameters.set("alpha", 0.05);
        parameters.set("structurePrior", 1.0);
        parameters.set("discretize", false);
        parameters.set("faithfulnessAssumed", false);
        parameters.set("symmetricFirstStep", true);
        parameters.set("maxPathLength", 0);
        parameters.set("completeRuleSetUsed", false);

        Statistics statistics = new Statistics();

        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());
        statistics.add(new ElapsedTime());

        Algorithms algorithms = new Algorithms();

//        algorithms.add(new IPcStable(new ConditionalGaussianLRT()));
//        algorithms.add(new woIPcStable(new ConditionalGaussianLRT()));

        algorithms.add(new IFges(new ConditionalGaussianBicScore()));
        algorithms.add(new woIFges(new ConditionalGaussianBicScore()));

//        algorithms.add(new IFci(new ConditionalGaussianLRT()));
//        algorithms.add(new woIFci(new ConditionalGaussianLRT()));

//        algorithms.add(new IGfci(new ConditionalGaussianLRT(), new ConditionalGaussianBicScore()));
//        algorithms.add(new woIGfci(new ConditionalGaussianLRT(), new ConditionalGaussianBicScore()));

        Simulations simulations = new Simulations();
        simulations.add(new CGISimulation(new RandomForward()));

        Comparison comparison = new Comparison();
        comparison.setShowAlgorithmIndices(true);
        comparison.setShowSimulationIndices(true);
        comparison.setSortByUtility(false);
        comparison.setShowUtilities(false);
        comparison.setParallelized(false);
        comparison.setSaveGraphs(true);

        for (Simulation simulation : simulations.getSimulations()) {
            comparison.saveToFiles("comparison", simulation, parameters);
        }

        comparison.compareFromSimulations("comparison", simulations, algorithms, statistics, parameters);

    }
}




