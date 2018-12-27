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
import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.independence.SemBicTest;
import edu.cmu.tetrad.algcomparison.simulation.LinearFisherModel;
import edu.cmu.tetrad.algcomparison.simulation.Simulations;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.util.Parameters;
import org.junit.Test;

import static java.lang.Math.abs;

/**
 * Pulling this test out for Madelyn.
 *
 * @author jdramsey
 */
public class ToyFaskBExample {

    @Test
    public void toy_positiveSkews() {
        testToyExample(true);
    }

    @Test
    public void toy_negativeSkews() {
        testToyExample(false);
    }

    @Test
    public void toy_mixedSkews() {
        testToyExample(true);
    }

    public void testToyExample(boolean positiveCoefs) {
        RandomGraph graph = new RandomForward();

        Simulations simulations = new Simulations();
        simulations.add(new LinearFisherModel(graph));

        Algorithms algorithms = new Algorithms();
        algorithms.add(new edu.cmu.tetrad.algcomparison.algorithm.multi.Fask_B(new SemBicTest()));

        Statistics statistics = new Statistics();
        statistics.add(new AdjacencyPrecision());
        statistics.add(new AdjacencyRecall());
        statistics.add(new ArrowheadPrecision());
        statistics.add(new ArrowheadRecall());

        Parameters parameters = new Parameters();

        parameters.set("numMeasures", 10);
        parameters.set("numLatents", 0);
        parameters.set("avgDegree", 2);
        parameters.set("maxDegree", 100);
        parameters.set("maxIndegree", 100);
        parameters.set("maxOutdegree", 100);
        parameters.set("connected", false);

//        parameters.add("coefLow");
//        parameters.add("coefHigh");
//        parameters.add("varLow");
//        parameters.add("varHigh");
        parameters.set("verbose", true);
        parameters.set("includePositiveCoefs", positiveCoefs);
        parameters.set("includeNegativeCoefs", !positiveCoefs);
        parameters.set("errorsNormal", false);
        parameters.set("betaLeftValue", 6);
        parameters.set("betaRightValue", 10);
        parameters.set("numRuns", 10);
//        parameters.add("percentDiscrete");
//        parameters.add("numCategories");
//        parameters.add("differentGraphs");
        parameters.set("sampleSize", 1000);
        parameters.set("intervalBetweenShocks", 20);
        parameters.set("intervalBetweenRecordings", 20);
        parameters.set("selfLoopCoef", 0.0);
        parameters.set("fisherEpsilon", 0.0001);
//        parameters.add("randomizeColumns");
//        parameters.add("measurementVariance");
//        parameters.add("saveLatentVars");

        parameters.set("depth", -1);
        parameters.set("skewEdgeAlpha", 0.0001);
        parameters.set("twoCycleAlpha", 0);

        parameters.set("useFasAdjacencies", true);
        parameters.set("useSkewAdjacencies", true);
        parameters.set("useMask", true);
        parameters.set("maskThreshold", 0.3);
        parameters.set("empirical", true);

        // Bootstrapping
//        parameters.add("numberResampling");
//        parameters.add("percentResampleSize");
//        parameters.add("resamplingWithReplacement");
//        parameters.add("resamplingEnsemble");
//        parameters.add("verbose");

        String type;

        if (positiveCoefs) {
            type = "_positive_coefs";
        } else {
            type = "_negative_coefs";
        }

        new Comparison().compareFromSimulations("/Users/user/tetrad/aatoyexample" + type, simulations, algorithms, statistics, parameters);
    }
}


