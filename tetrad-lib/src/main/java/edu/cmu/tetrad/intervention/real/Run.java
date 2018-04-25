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

package edu.cmu.tetrad.intervention.real;

import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Gfci;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.Rfci;
import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.*;
import edu.cmu.tetrad.algcomparison.independence.ChiSquare;
import edu.cmu.tetrad.algcomparison.independence.ConditionalGaussianLRT;
import edu.cmu.tetrad.algcomparison.score.BdeuScore;
import edu.cmu.tetrad.algcomparison.score.ConditionalGaussianBicScore;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * @author bandrews
 */

public class Run {
    public static void main(String... args) {
        Load load = new Load();
        String path = "tetrad-lib/src/main/java/edu/cmu/tetrad/Bryan/";
        List<DataSet> dataSets = load.load_data(path + "data", 5);
        edu.cmu.tetrad.intervention.real.Combine combine = new Combine(dataSets);
        DataSet dataSet = combine.combine_datasets();
        Knowledge knowledge = new Knowledge(dataSet, dataSets.get(0));
        Knowledge2 prior_knowledge = knowledge.get_knowledge();
        Parameters parameters = new Parameters();
        parameters.set("alpha", 0.01);
        parameters.set("structurePrior", 1);
//        Pc search = new Pc(new ChiSquare());
//        PcMax search = new PcMax(new ChiSquare(), false);
//        PcStable search = new PcStable(new ChiSquare());
//        CpcStable search = new CpcStable(new ChiSquare());
//        Rfci search = new Rfci(new ChiSquare());
//        Gfci search = new Gfci(new ConditionalGaussianLRT(), new ConditionalGaussianBicScore());
        Fges search = new Fges(new ConditionalGaussianBicScore());

        search.setKnowledge(prior_knowledge);
        Graph result = search.search(dataSet, parameters);
        try {
            File file;

            File dir = new File(path, "out");
            dir.mkdirs();
            file = new File(dir, "graph.txt");

            PrintStream out = new PrintStream(file);
            System.out.println("Printing graph to " + file.getAbsolutePath());
            out.println(result);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}




