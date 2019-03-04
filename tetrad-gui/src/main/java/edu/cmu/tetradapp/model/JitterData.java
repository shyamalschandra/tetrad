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

package edu.cmu.tetradapp.model;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradSerializableUtils;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.HashMap;
import java.util.List;

import static javax.swing.UIManager.put;

/**
 * Applies a logarithmic transform
 *
 * @author Joseph Ramsey
 */
public class JitterData extends DataWrapper {
    static final long serialVersionUID = 23L;

    //=============================CONSTRUCTORS==============================//

    public JitterData(DataWrapper wrapper, Parameters params) {
        DataModelList inList = wrapper.getDataModelList();
        DataModelList outList = new DataModelList();

        RealDistribution dist = new UniformRealDistribution(-.001, 0.001);

        for (DataModel model : inList) {
            if (!(model instanceof DataSet)) {
                throw new IllegalArgumentException("Not a data set: " + model.getName());
            }

            DataSet dataSet = ((DataSet) model).copy();
//            double[][] data = dataSet.getDoubleData().transpose().toArray();

            for (int c = 0; c < dataSet.getNumColumns(); c++) {
                Node node = dataSet.getVariable(c);

//                HashMap<Double, Integer> repeats = new HashMap<>();
//
//                for (int i = 0; i < dataSet.getNumRows(); i++) {
//                    double d = dataSet.getDouble(i, c);
//
//                    if (!repeats.containsKey(d)) {
//                        repeats.put(d, 0);
//                    }
//
//                    repeats.put(d, repeats.get(d) + 1);
//                }
//
//                double mean = StatUtils.mean(data[c]);
//                double sd = StatUtils.sd(data[c]);
//
//                for (int i = 0; i < dataSet.getNumRows(); i++) {
//                    if (repeats.get(dataSet.getDouble(i, c)) > 20) {
//                        dataSet.setDouble(i, c, new NormalDistribution(mean, sd).sample());
//                    }
//                }
//
                if (node instanceof ContinuousVariable) {
                    for (int r = 0; r < dataSet.getNumRows(); r++) {
                        dataSet.setDouble(r, c, dataSet.getDouble(r, c) + dist.sample());
                    }
                }
            }

            outList.add(dataSet);
        }

        setDataModel(outList);
        setSourceGraph(wrapper.getSourceGraph());

        LogDataUtils.logDataModelList("Jittered Data.", getDataModelList());

    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     *
     * @see TetradSerializableUtils
     */
    public static PcRunner serializableInstance() {
        return PcRunner.serializableInstance();
    }

}



