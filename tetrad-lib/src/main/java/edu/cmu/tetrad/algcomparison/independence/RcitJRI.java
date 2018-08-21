package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.search.IndTestKciMatlab;
import edu.cmu.tetrad.search.IndTestRcitJRI;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * Wrapper for Fisher Z test.
 *
 * @author jdramsey
 */
@TestOfIndependence(
        name = "RCIT (Java-R Interface)",
        command = "rcit-jri",
        dataType = DataType.Continuous
)
public class RcitJRI implements IndependenceWrapper {

    static final long serialVersionUID = 23L;


    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        final IndTestRcitJRI test = new IndTestRcitJRI(DataUtils.getContinuousDataSet(dataSet),
                parameters.getDouble("alpha"));
        test.setFastFDR(parameters.getBoolean("fastFDR"));
        return test;
    }

    @Override
    public String getDescription() {
        return "RCIT (Java-R Interface)";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();
        params.add("alpha");
        params.add("fastFDR");
        return params;
    }
}
