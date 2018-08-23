package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.search.IndTestFcitJRI;
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
        name = "FCIT (Python through Java-R Interface)",
        command = "fcit-jri",
        dataType = DataType.Continuous
)
public class FcitJRI implements IndependenceWrapper {

    static final long serialVersionUID = 23L;


    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        final IndTestFcitJRI test = new IndTestFcitJRI(DataUtils.getContinuousDataSet(dataSet),
                parameters.getDouble("alpha"));
        test.setFastFDR(parameters.getBoolean("fastFDR"));
        return test;
    }

    @Override
    public String getDescription() {
        return "FCIT (Python through Java-R Interface)";
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
