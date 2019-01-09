package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.search.IndTestRcitJRI;
import edu.cmu.tetrad.search.IndTestRcotJRI;
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
        name = "RCoT (Java-R Interface)",
        command = "rcot-jri",
        dataType = DataType.Continuous
)
public class RcotJRI implements IndependenceWrapper {

    static final long serialVersionUID = 23L;


    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        final IndTestRcotJRI test = new IndTestRcotJRI(DataUtils.getContinuousDataSet(dataSet),
                parameters.getDouble("alpha"));
        return test;
    }

    @Override
    public String getDescription() {
        return "RCoT (Java-R Interface)";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();
        params.add("alpha");
        return params;
    }
}
