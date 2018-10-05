package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.annotation.Gaussian;
import edu.cmu.tetrad.annotation.Linear;
import edu.cmu.tetrad.annotation.TestOfIndependence;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndTestFisherZSkew;
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
        name = "Fisher Z Skew Test",
        command = "fisher-z-skew",
        dataType = {DataType.Continuous, DataType.Covariance}
)
@Gaussian
@Linear
public class FisherZSkew implements IndependenceWrapper {

    static final long serialVersionUID = 23L;
    private double alpha = 0.001;

    @Override
    public IndependenceTest getTest(DataModel dataSet, Parameters parameters) {
        double alpha = parameters.getDouble("alpha");
        this.alpha = alpha;

        if (dataSet instanceof ICovarianceMatrix) {
            final IndTestFisherZSkew test = new IndTestFisherZSkew((ICovarianceMatrix) dataSet, alpha);
            test.setFastFDR(parameters.getBoolean("fastFDR"));
            return test;
        } else if (dataSet instanceof DataSet) {
            final IndTestFisherZSkew test = new IndTestFisherZSkew((DataSet) dataSet, alpha);
            test.setFastFDR(parameters.getBoolean("fastFDR"));
            return test;
        }

        throw new IllegalArgumentException("Expecting either a data set or a covariance matrix.");
    }

    @Override
    public String getDescription() {
        return "Fisher Z Skew test, alpha = " + alpha;
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
