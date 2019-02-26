package edu.cmu.tetrad.algcomparison.score;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndTestIndepRes;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.search.ScoredIndTest;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * Wrapper for Fisher Z test.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Score(
        name = "Independent Residuals Score",
        command = "indresscore",
        dataType = DataType.Continuous
)
public class IndTestIndepResScore implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;
    double alpha = 0.001;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;
        double alpha = parameters.getDouble("alpha");
        this.alpha = alpha;
        IndTestIndepRes test = new IndTestIndepRes((DataSet) dataSet, alpha);
        return new ScoredIndTest(test);
    }

    @Override
    public String getDescription() {
        return "Independent Residuals Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("alpha");
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return dataSet.getVariable(name);
    }
}
