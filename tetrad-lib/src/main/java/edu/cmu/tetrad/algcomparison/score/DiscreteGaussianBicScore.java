package edu.cmu.tetrad.algcomparison.score;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.ConditionalGaussianScore;
import edu.cmu.tetrad.search.DiscreteGaussianScore;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.util.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * Wrapper for the Discrete Gaussian score
 *
 * @author Bryan Andrews
 */
@edu.cmu.tetrad.annotation.Score(
        name = "Discrete Gaussian BIC Score",
        command = "disc-gauss-bic",
        dataType = DataType.Mixed
)
public class DiscreteGaussianBicScore implements ScoreWrapper {

    static final long serialVersionUID = 23L;
    private DataModel dataSet;

    @Override
    public Score getScore(DataModel dataSet, Parameters parameters) {
        this.dataSet = dataSet;
        final DiscreteGaussianScore discreteGaussianScore
                = new DiscreteGaussianScore(DataUtils.getMixedDataSet(dataSet), parameters.getDouble("structurePrior"));

        return discreteGaussianScore;
    }

    @Override
    public String getDescription() {
        return "Discrete Gaussian BIC Score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();

        parameters.add("structurePrior");
        return parameters;
    }

    @Override
    public Node getVariable(String name) {
        return  dataSet.getVariable(name);
    }
}
