package edu.cmu.tetrad.algcomparison.algorithm.multi;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.List;

/**
 * Wraps the IMaGES algorithm for continuous variables.
 * </p>
 * Requires that the parameter 'randomSelectionSize' be set to indicate how many
 * datasets should be taken at a time (randomly). This cannot given multiple values.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FASK_C",
        command = "fask_c",
        algoType = AlgType.forbid_latent_common_causes
)
public class Fask_C implements Algorithm, HasKnowledge, UsesScoreWrapper {
    static final long serialVersionUID = 23L;
    private ScoreWrapper score;
    private IKnowledge knowledge = new Knowledge2();

    public Fask_C() {

    }

    public Fask_C(ScoreWrapper score) {
        this.score = score;
    }

    private Graph getGraph(edu.cmu.tetrad.search.Fask_B search) {
        return search.search();
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (parameters.getInt("bootstrapSampleSize") < 1) {
            edu.cmu.tetrad.search.Fask_C search = new edu.cmu.tetrad.search.Fask_C((DataSet) dataSet, score.getScore((DataSet) dataSet, parameters));

            search.setDepth(parameters.getInt("depth"));
            search.setSkewEdgeAlpha(parameters.getDouble("skewEdgeAlpha"));
            search.setTwoCycleAlpha(parameters.getDouble("twoCycleAlpha"));
            search.setMaxIterations(parameters.getInt("maxIterations"));
            search.setVerbose(parameters.getBoolean("verbose"));
            search.setUseSkewAdjacencies(parameters.getBoolean("useSkewAdjacencies"));
            search.setUseFasAdjacencies(parameters.getBoolean("useFasAdjacencies"));

            search.setKnowledge(knowledge);
            return search.search();
        } else {
            Fask_C fask = new Fask_C(score);
            fask.setKnowledge(knowledge);

            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, fask, parameters.getInt("numberResampling"));
            search.setKnowledge(knowledge);

            search.setResamplingWithReplacement(parameters.getBoolean("resamplingWithReplacement"));

            ResamplingEdgeEnsemble edgeEnsemble = ResamplingEdgeEnsemble.Highest;
            switch (parameters.getInt("resamplingEnsemble", 1)) {
                case 0:
                    edgeEnsemble = ResamplingEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = ResamplingEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = ResamplingEdgeEnsemble.Majority;
            }
            search.setEdgeEnsemble(edgeEnsemble);
            search.setParameters(parameters);
            search.setVerbose(parameters.getBoolean("verbose"));
            return search.search();
        }
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return new EdgeListGraph(graph);
    }

    @Override
    public String getDescription() {
        return "FASK-C using " + score.getDescription();
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = score.getParameters();
        parameters.add("depth");
        parameters.add("skewEdgeAlpha");
        parameters.add("twoCycleAlpha");
        parameters.add("maxIterations");

        parameters.add("useFasAdjacencies");
        parameters.add("useSkewAdjacencies");

        // Bootstrapping
        parameters.add("numberResampling");
        parameters.add("percentResampleSize");
        parameters.add("resamplingWithReplacement");
        parameters.add("resamplingEnsemble");
        parameters.add("verbose");

        return parameters;
    }

    @Override
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    @Override
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    @Override
    public void setScoreWrapper(ScoreWrapper score) {
        this.score = score;
    }
}