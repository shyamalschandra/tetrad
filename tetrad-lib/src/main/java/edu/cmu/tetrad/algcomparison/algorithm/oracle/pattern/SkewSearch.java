package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.ArrayList;
import java.util.List;

@edu.cmu.tetrad.annotation.Algorithm(
        name = "SkewSearch",
        command = "skew-search",
        algoType = AlgType.forbid_latent_common_causes
)
public class SkewSearch implements Algorithm, HasKnowledge {
    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();

    public SkewSearch(IndependenceWrapper test) {
        this.test = test;
    }

    public SkewSearch(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    public SkewSearch() {
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (parameters.getInt("bootstrapSampleSize") < 1) {
            edu.cmu.tetrad.search.SkewSearch search = new edu.cmu.tetrad.search.SkewSearch((DataSet) dataSet);
            search.setDepth(parameters.getInt("depth"));
            search.setAlpha(parameters.getDouble("alpha"));
            search.setTwoCycleAlpha(parameters.getDouble("twoCycleAlpha"));
            search.setDelta(parameters.getDouble("faskDelta"));
            search.setKnowledge(knowledge);
            search.setVerbose(parameters.getBoolean("verbose"));
            return search.search();
        } else {
            SkewSearch algorithm = new SkewSearch(test);

            DataSet data = (DataSet) dataSet;

            GeneralResamplingTest search = new GeneralResamplingTest(data, algorithm, parameters.getInt("numberResampling"));
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
        return SearchGraphUtils.patternForDag(new EdgeListGraph(graph));
    }

    @Override
    public String getDescription() {
        return "SkewSearch";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
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
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("alpha");
        parameters.add("faskDelta2");
        parameters.add("verbose");
        // Bootstrapping
        parameters.add("bootstrapSampleSize");
        parameters.add("bootstrapEnsemble");
        return parameters;
    }
}