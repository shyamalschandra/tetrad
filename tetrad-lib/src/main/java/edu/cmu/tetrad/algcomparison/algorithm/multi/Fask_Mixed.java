package edu.cmu.tetrad.algcomparison.algorithm.multi;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.SachsUtils;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Experimental;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.GeneralMixedSearch;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.List;

/**
 * @author bandrews
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FASK_Mixed",
        command = "fask_mixed",
        algoType = AlgType.forbid_latent_common_causes
)
@Experimental
public class Fask_Mixed implements Algorithm, HasKnowledge, TakesIndependenceWrapper {
    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private IKnowledge knowledge = new Knowledge2();

    public Fask_Mixed() {

    }

    public Fask_Mixed(IndependenceWrapper test) {
        this.test = test;
    }

    private Graph getGraph(edu.cmu.tetrad.search.Fask_B search) {
        return search.search();
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {

        GeneralMixedSearch GMS = new GeneralMixedSearch((DataSet) dataSet);
        GMS.transform();
        DataSet tranformedDataset = GMS.getTrnasformedDataset();
        this.knowledge = GMS.getKnowledge();

        if (parameters.getInt("bootstrapSampleSize") < 1) {
            edu.cmu.tetrad.search.Fask_B search = new edu.cmu.tetrad.search.Fask_B(tranformedDataset, test.getTest(tranformedDataset, parameters));
            search.setDepth(parameters.getInt("depth"));
            search.setSkewEdgeAlpha(parameters.getDouble("skewEdgeAlpha"));
            search.setTwoCycleAlpha(parameters.getDouble("twoCycleAlpha"));
            search.setVerbose(parameters.getBoolean("verbose"));
            search.setUseSkewAdjacencies(parameters.getBoolean("useSkewAdjacencies"));
            search.setUseFasAdjacencies(parameters.getBoolean("useFasAdjacencies"));
            search.setUseMask(parameters.getBoolean("useMask"));
            search.setMaskThreshold(parameters.getDouble("maskThreshold"));
            search.setKnowledge(this.knowledge);

            return GMS.removeExogenous(getGraph(search));

        } else {
            Fask_Mixed fask = new Fask_Mixed(test);
            fask.setKnowledge(knowledge);
            GeneralResamplingTest search = new GeneralResamplingTest(tranformedDataset, fask, parameters.getInt("numberResampling"));
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

            return GMS.removeExogenous(search.search());
        }
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return new EdgeListGraph(graph);
    }

    @Override
    public String getDescription() {
        return "FASK-Mixed using " + test.getDescription();
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = test.getParameters();
        parameters.add("depth");
        parameters.add("skewEdgeAlpha");
        parameters.add("twoCycleAlpha");

        parameters.add("useFasAdjacencies");
        parameters.add("useSkewAdjacencies");
        parameters.add("useMask");
        parameters.add("maskThreshold");

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
    public void setIndependenceWrapper(IndependenceWrapper independenceWrapper) {
        this.test = independenceWrapper;
    }
}