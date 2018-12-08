package edu.cmu.tetrad.algcomparison.algorithm.oracle.pag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.CciTest;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.independence.KciMatlab;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.SachsUtils;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.List;

/**
 * FCI.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FCI",
        command = "fci",
        algoType = AlgType.allow_latent_common_causes
)
public class Fci_Sachs implements Algorithm, TakesInitialGraph, HasKnowledge, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();

    public Fci_Sachs() {
    }

    public Fci_Sachs(IndependenceWrapper test) {
        this.test = test;
    }

    public Fci_Sachs(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
    	if (parameters.getInt("bootstrapSampleSize") < 1) {
            if (algorithm != null) {
                initialGraph = algorithm.search(dataSet, parameters);
            }

            edu.cmu.tetrad.search.Fci search;

//            if (test.getClass() != KciMatlab.class && test.getClass() != CciTest.class) {
                search = new edu.cmu.tetrad.search.Fci(test.getTest(dataSet, parameters));
//            } else {
//                DataSet resampledDataSet = DataUtils.getResamplingDataset(((DataSet) dataSet), parameters.getInt("resampleSize"));
//                search = new edu.cmu.tetrad.search.Fci(test.getTest(resampledDataSet, parameters));
//            }

            search.setDepth(parameters.getInt("depth"));

            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            search.setKnowledge(knowledge);
            search.setMaxPathLength(parameters.getInt("maxPathLength"));
            search.setCompleteRuleSetUsed(parameters.getBoolean("completeRuleSetUsed"));
            search.setPossibleDsepSearchDone(parameters.getBoolean("possibleDsepDone"));
            search.setVerbose(parameters.getBoolean("verbose"));

//            if (initialGraph != null) {
//                search.setInitialGraph(initialGraph);
//            }

            Graph graph = search.search();

            return SU.pruneGraph(graph);
        } else {
            Fci_Sachs algorithm = new Fci_Sachs(test);
            //algorithm.setKnowledge(knowledge);
//          if (initialGraph != null) {
//      		algorithm.setInitialGraph(initialGraph);
//  		}
            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, algorithm, parameters.getInt("numberResampling"));
            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            search.setKnowledge(knowledge);
            search.setResampleSize(parameters.getInt("resampleSize"));
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
            return SU.pruneGraph(search.search());
        }
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return new DagToPag(new EdgeListGraph(graph)).convert();
    }

    public String getDescription() {
        return "FCI (Fast Causal Inference) using " + test.getDescription()
                + (algorithm != null ? " with initial graph from "
                        + algorithm.getDescription() : "");
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = test.getParameters();
        parameters.add("depth");
        parameters.add("maxPathLength");
        parameters.add("completeRuleSetUsed");
        parameters.add("possibleDsepDone");
        // Bootstrapping
        parameters.add("bootstrapSampleSize");
        parameters.add("bootstrapEnsemble");
        parameters.add("verbose");
        // Sachs
        parameters.add("forbidAmongInterventions");
        parameters.add("requiredEdgeKnowledge");
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
    public Graph getInitialGraph() {
        return initialGraph;
    }

    @Override
    public void setInitialGraph(Graph initialGraph) {
        this.initialGraph = initialGraph;
    }

    @Override
    public void setInitialGraph(Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    @Override
    public void setIndependenceWrapper(IndependenceWrapper test) {
        this.test = test;
    }

}
