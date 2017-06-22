package edu.cmu.tetrad.algcomparison.algorithm.intervention;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.util.Parameters;

import java.util.List;

/**
 * RFCI.
 *
 * @author jdramsey
 */
public class Rfci_I implements Algorithm, HasKnowledge {
    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private IKnowledge knowledge = new Knowledge2();

    public Rfci_I(IndependenceWrapper test) {
        this.test = test;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {

        //REMOVE EXTRA OBSERVATIONS

        CleanInterventions ci = new CleanInterventions();
        dataSet = ci.removeExtra(dataSet);

        //REMOVE EXTRA OBSERVATIONS

        // KNOWLEDGE ADDED

        Graph nodes = GraphUtils.emptyGraph(0);
        for (Node node : dataSet.getVariables()) {
            nodes.addNode(node);
        }
        InterventionalKnowledge k = new InterventionalKnowledge(nodes);
        this.knowledge = k.getKnowledge();

        // KNOWLEDGE ADDED

        edu.cmu.tetrad.search.Rfci search = new edu.cmu.tetrad.search.Rfci(test.getTest(dataSet, parameters));
        search.setKnowledge(knowledge);
        search.setDepth(parameters.getInt("depth"));
        search.setMaxPathLength(parameters.getInt("maxPathLength"));
        search.setCompleteRuleSetUsed(parameters.getBoolean("completeRuleSetUsed"));
        return search.search();
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return new DagToPag(new EdgeListGraph(graph)).convert();
    }

    public String getDescription() {
        return "RFCI Interventions (Really Fast Causal Inference) using " + test.getDescription();
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
}
