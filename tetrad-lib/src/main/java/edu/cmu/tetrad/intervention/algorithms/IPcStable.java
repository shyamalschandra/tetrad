package edu.cmu.tetrad.intervention.algorithms;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.intervention.CleanInterventions;
import edu.cmu.tetrad.intervention.InterventionalKnowledge;
import edu.cmu.tetrad.search.PcAll;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.bootstrap.BootstrapEdgeEnsemble;
import edu.pitt.dbmi.algo.bootstrap.GeneralBootstrapTest;

import java.util.List;

/**
 * PC-Stable.
 *
 * @author jdramsey
 */
public class IPcStable implements Algorithm, TakesInitialGraph, HasKnowledge, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();

    public IPcStable() {
    }

    public IPcStable(IndependenceWrapper test) {
        this.test = test;
    }

    public IPcStable(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (!parameters.getBoolean("bootstrapping")) {
            if (algorithm != null) {
//            	initialGraph = algorithm.search(dataSet, parameters);
            }

            //REMOVE EXTRA OBSERVATIONS
            CleanInterventions ci = new CleanInterventions();
            dataSet = ci.removeExtra(dataSet);

            //REMOVE INTERVENTIONAL CONTEXT
            dataSet = ci.removeInterventionContext(dataSet);

            //REMOVE CONTEXT
            if(parameters.getBoolean("removeContext")) {
                dataSet = ci.removeContext(dataSet);
            }

            //KNOWLEDGE ADDED
            Graph nodes = GraphUtils.emptyGraph(0);
            for (Node node : dataSet.getVariables()) {
                nodes.addNode(node);
            }
            InterventionalKnowledge k = new InterventionalKnowledge(nodes, parameters.getInt("hand"));
            this.knowledge = k.getKnowledge();

            edu.cmu.tetrad.search.PcAll search = new edu.cmu.tetrad.search.PcAll(test.getTest(dataSet, parameters), initialGraph);
            search.setDepth(parameters.getInt("depth"));
            search.setKnowledge(knowledge);
            search.setFasRule(PcAll.FasRule.FAS_STABLE);
            search.setColliderDiscovery(edu.cmu.tetrad.search.PcAll.ColliderDiscovery.FAS_SEPSETS);
            search.setConflictRule(PcAll.ConflictRule.PRIORITY);
            search.setVerbose(parameters.getBoolean("verbose"));
            return search.search();
        }else{
            IPcStable pcStable = new IPcStable(test, algorithm);

            pcStable.setKnowledge(knowledge);
            if (initialGraph != null) {
                pcStable.setInitialGraph(initialGraph);
            }
            DataSet data = (DataSet) dataSet;
            GeneralBootstrapTest search = new GeneralBootstrapTest(data, pcStable,
                    parameters.getInt("bootstrapSampleSize"));

            BootstrapEdgeEnsemble edgeEnsemble = BootstrapEdgeEnsemble.Highest;
            switch (parameters.getInt("bootstrapEnsemble", 1)) {
                case 0:
                    edgeEnsemble = BootstrapEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = BootstrapEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = BootstrapEdgeEnsemble.Majority;
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
        return "IPC-Stable (\"Peter and Clark\" Stable), Priority Rule, using " + test.getDescription() + (algorithm != null ? " with initial graph from " +
                algorithm.getDescription() : "");
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = test.getParameters();
        parameters.add("depth");
        parameters.add("verbose");
        parameters.add("hand");
        parameters.add("removeContext");
        // Bootstrapping
        parameters.add("bootstrapping");
        parameters.add("bootstrapSampleSize");
        parameters.add("bootstrapEnsemble");
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
    public void setIndependenceWrapper(IndependenceWrapper independenceWrapper) {
        this.test = independenceWrapper;
    }

}
