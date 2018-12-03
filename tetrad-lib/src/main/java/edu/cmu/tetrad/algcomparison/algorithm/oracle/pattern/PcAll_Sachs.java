package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import edu.cmu.tetrad.algcomparison.TrueGraphSetter;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.*;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.SachsUtils;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.csb.KCI;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.util.List;

/**
 * CPC.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "PC All",
        command = "pc-all",
        algoType = AlgType.forbid_latent_common_causes
)
public class PcAll_Sachs implements Algorithm, TakesInitialGraph, HasKnowledge, TakesIndependenceWrapper, TrueGraphSetter {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();
    private Graph trueGraph = null;

    public PcAll_Sachs() {
    }

    public PcAll_Sachs(IndependenceWrapper test) {
        this.test = test;
    }

    public PcAll_Sachs(IndependenceWrapper test, Algorithm algorithm) {
        this.test = test;
        this.algorithm = algorithm;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (parameters.getInt("numberResampling") < 1) {
            edu.cmu.tetrad.search.PcAll.ColliderDiscovery colliderDiscovery;

            switch (parameters.getInt("colliderDiscoveryRule")) {
                case 1:
                    colliderDiscovery = edu.cmu.tetrad.search.PcAll.ColliderDiscovery.FAS_SEPSETS;
                    break;
                case 2:
                    colliderDiscovery = edu.cmu.tetrad.search.PcAll.ColliderDiscovery.CONSERVATIVE;
                    break;
                case 3:
                    colliderDiscovery = edu.cmu.tetrad.search.PcAll.ColliderDiscovery.MAX_P;
                    break;
                 default:
                    throw new IllegalArgumentException("Not a choice: " + parameters.getInt("colliderDiscoveryRule"));
            }

            edu.cmu.tetrad.search.PcAll.ConflictRule conflictRule;

            switch (parameters.getInt("conflictRule")) {
                case 1:
                    conflictRule = edu.cmu.tetrad.search.PcAll.ConflictRule.OVERWRITE;
                    break;
                case 2:
                    conflictRule = edu.cmu.tetrad.search.PcAll.ConflictRule.BIDIRECTED;
                    break;
                case 3:
                    conflictRule = edu.cmu.tetrad.search.PcAll.ConflictRule.PRIORITY;
                    break;
                default:
                    throw new IllegalArgumentException("Not a choice.");
            }

            edu.cmu.tetrad.search.PcAll search;

//            if (test.getClass() != KciMatlab.class && test.getClass() != CciTest.class) {
                search = new edu.cmu.tetrad.search.PcAll(test.getTest(dataSet, parameters), initialGraph);
//            } else {
//                DataSet resampledDataSet = DataUtils.getResamplingDataset(((DataSet) dataSet), parameters.getInt("resampleSize"));
//                search = new edu.cmu.tetrad.search.PcAll(test.getTest(resampledDataSet, parameters), initialGraph);
//            }

            search.setDepth(parameters.getInt("depth"));

            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge();

            search.setKnowledge(knowledge);

            if (parameters.getInt("fasType") == 4) {
                search.setFasType(edu.cmu.tetrad.search.PcAll.FasType.Naive);
            } else if (parameters.getInt("fasType") == 5) {
                search.setFasType(edu.cmu.tetrad.search.PcAll.FasType.LiWang);
            } else if (parameters.getBoolean("stableFASFDR")) {
                search.setFasType(edu.cmu.tetrad.search.PcAll.FasType.StableFDR);
            } else {
                if (parameters.getBoolean("stableFAS")) {
                    search.setFasType(edu.cmu.tetrad.search.PcAll.FasType.STABLE);
                } else {
                    search.setFasType(edu.cmu.tetrad.search.PcAll.FasType.REGULAR);
                }
            }

            if (parameters.getBoolean("concurrentFAS")) {
                search.setConcurrent(edu.cmu.tetrad.search.PcAll.Concurrent.YES);
            } else {
                search.setConcurrent(edu.cmu.tetrad.search.PcAll.Concurrent.NO);
            }

            search.setColliderDiscovery(colliderDiscovery);
            search.setConflictRule(conflictRule);
            search.setUseHeuristic(parameters.getBoolean("useMaxPOrientationHeuristic"));
            search.setMaxPathLength(parameters.getInt("maxPOrientationMaxPathLength"));
            search.setVerbose(parameters.getBoolean("verbose"));

            search.setTrueGraph(trueGraph);

            return SU.pruneGraph(search.search());
        } else {
            PcAll_Sachs pcAll = new PcAll_Sachs(test, algorithm);

            if (initialGraph != null) {
                pcAll.setInitialGraph(initialGraph);
            }

            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, pcAll, parameters.getInt("numberResampling"));
            search.setParallelMode(false);

            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge();
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

    public void setTrueGraph(Graph trueGraph) {
        this.trueGraph = trueGraph;
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return SearchGraphUtils.patternForDag(new EdgeListGraph(graph));
    }

    @Override
    public String getDescription() {
        return "PC using " + test.getDescription() + (algorithm != null ? " with initial graph from "
                + algorithm.getDescription() : "");
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = test.getParameters();
        parameters.add("StableFDR");
        parameters.add("stableFAS");
        parameters.add("fasType");
        parameters.add("concurrentFAS");
        parameters.add("colliderDiscoveryRule");
        parameters.add("conflictRule");
        parameters.add("depth");
        parameters.add("useMaxPOrientationHeuristic");
        parameters.add("maxPOrientationMaxPathLength");
        // Resampling
        parameters.add("numberResampling");
        parameters.add("resampleSize");
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