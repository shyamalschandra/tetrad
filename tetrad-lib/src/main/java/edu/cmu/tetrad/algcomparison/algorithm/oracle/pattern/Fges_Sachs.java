package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.CciTest;
import edu.cmu.tetrad.algcomparison.independence.SemBicTest;
import edu.cmu.tetrad.algcomparison.score.CciScore;
import edu.cmu.tetrad.algcomparison.score.KciMatlabScore;
import edu.cmu.tetrad.algcomparison.score.KciScore;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.SachsUtils;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Experimental;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.SearchGraphUtils;
import edu.cmu.tetrad.search.SemBicScore;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * FGES (the heuristic version).
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "FGES",
        command = "fges",
        algoType = AlgType.forbid_latent_common_causes
)
@Experimental
public class Fges_Sachs implements Algorithm, TakesInitialGraph, HasKnowledge, UsesScoreWrapper {

    static final long serialVersionUID = 23L;

    private boolean compareToTrue = false;
    private ScoreWrapper score;
    private Algorithm algorithm = null;
    private Graph initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();

    public Fges_Sachs() {

    }

    public Fges_Sachs(ScoreWrapper score) {
        this.score = score;
        this.compareToTrue = false;
    }

    public Fges_Sachs(ScoreWrapper score, boolean compareToTrueGraph) {
        this.score = score;
        this.compareToTrue = compareToTrueGraph;
    }

    public Fges_Sachs(ScoreWrapper score, Algorithm algorithm) {
        this.score = score;
        this.algorithm = algorithm;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
    	if (parameters.getInt("numberResampling") < 1) {
            if (algorithm != null) {
//                initialGraph = algorithm.search(dataSet, parameters);
            }

            edu.cmu.tetrad.search.Fges search;
            search = new edu.cmu.tetrad.search.Fges(score.getScore(dataSet, parameters));


            search.setFaithfulnessAssumed(parameters.getBoolean("faithfulnessAssumed"));

            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            search.setKnowledge(knowledge);

            search.setVerbose(parameters.getBoolean("verbose"));
            search.setMaxDegree(parameters.getInt("maxDegree"));
            search.setSymmetricFirstStep(parameters.getBoolean("symmetricFirstStep"));

            Object obj = parameters.get("printStream");
            if (obj instanceof PrintStream) {
                search.setOut((PrintStream) obj);
            }

            if (initialGraph != null) {
                search.setInitialGraph(initialGraph);
            }

            return SU.pruneGraph(search.search());

        } else {
            Fges_Sachs fges = new Fges_Sachs(score, algorithm);

            //fges.setKnowledge(knowledge);
            DataSet data = (DataSet) dataSet;
            GeneralResamplingTest search = new GeneralResamplingTest(data, fges, parameters.getInt("numberResampling"));

            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            search.setKnowledge(knowledge);
            
            search.setPercentResampleSize(parameters.getInt("percentResampleSize"));
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
        if (compareToTrue) {
            return new EdgeListGraph(graph);
        } else {
            return SearchGraphUtils.patternForDag(new EdgeListGraph(graph));
        }
    }

    @Override
    public String getDescription() {
        return "FGES (Fast Greedy Equivalence Search) using " + score.getDescription();
    }

    @Override
    public DataType getDataType() {
        return score.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = score.getParameters();
        parameters.add("faithfulnessAssumed");
        parameters.add("symmetricFirstStep");
        parameters.add("maxDegree");
        parameters.add("verbose");
        // Resampling
        parameters.add("numberResampling");
        parameters.add("resampleSize");
        parameters.add("resamplingWithReplacement");
        parameters.add("resamplingEnsemble");
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

    public void setCompareToTrue(boolean compareToTrue) {
        this.compareToTrue = compareToTrue;
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
    public void setScoreWrapper(ScoreWrapper score) {
        this.score = score;
    }

}
