package edu.cmu.tetrad.algcomparison.algorithm.oracle.pag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.independence.CciTest;
import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.independence.KciMatlab;
import edu.cmu.tetrad.algcomparison.score.CciScore;
import edu.cmu.tetrad.algcomparison.score.KciMatlabScore;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.SachsUtils;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.Experimental;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.search.GFci;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.dbmi.algo.resampling.GeneralResamplingTest;
import edu.pitt.dbmi.algo.resampling.ResamplingEdgeEnsemble;

import java.io.PrintStream;
import java.util.List;

/**
 * GFCI.
 *
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "GFCI",
        command = "gfci",
        algoType = AlgType.allow_latent_common_causes
)
@Experimental
public class Gfci_Sachs implements Algorithm, HasKnowledge, UsesScoreWrapper, TakesIndependenceWrapper {

    static final long serialVersionUID = 23L;
    private IndependenceWrapper test;
    private ScoreWrapper score;
    private IKnowledge knowledge = new Knowledge2();

    public Gfci_Sachs() {
    }

    public Gfci_Sachs(IndependenceWrapper test, ScoreWrapper score) {
        this.test = test;
        this.score = score;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        if (parameters.getInt("numberResampling") < 1) {

            GFci search;
//            if (score.getClass() != KciMatlabScore.class && score.getClass() != CciScore.class && test.getClass() != KciMatlab.class && test.getClass() != CciTest.class) {
            search = new GFci(test.getTest(dataSet, parameters), score.getScore(dataSet, parameters));
//            } else {
//                DataModel resampledDataSet = DataUtils.getResamplingDataset(((DataSet) dataSet), parameters.getInt("resampleSize"));
//                search = new GFci(test.getTest(dataSet, parameters), score.getScore(resampledDataSet, parameters));
//            }

            search.setMaxDegree(parameters.getInt("maxDegree"));

            SachsUtils SU = new SachsUtils();
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            search.setKnowledge(knowledge);
            search.setVerbose(parameters.getBoolean("verbose"));
            search.setFaithfulnessAssumed(parameters.getBoolean("faithfulnessAssumed"));
            search.setMaxPathLength(parameters.getInt("maxPathLength"));
            search.setCompleteRuleSetUsed(parameters.getBoolean("completeRuleSetUsed"));

            Object obj = parameters.get("printStream");

            if (obj instanceof PrintStream) {
                search.setOut((PrintStream) obj);
            }

            return SU.pruneGraph(search.search());
        } else {
            Gfci_Sachs algorithm = new Gfci_Sachs(test, score);

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
        return new DagToPag(graph).convert();
    }

    @Override
    public String getDescription() {
        return "GFCI (Greedy Fast Causal Inference) using " + test.getDescription()
                + " and " + score.getDescription();
    }

    @Override
    public DataType getDataType() {
        return test.getDataType();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = test.getParameters();
        List<String> parameters1 = score.getParameters();

        for (String param : parameters1) {
            if (!parameters.contains(param)) {
                parameters.add(param);
            }
        }

        parameters.add("faithfulnessAssumed");
        parameters.add("maxDegree");
//        parameters.add("printStream");
        parameters.add("maxPathLength");
        parameters.add("completeRuleSetUsed");
        // Resampling
        parameters.add("numberResampling");
        parameters.add("resampleSize");
        parameters.add("resamplingWithReplacement");
        parameters.add("resamplingEnsemble");
        parameters.add("verbose");
        // Sachs
        parameters.add("forbidAmongInterventions");
        parameters.add("requiredEdgeKnowledge");
@Experimental
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            knowledge = SU.getKnowledge(
                    parameters.getBoolean("forbidAmongInterventions",true),
                    parameters.getBoolean("requiredEdgeKnowledge", false));
            search.setResampleSize(parameters.getInt("resampleSize"));
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
    public void setScoreWrapper(ScoreWrapper score) {
        this.score = score;
    }

    @Override
    public void setIndependenceWrapper(IndependenceWrapper test) {
        this.test = test;
    }

}
