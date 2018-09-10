package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.graph.SingleGraph;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.SemGraph;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.sem.TemplateExpander;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.RandomUtil;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jdramsey
 */
public class GeneralSemSimulationExample4 implements Simulation {
    static final long serialVersionUID = 23L;
    private RandomGraph randomGraph;
    private GeneralizedSemPm pm;
    private GeneralizedSemIm im;
    private List<DataSet> dataSets = new ArrayList<>();
    private List<DataSet> dataWithLatents = new ArrayList<>();
    private List<Graph> graphs = new ArrayList<>();
    private List<GeneralizedSemIm> ims = new ArrayList<>();

    private String[] nonlinearFunctions = {"sin(x)", "cos(x)", "(x)^2", "tanh(x)"};
    private String[] functions = {"sin(x)", "cos(x)", "(x)^2", "tanh(x)", "x"};
    private String variablesString = "TSUM(NEW(B) * $) + ERROR";
    private String errorString = "U(0, 1)^2";
    private String parametersString = "U(.4,.9)";

    public GeneralSemSimulationExample4(RandomGraph graph) {
        this.randomGraph = graph;
    }

    public GeneralSemSimulationExample4(GeneralizedSemPm pm) {
        SemGraph graph = pm.getGraph();
        graph.setShowErrorTerms(false);
        this.randomGraph = new SingleGraph(graph);
        this.pm = pm;
    }

    public GeneralSemSimulationExample4(GeneralizedSemIm im) {
        SemGraph graph = im.getSemPm().getGraph();
        graph.setShowErrorTerms(false);
        this.randomGraph = new SingleGraph(graph);
        this.im = im;
        this.ims = new ArrayList<>();
        ims.add(im);
        this.pm = im.getGeneralizedSemPm();
    }

    @Override
    public void createData(Parameters parameters) {
        errorString = parameters.getString("generalSemErrorTemplate");
        parametersString = parameters.getString("generalSemParameterTemplate");
        variablesString = parameters.getString("generalSemFunctionTemplateMeasured");

        Graph graph = randomGraph.createGraph(parameters);

        dataSets = new ArrayList<>();
        graphs = new ArrayList<>();
        ims = new ArrayList<>();

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            System.out.println("Simulating dataset #" + (i + 1));

            if (parameters.getBoolean("differentGraphs") && i > 0) {
                graph = randomGraph.createGraph(parameters);
            }

            graphs.add(graph);

            DataSet dataSet = simulate(graph, parameters);

            if (parameters.getBoolean("standardize")) {
                dataSet = DataUtils.standardizeData(dataSet);
            }

            double variance = parameters.getDouble("measurementVariance");

            if (variance > 0) {
                for (int k = 0; k < dataSet.getNumRows(); k++) {
                    for (int j = 0; j < dataSet.getNumColumns(); j++) {
                        double d = dataSet.getDouble(k, j);
                        double norm = RandomUtil.getInstance().nextNormal(0, Math.sqrt(variance));
                        dataSet.setDouble(k, j, d + norm);
                    }
                }
            }

            if (parameters.getBoolean("randomizeColumns")) {
                dataSet = DataUtils.reorderColumns(dataSet);
            }

            dataSet.setName("" + (i + 1));
            dataSets.add(DataUtils.restrictToMeasured(dataSet));
            dataWithLatents.add(dataSet);
        }
    }

    private synchronized DataSet simulate(Graph graph, Parameters parameters) {
        pm = getPm(graph, parameters);
        im = new GeneralizedSemIm(pm);
        ims.add(im);
        return im.simulateData(parameters.getInt("sampleSize"), true);

//        if (im == null) {
//            if (pm == null) {
//                pm = getPm(graph, parameters);// new GeneralizedSemPm(graph);
//                im = new GeneralizedSemIm(pm);
//                ims.add(im);
//                return im.simulateData(parameters.getInt("sampleSize"), true);
//            } else {
//                im = new GeneralizedSemIm(pm);
//                ims.add(im);
//                return im.simulateData(parameters.getInt("sampleSize"), true);
//            }
//        } else {
//            ims.add(im);
//            return im.simulateData(parameters.getInt("sampleSize"), true);
//        }
    }

    @Override
    public Graph getTrueGraph(int index) {
        return graphs.get(index);
    }

    @Override
    public int getNumDataModels() {
        return dataSets.size();
    }

    @Override
    public DataModel getDataModel(int index) {
        return dataSets.get(index);
    }

    @Override
    public DataModel getDataModelWithLatents(int index) {
        return null;
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    public String getDescription() {
        return "Nonlinear, non-Gaussian SEM simulation using " + randomGraph.getDescription();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();

        if (!(randomGraph instanceof SingleGraph)) {
            parameters.addAll(randomGraph.getParameters());
        }

        if (pm == null) {
            parameters.addAll(GeneralizedSemPm.getParameterNames());
        }

        parameters.add("numRuns");
        parameters.add("differentGraphs");
        parameters.add("sampleSize");

        return parameters;
    }

    private GeneralizedSemPm getPm(Graph graph, Parameters parameters) {
        GeneralizedSemPm pm = new GeneralizedSemPm(graph);

        List<Node> variablesNodes = pm.getVariableNodes();
        List<Node> errorNodes = pm.getErrorNodes();

        try {

            Map<Node, String> temp = new HashMap<>();

            for (Node x : variablesNodes) {
                List<Node> parents = graph.getParents(x);

                String formula = "";

                for (int i = 0; i < parents.size(); i++) {
                    Node p = parents.get(i);
                    formula += wrapRandom2(p.getName());
                    if (i < parents.size() - 1) formula += " + ";
                }

                temp.put(x, formula);

            }

            for (Node x : variablesNodes) {
                String formula = temp.get(x);
                if (!formula.isEmpty()) formula += " + ";
                formula += "U(-.5, .5)";
                formula = wrapRandom2(formula);
                pm.setNodeExpression(x, formula);
            }

            for (Node node : errorNodes) {
                String _template = TemplateExpander.getInstance().expandTemplate(
                        errorString, pm, node);
                pm.setNodeExpression(node, "U(-.5, .5)");
            }

            for (String parameter : pm.getParameters()) {
                pm.setParameterExpression(parameter, "U(.2, .7)");
            }

        } catch (ParseException e) {
            e.printStackTrace();
        }

        System.out.println(pm);

        return pm;
    }

    private String wrapRandom(String template) {
        String function = nonlinearFunctions[RandomUtil.getInstance().nextInt(nonlinearFunctions.length)];
        return function.replaceAll("x", template);
    }

    private String wrapRandom2(String template) {
        String function = functions[RandomUtil.getInstance().nextInt(functions.length)];
        return function.replaceAll("x", template);
    }

    public List<GeneralizedSemIm> getIms() {
        return ims;
    }
}
