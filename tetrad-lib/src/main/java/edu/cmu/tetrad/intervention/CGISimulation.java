package edu.cmu.tetrad.intervention;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.simulation.Simulation;
import edu.cmu.tetrad.bayes.BayesIm;
import edu.cmu.tetrad.bayes.BayesPm;
import edu.cmu.tetrad.bayes.MlBayesIm;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.ParamType;
import edu.cmu.tetrad.sem.Parameter;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.RandomUtil;
import org.apache.commons.lang3.RandomUtils;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.*;

/**
 * A simulation method based on the conditional Gaussian assumption which allows for interventions.
 * This method is an adaptation of jramsey's conditional Gaussian simulation.
 *
 * @author bandrews
 */
public class CGISimulation implements Simulation {
    static final long serialVersionUID = 23L;
    private RandomGraph randomGraph;
    private List<DataSet> dataSets = new ArrayList<>();
    private List<Graph> graphs = new ArrayList<>();
    private List<Graph> iGraphs = new ArrayList<>();
    private DataType dataType;
    private List<Node> shuffledOrder;
    private double varLow = 1;
    private double varHigh = 3;
    private double coefLow = 0.05;
    private double coefHigh = 1.5;
    private boolean coefSymmetric = true;
    private double meanLow = -1;
    private double meanHigh = 1;

    private int sampleSize;
    private int interventionSize;
    private int numInterventions;
    private double percentIDiscrete;
    private int minICategories;
    private int maxICategories;
    private int minEffected;
    private int maxEffected;
    private double minPotency;
    private double maxPotency;
    private double percentCInfluence;
    private double discreteCInfluence;
    private double continuousCInfluence;

    public CGISimulation(RandomGraph graph) {
        this.randomGraph = graph;
    }

    @Override
    public void createData(Parameters parameters) {

        // This method is largly code from Joe's CG Simulation.

        this.sampleSize = parameters.getInt("sampleSize");
        this.interventionSize = parameters.getInt("interventionSize");
        this.numInterventions = parameters.getInt("numInterventions");
        this.percentIDiscrete = parameters.getDouble("percentIDiscrete");
        this.minICategories = parameters.getInt("minICategories");
        this.maxICategories = parameters.getInt("maxICategories");
        this.minEffected = parameters.getInt("minEffected");
        this.maxEffected = parameters.getInt("maxEffected");
        this.minPotency = parameters.getDouble("minPotency");
        this.maxPotency = parameters.getDouble("maxPotency");
        this.percentCInfluence = parameters.getDouble("percentCInfluence");
        this.discreteCInfluence = parameters.getDouble("discreteCInfluence");
        this.continuousCInfluence = parameters.getDouble("continuousCInfluence");

        setVarLow(parameters.getDouble("varLow"));
        setVarHigh(parameters.getDouble("varHigh"));
        setCoefLow(parameters.getDouble("coefLow"));
        setCoefHigh(parameters.getDouble("coefHigh"));
        setCoefSymmetric(parameters.getBoolean("coefSymmetric"));
        setMeanLow(parameters.getDouble("meanLow"));
        setMeanHigh(parameters.getDouble("meanHigh"));

        this.sampleSize += this.numInterventions * this.interventionSize;

        double percentDiscrete = parameters.getDouble("percentDiscrete");

        boolean discrete = parameters.getString("dataType").equals("discrete");
        boolean continuous = parameters.getString("dataType").equals("continuous");

        if (discrete && percentDiscrete != 100.0) {
            throw new IllegalArgumentException("To simulate discrete data, 'percentDiscrete' must be set to 0.0.");
        } else if (continuous && percentDiscrete != 0.0) {
            throw new IllegalArgumentException("To simulate continuoue data, 'percentDiscrete' must be set to 100.0.");
        }

        if (discrete) this.dataType = DataType.Discrete;
        if (continuous) this.dataType = DataType.Continuous;

        this.shuffledOrder = null;

        Graph graph = randomGraph.createGraph(parameters);
        List<Intervention> interventions = createInterventions(graph);
        Graph iGraph = addInterventions(graph, interventions);

        dataSets = new ArrayList<>();
        graphs = new ArrayList<>();
        iGraphs = new ArrayList<>();

        System.out.println("");
        for (String p : parameters.getParametersNames()) {
            System.out.println(p + ": " + parameters.get(p));
        }
        System.out.println("");

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            System.out.println("Simulating dataset #" + (i + 1));

            if (parameters.getBoolean("differentGraphs") && i > 0) {
                graph = randomGraph.createGraph(parameters);
                interventions = createInterventions(graph);
                iGraph = addInterventions(graph, interventions);
            }

            graphs.add(graph);
            iGraphs.add(iGraph);

            DataSet dataSet = simulate(graph, interventions, parameters);

            // Add contextual influence

            Node context = new DiscreteVariable("C", this.numInterventions + 1);
            dataSet.addVariable(context);
            iGraph.addNode(context);
            for(int j = 0; j < this.sampleSize; j++) {
                dataSet.setInt(j, dataSet.getColumn(context), 0);
            }
            for(int j = 0; j < this.numInterventions; j++) {
                for(int k = 0; k < this.interventionSize; k++) {
                    int tmp_ind = this.sampleSize + j * this.interventionSize + k;
                    dataSet.setInt(tmp_ind, dataSet.getColumn(context), j+1);
                }
            }

            for(int j = 0; j < graph.getNumNodes(); j++) {
                for(int k = 0; k < this.numInterventions; k++) {
                    if (RandomUtils.nextDouble(0,1) < this.percentCInfluence) {
                        iGraph.addDirectedEdge(context, iGraph.getNode(dataSet.getVariable(j).getName()));
                        if (dataSet.getVariable(j) instanceof ContinuousVariable) {
                            double mod = 2*RandomUtils.nextDouble(0, this.continuousCInfluence) - this.continuousCInfluence;
                            for(int l = 0; l < this.interventionSize; l++) {
                                int tmp_ind = this.sampleSize + k * this.interventionSize + l;
                                dataSet.setDouble(tmp_ind, j, dataSet.getDouble(tmp_ind, j) + mod);
                            }
                        }
                    } else if(dataSet.getVariable(j) instanceof DiscreteVariable) {
                        double[] mod = new double[((DiscreteVariable) dataSet.getVariable(j)).getNumCategories()];
                        double norm = 0;
                        GammaDistribution g = new GammaDistribution(1, 1);
                        for(int l = 0; l < mod.length;  l++) {
                            mod[l] = g.sample();
                            norm += mod[l];
                        }
                        for(int l = 0; l < mod.length;  l++) {
                            mod[l] /= norm;
                        }
                        for(int l = 0; l < this.interventionSize; l++) {
                            if(RandomUtils.nextDouble(0,1) < this.discreteCInfluence) {
                                int tmp_ind = this.sampleSize + k * this.interventionSize + l;
                                double prob = RandomUtils.nextDouble(0,1);
                                for(int m = 0; m < mod.length; m++) {
                                    prob -= mod[m];
                                    if(prob <= 0) {
                                        dataSet.setInt(tmp_ind, j, m);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            List<Node> measured = new ArrayList<>();
            for (
                    Node node : iGraph.getNodes()) {
                if (node.getNodeType() == NodeType.MEASURED) {
                    measured.add(dataSet.getVariable(node.getName()));
                }
            }
            dataSet = dataSet.subsetColumns(measured);

            dataSet.setName("" + (i + 1));
            dataSets.add(dataSet);
        }
    }

    @Override
    public Graph getTrueGraph(int index) {
        return iGraphs.get(index);
    }

    @Override
    public DataModel getDataModel(int index) {
        return dataSets.get(index);
    }

    @Override
    public String getDescription() {
        return "Conditional Gaussian Intervention simulation using " + randomGraph.getDescription();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = randomGraph.getParameters();
        parameters.add("minCategories");
        parameters.add("maxCategories");
        parameters.add("percentDiscrete");
        parameters.add("numRuns");
        parameters.add("differentGraphs");
        parameters.add("sampleSize");
        parameters.add("varLow");
        parameters.add("varHigh");
        parameters.add("coefLow");
        parameters.add("coefHigh");
        parameters.add("coefSymmetric");
        parameters.add("meanLow");
        parameters.add("meanHigh");

        parameters.add("interventionSize");
        parameters.add("numInterventions");
        parameters.add("percentIDiscrete");
        parameters.add("minICategories");
        parameters.add("maxICategories");
        parameters.add("minEffected");
        parameters.add("maxEffected");
        parameters.add("minPotency");
        parameters.add("maxPotency");
        parameters.add("percentCInfluence");
        parameters.add("discreteCInfluence");
        parameters.add("continuousCInfluence");

        return parameters;
    }

    @Override
    public int getNumDataModels() {
        return dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return dataType;
    }

    private DataSet simulate(Graph G, List<Intervention> I, Parameters parameters) {

        // Simulate the data with interventions.

        HashMap<String, Integer> nd = new HashMap<>();

        List<Node> nodes = G.getNodes();

        Collections.shuffle(nodes);

        if (this.shuffledOrder == null) {
            List<Node> shuffledNodes = new ArrayList<>(nodes);
            Collections.shuffle(shuffledNodes);
            this.shuffledOrder = shuffledNodes;
        }

        // I assume that all interventions affect discrete variables.  Here the domain variables are assigned as
        // either discrete or continuous variables.  If they are assigned to be continuous, then I correct for that.

        for (int i = 0; i < nodes.size(); i++) {
            if (i < nodes.size() * parameters.getDouble("percentDiscrete") * 0.01) {
                final int minNumCategories = parameters.getInt("minCategories");
                final int maxNumCategories = parameters.getInt("maxCategories");
                final int value = pickNumCategories(minNumCategories, maxNumCategories);
                nd.put(shuffledOrder.get(i).getName(), value);
            } else {
                nd.put(shuffledOrder.get(i).getName(), 0);
                for (Intervention intervention : I) {
                    if (intervention.getEffected().contains(shuffledOrder.get(i).getName())) {
                        double epsilon = RandomUtils.nextDouble(this.varLow, this.varHigh);
                        if (intervention.isDiscrete()) {
                            intervention.changeEffect(shuffledOrder.get(i).getName(), this.meanLow, this.meanHigh, epsilon);
                        } else {
                            intervention.changeEffect(shuffledOrder.get(i).getName(), this.coefLow, this.coefHigh, epsilon);
                        }
                    }
                }
            }
        }

        // Generate the graphs.

        G = makeMixedGraph(G, nd);
        Graph IG = addInterventions(G, I);
        nodes = IG.getNodes();

        // Initialize the dataset.

        DataSet mixedData = new BoxDataSet(new MixedDataBox(nodes, this.sampleSize + this.numInterventions * this.interventionSize), nodes);

        // Set the intervened values in the dataset.

        for (Intervention intervention : I) {
            int col1 = mixedData.getColumn(IG.getNode("ID" + intervention.getName()));
            int col2 = mixedData.getColumn(IG.getNode("IV" + intervention.getName()));
            for (int i = 0; i < this.sampleSize + this.numInterventions * this.interventionSize; i++) {
                mixedData.setInt(i, col1, intervention.getDomain(i));
                if (!intervention.isDiscrete()) {
                    mixedData.setDouble(i, col2, intervention.getValue(i));
                }
            }
        }

        // Make the discrete versions of the continuous variables (this is for mixed variables simulation).

        List<Node> X = new ArrayList<>();
        List<Node> A = new ArrayList<>();

        for (Node node : G.getNodes()) {
            if (node instanceof ContinuousVariable) {
                X.add(node);
            } else {
                A.add(node);
            }
        }

        Graph AG = G.subgraph(A);
        Graph XG = G.subgraph(X);

        Map<ContinuousVariable, DiscreteVariable> erstatzNodes = new HashMap<>();
        Map<String, ContinuousVariable> erstatzNodesReverse = new HashMap<>();

        for (Node y : A) {
            for (Node x : G.getParents(y)) {
                if (x instanceof ContinuousVariable) {
                    DiscreteVariable ersatz = erstatzNodes.get(x);

                    if (ersatz == null) {
                        ersatz = new DiscreteVariable("Ersatz_" + x.getName(), RandomUtil.getInstance().nextInt(3)+2);
                        erstatzNodes.put((ContinuousVariable) x, ersatz);
                        erstatzNodesReverse.put(ersatz.getName(), (ContinuousVariable) x);
                        AG.addNode(ersatz);
                    }

                    AG.addDirectedEdge(ersatz, y);
                }
            }
        }

        // Parameterize the graph and get causal ordering.

        BayesPm bayesPm = new BayesPm(AG);
        BayesIm bayesIm = new MlBayesIm(bayesPm, MlBayesIm.RANDOM);
//        BayesIm bayesIm = new MlBayesIm_o(bayesPm, MlBayesIm_o.RANDOM);

        SemPm semPm = new SemPm(XG);

        Map<Combination, Double> paramValues = new HashMap<>();

        List<Node> tierOrdering = G.getCausalOrdering();

        int[] tiers = new int[tierOrdering.size()];

        for (int t = 0; t < tierOrdering.size(); t++) {
            tiers[t] = nodes.indexOf(tierOrdering.get(t));
        }

        Map<Integer, double[]> breakpointsMap = new HashMap<>();

        // In causal ordering simulate the data instances.  Also all this code is from Joe's CG Simulation.

        for (int mixedIndex : tiers) {
            for (int i = 0; i < this.sampleSize + this.numInterventions * this.interventionSize; i++) {
                if (nodes.get(mixedIndex) instanceof DiscreteVariable) {

                    int bayesIndex = bayesIm.getNodeIndex(nodes.get(mixedIndex));

                    int[] bayesParents = bayesIm.getParents(bayesIndex);
                    int[] parentValues = new int[bayesParents.length];

                    for (int k = 0; k < parentValues.length; k++) {
                        int bayesParentColumn = bayesParents[k];

                        Node bayesParent = bayesIm.getVariables().get(bayesParentColumn);
                        DiscreteVariable _parent = (DiscreteVariable) bayesParent;
                        int value;

                        ContinuousVariable orig = erstatzNodesReverse.get(_parent.getName());

                        if (orig != null) {
                            int mixedParentColumn = mixedData.getColumn(orig);
                            double d = mixedData.getDouble(i, mixedParentColumn);
                            double[] breakpoints = breakpointsMap.get(mixedParentColumn);

                            if (breakpoints == null) {
                                breakpoints = getBreakpoints(mixedData, _parent, mixedParentColumn);
                                breakpointsMap.put(mixedParentColumn, breakpoints);
                            }

                            value = breakpoints.length;

                            for (int j = 0; j < breakpoints.length; j++) {
                                if (d < breakpoints[j]) {
                                    value = j;
                                    break;
                                }
                            }
                        } else {
                            int mixedColumn = mixedData.getColumn(bayesParent);
                            value = mixedData.getInt(i, mixedColumn);
                        }

                        parentValues[k] = value;
                    }

                    int rowIndex = bayesIm.getRowIndex(bayesIndex, parentValues);
                    double sum = 0.0;

                    double r = RandomUtil.getInstance().nextDouble();
                    mixedData.setInt(i, mixedIndex, 0);

                    for (int k = 0; k < bayesIm.getNumColumns(bayesIndex); k++) {
                        double probability = bayesIm.getProbability(bayesIndex, rowIndex, k);
                        sum += probability;

                        if (sum >= r) {

                            // START ADDED FOR INTERVENTIONS

                            // If the variable and instance are both affected by an intervention and the
                            // intervention is potent enough, assign the intervention value.  This is where
                            // non-stationarity is introduced into the dataset.

                            String node = mixedData.getVariable(mixedIndex).getName();
                            for (Intervention intervention : I) {
                                if (intervention.getEffected().contains(node) && intervention.getDomain(i) > 0) {
                                    double prob = RandomUtils.nextDouble(0,1);
                                    if (prob <= intervention.getPotency()) {
                                        k = (int) intervention.getInterventionValue(node, i) - 1;
                                    }
                                    break;
                                }
                            }

                            // END ADDED FOR INTERVENTIONS

                            mixedData.setInt(i, mixedIndex, k);

                            break;
                        }
                    }
                } else {
                    Node y = nodes.get(mixedIndex);

                    Set<DiscreteVariable> discreteParents = new HashSet<>();
                    Set<ContinuousVariable> continuousParents = new HashSet<>();

                    for (Node node : G.getParents(y)) {
                        if (node instanceof DiscreteVariable) {
                            discreteParents.add((DiscreteVariable) node);
                        } else {
                            continuousParents.add((ContinuousVariable) node);
                        }
                    }

                    Parameter varParam = semPm.getParameter(y, y);
                    Parameter muParam = semPm.getMeanParameter(y);

                    Combination varComb = new Combination(varParam);
                    Combination muComb = new Combination(muParam);

                    for (DiscreteVariable v : discreteParents) {
                        varComb.addParamValue(v, mixedData.getInt(i, mixedData.getColumn(v)));
                        muComb.addParamValue(v, mixedData.getInt(i, mixedData.getColumn(v)));
                    }

                    double value = RandomUtil.getInstance().nextNormal(0, getParamValue(varComb, paramValues));

                    for (Node x : continuousParents) {
                        Parameter coefParam = semPm.getParameter(x, y);
                        Combination coefComb = new Combination(coefParam);

                        for (DiscreteVariable v : discreteParents) {
                            coefComb.addParamValue(v, mixedData.getInt(i, mixedData.getColumn(v)));
                        }

                        int parent = nodes.indexOf(x);
                        double parentValue = mixedData.getDouble(i, parent);
                        double parentCoef = getParamValue(coefComb, paramValues);
                        value += parentValue * parentCoef;
                    }

                    value += getParamValue(muComb, paramValues);

                    // START ADDED FOR INTERVENTIONS

                    // If the variable and instance are both affected by an intervention and the
                    // intervention is potent enough, assign the intervention value.  This is where
                    // non-stationarity is introduced into the dataset.

                    String node = mixedData.getVariable(mixedIndex).getName();
                    for (Intervention intervention : I) {
                        if (intervention.getEffected().contains(node) && intervention.getDomain(i) > 0) {
                            double prob = RandomUtils.nextDouble(0,1);
                            if (prob <= intervention.getPotency()) {
                                value = intervention.getInterventionValue(node, i);
                            }
                            break;
                        }
                    }

                    // END ADDED FOR INTERVENTIONS

                    mixedData.setDouble(i, mixedIndex, value);

                }
            }
        }

        return mixedData;
    }

    private double[] getBreakpoints(DataSet mixedData, DiscreteVariable _parent, int mixedParentColumn) {
        double[] data = new double[mixedData.getNumRows()];

        for (int r = 0; r < mixedData.getNumRows(); r++) {
            data[r] = mixedData.getDouble(r, mixedParentColumn);
        }

        return Discretizer.getEqualFrequencyBreakPoints(data, _parent.getNumCategories());
    }

    private Double getParamValue(Combination values, Map<Combination, Double> map) {
        Double d = map.get(values);

        if (d == null) {
            Parameter parameter = values.getParameter();

            if (parameter.getType() == ParamType.VAR) {
                d = RandomUtil.getInstance().nextUniform(varLow, varHigh);
                map.put(values, d);
            } else if (parameter.getType() == ParamType.COEF) {
                double min = coefLow;
                double max = coefHigh;
                double value = RandomUtil.getInstance().nextUniform(min, max);
                d = RandomUtil.getInstance().nextUniform(0, 1) < 0.5 && coefSymmetric ? -value : value;
                map.put(values, d);
            } else if (parameter.getType() == ParamType.MEAN) {
                d = RandomUtil.getInstance().nextUniform(meanLow, meanHigh);
                map.put(values, d);
            }
        }

        return d;
    }

    public void setVarLow(double varLow) {
        this.varLow = varLow;
    }

    public void setVarHigh(double varHigh) {
        this.varHigh = varHigh;
    }

    public void setCoefLow(double coefLow) {
        this.coefLow = coefLow;
    }

    public void setCoefHigh(double coefHigh) {
        this.coefHigh = coefHigh;
    }

    public void setCoefSymmetric(boolean coefSymmetric) {
        this.coefSymmetric = coefSymmetric;
    }

    public void setMeanLow(double meanLow) {
        this.meanLow = meanLow;
    }

    public void setMeanHigh(double meanHigh) {
        this.meanHigh = meanHigh;
    }

    private class Combination {
        private Parameter parameter;
        private Set<VariableValues> paramValues;

        public Combination(Parameter parameter) {
            this.parameter = parameter;
            this.paramValues = new HashSet<>();
        }

        public void addParamValue(DiscreteVariable variable, int value) {
            this.paramValues.add(new VariableValues(variable, value));
        }

        public int hashCode() {
            return parameter.hashCode() + paramValues.hashCode();
        }

        public boolean equals(Object o) {
            if (o == this) return true;
            if (!(o instanceof Combination)) return false;
            Combination v = (Combination) o;
            return v.parameter == this.parameter && v.paramValues.equals(this.paramValues);
        }

        public Parameter getParameter() {
            return parameter;
        }
    }

    private class VariableValues {
        private DiscreteVariable variable;
        private int value;

        public VariableValues(DiscreteVariable variable, int value) {
            this.variable = variable;
            this.value = value;
        }

        public DiscreteVariable getVariable() {
            return variable;
        }

        public int getValue() {
            return value;
        }

        public int hashCode() {
            return variable.hashCode() + value;
        }

        public boolean equals(Object o) {
            if (o == this) return true;
            if (!(o instanceof VariableValues)) return false;
            VariableValues v = (VariableValues) o;
            return v.variable.equals(this.variable) && v.value == this.value;
        }
    }

    private static Graph makeMixedGraph(Graph g, Map<String, Integer> m) {
        List<Node> nodes = g.getNodes();
        for (int i = 0; i < nodes.size(); i++) {
            Node n = nodes.get(i);
            int nL = m.get(n.getName());
            if (nL > 0) {
                Node nNew = new DiscreteVariable(n.getName(), nL);
                nodes.set(i, nNew);
            } else {
                Node nNew = new ContinuousVariable(n.getName());
                nodes.set(i, nNew);
            }
        }

        Graph outG = new EdgeListGraph(nodes);

        for (Edge e : g.getEdges()) {
            Node n1 = e.getNode1();
            Node n2 = e.getNode2();
            Edge eNew = new Edge(outG.getNode(n1.getName()), outG.getNode(n2.getName()), e.getEndpoint1(), e.getEndpoint2());
            outG.addEdge(eNew);
        }

        return outG;
    }

    private int pickNumCategories(int min, int max) {
        return RandomUtils.nextInt(min, max + 1);
    }

    private List<Intervention> createInterventions(Graph G) {

        // Given a graph, this methods generates a list of interventions on the variables in the input graph.

        List<Intervention> interventions = new ArrayList<>();
        List<Node> nodes = G.getNodes();
        Collections.shuffle(nodes);

        int index = 0;

        for (int i = 0; i < this.numInterventions; i++) {
            Intervention I;
            int numEffected = RandomUtils.nextInt(this.minEffected, this.maxEffected + 1);
            double potency = RandomUtils.nextDouble(this.minPotency, this.maxPotency);

            if (i < this.numInterventions * 0.01 * this.percentIDiscrete) {
                int numICategories = RandomUtils.nextInt(this.minICategories, this.maxICategories + 1);
                I = new Intervention(true, potency, this.sampleSize + this.interventionSize * this.numInterventions, this.interventionSize, this.sampleSize + i * this.interventionSize, numICategories, 0, 0);
                interventions.add(I);
            } else {
                double mean = RandomUtils.nextDouble(this.meanLow, this.meanHigh);
                double var = RandomUtils.nextDouble(this.varLow, this.varHigh);
                I = new Intervention(false, potency, this.sampleSize + this.interventionSize * this.numInterventions, this.interventionSize, this.sampleSize + i * this.interventionSize, 1, mean, var);
                interventions.add(I);
            }

            int j = 0;
            while (j < numEffected) {
                if (nodes.get(index).getNodeType() == NodeType.MEASURED) {
                    String effected = nodes.get(index).getName();
                    I.addEffect(effected);
                    j++;
                }
                index++;
            }
        }
        return interventions;
    }


    private Graph addInterventions(Graph G, List<Intervention> I) {

        // Given a list of interventions and a graph, this method generates a graph based on the input graph which includes the interventions in the list of interventions.

        GraphUtils GU = new GraphUtils();
        Graph IG = GU.emptyGraph(0);
        for (Node node : G.getNodes()) {
            IG.addNode(node);
        }
        for (Edge edge : G.getEdges()) {
            IG.addEdge(edge);
        }
        for (Intervention intervention : I) {
            if (!intervention.isDiscrete()) {
                ContinuousVariable iv = new ContinuousVariable("IV" + intervention.getName());
                IG.addNode(iv);
                for (String node : intervention.getEffected()) {
                    Edge edge = new Edge(iv, G.getNode(node), Endpoint.TAIL, Endpoint.ARROW);
                    IG.addEdge(edge);
                }
            }
            DiscreteVariable id = new DiscreteVariable("ID" + intervention.getName());
            IG.addNode(id);
            for (String node : intervention.getEffected()) {
                Edge edge = new Edge(id, G.getNode(node), Endpoint.TAIL, Endpoint.ARROW);
                IG.addEdge(edge);
            }
        }
        return IG;
    }

}
