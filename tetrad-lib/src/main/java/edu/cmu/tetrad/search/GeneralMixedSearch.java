package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * ALL (or no) VARIABLES MUST BE IN A TIER
 */

public class GeneralMixedSearch {

    // The original dataset.
    private DataSet dataSet;

    // The transformed dataset.
    private DataSet transformedDataSet;

    // Knowledge of forbidden and required edges.
    private IKnowledge knowledge;

    // The list of original variables.
    private List<Node> variables;

    // The list of transformed variables.
    private List<Node> transformedVariables;

    // The list of discrete variables.
    private List<Node> exogenousVariables;

    // The list of continuous variables.
    private List<Node> domainVariables;

    // Constructor
    public GeneralMixedSearch(DataSet dataSet) {
        this.dataSet = dataSet;
        this.variables = dataSet.getVariables();
        this.knowledge = new Knowledge2();
    }

    // Transform
    public void transform() {
        exogenousKnowledge();
        transformDiscreteToContinuous();
    }

    // Get tranformed dataset
    public DataSet getTrnasformedDataset() {
        return this.transformedDataSet;
    }

    // Get tranformed variables
    public List<Node> getTrnasformedVariables() {
        return this.transformedVariables;
    }

    // Get exogenous variables
    public List<Node> getExogenousVariables() {
        return this.exogenousVariables;
    }

    // Get domain variables
    public List<Node> getDomainVariables() {
        return this.domainVariables;
    }

    // Get knowledge
    public IKnowledge getKnowledge() {
        return this.knowledge;
    }

    // Set knowledge
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    // Create exogenous knowledge
    private void exogenousKnowledge() {

        // get the variable names
        List<String> discreteNames = new ArrayList<>();
        List<String> continuousNames = new ArrayList<>();
        for (Node variable : this.variables) {
            if (variable instanceof DiscreteVariable) {
                discreteNames.add(variable.getName());
            } else if (variable instanceof ContinuousVariable) {
                continuousNames.add(variable.getName());
            }
        }

        // find the min and max tiers
        int minTier = -1;
        int maxTier = 0;
        for (Node variable : this.variables) {
            int tier = this.knowledge.isInWhichTier(variable);
            if (minTier == -1 && tier > -1) {
                minTier = tier;
            } else if (tier < minTier) {
                minTier = tier;
            }
            if (tier > maxTier) {
                maxTier = tier;
            }
        }

        // find exogenous variables are remove them from current knowledge
        List<String> tier0 = new ArrayList<>();
        if (minTier > -1) {
            List<String> tier = this.knowledge.getTier(minTier);
            boolean flag = false;
            if (this.knowledge.isTierForbiddenWithin(minTier)) {
                flag = true;
            } else {
                for (String node : discreteNames) {
                    if (tier.contains(node)) {
                        flag = true;
                        break;
                    }
                }
            }
            if (flag) {
                for (String node1 : tier) {
                    tier0.add(node1);
                    this.knowledge.removeFromTiers(node1);
                    for (String node2 : continuousNames) {
                        this.knowledge.removeRequired(node2, node1);
                    }
                }
            }
        } else if (!discreteNames.isEmpty()) {
            for (String node1 : discreteNames) {
                tier0.add(node1);
                this.knowledge.removeFromTiers(node1);
                for (String node2 : continuousNames) {
                    this.knowledge.removeRequired(node2, node1);
                }
            }
            for (String node : continuousNames) {
                this.knowledge.addToTier(0, node);
            }
        }

        // shift tiers
        if (!tier0.isEmpty()) {
            for (int i = maxTier + 1; i > 0; i -= 1) {
                boolean forbidden = this.knowledge.isTierForbiddenWithin(i - 1);
                this.knowledge.setTier(i, this.knowledge.getTier(i - 1));
                this.knowledge.setTier(i - 1, new ArrayList<>());
                this.knowledge.setTierForbiddenWithin(i, forbidden);
            }
            for (String node : tier0) {
                this.knowledge.addToTier(0, node);
            }
            this.knowledge.setTierForbiddenWithin(0, true);
        }

        // assign variables
        this.exogenousVariables = new ArrayList<>();
        this.domainVariables = new ArrayList<>();
        for (Node variable : this.variables) {
            if (tier0.contains(variable.getName())) {
                this.exogenousVariables.add(variable);
            } else {
                this.domainVariables.add(variable);
            }
        }
    }

    // Transform a mixed dataset into a continuous dataset.
    private void transformDiscreteToContinuous(){

        int N = this.dataSet.getNumRows();
        Map map = new HashMap<>();

        List<Node> A = new ArrayList<>();
        List<double[]> B = new ArrayList<>();

        int i = 0;
        int i_ = 0;
        while (i_ < this.variables.size()) {

            Node v = this.variables.get(i_);

            List<String> required = new ArrayList<>();
            for (Node domain : this.domainVariables) {
                if (this.knowledge.isRequired(v.getName(), domain.getName())) {
                    required.add(domain.getName());
                }
            }

            if (v instanceof DiscreteVariable) {

                Map<String, Integer> keys = new HashMap<>();
                for (int j = 0; j < N; j++) {
                    String key = v.getName().concat("_");
                    key = key.concat(Integer.toString(this.dataSet.getInt(j, i_)));
                    if (!keys.containsKey(key)) {
                        keys.put(key, i);
                        Node v_ = new ContinuousVariable(key);
                        A.add(v_);
                        B.add(new double[N]);
                        i++;
                        this.exogenousVariables.add(v_);
                        for (String other : required) {
                            this.knowledge.setRequired(key, other);
                        }
                    }
                    B.get(keys.get(key))[j] = 1;
                }

                /*
                 * Remove a degenerate dimension.
                 */
                i--;
                keys.remove(A.get(i).getName());
                A.remove(i);
                B.remove(i);

                List<Integer> index = new ArrayList<>();
                index.addAll(keys.values());
                map.put(i_, index);
                i_++;

            } else {

                A.add(v);
                double[] b = new double[N];
                for (int j = 0; j < N; j++) {
                    b[j] = this.dataSet.getDouble(j,i_);
                }
                B.add(org.apache.commons.math3.stat.StatUtils.normalize(b));

                List<Integer> index = new ArrayList<>();
                index.add(i);
                map.put(i_, index);
                i++;
                i_++;

            }
        }

        double[][] B_ = new double[N][B.size()];
        for (int j = 0; j < B.size(); j++) {
            for (int k = 0; k < N; k++) {
                B_[k][j] = B.get(j)[k];
            }
        }

        this.transformedVariables = A;
        List<Node> toRemove = new ArrayList<>();
        for (Node node : this.exogenousVariables) {
            if (!A.contains(node)) {
                toRemove.add(node);
                this.knowledge.removeVariable(node.getName());
            } else if(!this.variables.contains(node)) {
                this.knowledge.addToTier(0, node.getName());
            }
        }
        this.exogenousVariables.removeAll(toRemove);

        RealMatrix D = new BlockRealMatrix(B_);
        DataSet dataSet = new BoxDataSet(new DoubleDataBox(D.getData()), A);
        this.transformedDataSet = DataUtils.center(dataSet);

    }

    // Remove the exogenous variables from a graph.
    public Graph removeExogenous(Graph graph) {
        for (Node node : this.exogenousVariables) {
            graph.removeNode(node);
        }
        return graph;
    }

}
