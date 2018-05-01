package edu.cmu.tetrad.intervention.real;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.util.dist.Discrete;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by bandrews on 9/14/17.
 */

public class Combine {

    private List<DataSet> dataSets;
    private List<Integer> interventions = new ArrayList<>();
    private List<Integer> domains = new ArrayList<>();
    private List<Node> variables = new ArrayList<>();
    private List<Integer> dataset_index = new ArrayList<>();

    public Combine(List<DataSet> dataSets) {
        this.dataSets = dataSets;
        List<Node> vars = dataSets.get(0).getVariables();
        for (Node var : vars) {
            if (var.getName().startsWith("I_")) {
                this.interventions.add(dataSets.get(0).getColumn(var));
            } else {
                this.domains.add(dataSets.get(0).getColumn(var));
            }
        }
    }

    public DataSet combine_datasets() {

        for (Integer intervention : this.interventions) {
            Node new_context = this.dataSets.get(1).getVariable(intervention);
            boolean non_deterministic = true;
            for (Node variable : this.variables) {
                if (check_deterministic(new_context, variable, 1, 1)) {
                    variable.setName(variable.getName() + ":IC" + new_context.getName().substring(1));
                    non_deterministic = false;
                }
                if (!non_deterministic) {
                    break;
                }
            }
            if (non_deterministic) {
                new_context.setName("IC" + new_context.getName().substring(1));
                this.variables.add(new_context);
                this.dataset_index.add(1);
            }
        }

        for (Integer intervention : this.interventions) {
            Node new_value = this.dataSets.get(0).getVariable(intervention);
            boolean non_deterministic = true;
            for (Node variable : this.variables) {
                int j = 0;
                if (variable.getName().startsWith("IC")) {
                    j = 1;
                }
                if (check_deterministic(new_value, variable, 0, j)) {
                    variable.setName(variable.getName() + ":IV" + new_value.getName().substring(1));
                    non_deterministic = false;
                }
                if (!non_deterministic) {
                    break;
                }
            }
            if (non_deterministic) {
                new_value.setName("IV" + new_value.getName().substring(1));
                this.variables.add(new_value);
                this.dataset_index.add(0);
            }
        }

        for (Integer domain : this.domains) {
            Node new_variable = this.dataSets.get(0).getVariable(domain);
            boolean non_deterministic = true;
            for (Node variable : this.variables) {
                int j = 0;
                if (variable.getName().startsWith("IC")) {
                    j = 1;
                }
                if (check_deterministic(new_variable, variable, 0, j)) {
                    variable.setName(variable.getName() + ":" + new_variable.getName());
                    non_deterministic = false;
                }
                if (!non_deterministic) {
                    break;
                }
            }
            if (non_deterministic) {
                new_variable.setName(new_variable.getName());
                this.variables.add(new_variable);
                this.dataset_index.add(0);
            }
        }

        DataSet dataSet;
        if (this.dataSets.get(0).isContinuous() && this.dataSets.get(1).isContinuous()) {
            dataSet = new BoxDataSet(new DoubleDataBox(this.dataSets.get(0).getNumRows(), this.variables.size()), this.variables);
        } else if (this.dataSets.get(0).isDiscrete() && this.dataSets.get(1).isDiscrete()) {
            dataSet = new BoxDataSet(new IntDataBox(this.dataSets.get(0).getNumRows(), this.variables.size()), this.variables);
        } else {
            dataSet = new BoxDataSet(new MixedDataBox(this.variables, this.dataSets.get(0).getNumRows()), this.variables);
        }
        for (int i = 0; i < this.dataSets.get(0).getNumRows(); i++) {
            for (int j = 0; j < this.variables.size(); j++) {
                if (this.variables.get(j) instanceof DiscreteVariable) {
                    dataSet.setInt(i, j, this.dataSets.get(this.dataset_index.get(j)).getInt(i, this.dataSets.get(this.dataset_index.get(j)).getColumn(this.variables.get(j))));
                } else {
                    dataSet.setDouble(i, j, this.dataSets.get(this.dataset_index.get(j)).getDouble(i, this.dataSets.get(this.dataset_index.get(j)).getColumn(this.variables.get(j))));
                }
            }
        }
        return dataSet;
    }

    private boolean check_deterministic(Node a, Node b, int i, int j) {

        boolean deterministic = true;

        int ac = this.dataSets.get(i).getColumn(a);
        int bc = this.dataSets.get(j).getColumn(b);

        HashMap<Double, Double> mapping_forward = new HashMap<>();
        HashMap<Double, Double> mapping_reverse = new HashMap<>();

        if (a instanceof DiscreteVariable && b instanceof DiscreteVariable) {
            for (int k = 0; k < this.dataSets.get(0).getNumRows(); k++) {
                double v1 = this.dataSets.get(i).getInt(k,ac);
                double v2 = this.dataSets.get(j).getInt(k,bc);
                if (mapping_forward.containsKey(v1)) {
                    if (mapping_forward.get(v1) != v2) {
                        deterministic = false;
                        break;
                    }
                } else {
                    mapping_forward.put(v1, v2);
                }
                if (mapping_reverse.containsKey(v2)) {
                    if (mapping_reverse.get(v2) != v1) {
                        deterministic = false;
                        break;
                    }
                } else {
                    mapping_reverse.put(v2, v1);
                }
            }

        } else {
            deterministic = false;
        }

        return deterministic;
    }
}