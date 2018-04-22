package edu.cmu.tetrad.intervention;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by bandrews on 6/21/17.
 */
public class CleanInterventions {

    public CleanInterventions() {
    }

    public Graph removeEdges(Graph I) {
        Graph G = clone(I);
        for (Edge edge : G.getEdges()) {
            if (edge.getNode1().getName().startsWith("I") || edge.getNode2().getName().startsWith("I")) {
                G.removeEdge(edge);
            }
            if (edge.getNode1().getName().startsWith("C") || edge.getNode2().getName().startsWith("C")) {
                G.removeEdge(edge);
            }
        }
        return G;
    }

    public Graph removeNodes(Graph I) {
        Graph G = clone(I);
        for (Node node : G.getNodes()) {
            if (node.getName().startsWith("I")) {
                G.removeNode(node);
            }
            if (node.getName().startsWith("C")) {
                G.removeNode(node);
            }
        }
        return G;
    }

    private Graph clone(Graph G1) {
        Graph G2 = GraphUtils.emptyGraph(0);
        for (Node node : G1.getNodes()) {
            G2.addNode(node);
        }
        for (Edge edge : G1.getEdges()) {
            G2.addEdge(edge);
        }
        return G2;
    }

    public DataModel removeVars(DataModel I) {
        DataSet D = (DataSet) I.copy();
        for (Node col : D.getVariables()) {
            if (col.getName().startsWith("I")) {
                D.removeColumn(col);
            }
        }
        return D;
    }

    public DataModel removeInterventionContext(DataModel I) {
        DataSet D = (DataSet) I.copy();
        for (Node col : D.getVariables()) {
            if (col.getName().startsWith("IV")) {
                D.removeColumn(D.getVariable(col.getName().replaceFirst("V","D")));
            }
        }
        return D;
    }

    public DataModel removeInterventions(DataModel I) {
        DataSet D = (DataSet) I.copy();
        for (Node col : D.getVariables()) {
            if (col.getName().startsWith("I")) {
                D.removeColumn(D.getVariable(col.getName()));
                D.removeColumn(D.getVariable(col.getName()));
            }
        }
        return D;
    }

    public DataModel removeContext(DataModel I) {
        DataSet D = (DataSet) I.copy();
        for (Node col : D.getVariables()) {
            if (col.getName().startsWith("C")) {
                D.removeColumn(col);
            }
        }
        return D;
    }

    public DataModel removeRows(DataModel I) {
        DataSet D = (DataSet) I.copy();
        List<Integer> selectedRows = new ArrayList<>();
        for (int i = 0; i < D.getNumRows(); i++) {
            boolean rowI = false;
            for (Node node : D.getVariables()) {
                if (node.getName().startsWith("ID") && D.getInt(i,D.getColumn(node)) > 0) {
                    rowI = true;
                    break;
                }
            }
            if (!rowI) {
                selectedRows.add(i);
            }
        }

        DataSet newD = new BoxDataSet(new MixedDataBox(D.getVariables(), selectedRows.size()), D.getVariables());

        for (int row = 0; row < selectedRows.size(); row++) {
            for (int col = 0; col < D.getNumColumns(); col++) {
                if (D.getVariable(col) instanceof DiscreteVariable) {
                    newD.setInt(row, col, D.getInt(selectedRows.get(row), col));
                } else {
                    newD.setDouble(row, col, D.getDouble(selectedRows.get(row), col));
                }
            }
        }

        return newD;
    }

    public DataModel removeExtra(DataModel I) {
        DataSet D = (DataSet) I.copy();
        List<Integer> observedRows = new ArrayList<>();
        List<Integer> intervenedRows = new ArrayList<>();
        for (int i = 0; i < D.getNumRows(); i++) {
            boolean rowI = false;
            for (Node node : D.getVariables()) {
                if (node.getName().startsWith("ID") && D.getInt(i,D.getColumn(node)) > 0) {
                    intervenedRows.add(i);
                    rowI = true;
                    break;
                }
            }
            if (!rowI) {
                observedRows.add(i);
            }
        }

        DataSet newD = new BoxDataSet(new MixedDataBox(D.getVariables(), observedRows.size()), D.getVariables());

        int numObv = observedRows.size() - intervenedRows.size();

        for (int row = 0; row < numObv; row++) {
            for (int col = 0; col < D.getNumColumns(); col++) {
                if (D.getVariable(col) instanceof DiscreteVariable) {
                    newD.setInt(row, col, D.getInt(observedRows.get(row), col));
                } else {
                    newD.setDouble(row, col, D.getDouble(observedRows.get(row), col));
                }
            }
        }

        for (int row = 0; row < intervenedRows.size(); row++) {
            for (int col = 0; col < D.getNumColumns(); col++) {
                if (D.getVariable(col) instanceof DiscreteVariable) {
                    newD.setInt(numObv + row, col, D.getInt(intervenedRows.get(row), col));
                } else {
                    newD.setDouble(numObv + row, col, D.getDouble(intervenedRows.get(row), col));
                }
            }
        }

        return newD;
    }

}
