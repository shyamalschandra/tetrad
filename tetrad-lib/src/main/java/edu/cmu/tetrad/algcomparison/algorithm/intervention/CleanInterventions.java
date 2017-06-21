package edu.cmu.tetrad.algcomparison.algorithm.intervention;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

/**
 * Created by bandrews on 6/21/17.
 */
public class CleanInterventions {

    public CleanInterventions() {
    }

    public Graph removeEdges(Graph I) {
        Graph G = GraphUtils.emptyGraph(0);
        for (Edge edge : G.getEdges()) {
            if (edge.getNode1().getName().startsWith("I") || edge.getNode2().getName().startsWith("I")) {
                G.removeEdge(edge);
            }
        }
        return G;
    }

    public Graph removeNodes(Graph I) {
        Graph G = GraphUtils.emptyGraph(0);
        for (Node node : G.getNodes()) {
            if (node.getName().startsWith("I")) {
                G.removeNode(node);
            }
        }
        return G;
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

}
