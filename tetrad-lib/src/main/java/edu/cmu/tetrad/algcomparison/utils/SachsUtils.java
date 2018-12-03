package edu.cmu.tetrad.algcomparison.utils;

import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.util.ArrayList;
import java.util.List;

public class SachsUtils {

    private List<String> interventions;

    private List<String> domains;

    public SachsUtils() {
        interventions = new ArrayList<>();
        interventions.add("cd3_cd28");
        interventions.add("icam2");
        interventions.add("aktinhib");
        interventions.add("g0076");
        interventions.add("psitect");
        interventions.add("u0126");
        interventions.add("ly");
        interventions.add("pma");
        interventions.add("b2camp");

        domains = new ArrayList<>();
        domains.add("raf");
        domains.add("mek");
        domains.add("plc");
        domains.add("pip2");
        domains.add("pip3");
        domains.add("erk");
        domains.add("akt");
        domains.add("pka");
        domains.add("pkc");
        domains.add("p38");
        domains.add("jnk");
    }

    public Knowledge2 getKnowledge(){
        Knowledge2 knowledge = new Knowledge2();
        for (String intervention : interventions) {
            knowledge.addToTier(0, intervention);
        }
        for (String domain : domains) {
            knowledge.addToTier(1, domain);
        }
        knowledge.setTierForbiddenWithin(0, true);

        return knowledge;
    }

    public Graph pruneGraph(Graph graph) {
        for (
                Edge edge : graph.getEdges()) {
            if (interventions.contains(edge.getNode1().getName())) {
                graph.removeEdge(edge);
            } else if (interventions.contains(edge.getNode2().getName())) {
                graph.removeEdge(edge);
            }
        }

        for (
                Node node : graph.getNodes()) {
            if (interventions.contains(node.getName())) {
                graph.removeNode(node);
            }
        }

        return graph;
    }
}