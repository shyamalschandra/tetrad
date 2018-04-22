package edu.cmu.tetrad.intervention;

import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by bandrews on 6/19/17.
 */
public class InterventionalKnowledge {

    private List<Node> measured;
    private List<Node> IDs;
    private List<Node> IVs;
    private List<Node> C;
    private Knowledge2 knowledge;

    public InterventionalKnowledge(Graph graph, int hand) {
        this.measured = new ArrayList<>();
        this.IDs = new ArrayList<>();
        this.IVs = new ArrayList<>();
        this.C = new ArrayList<>();
        this.knowledge = new Knowledge2();
        for (Node node : graph.getNodes()) {
            if (node.getName().startsWith("ID")) {
                this.IDs.add(node);
            } else if (node.getName().startsWith("IV")) {
                this.IVs.add(node);
            } else if (node.getName().startsWith("C")) {
                this.C.add(node);
            } else if (node.getNodeType() == NodeType.MEASURED) {
                this.measured.add(node);
            }
        }
        defineKnowledge(hand);
    }

    public Knowledge2 getKnowledge(){
        return this.knowledge;
    }

    private void defineKnowledge(int hand) {

        // Tiers
        for (Node node : this.measured) {
            this.knowledge.addToTier(1, node.getName());
        }
        for (Node node : this.IDs) {
            this.knowledge.addToTier(0, node.getName());
        }
        for (Node node : this.IVs) {
            this.knowledge.addToTier(0, node.getName());
        }
        for (Node node : this.C) {
            this.knowledge.addToTier(0, node.getName());
        }

        this.knowledge.setTierForbiddenWithin(0, true);

        if(hand > 0) {
            // Require
            for (Node node : this.IDs) {
                for (String intervened : getIntervened(node.getName())) {
                    this.knowledge.setRequired(node.getName(), intervened);
                }
            }
            for (Node node : this.IVs) {
                for (String intervened : getIntervened(node.getName())) {
                    this.knowledge.setRequired(node.getName(), intervened);
                }
            }
        }

        if(hand > 1) {
            // Forbidden
            for (Node node : this.IDs) {
                List<String> intervened = getIntervened(node.getName());
                for (Node measuredNode : this.measured) {
                    if (!intervened.contains(measuredNode.getName())) {
                        this.knowledge.setForbidden(node.getName(), measuredNode.getName());
                    }
                }
            }
            for (Node node : this.IVs) {
                List<String> intervened = getIntervened(node.getName());
                for (Node measuredNode : this.measured) {
                    if (!intervened.contains(measuredNode.getName())) {
                        this.knowledge.setForbidden(node.getName(), measuredNode.getName());
                    }
                }
            }
        }
    }

    private List<String> getIntervened(String intervention) {
        List<String> interventionNodes = new ArrayList<>();
        for (String intervened : intervention.split("_")) {
            interventionNodes.add(intervened);
        }
        interventionNodes.remove(0);
        return interventionNodes;
    }

}
