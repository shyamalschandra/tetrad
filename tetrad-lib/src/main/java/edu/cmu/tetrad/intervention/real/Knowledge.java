package edu.cmu.tetrad.intervention.real;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.Knowledge2;

import java.util.*;

/**
 * Created by bandrews on 9/14/17.
 */

public class Knowledge {

    private DataSet meta_data;
    private DataSet dataSet;

    public Knowledge(DataSet dataSet, DataSet meta_data) {
        this.dataSet = dataSet;
        this.meta_data = meta_data;
    }

    public Knowledge2 get_knowledge() {

        Knowledge2 knowledge = new Knowledge2(this.dataSet.getVariableNames());

        List<String> tier0 = new ArrayList<>();
        List<String> tier1 = new ArrayList<>();
        List<String> tier2 = new ArrayList<>();
        for (String var : this.dataSet.getVariableNames()) {
            if (var.startsWith("IC")) {
                tier0.add(var);
            } else if (var.startsWith("IV")) {
                tier1.add(var);
            } else {
                tier2.add(var);
            }
        }
        knowledge.setTier(0, tier0);
        knowledge.setTier(1, tier1);
        knowledge.setTier(2, tier2);
//        knowledge.setTierForbiddenWithin(0, true);
//        knowledge.setTierForbiddenWithin(1, true);
//
//        for (String a : knowledge.getTier(0)) {
//            for (String b : knowledge.getTier(1)) {
//                knowledge.setForbidden(a, b);
//            }
//        }

        List<String> nodes = new ArrayList<>();
        for (String a : knowledge.getTier(0)) {
            nodes.add(a);
        }
        for (String b : knowledge.getTier(1)) {
            nodes.add(b);
        }
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i+1; j < nodes.size(); j++) {
                knowledge.setRequired(nodes.get(i), nodes.get(j));
            }
        }

        for (String as : knowledge.getTier(0)) {
            for (String a : as.split(":")) {
                for (String bs : knowledge.getTier(1)) {
                    for (String b : bs.split(":")) {
                        if(a.substring(2) == b.substring(2)) {
                            knowledge.removeForbidden(as, bs);
                            knowledge.setRequired(as, bs);
                        }
                    }
                }
            }
        }

        HashMap<String, Set<String>> required = new HashMap<>();
        HashMap<String, Set<String>> forbidden = new HashMap<>();
        for (int i = 0; i < this.meta_data.getNumRows(); i++) {
            List<String> Is = new ArrayList<>();
            List<String> Drs = new ArrayList<>();
            List<String> Dfs = new ArrayList<>();
            for (int j = 0; j < this.meta_data.getNumColumns(); j++) {
                if (this.meta_data.getInt(i, j) != 0) {
                    if (this.meta_data.getVariable(j).getName().startsWith("I_")) {
                        Is.add(this.meta_data.getVariable(j).getName());
                    } else {
                        if (this.meta_data.getInt(i, j) == 1) {
                            Drs.add(this.meta_data.getVariable(j).getName());
                        } else {
                            Dfs.add(this.meta_data.getVariable(j).getName());
                        }
                    }
                }
            }
            for (String I : Is) {
                for (String Dr : Drs) {
                    if (!required.containsKey(I)) {
                        required.put(I, new HashSet<String>());
                    }
                    required.get(I).add(Dr);
                }
                for (String Df : Dfs) {
                    if (!forbidden.containsKey(I)) {
                        forbidden.put(I, new HashSet<String>());
                    }
                    forbidden.get(I).add(Df);

                }
            }
        }

        for (String as : knowledge.getTier(1)) {
            for (String a : as.split(":")) {
                for (String bs : knowledge.getTier(2)) {
                    for (String b : bs.split(":")) {
                        if (required.get("I" + a.substring(2)).contains(b)) {
                            knowledge.setRequired(as, bs);
                        } else if (forbidden.get("I" + a.substring(2)).contains(b)) {
                            knowledge.setForbidden(as, bs);
                        }
                    }
                }
            }
        }

        return knowledge;
    }
}