package edu.cmu.tetrad.intervention.real;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.Knowledge2;

import java.util.*;

/**
 * Created by bandrews on 9/14/17.
 */

public class Knowledge {

    private DataSet dataSet;

    public Knowledge(DataSet dataSet) {
        this.dataSet = dataSet;
    }

    public Knowledge2 get_knowledge() {

        Knowledge2 knowledge = new Knowledge2(this.dataSet.getVariableNames());

        List<String> tier0 = new ArrayList<>();
        List<String> tier1 = new ArrayList<>();
        for (String var : this.dataSet.getVariableNames()) {
            if (var.startsWith("I")) {
                tier0.add(var);
            } else {
                tier1.add(var);
            }
        }
        knowledge.setTier(0, tier0);
        knowledge.setTier(1, tier1);
        knowledge.setTierForbiddenWithin(0, true);

        return knowledge;
    }
}