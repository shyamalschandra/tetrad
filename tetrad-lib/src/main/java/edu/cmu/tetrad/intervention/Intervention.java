package edu.cmu.tetrad.intervention;

import org.apache.commons.lang3.RandomUtils;

import java.util.*;

/**
 * Created by bandrews on 6/20/17.
 * This class hold information regarding interventions for the conditional
 * Gaussian simulation CGISimulation.
 */

// KNOWN ISSUE, ERRORS OCCUR WHEN NUMBER OF INTERVENTIONS CAN BE GREATER THAN NUMBER OF CATEGORIES //

public class Intervention {

    private List<Integer> domains; // This list hold information regarding what domain the data at each instance is generated from (observed / intervention #).
    private List<Double> interventionValues;  // This list holds the value of the intervention at each instance.
    private HashMap<String, List<Double>> effectValues = new HashMap<>(); // Maps a domain variable's name to the instance values it takes on given the intervention.
    private boolean isDiscrete; // Is the intervention discrete.
    private double potency; // How potent is the intervention (soft / hard).
    private int samplesSize; // Observational sample size.
    private int interventionSize; // Sample size of each intervention.
    private String name = ""; // ID of intervention.
    private Random rng = new Random(); // A random number generator.


    public Intervention(boolean isDiscrete, double potency, int samplesSize, int interventionSize, int startIndex, int numICategories, double mean, double var) {
        this.isDiscrete = isDiscrete;
        this.potency = potency;
        this.samplesSize = samplesSize;
        this.interventionSize = interventionSize;
        this.domains = new ArrayList<>();
        this.interventionValues = new ArrayList<>();
        for (int i = 0; i < startIndex; i++) {
            this.domains.add(0);
            if (isDiscrete) {
                int category = this.rng.nextInt(numICategories) + 1;
                this.interventionValues.add((double) category);
            } else {
                double value = mean + rng.nextGaussian() * Math.sqrt(var);
                this.interventionValues.add(value);
            }
        }
        if (isDiscrete) {
            for (int i = 0; i < this.interventionSize; i++) {
                int category = this.rng.nextInt(numICategories) + 1;
                this.domains.add(category);
                this.interventionValues.add((double) category);
            }
        } else {
            for (int i = 0; i < this.interventionSize; i++) {
                double value = mean + rng.nextGaussian() * Math.sqrt(var);
                this.domains.add(1);
                this.interventionValues.add(value);
            }
        }
        while (this.domains.size() < this.samplesSize) {
            this.domains.add(0);
            if (isDiscrete) {
                int category = this.rng.nextInt(numICategories) + 1;
                this.interventionValues.add((double) category);
            } else {
                double value = mean + rng.nextGaussian() * Math.sqrt(var);
                this.interventionValues.add(value);
            }
        }

    }

    public void changeEffect(String effect, double low, double high, double epsilon) {

        // addEffect assume that the domain variables targeted by the intervention is discrete; this method
        // corrects for the case where that assumption is incorrect and the domain variable is continuous.

        List<Double> values = new ArrayList<>();
        if (this.isDiscrete) {
            HashMap<Integer, Double> map = new HashMap<>();
            for (Integer value : this.domains) {
                if (!map.containsKey(value)) {
                    map.put(value, RandomUtils.nextDouble(low, high));
                }
                values.add(map.get(value) + this.rng.nextGaussian() * Math.sqrt(epsilon));
            }
        } else {
            double beta = RandomUtils.nextDouble(low, high) * (2*RandomUtils.nextInt(0,2)-1);
            for (int i = 0; i < this.samplesSize; i++) {
                values.add(i, beta * this.interventionValues.get(i) + this.rng.nextGaussian() * Math.sqrt(epsilon));
            }
        }
        this.effectValues.remove(effect);
        this.effectValues.put(effect, values);
    }

    public void addEffect(String effect) {

        // Add a domain variables affected by the intervention.

        this.name += "_" + effect;
        if (this.isDiscrete) {
            this.effectValues.put(effect, this.interventionValues);
        } else {
            double[] sorted = new double[this.samplesSize];
            for (int i = 0; i < this.samplesSize; i++) {
                sorted[i] = this.interventionValues.get(i);
            }
            Arrays.sort(sorted);
            double a = sorted[this.samplesSize * 1/2];
            List<Double> effectValue = new ArrayList<>();
            for (Double value : this.interventionValues) {
                if (value < a) {
                    effectValue.add(1.0);
                } else {
                    effectValue.add(2.0);
                }
            }
            this.effectValues.put(effect, effectValue);
        }
    }

    // Getters and setters

    public boolean isDiscrete() {
        return isDiscrete;
    }

    public Set<String> getEffected() {
        return this.effectValues.keySet();
    }

    public double getPotency() {
        return potency;
    }

    public String getName() {
        return name;
    }

    public int getDomain(int index) {
        return domains.get(index);
    }

    public double getValue(int index) {
        return interventionValues.get(index);
    }

    public double getInterventionValue(String effected, int index) {
        return effectValues.get(effected).get(index);
    }
}
