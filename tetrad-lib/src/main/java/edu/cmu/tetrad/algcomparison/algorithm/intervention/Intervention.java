package edu.cmu.tetrad.algcomparison.algorithm.intervention;

import org.apache.commons.lang3.RandomUtils;

import java.util.*;

/**
 * Created by bandrews on 6/20/17.
 */

// KNOWN ISSUE, ERRORS OCCUR WHEN NUMBER OF INTERVENTIONS CAN BE GREATER THAN NUMBER OF CATEGORIES //

public class Intervention {

    private List<Integer> domains;
    private List<Double> interventionValues;
    private HashMap<String, List<Double>> effectValues = new HashMap<>();
    private boolean isDiscrete;
    private double potency;
    private int samplesSize;
    private int interventionSize;
    private String name = "";
    private Random rng = new Random();


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
            double beta = RandomUtils.nextDouble(low, high);
            for (int i = 0; i < this.samplesSize; i++) {
                values.add(i, beta * this.interventionValues.get(i) + this.rng.nextGaussian() * Math.sqrt(epsilon));
            }
        }
        this.effectValues.remove(effect);
        this.effectValues.put(effect, values);
    }

    public void addEffect(String effect) {
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
            //double a = sorted[this.samplesSize * 1/3];
            //double b = sorted[this.samplesSize * 2/3];
            List<Double> effectValue = new ArrayList<>();
            for (Double value : this.interventionValues) {
                if (value < a) {
                    effectValue.add(1.0);
                //} else if(value < b) {
                //    effectValue.add(1.0);
                } else {
                    effectValue.add(2.0);
                }
            }
            this.effectValues.put(effect, effectValue);
        }
    }

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
