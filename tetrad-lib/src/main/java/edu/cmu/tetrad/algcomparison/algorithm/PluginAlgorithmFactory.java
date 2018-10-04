package edu.cmu.tetrad.algcomparison.algorithm;

import java.lang.reflect.AnnotatedType;

import edu.cmu.tetrad.algcomparison.independence.IndependenceWrapper;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.plugin.algorithm.PluginAlgorithmInitialGraphKnowledgeScoreWrapper;

/**
 * Oct 3, 2018 4:40:32 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
public class PluginAlgorithmFactory {

	public static Algorithm create(Class<?> algoClass, IndependenceWrapper test, ScoreWrapper score)
			throws IllegalAccessException, InstantiationException {
        if (algoClass == null) {
            throw new IllegalArgumentException("Algorithm class cannot be null.");
        }

		AnnotatedType[] annotatedTypes = algoClass.getAnnotatedInterfaces();
		boolean scoreRequired = false;
		boolean testRequired = false;
		if (annotatedTypes != null) {
			for (int i = 0; i < annotatedTypes.length; i++) {
				AnnotatedType annotatedType = annotatedTypes[i];
				String strAnnotatedType = annotatedType.getType().toString();

				if (strAnnotatedType.contains("edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper")) {
					scoreRequired = true;
				}
				if (strAnnotatedType.contains("edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper")) {
					testRequired = true;
				}
			}
		}
		
		if (testRequired && test == null) {
            throw new IllegalArgumentException("Test of independence is required.");
        }

        if (scoreRequired && score == null) {
            throw new IllegalArgumentException("Score is required.");
        }
        
        Algorithm algorithm =
				new PluginAlgorithmInitialGraphKnowledgeScoreWrapper(algoClass);
        if (testRequired) {
            ((TakesIndependenceWrapper) algorithm).setIndependenceWrapper(test);
        }
        if (scoreRequired) {
            ((UsesScoreWrapper) algorithm).setScoreWrapper(score);
        }
		
		return algorithm;
	}

	public static Algorithm create(Class<?> algoClass, Class<? extends IndependenceWrapper> indTestClass, Class<? extends ScoreWrapper> scoreClass)
			throws IllegalAccessException, InstantiationException {
		IndependenceWrapper test = (indTestClass == null) ? null : indTestClass.newInstance();
        ScoreWrapper score = (scoreClass == null) ? null : scoreClass.newInstance();
        
        return create(algoClass, test, score);
	}
	
}
