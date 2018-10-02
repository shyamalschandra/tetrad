package edu.cmu.tetrad.plugin.algorithm;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.util.Parameters;

/**
 * Oct 1, 2018 5:34:22 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
public class PluginAlgorithmScoreWrapper extends PluginAlgorithmWrapper implements ScoreWrapper {

	private static final long serialVersionUID = 1L;
	
	private final Class<?> extensionClass;
	private final Object algorithm;

	/**
	 * @param clazz
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 */
	public PluginAlgorithmScoreWrapper(Class<?> extensionClass) throws InstantiationException, IllegalAccessException {
		super(extensionClass);
		this.extensionClass = extensionClass;
		this.algorithm = super.getAlgorithm();
	}

	@Override
	public Score getScore(DataModel dataSet, Parameters parameters) {
		Score result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getScore", DataModel.class, Parameters.class);
			result = (Score) method.invoke(algorithm, dataSet, parameters);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
		return result;
	}

	@Override
	public Node getVariable(String name) {
		Node result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getVariable", String.class);
			result = (Node) method.invoke(algorithm, name);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
		return result;
	}

}
