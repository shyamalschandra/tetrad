/**
 * 
 */
package edu.cmu.tetrad.plugin.algorithm;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.util.Parameters;

/**
 * Oct 1, 2018 5:56:14 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
public class PluginAlgorithmKnowledgeScoreWrapper extends PluginAlgorithmWrapper implements HasKnowledge, ScoreWrapper {

	private static final long serialVersionUID = 1L;
	
	private final Class<?> extensionClass;
	private final Object algorithm;

	/**
	 * @param extensionClass
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 */
	public PluginAlgorithmKnowledgeScoreWrapper(Class<?> extensionClass)
			throws InstantiationException, IllegalAccessException {
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

	@Override
	public IKnowledge getKnowledge() {
		IKnowledge result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getKnowledge");
			result = (IKnowledge) method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
		return result;
	}

	@Override
	public void setKnowledge(IKnowledge knowledge) {
		try {
			Method method = extensionClass.getDeclaredMethod("setKnowledge", IKnowledge.class);
			method.invoke(algorithm, knowledge);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
	}

}
