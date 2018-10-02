package edu.cmu.tetrad.plugin.algorithm;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.List;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.Parameters;

/**
 * Oct 1, 2018 4:51:25 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
public class PluginAlgorithmWrapper implements Algorithm {

	private static final long serialVersionUID = 1L;
	
	private final Class<?> extensionClass;
	private final Object algorithm;

	public PluginAlgorithmWrapper(Class<?> extensionClass) throws InstantiationException, IllegalAccessException {
		this.extensionClass = extensionClass;
		this.algorithm = extensionClass.newInstance();
	}

	@Override
	public Graph search(DataModel dataSet, Parameters parameters) {
		Graph result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("search", DataModel.class, Parameters.class);
			result = (Graph) method.invoke(algorithm, dataSet, parameters);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	@Override
	public Graph getComparisonGraph(Graph graph) {
		Graph result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getComparisonGraph", Graph.class);
			result = (Graph) method.invoke(algorithm, graph);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	@Override
	public String getDescription() {
		String result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getDescription");
			result = (String) method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	@Override
	public DataType getDataType() {
		DataType result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getDataType");
			result = (DataType) method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	@Override
	public List<String> getParameters() {
		List<String> result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getParameters");
			result = (List<String>) method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	public Class<?> getExtensionClass() {
		return extensionClass;
	}

	protected Object getAlgorithm() {
		return algorithm;
	}
}
