/**
 * 
 */
package edu.cmu.tetrad.plugin.algorithm;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.List;

import com.google.gson.Gson;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.plugin.PluginExtension;
import edu.cmu.tetrad.util.Parameters;

/**
 * Oct 1, 2018 6:11:57 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
public class PluginAlgorithmInitialGraphKnowledgeScoreWrapper implements PluginExtension, Algorithm, TakesInitialGraph, HasKnowledge, UsesScoreWrapper {

	private static final long serialVersionUID = 1L;
	
	private final Class<?> extensionClass;
	private final Object algorithm;
	
	private Gson gson = new Gson();

	/**
	 * @param extensionClass
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 */
	public PluginAlgorithmInitialGraphKnowledgeScoreWrapper(Class<?> extensionClass)
			throws InstantiationException, IllegalAccessException {
		this.extensionClass = extensionClass;
		this.algorithm = extensionClass.newInstance();
		Method[] mothods = extensionClass.getDeclaredMethods();
		if (mothods != null) {
			for (int i = 0; i < mothods.length; i++) {
				Method method = mothods[i];
				String methodName = method.getName();
				if (methodName.equalsIgnoreCase("getAlgorithmDescriptions")) {
					try {
						String description = (String) method.invoke(algorithm);
						System.out.println("getAlgorithmDescriptions: " + description);
					} catch (IllegalAccessException | IllegalArgumentException
							| InvocationTargetException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
				}
			}
		}
	}

	@Override
	public Graph search(DataModel dataSet, Parameters parameters) {
		Graph result = null;
		try {
			String dataSetJson = gson.toJson(dataSet);
			String parametersJson = gson.toJson(parameters);
			System.out.println("dataSetJson: ");
			System.out.println(dataSetJson);
			Method method = extensionClass.getDeclaredMethod("search", String.class, String.class);
			result = (Graph) method.invoke(algorithm, dataSetJson, parametersJson);
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
			Object obj = method.invoke(algorithm);
			result = DataType.valueOf(obj.toString());
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	@Override
	public List<String> getParameters() {
		List<String> result = null;
		/*try {
			Method method = extensionClass.getDeclaredMethod("getParameters");
			result = (List<String>) method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}*/
		
		Method[] mothods = extensionClass.getDeclaredMethods();
		if (mothods != null) {
			for (int i = 0; i < mothods.length; i++) {
				Method method = mothods[i];
				String methodName = method.getName();
				if (methodName.equalsIgnoreCase("getParameters")) {
					try {
						result = (List<String>)method.invoke(algorithm);
					} catch (IllegalAccessException | IllegalArgumentException
							| InvocationTargetException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					break;
				}
			}
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
			Method method = extensionClass.getDeclaredMethod("setKnowledge", String.class);
			String json = gson.toJson(knowledge);
			method.invoke(algorithm, json);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
		
		/*Method[] mothods = extensionClass.getDeclaredMethods();
		if (mothods != null) {
			for (int i = 0; i < mothods.length; i++) {
				Method method = mothods[i];
				String methodName = method.getName();
				if (methodName.equalsIgnoreCase("setKnowledge")) {
					Class<?>[] paramTypeClasses = method.getParameterTypes();
					
					if(paramTypeClasses != null) {
						for(int j=0;j<paramTypeClasses.length;j++) {
							Class<?> paramTypeClass = paramTypeClasses[j];
							System.out.println("setKnowledge: paramTypeClass: " + paramTypeClass.getCanonicalName());
							if(paramTypeClass.getCanonicalName().contains("java.lang.String")) {
								try {
									
									String json = gson.toJson(knowledge);
									method.invoke(algorithm, json);
								} catch (IllegalAccessException | IllegalArgumentException
										| InvocationTargetException e1) {
									e1.printStackTrace();
								}
								break;
							}
						}
						
					}
					
				}
			}
		}*/
	}

	@Override
	public Graph getInitialGraph() {
		Graph result = null;
		try {
			Method method = extensionClass.getDeclaredMethod("getInitialGraph");
			result = (Graph) method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}

		return result;
	}

	@Override
	public void setInitialGraph(Graph initialGraph) {
		try {
			Method method = extensionClass.getDeclaredMethod("setInitialGraph",Graph.class);
			method.invoke(algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void setInitialGraph(Algorithm algorithm) {
		try {
			Method method = extensionClass.getDeclaredMethod("setInitialGraph",Algorithm.class);
			method.invoke(this.algorithm, algorithm);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException
				| SecurityException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void setScoreWrapper(ScoreWrapper score) {
		Method[] mothods = extensionClass.getDeclaredMethods();
		if (mothods != null) {
			for (int i = 0; i < mothods.length; i++) {
				Method method = mothods[i];
				String methodName = method.getName();
				if (methodName.equalsIgnoreCase("setScoreWrapper")) {
					Class<?>[] paramTypeClasses = method.getParameterTypes();
					
					if(paramTypeClasses != null) {
						for(int j=0;j<paramTypeClasses.length;j++) {
							Class<?> paramTypeClass = paramTypeClasses[j];
							if(paramTypeClass.getCanonicalName().contains("java.lang.String")) {
								try {
									method.invoke(algorithm, score.getClass().getCanonicalName());
								} catch (IllegalAccessException | IllegalArgumentException
										| InvocationTargetException e1) {
									e1.printStackTrace();
								}
								break;
							}
						}
						
					}
					
				}
			}
		}
		
	}

	@Override
	public Class<?> getExtensionClass() {
		return extensionClass;
	}
}
