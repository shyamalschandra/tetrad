package edu.cmu.tetrad.util;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import com.google.gson.Gson;

import edu.cmu.tetrad.data.BoxDataSet;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.data.VerticalIntDataBox;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraphSingleConnections;
import edu.cmu.tetrad.graph.EdgeTypeProbability;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphNode;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.graph.Triple;
import edu.cmu.tetrad.graph.EdgeTypeProbability.EdgeType;

/**
 * 
 * Dec 9, 2016 5:43:47 PM
 * 
 * @author Chirayu (Kong) Wongchokprasitti, PhD
 * 
 */
public class JsonUtils {

	public static Graph parseJSONObjectToTetradGraph(String jsonResponse) {
		return parseJSONObjectToTetradGraph(new JSONObject(jsonResponse));
	}

	public static Graph parseJSONObjectToTetradGraph(JSONObject jObj) {
		if(!jObj.isNull("graph")) {
			return parseJSONObjectToTetradGraph(jObj.getJSONObject("graph"));
		}
		
		// Node
		List<Node> nodes = parseJSONArrayToTetradNodes(jObj.getJSONArray("nodes"));
		EdgeListGraphSingleConnections graph = new EdgeListGraphSingleConnections(nodes);

		// Edge
		Set<Edge> edges = parseJSONArrayToTetradEdges(graph, jObj.getJSONArray("edgesSet"));
		for (Edge edge : edges) {
			graph.addEdge(edge);
		}

		// ambiguousTriples
		Set<Triple> ambiguousTriples = parseJSONArrayToTetradTriples(jObj.getJSONArray("ambiguousTriples"));
		for (Triple triple : ambiguousTriples) {
			graph.addAmbiguousTriple(triple.getX(), triple.getY(), triple.getZ());
		}

		// underLineTriples
		Set<Triple> underLineTriples = parseJSONArrayToTetradTriples(jObj.getJSONArray("underLineTriples"));
		for (Triple triple : underLineTriples) {
			graph.addUnderlineTriple(triple.getX(), triple.getY(), triple.getZ());
		}

		// dottedUnderLineTriples
		Set<Triple> dottedUnderLineTriples = parseJSONArrayToTetradTriples(jObj.getJSONArray("dottedUnderLineTriples"));
		for (Triple triple : dottedUnderLineTriples) {
			graph.addDottedUnderlineTriple(triple.getX(), triple.getY(), triple.getZ());
		}

		// stuffRemovedSinceLastTripleAccess
		boolean stuffRemovedSinceLastTripleAccess = jObj.getBoolean("stuffRemovedSinceLastTripleAccess");
		graph.setStuffRemovedSinceLastTripleAccess(stuffRemovedSinceLastTripleAccess);

		// highlightedEdges
		Set<Edge> highlightedEdges = parseJSONArrayToTetradEdges(graph, jObj.getJSONArray("highlightedEdges"));
		for (Edge edge : highlightedEdges) {
			graph.setHighlighted(edge, true);
		}

		return graph;
	}

	public static Set<Triple> parseJSONArrayToTetradTriples(JSONArray jArray) {
		Set<Triple> triples = new HashSet<>();

		for (int i = 0; i < jArray.length(); i++) {
			Triple triple = parseJSONArrayToTetradTriple(jArray.getJSONObject(i));
			triples.add(triple);
		}

		return triples;
	}

	public static Triple parseJSONArrayToTetradTriple(JSONObject jObj) {
		Node x = parseJSONObjectToTetradNode(jObj.getJSONObject("x"));
		Node y = parseJSONObjectToTetradNode(jObj.getJSONObject("y"));
		Node z = parseJSONObjectToTetradNode(jObj.getJSONObject("z"));

		return new Triple(x, y, z);
	}

	public static Set<Edge> parseJSONArrayToTetradEdges(Graph graph, JSONArray jArray) {
		Set<Edge> edges = new HashSet<>();

		for (int i = 0; i < jArray.length(); i++) {
			Edge edge = parseJSONObjectToTetradEdge(graph, jArray.getJSONObject(i));
			edges.add(edge);
		}

		return edges;
	}

	public static Edge parseJSONObjectToTetradEdge(Graph graph, JSONObject jObj) {
		Node node1 = graph.getNode(jObj.getJSONObject("node1").getString("name"));
		Node node2 = graph.getNode(jObj.getJSONObject("node2").getString("name"));
		Endpoint endpoint1 = Endpoint.TYPES[jObj.getJSONObject("endpoint1").getInt("ordinal")];
		Endpoint endpoint2 = Endpoint.TYPES[jObj.getJSONObject("endpoint2").getInt("ordinal")];
		Edge edge = new Edge(node1, node2, endpoint1, endpoint2);
		
		try {
		    // properties
		    JSONArray jArray = jObj.getJSONArray("properties");
		    if(jArray != null){
		    	for (int i = 0; i < jArray.length(); i++) {
		    		edge.addProperty(parseJSONObjectToEdgeProperty(jArray.getString(i)));
		    	}
		    }
		} catch (JSONException e) {
		    // TODO Auto-generated catch block
		    //e.printStackTrace();
		}
		
		try {
		    // properties
		    JSONArray jArray = jObj.getJSONArray("edgeTypeProbabilities");
		    if(jArray != null){
		    	for (int i = 0; i < jArray.length(); i++) {
		    		edge.addEdgeTypeProbability(parseJSONObjectToEdgeTypeProperty(jArray.getJSONObject(i)));
		    	}
		    }
		} catch (JSONException e) {
		    // TODO Auto-generated catch block
		    //e.printStackTrace();
		}
		
		return edge;
	}
	
	public static EdgeTypeProbability parseJSONObjectToEdgeTypeProperty(JSONObject jObj){
		String _edgeType = jObj.getString("edgeType");
		EdgeType edgeType = EdgeType.nil;
		switch(_edgeType){
			case "ta" : edgeType = EdgeType.ta; break;
			case "at" : edgeType = EdgeType.at; break;
			case "ca" : edgeType = EdgeType.ca; break;
			case "ac" : edgeType = EdgeType.ac; break;
			case "cc" : edgeType = EdgeType.cc; break;
			case "aa" : edgeType = EdgeType.aa; break;
			case "tt" : edgeType = EdgeType.tt; break;
		}
		double probability = jObj.getDouble("probability");
		EdgeTypeProbability edgeTypeProbability = new EdgeTypeProbability(edgeType, probability);
		
		try {
		    // properties
		    JSONArray jArray = jObj.getJSONArray("properties");
		    if(jArray != null){
		    	for (int i = 0; i < jArray.length(); i++) {
		    		edgeTypeProbability.addProperty(parseJSONObjectToEdgeProperty(jArray.getString(i)));
		    	}
		    }
		} catch (JSONException e) {
		    // TODO Auto-generated catch block
		    //e.printStackTrace();
		}
		return edgeTypeProbability;
	}

	public static Edge.Property parseJSONObjectToEdgeProperty(String prop){
		if(prop.equalsIgnoreCase("dd")){
			return Edge.Property.dd;
		}
		if(prop.equalsIgnoreCase("nl")){
			return Edge.Property.nl;
		}
		if(prop.equalsIgnoreCase("pd")){
			return Edge.Property.pd;
		}
		if(prop.equalsIgnoreCase("pl")){
			return Edge.Property.pl;
		}
		return null;
	}
	
	public static List<Node> parseJSONArrayToTetradNodes(JSONArray jArray) {
		List<Node> nodes = new ArrayList<>();

		for (int i = 0; i < jArray.length(); i++) {
			Node node = parseJSONObjectToTetradNode(jArray.getJSONObject(i));
			nodes.add(node);
		}

		return nodes;
	}

	public static Node parseJSONObjectToTetradNode(JSONObject jObj) {
		JSONObject nodeType = jObj.getJSONObject("nodeType");
		int ordinal = nodeType.getInt("ordinal");
		int centerX = jObj.getInt("centerX");
		int centerY = jObj.getInt("centerY");
		String name = jObj.getString("name");

		GraphNode graphNode = new GraphNode(name);
		graphNode.setNodeType(NodeType.TYPES[ordinal]);
		graphNode.setCenter(centerX, centerY);

		return graphNode;
	}
	
	public static DataSet parseJSONObjectToDataSet(String dataSetJson, DataType dataType) {
		return parseJSONObjectToDataSet(new JSONObject(dataSetJson), dataType);
	}

	public static DataSet parseJSONObjectToDataSet(JSONObject jObj, DataType dataType) {
		//Name
		String name = jObj.getString("name");
		
		//Variables
		List<Node> nodes = parseJSONArrayToTetradNodes(jObj.getJSONArray("variables"));
		int columns = nodes.size();
		System.out.println("columns: " + columns);
		
		//DataBox
		BoxDataSet boxDataSet = null;
		int rows = jObj.getJSONObject("dataBox").getJSONArray("data").getJSONArray(0).length();
		System.out.println("rows: " + rows);
		
		switch(dataType) {
		case Continuous:
			
			
			break;
		case Discrete:
			VerticalIntDataBox dataBox = new VerticalIntDataBox(rows, columns);
			List<Node> node_list = new ArrayList<>();
			for(Node node : nodes) {
				Node discreteVar = new DiscreteVariable(node.getName());
				node_list.add(discreteVar);
			}
			for(int column=0;column<columns;column++) {
				JSONArray jArray = jObj.getJSONObject("dataBox").getJSONArray("data").getJSONArray(column);
				for(int row=0;row<jArray.length();row++) {
					int value = jArray.getInt(row);
					dataBox.set(row, column, value);
				}
			}
			boxDataSet = new BoxDataSet(dataBox, node_list);
			
			break;
			
		case Mixed:
			
			break;
			
		
		default:
			break;
		}
		
		boxDataSet.setName(name);
		
		//Selection
		List<Node> selections = parseJSONArrayToTetradNodes(jObj.getJSONArray("variables"));
		if(!selections.isEmpty()) {
			for(Node node : selections) {
				boxDataSet.setSelected(node, true);
			}
		}
		
		//CaseIds
		
		//multipliers
		
		//knowledge
		String knowledgeJson = jObj.getJSONObject("knowledge").toString();
		Gson gson = new Gson();
		IKnowledge iknowledge = gson.fromJson(knowledgeJson, Knowledge2.class);
		boxDataSet.setKnowledge(iknowledge);
		
		//outputDelimiter
		
		return boxDataSet;
	}
	
}
