/**
 * 
 */
package edu.cmu.tetradapp.plugin;

import java.util.List;

import org.springframework.beans.factory.annotation.Autowired;

import edu.pitt.dbmi.tetrad.plugin.api.PluginAlgorithm;

/**
 * Sep 20, 2018 5:25:05 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
public class PluginAlgorithms {

	@Autowired
	private List<PluginAlgorithm> pluginAlgorithms;
	
	public void printPluginAlgorithms() {
		System.out.println(String.format("Found %d plugin algorithms for extension point '%s'", pluginAlgorithms.size(), PluginAlgorithm.class.getName()));
		for(PluginAlgorithm algorithm : pluginAlgorithms) {
			System.out.println("Plugin Algorithm: " + algorithm.getClass().getName()); 
		}
	}
	
}
