/**
 * 
 */
package edu.cmu.tetradapp.plugin;

import org.pf4j.spring.SpringPluginManager;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

/**
 * Sep 20, 2018 5:23:17 PM
 *
 * @author Chirayu Kong Wongchokprasitti, PhD (chw20@pitt.edu)
 *
 */
@Configuration
public class SpringConfiguration {

	@Bean
    public SpringPluginManager pluginManager() {
        return new SpringPluginManager();
    }
	
	
}
