/*
 * Copyright (C) 2017 University of Pittsburgh.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package edu.cmu.tetradapp.ui.model;

import edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper;
import edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper;
import edu.cmu.tetrad.annotation.Algorithm;
import edu.cmu.tetrad.annotation.AnnotatedClass;
import edu.cmu.tetrad.util.AlgorithmDescriptions;
import java.io.Serializable;
import java.lang.annotation.Annotation;
import java.lang.reflect.AnnotatedType;

/**
 *
 * Nov 30, 2017 4:41:37 PM
 *
 * @author Kevin V. Bui (kvb2@pitt.edu)
 */
public class AlgorithmModel implements Serializable, Comparable<AlgorithmModel> {

    private static final long serialVersionUID = 8599854464475682558L;

    private final AnnotatedClass<Algorithm> algorithm;
    private final boolean requiredScore;
    private final boolean requiredTest;
    private final String description;

    public AlgorithmModel(AnnotatedClass<Algorithm> algorithm) {
        if (algorithm == null) {
            throw new IllegalArgumentException("Algorithm annotation cannot be null.");
        }

        this.algorithm = algorithm;
        if(algorithm.getAnnotation() == null) {
        	AnnotatedType[] annotatedTypes = algorithm.getClazz().getAnnotatedInterfaces();
    		boolean requiredScore = false;
    		boolean requiredTest = false;
        	if(annotatedTypes != null) {
				for(int i=0;i<annotatedTypes.length;i++) {
					AnnotatedType annotatedType = annotatedTypes[i];
					String strAnnotatedType = annotatedType.getType().toString();
					
					if(strAnnotatedType.contains("edu.cmu.tetrad.algcomparison.utils.UsesScoreWrapper")) {
						requiredScore = true;
					}
					if(strAnnotatedType.contains("edu.cmu.tetrad.algcomparison.utils.TakesIndependenceWrapper")) {
						requiredTest = true;
					}
				}
			}
        	
            this.requiredScore = requiredScore;
            this.requiredTest = requiredTest;
        }else {
            this.requiredScore = UsesScoreWrapper.class.isAssignableFrom(algorithm.getClazz());
            this.requiredTest = TakesIndependenceWrapper.class.isAssignableFrom(algorithm.getClazz());
        }
        String description = "";
        try {
			description = AlgorithmDescriptions.getInstance().get(algorithm.getAnnotation().command());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		}
        this.description = description;
    }

    @Override
    public int compareTo(AlgorithmModel other) {
    	String name = "";
    	try {
			name = algorithm.getAnnotation().name();
		} catch (Exception e) {
			Annotation[] annotations = algorithm.getClazz().getDeclaredAnnotations();
			if(annotations != null) {
				for(int i=0;i<annotations.length;i++) {
					Annotation annotation = annotations[i];
					String annotationType = annotation.toString();
					
					if(annotationType.contains("edu.cmu.tetrad.annotation.Algorithm")) {
						int left_paraphase = annotationType.indexOf("(");
						if(left_paraphase > -1) {
							annotationType = annotationType.substring(left_paraphase + 1);
						}
						int right_paraphase = annotationType.indexOf(")");
						if(right_paraphase > -1) {
							annotationType = annotationType.substring(0, right_paraphase);
						}
						String[] tags = annotationType.split(",");
						if(tags != null) {
							for(int j=0;j<tags.length;j++) {
								String[] tokens = tags[j].split("=");
								String tag_name = tokens[0].trim();
								String value = tokens[1].trim();
								
								if(tag_name.equalsIgnoreCase("name")) {
									name = value;
									break;
								}
							}
						}
						break;
					}
				}
			}
		}
    	String otherName = "";
    	try {
			otherName = other.algorithm.getAnnotation().name();
		} catch (Exception e) {
			Annotation[] annotations = other.algorithm.getClazz().getDeclaredAnnotations();
			if(annotations != null) {
				for(int i=0;i<annotations.length;i++) {
					Annotation annotation = annotations[i];
					String annotationType = annotation.toString();
					
					if(annotationType.contains("edu.cmu.tetrad.annotation.Algorithm")) {
						int left_paraphase = annotationType.indexOf("(");
						if(left_paraphase > -1) {
							annotationType = annotationType.substring(left_paraphase + 1);
						}
						int right_paraphase = annotationType.indexOf(")");
						if(right_paraphase > -1) {
							annotationType = annotationType.substring(0, right_paraphase);
						}
						String[] tags = annotationType.split(",");
						if(tags != null) {
							for(int j=0;j<tags.length;j++) {
								String[] tokens = tags[j].split("=");
								String tag_name = tokens[0].trim();
								String value = tokens[1].trim();
								
								if(tag_name.equalsIgnoreCase("name")) {
									name = value;
									break;
								}
							}
						}
						break;
					}
				}
			}
		} 
        return name.compareTo(otherName);
    }

    @Override
    public String toString() {
    	String name = "";
    	try {
			name = algorithm.getAnnotation().name();
		} catch (Exception e) {
			Annotation[] annotations = algorithm.getClazz().getDeclaredAnnotations();
			if(annotations != null) {
				for(int i=0;i<annotations.length;i++) {
					Annotation annotation = annotations[i];
					String annotationType = annotation.toString();
					
					if(annotationType.contains("edu.cmu.tetrad.annotation.Algorithm")) {
						int left_paraphase = annotationType.indexOf("(");
						if(left_paraphase > -1) {
							annotationType = annotationType.substring(left_paraphase + 1);
						}
						int right_paraphase = annotationType.indexOf(")");
						if(right_paraphase > -1) {
							annotationType = annotationType.substring(0, right_paraphase);
						}
						String[] tags = annotationType.split(",");
						if(tags != null) {
							for(int j=0;j<tags.length;j++) {
								String[] tokens = tags[j].split("=");
								String tag_name = tokens[0].trim();
								String value = tokens[1].trim();
								
								if(tag_name.equalsIgnoreCase("name")) {
									name = value;
									break;
								}
							}
						}
						break;
					}
				}
			}
		}
        return name;
    }

    public AnnotatedClass<Algorithm> getAlgorithm() {
        return algorithm;
    }

    public boolean isRequiredScore() {
        return requiredScore;
    }

    public boolean isRequiredTest() {
        return requiredTest;
    }

    public String getDescription() {
        return description;
    }

}
