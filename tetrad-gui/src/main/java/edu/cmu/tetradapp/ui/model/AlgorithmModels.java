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

import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.annotation.AlgorithmAnnotations;
import edu.cmu.tetrad.data.DataType;
import java.util.Collections;
import java.util.EnumMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 *
 * Nov 30, 2017 4:20:43 PM
 *
 * @author Kevin V. Bui (kvb2@pitt.edu)
 */
public class AlgorithmModels {

    private static final AlgorithmModels INSTANCE = new AlgorithmModels();

    private final List<AlgorithmModel> models;
    private final Map<AlgType, List<AlgorithmModel>> modelMap;

    private AlgorithmModels() {
        AlgorithmAnnotations algoAnno = AlgorithmAnnotations.getInstance();
        List<AlgorithmModel> list = algoAnno.filterOutExperimental(algoAnno.getAnnotatedClasses()).stream()
                .map(e -> new AlgorithmModel(e))
                .sorted()
                .collect(Collectors.toList());
        models = Collections.unmodifiableList(list);

        Map<AlgType, List<AlgorithmModel>> map = new EnumMap<>(AlgType.class);

        // initialize enum map
        for (AlgType algType : AlgType.values()) {
            map.put(algType, new LinkedList<>());
        }

        // group by datatype
        models.stream().forEach(e -> {
            map.get(e.getAlgorithm().getAnnotation().algoType()).add(e);
        });

        // make it unmodifiable
        map.forEach((k, v) -> map.put(k, Collections.unmodifiableList(v)));
        modelMap = Collections.unmodifiableMap(map);
    }

    public static AlgorithmModels getInstance() {
        return INSTANCE;
    }

    private List<AlgorithmModel> filterInclusivelyByAllOrSpecificDataType(List<AlgorithmModel> algorithmModels, DataType dataType, boolean multiDataSetAlgorithm) {
        AlgorithmAnnotations algoAnno = AlgorithmAnnotations.getInstance();

        return (dataType == DataType.All)
                ? algorithmModels
                : algorithmModels.stream()
                        .filter(e -> {
                            return multiDataSetAlgorithm
                                    ? algoAnno.acceptMultipleDataset(e.getAlgorithm().getClazz())
                                    : true;
                        })
                        .filter(e -> {
                            for (DataType dt : e.getAlgorithm().getAnnotation().dataType()) {
                                if (dt == DataType.All || dt == dataType) {
                                    return true;
                                }
                            }
                            return false;
                        })
                        .collect(Collectors.toList());
    }

    public List<AlgorithmModel> getModels(DataType dataType, boolean multiDataSetAlgorithm) {
        return filterInclusivelyByAllOrSpecificDataType(models, dataType, multiDataSetAlgorithm);
    }

    public List<AlgorithmModel> getModels(AlgType algType, DataType dataType, boolean multiDataSetAlgorithm) {
        return modelMap.containsKey(algType)
                ? filterInclusivelyByAllOrSpecificDataType(modelMap.get(algType), dataType, multiDataSetAlgorithm)
                : Collections.EMPTY_LIST;
    }

}
