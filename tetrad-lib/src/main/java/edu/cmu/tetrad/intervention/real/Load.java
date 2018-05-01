package edu.cmu.tetrad.intervention.real;

import edu.cmu.tetrad.data.DataReader;
import edu.cmu.tetrad.data.DataSet;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by bandrews on 9/14/17.
 */

public class Load {

    public List<DataSet> load_data(String path, int max_discrete) {

        path = new File(path).getAbsolutePath();
        List<DataSet> dataSets = new ArrayList<>();

        try {
            File data_file = new File(path + "/data_interventions_continuous.txt");
            File meta_file = new File(path + "/meta_interventions.txt");
            DataReader reader = new DataReader();
            reader.setVariablesSupplied(true);
            reader.setMaxIntegralDiscrete(max_discrete);
            dataSets.add(reader.parseTabular(data_file));
            dataSets.add(reader.parseTabular(meta_file));
        } catch(IOException e) {
            e.printStackTrace();
        }

        return dataSets;

    }
}