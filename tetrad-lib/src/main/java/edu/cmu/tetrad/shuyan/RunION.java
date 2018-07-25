package edu.cmu.tetrad.shuyan;


import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.IonJoeModifications;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by user on 3/26/18.
 */
public class RunION {

    public RunION(){

    }

    public List<Graph> run(List<Graph> pags){

        IonJoeModifications2 ion = new IonJoeModifications2(pags);
        return ion.search();
    }

    public static void main(String...args){
        Graph pag1 = GraphUtils.loadGraphTxt(new File("C:\\Users\\wslxx\\OneDrive\\Documents\\compareION_Graph\\xya.txt"));
        Graph pag2 = GraphUtils.loadGraphTxt(new File("C:\\Users\\wslxx\\OneDrive\\Documents\\compareION_Graph\\xyb.txt"));
        Graph pag3 = GraphUtils.loadGraphTxt(new File("C:\\Users\\wslxx\\OneDrive\\Documents\\compareION_Graph\\xyc.txt"));

        ArrayList<Graph> pags = new ArrayList<>();
        pags.add(pag1);
        pags.add(pag2);
        pags.add(pag3);

        RunION pedro = new RunION();

        List<Graph> results = pedro.run(pags);

        for(int i = 0; i < results.size(); i++){
            System.out.println(results.get(i).toString());
            //System.out.println("finish");
            System.out.println();
            System.out.println();
            if (i+1 == results.size()) System.out.println(results.size() + " possible graphs");
        }
    }
}
