package edu.cmu.tetrad.search;

import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

public class TextConsole implements RMainLoopCallbacks {
    @Override
    public void rWriteConsole(Rengine rengine, String s, int i) {

    }

    @Override
    public void rBusy(Rengine rengine, int i) {

    }

    @Override
    public String rReadConsole(Rengine rengine, String s, int i) {
        return null;
    }

    @Override
    public void rShowMessage(Rengine rengine, String s) {

    }

    @Override
    public String rChooseFile(Rengine rengine, int i) {
        return null;
    }

    @Override
    public void rFlushConsole(Rengine rengine) {

    }

    @Override
    public void rSaveHistory(Rengine rengine, String s) {

    }

    @Override
    public void rLoadHistory(Rengine rengine, String s) {

    }
}
