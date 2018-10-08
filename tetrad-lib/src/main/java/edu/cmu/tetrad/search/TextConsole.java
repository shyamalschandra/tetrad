//package edu.cmu.tetrad.search;
//
//import java.awt.FileDialog;
//import java.awt.Frame;
//import java.io.BufferedReader;
//import java.io.InputStreamReader;
//import org.rosuda.JRI.RMainLoopCallbacks;
//import org.rosuda.JRI.Rengine;
//
//public class TextConsole implements RMainLoopCallbacks {
//    public TextConsole() {
//    }
//
//    public void rWriteConsole(Rengine var1, String var2, int var3) {
//        System.out.print(var2);
//    }
//
//    public void rBusy(Rengine var1, int var2) {
//        System.out.println("rBusy(" + var2 + ")");
//    }
//
//    public String rReadConsole(Rengine var1, String var2, int var3) {
//        System.out.print(var2);
//
//        try {
//            BufferedReader var4 = new BufferedReader(new InputStreamReader(System.in));
//            String var5 = var4.readLine();
//            return var5 != null && var5.length() != 0 ? var5 + "\n" : var5;
//        } catch (Exception var6) {
//            System.out.println("jriReadConsole exception: " + var6.getMessage());
//            return null;
//        }
//    }
//
//    public void rShowMessage(Rengine var1, String var2) {
//        System.out.println("rShowMessage \"" + var2 + "\"");
//    }
//
//    public String rChooseFile(Rengine var1, int var2) {
//        FileDialog var3 = new FileDialog(new Frame(), var2 == 0 ? "Select a file" : "Select a new file", var2 == 0 ? 0 : 1);
//        var3.show();
//        String var4 = null;
//        if (var3.getDirectory() != null) {
//            var4 = var3.getDirectory();
//        }
//
//        if (var3.getFile() != null) {
//            var4 = var4 == null ? var3.getFile() : var4 + var3.getFile();
//        }
//
//        return var4;
//    }
//
//    public void rFlushConsole(Rengine var1) {
//    }
//
//    public void rLoadHistory(Rengine var1, String var2) {
//    }
//
//    public void rSaveHistory(Rengine var1, String var2) {
//    }
//}
