package org.orangepalantir.genericpca;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.EventQueue;
import java.awt.FileDialog;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by melkor on 10/28/17.
 */
public class DisplayInverseFFT {
    int length;

    public void displayData(List<double[]> vectors){
        length = (int)Math.sqrt(vectors.get(0).length);
        ImageStack original = new ImageStack(length, length);
        ImageStack inverse = new ImageStack(length, length);
        for(double[] vector: vectors){
            ImageProcessor proc = new FloatProcessor(length, length);
            for(int i = 0; i<length; i++){
                for(int j = 0; j<length; j++){
                    proc.setf(i, j, (float)vector[i*length + j]);
                }
            }
            original.addSlice(proc);
            FHT fht = new FHT(proc, false);
            fht.inverseTransform();
            inverse.addSlice(fht);
        }
        ImageJ.main(new String[]{});
        new ImagePlus("originals", original).show();
        new ImagePlus("inverse", inverse).show();
    }
    public static void main(String[] args){
        DisplayInverseFFT display = new DisplayInverseFFT();
        FileDialog log = new FileDialog((JFrame)null, "Choose vector data");
        log.setMode(FileDialog.LOAD);
        log.setVisible(true);
        String f = log.getFile();
        if(f==null) return;
        String d = log.getDirectory();
        try {

            display.displayData(SimpleReader.loadData(Paths.get(d,f)));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
