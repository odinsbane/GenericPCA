package org.orangepalantir.genericpca;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Reads in a text file where each line represents a double[], trains and performs an analysis dump.
 *
 * Created on 8/2/17.
 */
public class SimpleReader {
    public static void main(String[] args) throws IOException {
            String s;
        List<double[]> inputs =loadData(Paths.get(args[0]));
            Trainer trainer = new Trainer(inputs);
            trainer.calculateEigenVectors();
            new AnalysisDump(trainer).dump();
            System.out.println(inputs.size());

    }

    static List<double[]> loadData(Path input) throws IOException {
        List<double[]> inputs = new ArrayList<>();
        try (BufferedReader reads = Files.newBufferedReader(input)) {
            String s;
            while((s=reads.readLine())!=null){
                if(s.length()==0||s.charAt(0)=='#') continue;
                double[] vector = Arrays.stream(s.split("\t")).mapToDouble(Double::valueOf).toArray();
                double[] vector2 = new double[vector.length/2];
                for(int i = 0; i<vector2.length; i++){
                    vector2[i] = vector[2*i];
                }
                vector = vector2;
                inputs.add(vector);
            }

        }
        return inputs;
    }

    public void dumpCoefficients(Trainer trainer, List<double[]> samples, Path out){

    }
}
