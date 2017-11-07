package org.orangepalantir.genericpca;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by msmith on 01.08.17.
 */
public class IndexedCoefficient{
    final int i;
    final double m;
    public IndexedCoefficient(int index,double coefficient){
        this.i = index;
        this.m = coefficient;
    }

    public double getMagnitude(){
        return m*m;
    }
    public double getCoefficient(){
        return m;
    }
    public static List<List<IndexedCoefficient>> readCoefficients(Path path) throws IOException {
        List<List<IndexedCoefficient>> in = new ArrayList<>();
        try(BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)){
            String s;
            List<IndexedCoefficient> current = new ArrayList<>();
            while((s=reader.readLine())!=null){
                if(s.length()==0){
                    in.add(current);
                    current = new ArrayList<>();
                    continue;
                }
                if(s.charAt(0)=='#'){
                    continue;
                }
                String[] tokens = s.split("\\s");
                int index = Integer.parseInt(tokens[0]);
                double value = Double.parseDouble(tokens[1]);
                current.add(new IndexedCoefficient(index, value));
            }
            if(current.size()>0) {
                in.add(current);
            }
        }
        return in;
    }
}
