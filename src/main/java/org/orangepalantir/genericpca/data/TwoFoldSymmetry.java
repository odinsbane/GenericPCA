package org.orangepalantir.genericpca.data;

import org.orangepalantir.genericpca.SimpleReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Input shapes are 1D arrays representing square images. They will be quartered and averaged, as though there is
 * 2 fold symmetry.
 *
 * Created on 11/10/17.
 */
public class TwoFoldSymmetry {
    public static List<double[]> symmetrize(double[] data){
        int l = (int)Math.sqrt(data.length);
        int n = data.length/4;
        int l2 = l/2;
        System.out.printf("From %dx%d to %dx%d\n", l, l, l2, l2);
        double[] ave0 = new double[n];
        double[] ave1 = new double[n];
        double[] ave2 = new double[n];
        double[] ave3 = new double[n];

        for(int i = 0; i<l2; i++){
            for(int j = 0; j<l2; j++){
                ave0[j + i*l2] = data[j + i*l];
                ave1[j + i*l2] = data[j + (l - i - 1)*l];
                ave2[j + i*l2] = data[l-j-1 + i*l];
                ave3[j + i*l2] = data[l-j-1 + (l - i -1)*l];
            }
        }
        return Arrays.asList(ave0, ave1, ave2, ave3);
    }

    public static void writeOut(List<double[]> data, Path path){
        try(BufferedWriter construct = Files.newBufferedWriter(path, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
        ) {
            for(double[] vector: data){

                for(int i = 0; i<vector.length; i++){
                    construct.write(Double.toHexString(vector[i]));
                    if(i<vector.length-1){
                        construct.write('\t');
                    } else{
                        construct.write('\n');
                    }
                }
            }
        } catch (IOException e) {
        }
    }
    public static void main(String[] args) throws IOException {
        Path input = Paths.get(args[0]);
        Path output = Paths.get(args[1]);
        List<double[]> inputs = SimpleReader.loadData(input);
        List<double[]> outputs = inputs.stream().map(TwoFoldSymmetry::symmetrize).flatMap(List::stream).collect(Collectors.toList());
        writeOut(outputs, output);

    }
}
