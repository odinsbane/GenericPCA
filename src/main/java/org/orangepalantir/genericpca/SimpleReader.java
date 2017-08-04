package org.orangepalantir.genericpca;

import org.orangepalantir.genericpca.display.TwoDHeatMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Reads in a text file where each line represents a double[], trains and performs an analysis dump.
 *
 * Created on 8/2/17.
 */
public class SimpleReader {
    public static void main(String[] args) throws IOException {
        List<double[]> inputs =loadSampledData(Paths.get(args[0]));
        Trainer trainer = new Trainer(inputs);
        trainer.calculateEigenVectors();
        new AnalysisDump(trainer).dump();
        List<double[]> values = new ArrayList<>(trainer.eigenVectors.size() + 1);
        values.add(trainer.average);
        values.addAll(trainer.eigenVectors);

        for(int i = 1; i<args.length; i++){
            Path path = Paths.get(args[i]).toAbsolutePath();
            List<double[]> samples = loadSampledData(path);
            Path dir = path.getParent();
            String name = path.getFileName().toString();
            int j = name.lastIndexOf('.');
            j = j>0?j:name.length();

            String newName = name.substring(0, j) + "-coef.txt";
            System.out.println(name + " analysis is written to: " + newName);
            bestCoefficients(trainer, samples, dir.resolve(newName));


        }



    }

    /**
     * The input files are i00, std00, i01, std01, ... this method skips the std data.
     * @param input
     * @return
     * @throws IOException
     */
    public static List<double[]> loadSampledData(Path input) throws IOException {
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

    public static List<double[]> loadData(Path input) throws IOException {
        List<double[]> inputs = new ArrayList<>();
        try (BufferedReader reads = Files.newBufferedReader(input)) {
            String s;
            while((s=reads.readLine())!=null){
                if(s.length()==0||s.charAt(0)=='#') continue;
                double[] vector = Arrays.stream(s.split("\t")).mapToDouble(Double::valueOf).toArray();
                inputs.add(vector);
            }

        }
        return inputs;
    }

    public static void dumpCoefficients(Trainer trainer, List<double[]> samples, Path out){
        try(BufferedWriter writer = Files.newBufferedWriter(out, StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)){
            for(double[] vector: samples) {
                List<IndexedCoefficient> coefficients = trainer.getCoefficients(vector);

                for (int k = 0; k < coefficients.size(); k++) {
                    IndexedCoefficient ic = coefficients.get(coefficients.size() - k - 1);
                    writer.write(String.format("%f\t%f\n", 1.0 * ic.i, ic.getCoefficient()));
                }
                writer.write("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void bestCoefficients(Trainer trainer, List<double[]> samples, Path out){
        try(BufferedWriter writer = Files.newBufferedWriter(out, StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)){
            double[] average = new double[samples.get(0).length];
            for(double[] vector: samples) {
                List<IndexedCoefficient> coefficients = trainer.getCoefficients(vector);

                for (int k = 0; k < coefficients.size(); k++) {

                    IndexedCoefficient ic = coefficients.get(coefficients.size() - k - 1);
                    average[k] += ic.getCoefficient();
                }

                coefficients.sort(Comparator.comparingDouble(IndexedCoefficient::getMagnitude));
                int tasl = coefficients.size() -1;
                for(int i = 0; i<3; i++){
                    IndexedCoefficient ic = coefficients.get(tasl - i);
                    writer.write(String.format("%d\t%f\n", ic.i, ic.getCoefficient() ));
                }
                writer.write('\n');

            }

            /*
            double n = samples.size();
            for(int i = 0; i<average.length; i++){
                writer.write(String.format("%d\t%f\n", i, average[i]/n));
            }
            writer.write('\n');
            */
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
