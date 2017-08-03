package org.orangepalantir.genericpca;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Comparator;
import java.util.List;

/**
 * Created by msmith on 01.08.17.
 */
public class AnalysisDump {
    Path con;
    Path coef;
    Path eig;
    int terms = 5;
    Trainer trainer;
    public AnalysisDump(Trainer trainer){
        con = Paths.get("constructs.txt");
        coef = Paths.get("coefficients.txt");
        eig = Paths.get("eigens.txt");
        this.trainer = trainer;
    }

    public void dump() throws IOException {
        BufferedWriter eigens = Files.newBufferedWriter(eig, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
        BufferedWriter construct = Files.newBufferedWriter(con, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
        BufferedWriter constants = Files.newBufferedWriter(coef, StandardCharsets.UTF_8, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);

        double[] vector = trainer.average;
        List<double[]> shapes = trainer.featurePoints;



        for(int j = 0; j<shapes.size(); j++){
            double[] input = shapes.get(j);
            List<IndexedCoefficient> coefficients = trainer.getCoefficients(input);

            for(int k = 0; k<coefficients.size(); k++){
                IndexedCoefficient ic = coefficients.get(coefficients.size() - k - 1);
                constants.write(String.format("%f\t%f\n", 1.0*ic.i, ic.getCoefficient()));
            }

            coefficients.sort(Comparator.comparingDouble(IndexedCoefficient::getMagnitude));
            double[] v   = new double[vector.length];

            for(int k = 0; k<coefficients.size(); k++){
                IndexedCoefficient ic = coefficients.get(coefficients.size() - k - 1);

                if(k<terms) {
                    double[] ev = trainer.getVector(ic.i);

                    Trainer.add(v, 1, ev, ic.getCoefficient(), v);
                }
            }
            constants.write("\n");
            for(int i = 0; i<v.length/2; i++){
                construct.write(
                        String.format("%f\t%f\t%f\t%f\t%f\t%f\n",
                                vector[2*i], vector[2*i+1],
                                v[2*i], v[2*i+1],
                                input[2*i], input[2*i+1]));
            }
            construct.write("\n");

        }

        for(int j = 0; j<vector.length; j++){

            double[] ev = trainer.getVector(j);
            for(int i = 0; i<ev.length; i++){
                eigens.write(String.format("%f", ev[i]));
                if(i<ev.length-1){
                    eigens.write('\t');
                } else{
                    eigens.write('\n');
                }
            }
        }

        eigens.close();
        construct.close();
        constants.close();
    }
}
