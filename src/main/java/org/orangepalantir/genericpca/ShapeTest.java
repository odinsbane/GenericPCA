package org.orangepalantir.genericpca;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

/**
 *
 * For creating PC of 2D shapes.
 *
 * Created by msmith on 31.07.17.
 */
public class ShapeTest {
    static int count = 25;
    static Random random = new Random(1);
    public static void main(String[] args) throws IOException {
        List<double[]> shapes = new ArrayList<>();

        for(int i =0; i<10; i++){
            double[] start = circle();
            shiftX(0.1, start);
            noise(0.2, start);
            shapes.add(start);
        }

        for(int i =0; i<10; i++){
            double[] start = circle();
            shiftX(-0.1, start);
            noise(0.2, start);
            shapes.add(start);
        }

        for(int i = 0; i<10; i++){
            double[] start = circle();
            radial(4, 0.1, start);
            noise(0.2, start);
            shapes.add(start);
        }

        Trainer trainer = new Trainer(shapes);
        trainer.calculateEigenVectors();

        new AnalysisDump(trainer).dump();


    }

    /**
     * a + mag*b
     * @param a
     * @param b
     * @param mag
     * @return
     */
    static void add(double[] a, double magA,  double[] b, double magB, double[] target){
        for(int i = 0; i<a.length; i++){
            target[i] = a[i]*magA + b[i]*magB;
        }
    }

    static void noise(double M, double[] start){
        for(int i = 0;i<start.length; i++){
            start[i] += random.nextGaussian()*M;
        }
    }

    static void radial(int N, double M, double[] start){

        double dtheta = 2*N*Math.PI/count;
        for(int i = 0; i<count; i++){
            start[2*i] += M*Math.sin(dtheta*i);
            start[2*i+1] += M*Math.cos(dtheta*i);
        }
    }

    static void shiftX(double M, double[] start){
        double dtheta = 2*Math.PI/count;
        for(int i = 0; i<count; i++){
            start[2*i] += M*Math.sin(dtheta*i);
        }
    }

    static double[] circle(){
        double dtheta = 2*Math.PI/count;

        double[] val = new double[2*count];

        for(int i = 0;i<count; i++){
            val[2*i] = Math.sin(dtheta*i);
            val[2*i + 1] = Math.cos(dtheta*i);
        }
        return val;
    }

}
