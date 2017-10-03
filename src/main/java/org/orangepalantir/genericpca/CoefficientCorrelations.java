package org.orangepalantir.genericpca;

import lightgraph.DataSet;
import lightgraph.Graph;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Created by melkor on 10/3/17.
 */



public class CoefficientCorrelations {
    List<double[]> coefficients = new ArrayList<>();
    public void setInput(List<List<IndexedCoefficient>> input){
        for(List<IndexedCoefficient> coefs: input){
            double[] column = new double[coefs.size()];
            for(IndexedCoefficient coef: coefs){
                column[coef.i] = coef.getCoefficient();
            }
            coefficients.add(column);
        }
    }

    public void showCorrelations(int terms){
        double[] aves = new double[terms];
        double[] sqs = new double[terms];
        int top = coefficients.get(0).length-1;

        for(double[] set: coefficients){
            System.out.println(set.length);
            for(int i = 0; i<terms; i++){
                int dex = top-i;
                aves[i] += set[dex];
                sqs[i] += set[dex]*set[dex];
            }
        }
        int n = coefficients.size();
        for(int i = 0; i<terms; i++){
            aves[i] = aves[i]/n;
            sqs[i] = Math.sqrt(sqs[i]/n - aves[i]*aves[i]);
        }

        List<double[]> correlations = new ArrayList<>();

        for(int i = 0; i<terms; i++){
            double[] values = new double[terms];
            for(double[] set: coefficients){

                for(int j = 0; j<terms; j++){
                    values[j] += (set[top-i] - aves[i])*(set[top-j] - aves[j]);
                }

            }
            for(int j = 0; j<terms; j++){
                if(i==j){
                    values[j] = 0;
                } else{
                    values[j] = values[j]/(sqs[i]*sqs[j]*n);
                }
            }
            correlations.add(values);
        }

        double[] xaxis = new double[terms];
        for(int i = 0; i<terms; i++){
            xaxis[i] = top - i;
        }

        Graph graph = new Graph();
        int index = 0;
        for(double[] row: correlations){
            DataSet ds = graph.addData(xaxis, row);
            ds.setLabel("" + index);
            index++;
        }

        graph.show(true);

    }

    public static void main(String[] args) throws IOException {
        List<List<IndexedCoefficient>> coefficients = IndexedCoefficient.readCoefficients(Paths.get(args[0]));

        CoefficientCorrelations cc = new CoefficientCorrelations();
        cc.setInput(coefficients);
        cc.showCorrelations(100);
    }
}
