package org.orangepalantir.genericpca;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Generic method for taking a collection of vectors (double[]) and finding the principle components.
 *
 * Created on 27.07.17.
 */
public class Trainer {
    List<double[]> featurePoints = new ArrayList<>();
    double[] average;
    final int N;

    EigenvalueDecomposition eigens;
    List<double[]> eigenVectors = new ArrayList<>();

    public Trainer(Collection<double[]> featurePoints){
        this.featurePoints.addAll(featurePoints);
        N = this.featurePoints.get(0).length;
        average = new double[N];
    }

    public void calculateEigenVectors(){
        calculateAverage();
        List<double[]> deltas = calculateDeviation();
        double[][] cov = getCovarianceMatrix(deltas);
        eigens = new Matrix(cov).eig();
        double[][] rawEigenVectors = eigens.getV().getArray();

        for(int v = 0; v<N; v++){
            double[] vector = new double[N];
            double sum = 0;
            for(int i = 0; i<N; i++){
                double value = rawEigenVectors[i][v];
                vector[i] = value;
                sum += value*value;
            }
            sum = 1.0/Math.sqrt(sum);
            for(int i = 0; i<N; i++){
                vector[i] = vector[i]*sum;
            }

            eigenVectors.add(vector);
        }
    }

    public double[][] getCovarianceMatrix(List<double[]> deltas){
        double[][] covariance = new double[N][N];
        for(double[] delta: deltas){
            for(int i = 0; i<N; i++){
                for(int j = 0; j<N; j++){
                    covariance[i][j] += delta[i]*delta[j];
                }
            }
        }
        return covariance;
    }
    public List<double[]> calculateDeviation(){
        List<double[]> deviations = new ArrayList<>();
        for(double[] d: featurePoints){
            double[] delta = new double[N];
            for(int i = 0; i<N; i++){
                delta[i] = d[i] - average[i];
            }
            deviations.add(delta);
        }

        return deviations;
    }

    public double[] calculatDelta(double[] vector){
        double[] delta = new double[N];
        for(int i = 0; i<N; i++){
            delta[i] = vector[i] - average[i];
        }
        return delta;
    }

    public void calculateAverage(){
        for(int i = 0; i<average.length; i++) average[i] = 0;
        for(double[] vector: featurePoints){
            for(int i = 0; i<N; i++){
                average[i] += vector[i];
            }
        }
        double f = 1.0/featurePoints.size();
        for(int i = 0; i<average.length; i++) average[i] = average[i]*f;

    }

    List<IndexedCoefficient> getCoefficients(double[] vector){
        List<IndexedCoefficient> coefficients = new ArrayList<>(N);
        double[] delta = calculatDelta(vector);
        for(int i = 0; i<eigenVectors.size(); i++){
            double dot = 0;
            double[] ev = eigenVectors.get(i);
            for(int j =0; j<delta.length; j++){

                dot += ev[j]*delta[j];

            }
            coefficients.add(new IndexedCoefficient(i, dot));

        }

        return coefficients;
    }

    public double[] getVector(int index){
        return eigenVectors.get(index);
    }


    /**
     * Sums two vectors. places the result in target.
     *
     *    a[i]*magA + b[i]*magB
     *
     * @param a
     * @param magA
     * @param b
     * @param magB
     * @param target
     */
    static void add(double[] a, double magA,  double[] b, double magB, double[] target){
        for(int i = 0; i<a.length; i++){
            target[i] = a[i]*magA + b[i]*magB;
        }
    }

}

