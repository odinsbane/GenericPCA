package org.orangepalantir.genericpca;

/**
 * Created by msmith on 01.08.17.
 */
class IndexedCoefficient{
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
}
