Generica Principle Component Analysis in Java.
=====================

Java PCA library
-----------------------------------

This is an extremely small library for performing PCA.

Example
-------

Start with a vector.

    double[] v = {x0, y0, x1, y1, x2, y2, ...};

We have a list of vectors, and we want to perform pca we create the covariance matrix.

   double[][] A;

Each element is defined as

   A_ij = Sum_n[ (vn_i * vn_j) ]

Then we find the eigen vectors and express the input vectors as sums of the eigen vectors.