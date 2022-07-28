using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

public class DensityDistribution : MonoBehaviour
{
    // default is linear
    private bool stepwise = false;

    private float[] proportions;
    public float[] m = {0,0}, c = {0,0};
    public float midpoint;
    
    // sets mode to quadratic and sets proportions
    public void set_stepwise(float[] proportions) {
        this.stepwise = true;
        this.proportions = proportions;
    }

    /*
    *   Finds m and c for a given set of points, and saves to private attributes based on i
    */
    private void solve(float x1, float y1, float x2, float y2, int i) {
        // calculate first series of equations
        Matrix<double> X = DenseMatrix.OfArray(new double[,] { {x1, 1}, {x2, 1} });
        Matrix<double> y = DenseMatrix.OfArray(new double[,] { {y1}, {y2} });

        Matrix<double> sol = X.Inverse() * y;
        m[i] = (float) sol[0,0]; c[i] = (float) sol[1,0]; 
    }

    // calculate coefficients
    public void calculate(float xmin, float xmax) {
        if (stepwise) {
            float[] pos = {proportions[0] * (xmax - xmin) + xmin, 
                       proportions[1] * (xmax - xmin) + xmin};

            // save midpoint
            midpoint = pos[0];

            // solve
            solve(xmin, xmin, pos[0], pos[1], 0);
            solve(pos[0], pos[1], xmax, xmax, 1);
        }
    }

    // convert
    public float convert(float rand) {
        if (!stepwise) return rand;
        else {
            int ind = (rand <= midpoint) ? 0 : 1;
            return m[ind] * rand + c[ind];
        } 
    }
}