package com.sachkomaxim.dcs.lab1add;

import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;

public class LocalizedWorker implements Runnable {
    private final int N;
    private final int start, end, localRange;

    private final double[][] localMC, localMZ, localMT, localME;
    private final double[] localB, localD, localE;
    private final double[] wLocalResults;

    private double[] localAPart;
    private double[][] localMGPart;

    private final CyclicBarrier barrier;

    private final int threadNumber;

    LocalizedWorker(int N, int start, int end, Data data, double[] wLocalResults,
                    CyclicBarrier barrier, int threadNumber) {
        this.N = N;
        this.start = start;
        this.end = end;
        this.threadNumber = threadNumber;
        this.localRange = end - start;

        this.localB = deepCopyVector(data.B);
        this.localD = deepCopyVector(data.D);
        this.localE = deepCopyVector(data.E);

        this.localMC = deepCopyMatrix(data.MC);
        this.localMZ = deepCopyMatrix(data.MZ);
        this.localMT = deepCopyMatrix(data.MT);
        this.localME = deepCopyMatrix(data.ME);

        this.wLocalResults = wLocalResults;
        this.barrier = barrier;

        this.localAPart = new double[localRange];
        this.localMGPart = new double[localRange][N];
    }

    @Override
    public void run() {
        System.out.printf("Thread %d has started its computation%n", threadNumber);

        localAPart = calculateLocalAPart();
        localMGPart = calculateLocalMGPart();

        System.out.printf("Thread %d has finished its computation%n", threadNumber);
    }

    private double[] calculateLocalAPart() {
        double[] QPart = partialMultiplyVectorAndMatrix(localB, localMC);
        double[] RPart = partialMultiplyVectorAndMatrix(localD, localMZ);
        double[] SPart = partialMultiplyVectorAndMatrix(localE, localMT);
        return partialAddOfThreeOrMoreVectors(QPart, RPart, SPart);
    }

    private double[][] calculateLocalMGPart() {
        wLocalResults[threadNumber - 1] = partialFindMinWlocalInVector(localD, localB);

        try {
            barrier.await();
        } catch (InterruptedException | BrokenBarrierException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }

        double wGlobal = Double.MAX_VALUE;
        for (double wLocalResult : wLocalResults) {
            wGlobal = Math.min(wGlobal, wLocalResult);
        }

        double[][] MCMTwPart = partialCalculateMCMTw(wGlobal, localMC, localMT);
        double[][] MZMEPart = partialMultiplyTwoMatrices(localMZ, localME);
        return partialSubtractOfTwoMatrices(MCMTwPart, MZMEPart);
    }

    private double[] partialAddOfThreeOrMoreVectors(double[] vector1, double[]... vectors) {
        double[] result = new double[localRange];

        for (int i = 0; i < localRange; i++) {
            double sum = vector1[i];
            double c = 0.0;
            for (double[] vector : vectors) {
                if (vector1.length != vector.length) {
                    throw new IllegalArgumentException("Vectors must have the same length");
                }

                double y = vector[i] - c;
                double t = sum + y;
                c = (t - sum) - y;
                sum = t;
            }
            result[i] = sum;
        }
        return result;
    }

    private double[][] partialSubtractOfTwoMatrices(double[][] matrix1, double[][] matrix2) {
        if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length) {
            throw new IllegalArgumentException("Matrices must have the same length");
        }

        double[][] result = new double[localRange][N];

        for (int i = 0; i < localRange; i++) {
            for (int j = 0; j < N; j++) {
                result[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }
        return result;
    }

    private double[] partialMultiplyVectorAndMatrix(double[] vector, double[][] matrix) {
        if (vector.length != matrix.length) {
            throw new IllegalArgumentException("Matrix must have the same number of rows as length of Vector");
        }

        double[] result = new double[localRange];

        for (int i = start; i < end; i++) {
            double sum = 0.0;
            double c = 0.0;
            for (int j = 0; j < N; j++) {
                double term = vector[j] * matrix[j][i];
                double y = term - c;
                double t = sum + y;
                c = (t - sum) - y;
                sum = t;
            }
            result[i - start] = sum;
        }
        return result;
    }

    private double[][] partialMultiplyTwoMatrices(double[][] matrix1, double[][] matrix2) {
        if (matrix1[0].length != matrix2.length) {
            throw new IllegalArgumentException("Matrix2 must have the same number of rows as number of columns in Matrix1");
        }

        double[][] result = new double[localRange][N];

        for (int i = start; i < end; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                double c = 0.0;
                for (int k = 0; k < N; k++) { // Спільна розмірність
                    double term = matrix1[i][k] * matrix2[k][j];
                    double y = term - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                result[i - start][j] = sum;
            }
        }
        return result;
    }

    public double partialFindMinWlocalInVector (double[] vector1, double[] vector2) {
        if (vector1.length != vector2.length) {
            throw new IllegalArgumentException("Vectors must have the same length");
        }

        double min = Double.MAX_VALUE;

        for (int i = start; i < end; i++) {
            double sum = vector1[i] + vector2[i];

            min = Math.min(min, sum);
        }

        return min;
    }

    private double[][] partialCalculateMCMTw(double w, double[][] MC, double[][] MT) {
        if (MC[0].length != MT.length) {
            throw new IllegalArgumentException("MC must have the same number of rows as number of columns in MT");
        }

        double[][] result = new double[localRange][N];

        for (int i = start; i < end; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                double c = 0.0;
                for (int k = 0; k < N; k++) {
                    double term = MC[i][k] * MT[k][j];
                    double y = term - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                result[i - start][j] = w * sum;
            }
        }
        return result;
    }

    private static double[][] deepCopyMatrix(double[][] src) {
        int n = src.length;
        double[][] dst = new double[n][];
        for (int i = 0; i < n; i++) {
            dst[i] = new double[src[i].length];
            System.arraycopy(src[i], 0, dst[i], 0, src[i].length);
        }
        return dst;
    }

    private static double[] deepCopyVector(double[] src) {
        double[] dst = new double[src.length];
        System.arraycopy(src, 0, dst, 0, src.length);
        return dst;
    }

    public double[] getLocalAResultPart() {
        return localAPart;
    }

    public double[][] getLocalMGResultPart() {
        return localMGPart;
    }
}
