package com.sachkomaxim.dcs.lab1add;

import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.locks.ReentrantLock;

public class SharedWorker implements Runnable {
    private final int N;
    private final int start, end;

    private final SharedContent ctn;

    private final CyclicBarrier barrier;

    private final int threadNumber;

    SharedWorker(int N, int start, int end, CyclicBarrier barrier, int threadNumber, SharedContent ctn) {
        this.N = N;
        this.start = start;
        this.end = end;
        this.threadNumber = threadNumber;

        this.barrier = barrier;

        this.ctn = ctn;
    }

    @Override
    public void run() {
        System.out.printf("Thread %d has started its computation%n", threadNumber);

        calculateAResultPart();
        calculateMGResultPart();

        System.out.printf("Thread %d has finished its computation%n", threadNumber);
    }

    private void calculateAResultPart() {
        partialMultiplyVectorAndMatrix(ctn.B, ctn.MC, ctn.Q, ctn.QLock);
        partialMultiplyVectorAndMatrix(ctn.D, ctn.MZ, ctn.R, ctn.RLock);
        partialMultiplyVectorAndMatrix(ctn.E, ctn.MT, ctn.S, ctn.SLock);
        partialAddOfThreeOrMoreVectors(ctn.Q, ctn.R, ctn.S);
    }

    private void calculateMGResultPart() {
        double wLocal = partialFindMinWlocalInVector(ctn.D, ctn.B);

        ctn.wLock.lock();
        try {
            ctn.w = Math.min(ctn.w, wLocal);
        }
        finally {
            ctn.wLock.unlock();
        }

        try {
            barrier.await();
        } catch (InterruptedException | BrokenBarrierException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }

        partialCalculateMCMTw(ctn.w, ctn.MC, ctn.MT);
        partialMultiplyTwoMatrices(ctn.MZ, ctn.ME);
        partialSubtractOfTwoMatrices(ctn.MCMTw, ctn.MZME);
    }

    private void partialAddOfThreeOrMoreVectors(double[] vector1, double[]... vectors) {
        for (int i = start; i < end; i++) {
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
            ctn.ALock.lock();
            try{
                ctn.A[i] = sum;
            }
            finally {
                ctn.ALock.unlock();
            }
        }
    }

    private void partialSubtractOfTwoMatrices(double[][] matrix1, double[][] matrix2) {
        if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length) {
            throw new IllegalArgumentException("Matrices must have the same length");
        }

        for (int i = start; i < end; i++) {
            for (int j = 0; j < N; j++) {
                ctn.MGLock.lock();
                try {
                    ctn.MG[i][j] = matrix1[i][j] - matrix2[i][j];
                }
                finally {
                    ctn.MGLock.unlock();
                }
            }
        }
    }

    private void partialMultiplyVectorAndMatrix(double[] vector, double[][] matrix, double[] result, ReentrantLock MultVectAndMtrxLock) {
        if (vector.length != matrix.length) {
            throw new IllegalArgumentException("Matrix must have the same number of rows as length of Vector");
        }

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
            MultVectAndMtrxLock.lock();
            try {
                result[i] = sum;
            }
            finally {
                MultVectAndMtrxLock.unlock();
            }
        }
    }

    private void partialMultiplyTwoMatrices(double[][] matrix1, double[][] matrix2) {
        if (matrix1[0].length != matrix2.length) {
            throw new IllegalArgumentException("Matrix2 must have the same number of rows as number of columns in Matrix1");
        }

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
                ctn.MZMELock.lock();
                try {
                    ctn.MZME[i][j] = sum;
                }
                finally {
                    ctn.MZMELock.unlock();
                }
            }
        }
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

    private void partialCalculateMCMTw(double w, double[][] MC, double[][] MT) {
        if (MC[0].length != MT.length) {
            throw new IllegalArgumentException("MC must have the same number of rows as number of columns in MT");
        }

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
                ctn.MCMTwLock.lock();
                try {
                    ctn.MCMTw[i][j] = w * sum;
                }
                finally {
                    ctn.MCMTwLock.unlock();
                }
            }
        }
    }
}
