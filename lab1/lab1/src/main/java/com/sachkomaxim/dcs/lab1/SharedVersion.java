package com.sachkomaxim.dcs.lab1;

import java.io.*;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;

public class SharedVersion {
    private static final int N = Data.N;
    private static final int THREADS = 4;

    private static final String OUTPUT_FILE = "output_shared.txt";

    private static double[][] MC, MZ, MT, ME, MG, MCMTw, MZME;
    private static double[] B, D, E, A, Q, R, S;
    private static double w;

    private static final Object wLock = new Object();
    private static final Object ALock = new Object();
    private static final Object MGLock = new Object();
    private static final Object QLock = new Object();
    private static final Object RLock = new Object();
    private static final Object SLock = new Object();
    private static final Object MZMELock = new Object();
    private static final Object MCMTwLock = new Object();

    static class ComputingRunnable implements Runnable {
        private final int start, stop;
        private final CyclicBarrier barrier;
        private final int threadNumber;

        ComputingRunnable(int start, int stop, CyclicBarrier barrier, int threadNumber) {
            this.start = start;
            this.stop = stop;
            this.threadNumber = threadNumber;

            this.barrier = barrier;
        }

        @Override
        public void run() {
            System.out.printf("Thread %d has started its computation%n", threadNumber);

            calculateAResultPart();
            calculateMGResultPart();

            System.out.printf("Thread %d has finished its computation%n", threadNumber);
        }

        private void calculateAResultPart() {
            partialMultiplyVectorAndMatrix(B, MC, Q, QLock);
            partialMultiplyVectorAndMatrix(D, MZ, R, RLock);
            partialMultiplyVectorAndMatrix(E, MT, S, SLock);
            partialAddOfThreeOrMoreVectors(Q, R, S);
        }

        private void calculateMGResultPart() {
            double wLocal = partialFindMinWlocalInVector(D, B);

            synchronized (wLock){
                w = Math.min(w, wLocal);
            }

            try {
                barrier.await();
            } catch (InterruptedException | BrokenBarrierException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException(e);
            }

            partialCalculateMCMTw(w, MC, MT);
            partialMultiplyTwoMatrices(MZ, ME);
            partialSubtractOfTwoMatrices(MCMTw, MZME);
        }

        private void partialAddOfThreeOrMoreVectors(double[] vector1, double[]... vectors) {
            for (int i = start; i < stop; i++) {
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
                synchronized (ALock) {
                    A[i] = sum;
                }
            }
        }

        private void partialSubtractOfTwoMatrices(double[][] matrix1, double[][] matrix2) {
            if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length) {
                throw new IllegalArgumentException("Matrices must have the same length");
            }

            for (int i = start; i < stop; i++) {
                for (int j = 0; j < N; j++) {
                    synchronized (MGLock) {
                        MG[i][j] = matrix1[i][j] - matrix2[i][j];
                    }
                }
            };
        }

        private void partialMultiplyVectorAndMatrix(double[] vector, double[][] matrix, double[] result, Object MultVectAndMtrxLock) {
            if (vector.length != matrix.length) {
                throw new IllegalArgumentException("Matrix must have the same number of rows as length of Vector");
            }

            for (int i = start; i < stop; i++) {
                double sum = 0.0;
                double c = 0.0;
                for (int j = 0; j < N; j++) {
                    double term = vector[j] * matrix[j][i];
                    double y = term - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                synchronized (MultVectAndMtrxLock) {
                    result[i] = sum;
                }
            }
        }

        private void partialMultiplyTwoMatrices(double[][] matrix1, double[][] matrix2) {
            if (matrix1[0].length != matrix2.length) {
                throw new IllegalArgumentException("Matrix2 must have the same number of rows as number of columns in Matrix1");
            }

            for (int i = start; i < stop; i++) {
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
                    synchronized (MZMELock) {
                        MZME[i][j] = sum;
                    }
                }
            }
        }

        public double partialFindMinWlocalInVector (double[] vector1, double[] vector2) {
            if (vector1.length != vector2.length) {
                throw new IllegalArgumentException("Vectors must have the same length");
            }

            double min = Double.MAX_VALUE;

            for (int i = start; i < stop; i++) {
                double sum = vector1[i] + vector2[i];

                min = Math.min(min, sum);
            }

            return min;
        }

        private void partialCalculateMCMTw(double w, double[][] MC, double[][] MT) {
            if (MC[0].length != MT.length) {
                throw new IllegalArgumentException("MC must have the same number of rows as number of columns in MT");
            }

            for (int i = start; i < stop; i++) {
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
                    synchronized (MCMTwLock) {
                        MCMTw[i][j] = w * sum;
                    }
                }
            }
        }
    }

    public static void main(String[] args) {
        try {
            System.out.println("Starting computation with " + THREADS + " threads");

            long startTime = System.currentTimeMillis();
            Data data = loadData();

            initGlobals(data);

            Thread[] threads = new Thread[THREADS];
            ComputingRunnable[] workers = new ComputingRunnable[THREADS];
            CyclicBarrier barrier = new CyclicBarrier(THREADS);
            int rowsPerThread = N / THREADS;
            int startRow = 0;

            for (int i = 0; i < THREADS; i++) {
                int endRow = (i == THREADS - 1) ? N : startRow + rowsPerThread;
                workers[i] = new ComputingRunnable(startRow, endRow, barrier, i + 1);
                threads[i] = new Thread(workers[i]);
                threads[i].start();
                startRow = endRow;
            }

            for (int i = 0; i < THREADS; i++) {
                threads[i].join();
            }

            saveAndPrintResults(A, MG);
            System.out.println("Results calculated and written to " + OUTPUT_FILE);

            double executionTime = (System.currentTimeMillis() - startTime) / 1000.0;
            System.out.printf("%nComputation completed in %.3f seconds%n", executionTime);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void initGlobals(Data data) {
        MC = data.MC;
        MZ = data.MZ;
        MT = data.MT;
        ME = data.ME;
        B = data.B;
        D = data.D;
        E = data.E;

        A = new double[N];
        MG = new double[N][N];

        w = Double.MAX_VALUE;
        Q = new double[N];
        R = new double[N];
        S = new double[N];
        MCMTw = new double[N][N];
        MZME = new double[N][N];
    }

    private static void saveAndPrintResults(double[] A, double[][] MG) {
        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(OUTPUT_FILE)))) {
            if (N <= 10) System.out.println("\nVector A:");
            writer.println("\nVector A:");
            for (int i = 0; i < N; i++) {
                String value = String.format("%f ", A[i]);
                if (N <= 10) System.out.print(value);
                writer.print(value);
            }
            if (N <= 10) System.out.println();
            writer.println();
            writer.println();

            if (N <= 10) System.out.println("\nMatrix MG:");
            writer.println("Matrix MG:");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    String value = String.format("%f ", MG[i][j]);
                    if (N <= 10) System.out.print(value);
                    writer.print(value);
                }
                if (N <= 10) System.out.println();
                writer.println();
            }
            if (N <= 10) System.out.println();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Data loadData() {
        Data data = null;
        try (ObjectInputStream ois = new ObjectInputStream(
                new BufferedInputStream(new FileInputStream(Data.INPUT_FILE)))) {
            data = (Data) ois.readObject();
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
        return data;
    }
}
