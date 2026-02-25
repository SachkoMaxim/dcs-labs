package com.sachkomaxim.dcs.lab1;

import java.io.*;
import java.util.Arrays;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.BrokenBarrierException;

public class LocalizedVersion {
    private static final int N = Data.N;
    private static final int THREADS = 4;

    private static final String OUTPUT_FILE = "output_localized.txt";

    static class ComputingRunnable implements Runnable {
        private final int start, stop, localRange;
        private final double[][] localMC, localMZ, localMT, localME;
        private final double[] localB, localD, localE;
        private final double[] wLocalResults;
        private final CyclicBarrier barrier;
        private double[] localAResultPart;
        private double[][] localMGResultPart;
        private final int threadNumber;

        ComputingRunnable(int start, int stop, Data data, double[] wLocalResults, CyclicBarrier barrier, int threadNumber) {
            this.start = start;
            this.stop = stop;
            this.threadNumber = threadNumber;
            this.localRange = stop - start;

            this.localB = Arrays.copyOf(data.B, data.B.length);
            this.localD = Arrays.copyOf(data.D, data.D.length);
            this.localE = Arrays.copyOf(data.E, data.E.length);

            this.localMC = new double[N][N];
            this.localMZ = new double[N][N];
            this.localMT = new double[N][N];
            this.localME = new double[N][N];

            for (int i = 0; i < N; i++) {
                this.localMC[i] = Arrays.copyOf(data.MC[i], data.MC[i].length);
                this.localMZ[i] = Arrays.copyOf(data.MZ[i], data.MZ[i].length);
                this.localMT[i] = Arrays.copyOf(data.MT[i], data.MT[i].length);
                this.localME[i] = Arrays.copyOf(data.ME[i], data.ME[i].length);
            }

            this.wLocalResults = wLocalResults;
            this.barrier = barrier;

            this.localAResultPart = new double[localRange];
            this.localMGResultPart = new double[localRange][N];
        }

        @Override
        public void run() {
            System.out.printf("Thread %d has started its computation%n", threadNumber);

            localAResultPart = calculateLocalAResultPart();
            localMGResultPart = calculateLocalMGResultPart();

            System.out.printf("Thread %d has finished its computation%n", threadNumber);
        }

        private double[] calculateLocalAResultPart() {
            double[] QPart = partialMultiplyVectorAndMatrix(localB, localMC);
            double[] RPart = partialMultiplyVectorAndMatrix(localD, localMZ);
            double[] SPart = partialMultiplyVectorAndMatrix(localE, localMT);
            return partialAddOfThreeOrMoreVectors(QPart, RPart, SPart);
        }

        private double[][] calculateLocalMGResultPart() {
            wLocalResults[threadNumber - 1] = partialFindMinWlocalInVector(localD, localB);

            try {
                barrier.await();
            } catch (InterruptedException | BrokenBarrierException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException(e);
            }

            double wGlobal = Double.MAX_VALUE;
            for (int i = 0; i < THREADS; i++) {
                wGlobal = Math.min(wGlobal, wLocalResults[i]);
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
                result[i - start] = sum;
            }
            return result;
        }

        private double[][] partialMultiplyTwoMatrices(double[][] matrix1, double[][] matrix2) {
            if (matrix1[0].length != matrix2.length) {
                throw new IllegalArgumentException("Matrix2 must have the same number of rows as number of columns in Matrix1");
            }

            double[][] result = new double[localRange][N];

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

            for (int i = start; i < stop; i++) {
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
                    result[i - start][j] = w * sum;
                }
            }
            return result;
        }

        public double[] getLocalAResultPart() {
            return localAResultPart;
        }

        public double[][] getLocalMGResultPart() {
            return localMGResultPart;
        }
    }

    public static void main(String[] args) {
        try {
            System.out.println("Starting computation with " + THREADS + " threads");

            long startTime = System.currentTimeMillis();
            Data data = loadData();

            Thread[] threads = new Thread[THREADS];
            ComputingRunnable[] workers = new ComputingRunnable[THREADS];
            CyclicBarrier barrier = new CyclicBarrier(THREADS);
            double[] wLocalResults = new double[THREADS];
            int rowsPerThread = N / THREADS;
            int startRow = 0;

            for (int i = 0; i < THREADS; i++) {
                int endRow = (i == THREADS - 1) ? N : startRow + rowsPerThread;
                workers[i] = new ComputingRunnable(startRow, endRow, data, wLocalResults, barrier, i + 1);
                threads[i] = new Thread(workers[i]);
                threads[i].start();
                startRow = endRow;
            }

            double[] finalA = new double[N];
            double[][] finalMG = new double[N][N];

            for (int i = 0; i < THREADS; i++) {
                threads[i].join();
                System.arraycopy(workers[i].getLocalAResultPart(), 0, finalA,
                        i * rowsPerThread, workers[i].getLocalAResultPart().length);
                for (int j = 0; j < workers[i].getLocalMGResultPart().length; j++) {
                    System.arraycopy(workers[i].getLocalMGResultPart()[j], 0,
                            finalMG[i * rowsPerThread + j], 0, N);
                }
            }

            saveAndPrintResults(finalA, finalMG);
            System.out.println("Results calculated and written to " + OUTPUT_FILE);

            double executionTime = (System.currentTimeMillis() - startTime) / 1000.0;
            System.out.printf("%nComputation completed in %.3f seconds%n", executionTime);

        } catch (Exception e) {
            e.printStackTrace();
        }
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
