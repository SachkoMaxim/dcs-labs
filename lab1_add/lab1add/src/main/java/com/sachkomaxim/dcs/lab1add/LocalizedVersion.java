package com.sachkomaxim.dcs.lab1add;

import java.util.concurrent.CyclicBarrier;

public class LocalizedVersion extends ProgramVersion {
    public static final String OUTPUT_FILE = "output_localized.txt";

    public static void main(String[] args) {
        try {
            System.out.println("Starting computation with " + THREADS + " threads");

            long startTime = System.currentTimeMillis();
            Data data = loadData();

            Thread[] threads = new Thread[THREADS];
            LocalizedWorker[] workers = new LocalizedWorker[THREADS];
            CyclicBarrier barrier = new CyclicBarrier(THREADS);
            double[] wLocalResults = new double[THREADS];
            int rowsPerThread = N / THREADS;
            int startRow = 0;

            for (int i = 0; i < THREADS; i++) {
                int endRow = (i == THREADS - 1) ? N : startRow + rowsPerThread;
                workers[i] = new LocalizedWorker(N, startRow, endRow, data, wLocalResults, barrier, i + 1);
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

            saveAndPrintResults(finalA, finalMG, OUTPUT_FILE);
            System.out.println("Results calculated and written to " + OUTPUT_FILE);

            double executionTime = (System.currentTimeMillis() - startTime) / 1000.0;
            System.out.printf("%nComputation completed in %.3f seconds%n", executionTime);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
