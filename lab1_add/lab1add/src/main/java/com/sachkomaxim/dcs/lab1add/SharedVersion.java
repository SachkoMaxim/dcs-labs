package com.sachkomaxim.dcs.lab1add;

import java.util.concurrent.CyclicBarrier;

public class SharedVersion extends ProgramVersion {
    public static final String OUTPUT_FILE = "output_shared.txt";

    public static void main(String[] args) {
        try {
            System.out.println("Starting computation with " + THREADS + " threads");

            long startTime = System.currentTimeMillis();
            Data data = loadData();
            SharedContent ctn = new SharedContent(N, data);

            Thread[] threads = new Thread[THREADS];
            SharedWorker[] workers = new SharedWorker[THREADS];
            CyclicBarrier barrier = new CyclicBarrier(THREADS);
            int rowsPerThread = N / THREADS;
            int startRow = 0;

            for (int i = 0; i < THREADS; i++) {
                int endRow = (i == THREADS - 1) ? N : startRow + rowsPerThread;
                workers[i] = new SharedWorker(N, startRow, endRow, barrier, i + 1, ctn);
                threads[i] = new Thread(workers[i]);
                threads[i].start();
                startRow = endRow;
            }

            for (int i = 0; i < THREADS; i++) {
                threads[i].join();
            }

            saveAndPrintResults(ctn.A, ctn.MG, OUTPUT_FILE);
            System.out.println("Results calculated and written to " + OUTPUT_FILE);

            double executionTime = (System.currentTimeMillis() - startTime) / 1000.0;
            System.out.printf("%nComputation completed in %.3f seconds%n", executionTime);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
