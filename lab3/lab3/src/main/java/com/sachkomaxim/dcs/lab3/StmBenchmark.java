package com.sachkomaxim.dcs.lab3;

import org.multiverse.api.StmUtils;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class StmBenchmark {
    private static int THREADS;
    private static int ITERATIONS_FOR_THREAD;

    public static void main(int threads, int iterations, int startRow) {
        THREADS = threads;
        ITERATIONS_FOR_THREAD = iterations;

        long start = System.nanoTime();
        executeUpdateQueriesWithStm(startRow);
        long end = System.nanoTime();
        System.out.println("Queries execution time with STM: " + (end - start) / 1000000F + " ms");
    }

    public static void executeUpdateQueriesWithStm(int startRow) {
        ExecutorService executor = Executors.newFixedThreadPool(THREADS);
        try {
            for (int t = 0; t < THREADS; t++) {
                executor.execute(() -> {
                    try (Connection conn = DriverManager.getConnection(Main.DB_URL, Main.DB_USER, Main.DB_PASSWORD)) {
                        PreparedStatement pstmt = conn.prepareStatement(Main.UPDATE_QUERY);

                        for (int i = 0; i < ITERATIONS_FOR_THREAD; i++) {
                            int targetId = startRow + (i % 100);

                            StmUtils.atomic(() -> {
                                try {
                                    pstmt.setInt(1, (int)(Math.random() * 1350 + 25));
                                    pstmt.setInt(2, targetId);
                                    pstmt.executeUpdate();
                                } catch (SQLException e) {
                                    throw new RuntimeException(e);
                                }
                            });
                        }
                    } catch (SQLException e) {
                        throw new RuntimeException(e);
                    }
                });
            }
        } finally{
            executor.shutdown();
            try {
                executor.awaitTermination(2, java.util.concurrent.TimeUnit.MINUTES);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}
