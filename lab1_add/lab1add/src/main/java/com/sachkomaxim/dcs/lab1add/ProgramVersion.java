package com.sachkomaxim.dcs.lab1add;

import java.io.*;

abstract class ProgramVersion {
    public static final int N = Data.N;
    public static final int THREADS = 4;

    public static void saveAndPrintResults(double[] A, double[][] MG, String outputFile) {
        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)))) {
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

    public static Data loadData() {
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
