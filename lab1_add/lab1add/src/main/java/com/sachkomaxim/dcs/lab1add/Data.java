package com.sachkomaxim.dcs.lab1add;

import java.io.*;
import java.util.Random;

public class Data implements Serializable {
    public static final String INPUT_FILE = "input.bin";

    public static final int N = 4;

    public double[][] MC;
    public double[] B;
    public double[][] MZ;
    public double[] D;
    public double[][] MT;
    public double[] E;
    public double[][] ME;

    public Data() {
        MC = new double[N][N];
        B = new double[N];
        MZ = new double[N][N];
        D = new double[N];
        MT = new double[N][N];
        E = new double[N];
        ME = new double[N][N];
    }

    public static void generateContent() {
        Data data = new Data();
        Random random = new Random();

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                data.MC[i][j] = getRandomDouble(random);
                data.MZ[i][j] = getRandomDouble(random);
                data.MT[i][j] = getRandomDouble(random);
                data.ME[i][j] = getRandomDouble(random);
            }
        }

        for (int i = 0; i < N; i++) {
            data.B[i] = getRandomDouble(random);
            data.D[i] = getRandomDouble(random);
            data.E[i] = getRandomDouble(random);
        }

        saveResults(data);

        try (ObjectOutputStream oos = new ObjectOutputStream(
                new BufferedOutputStream(new FileOutputStream("input.bin")))) {
            oos.writeObject(data);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static double getRandomDouble(Random random) {
        int expnt = random.nextInt(7) - 2;
        double coeff = 1 + random.nextDouble() * 9;
        return coeff * Math.pow(10, expnt);
    }

    private static void saveResults(Data data) {
        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("input.txt")))) {
            writer.println("\nVector B:");
            for (int i = 0; i < N; i++) {
                String value = String.format("%f ", data.B[i]);
                writer.print(value);
            }
            writer.println();

            writer.println("\nVector D:");
            for (int i = 0; i < N; i++) {
                String value = String.format("%f ", data.D[i]);
                writer.print(value);
            }
            writer.println();

            writer.println("\nVector E:");
            for (int i = 0; i < N; i++) {
                String value = String.format("%f ", data.E[i]);
                writer.print(value);
            }
            writer.println();
            writer.println();

            writer.println("Matrix MC:");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    String value = String.format("%f ", data.MC[i][j]);
                    writer.print(value);
                }
                writer.println();
            }
            writer.println();

            writer.println("Matrix MZ:");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    String value = String.format("%f ", data.MZ[i][j]);
                    writer.print(value);
                }
                writer.println();
            }
            writer.println();

            writer.println("Matrix MT:");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    String value = String.format("%f ", data.MT[i][j]);
                    writer.print(value);
                }
                writer.println();
            }
            writer.println();

            writer.println("Matrix ME:");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    String value = String.format("%f ", data.ME[i][j]);
                    writer.print(value);
                }
                writer.println();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
