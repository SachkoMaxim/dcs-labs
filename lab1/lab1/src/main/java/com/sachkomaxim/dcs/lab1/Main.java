package com.sachkomaxim.dcs.lab1;

public class Main {
    public static void main(String[] args) {
        Data.generateContent();

        System.out.println("Starting localized variant of the program...");
        System.out.println("------------------------------------------------");
        LocalizedVersion.main(args);

        System.out.println();
        System.out.println("================================================");
        System.out.println("================================================");
        System.out.println();

        System.out.println("Starting shared variant of the program...");
        System.out.println("------------------------------------------------");
        SharedVersion.main(args);
    }
}
