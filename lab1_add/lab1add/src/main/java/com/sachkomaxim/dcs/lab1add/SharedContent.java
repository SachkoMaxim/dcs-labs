package com.sachkomaxim.dcs.lab1add;

import java.util.concurrent.locks.ReentrantLock;

public class SharedContent {
    // Дані
    public double[][] MC, MZ, MT, ME, MG, MCMTw, MZME;
    public double[] B, D, E, A, Q, R, S;
    public double w;
    public final int N;

    // Об'єкти синхронізації (Locks)
    public static final ReentrantLock wLock = new ReentrantLock();
    public final ReentrantLock ALock = new ReentrantLock();
    public final ReentrantLock MGLock = new ReentrantLock();
    public final ReentrantLock QLock = new ReentrantLock();
    public final ReentrantLock RLock = new ReentrantLock();
    public final ReentrantLock SLock = new ReentrantLock();
    public final ReentrantLock MZMELock = new ReentrantLock();
    public final ReentrantLock MCMTwLock = new ReentrantLock();

    public SharedContent(int n, Data data) {
        this.N = n;

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
}
