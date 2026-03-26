package com.sachkomaxim.dcs.lab3;

import java.sql.*;
import java.util.Random;

public class Main {
    public static final String DB_URL = "jdbc:h2:mem:testdb;DB_CLOSE_DELAY=-1";
    public static final String DB_USER = "user";
    public static final String DB_PASSWORD = "userPass";

    private static final int THREADS = 8;
    private static final int ITERATIONS_FOR_THREAD = 57;

    static final String CREATE_TABLE_QUERY = """
            CREATE TABLE IF NOT EXISTS lego(
                   id INT PRIMARY KEY,
                   name VARCHAR(255) NOT NULL,
                   theme VARCHAR(255) NOT NULL,
                   pieces INT NOT NULL,
                   cost INT NOT NULL
            );""";

    static final String INSERT_VALUES_QUERY = """
            INSERT INTO lego (id, name, theme, pieces, cost)
            SELECT 
                X,
                'Set ' || X,
                ARRAY_GET(
                    ARRAY['Star Wars', 'Marvel', 'Minecraft', 'Technic', 'Ideas', 'Ninjago', 'City'],\s
                    CAST(FLOOR(RAND() * 7) + 1 AS INT)
                ),
                CAST(RAND() * 151 + 20 AS INT),
                CAST(RAND() * 1411 + 230 AS INT)
            FROM SYSTEM_RANGE(1, 300);""";

    public static final String UPDATE_QUERY = "UPDATE lego SET cost = ? WHERE id = ?";

    public static void main(String[] args) {
        initializeTableWithValues();

        Random random = new Random();
        int startRow = random.nextInt(153) + 1;

        System.out.println("\nStarting HTM test...");
        HtmBenchmark.main(THREADS, ITERATIONS_FOR_THREAD, startRow);

        System.out.println("\nStarting STM test...");
        StmBenchmark.main(THREADS, ITERATIONS_FOR_THREAD, startRow);
    }

    public static void initializeTableWithValues() {
        try (Connection conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASSWORD)) {
            Statement stmt = conn.createStatement();

            stmt.executeUpdate(CREATE_TABLE_QUERY);
            System.out.println("Created the table");

            stmt.executeUpdate(INSERT_VALUES_QUERY);
            System.out.println("Inserted 300 rows into the table");

            System.out.println("H2 Database initialized with 300 rows.");
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }
}
