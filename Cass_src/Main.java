package org.cass;

import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Main {
    public static void main(String[] args) throws IOException, ClassNotFoundException, NoSuchMethodException {
        String dir = "/media/j/MY WORK/TUD AM/PhyloTrees/cass_data/restricted_grass_clusters";
        // Read all files in the directory
        var filenames = Stream.of(Objects.requireNonNull(new File(dir).listFiles()))
                .filter(file -> !file.isDirectory())
                .map(File::getName)
                .collect(Collectors.toSet());

        for (String fname : filenames) {
            System.out.println("Handling " + fname);
            Class<?> cass = Class.forName("org.cass.CassAlgorithm");
            Method main = cass.getDeclaredMethod("main", String[].class);
            var cassArgs = new String[]{fname};
            try (ByteArrayOutputStream out = new ByteArrayOutputStream();
                 PrintStream ps = new PrintStream(out)) {
                System.setOut(ps);
                main.invoke(main, cassArgs);
                System.out.print(out);
            } catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            } finally {
                // Reset to the console
                System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
            }

            break;
        }
    }
}