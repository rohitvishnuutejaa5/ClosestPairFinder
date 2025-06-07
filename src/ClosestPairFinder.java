// Implementation #2: ClosestPairFinder.java
import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.util.List;
import java.util.ArrayList;

class Coordinate2D {
    double xPos, yPos;
    int id;
    
    public Coordinate2D(double xPos, double yPos, int id) {
        this.xPos = xPos;
        this.yPos = yPos;
        this.id = id;
    }
    
    public double distanceTo(Coordinate2D other) {
        return Math.sqrt(Math.pow(this.xPos - other.xPos, 2) + Math.pow(this.yPos - other.yPos, 2));
    }
    
    @Override
    public String toString() {
        return String.format("Point_%d(%.2f, %.2f)", id, xPos, yPos);
    }
}

class PairResult {
    Coordinate2D first, second;
    double minDistance;
    
    public PairResult(Coordinate2D first, Coordinate2D second, double minDistance) {
        this.first = first;
        this.second = second;
        this.minDistance = minDistance;
    }
}

class ChartGenerator {
    static final int CHART_WIDTH = 800;
    static final int CHART_HEIGHT = 600;
    static final int MARGIN = 60;
    
    public static void generateRuntimeComparisonChart(double[] alg1Times, double[] alg2Times, int[] nValues, String outputPath) {
        BufferedImage chart = new BufferedImage(CHART_WIDTH, CHART_HEIGHT, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = chart.createGraphics();
        
        // Set up chart
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, CHART_WIDTH, CHART_HEIGHT);
        
        // Find maximum time for scaling
        double maxTime = 0;
        for (int i = 0; i < alg1Times.length; i++) {
            maxTime = Math.max(maxTime, Math.max(alg1Times[i], alg2Times[i]));
        }
        
        // Add padding to max value
        maxTime *= 1.1;
        
        // Draw axes
        g2d.setColor(Color.BLACK);
        g2d.drawLine(MARGIN, CHART_HEIGHT - MARGIN, CHART_WIDTH - MARGIN, CHART_HEIGHT - MARGIN); // x-axis
        g2d.drawLine(MARGIN, MARGIN, MARGIN, CHART_HEIGHT - MARGIN); // y-axis
        
        // Draw grid lines and labels for y-axis
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 10));
        int numYGrids = 10;
        for (int i = 0; i <= numYGrids; i++) {
            int y = CHART_HEIGHT - MARGIN - (i * (CHART_HEIGHT - 2 * MARGIN) / numYGrids);
            g2d.setColor(new Color(220, 220, 220));
            g2d.drawLine(MARGIN, y, CHART_WIDTH - MARGIN, y);
            g2d.setColor(Color.BLACK);
            String label = String.format("%.2f", (i * maxTime / numYGrids));
            g2d.drawString(label, 5, y + 4);
        }
        
        // Plot data points and lines
        int dataWidth = CHART_WIDTH - 2 * MARGIN;
        int dataHeight = CHART_HEIGHT - 2 * MARGIN;
        
        // Draw X-axis labels and vertical grid lines
        for (int i = 0; i < nValues.length; i++) {
            int x = MARGIN + (i * dataWidth / (nValues.length - 1));
            g2d.setColor(new Color(220, 220, 220));
            g2d.drawLine(x, MARGIN, x, CHART_HEIGHT - MARGIN);
            g2d.setColor(Color.BLACK);
            g2d.drawString(String.valueOf(nValues[i]), x - 15, CHART_HEIGHT - MARGIN + 15);
        }
        
        // Plot ALG1 data points and line
        g2d.setColor(Color.RED);
        int[] alg1XPoints = new int[nValues.length];
        int[] alg1YPoints = new int[nValues.length];
        
        for (int i = 0; i < nValues.length; i++) {
            alg1XPoints[i] = MARGIN + (i * dataWidth / (nValues.length - 1));
            alg1YPoints[i] = CHART_HEIGHT - MARGIN - (int)(alg1Times[i] * dataHeight / maxTime);
            g2d.fillOval(alg1XPoints[i] - 3, alg1YPoints[i] - 3, 6, 6);
        }
        
        for (int i = 0; i < nValues.length - 1; i++) {
            g2d.drawLine(alg1XPoints[i], alg1YPoints[i], alg1XPoints[i+1], alg1YPoints[i+1]);
        }
        
        // Plot ALG2 data points and line
        g2d.setColor(Color.BLUE);
        int[] alg2XPoints = new int[nValues.length];
        int[] alg2YPoints = new int[nValues.length];
        
        for (int i = 0; i < nValues.length; i++) {
            alg2XPoints[i] = MARGIN + (i * dataWidth / (nValues.length - 1));
            alg2YPoints[i] = CHART_HEIGHT - MARGIN - (int)(alg2Times[i] * dataHeight / maxTime);
            g2d.fillOval(alg2XPoints[i] - 3, alg2YPoints[i] - 3, 6, 6);
        }
        
        for (int i = 0; i < nValues.length - 1; i++) {
            g2d.drawLine(alg2XPoints[i], alg2YPoints[i], alg2XPoints[i+1], alg2YPoints[i+1]);
        }
        
        // Add title and labels
        g2d.setColor(Color.BLACK);
        g2d.setFont(new Font("SansSerif", Font.BOLD, 16));
        g2d.drawString("Empirical Running Time Comparison", CHART_WIDTH / 2 - 130, 30);
        
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 12));
        g2d.drawString("Input Size (n)", CHART_WIDTH / 2 - 40, CHART_HEIGHT - 15);
        
        // Add legend
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 10));
        g2d.setColor(Color.RED);
        g2d.fillRect(CHART_WIDTH - 150, 40, 10, 10);
        g2d.setColor(Color.BLACK);
        g2d.drawString("Brute Force (ALG1)", CHART_WIDTH - 135, 50);
        
        g2d.setColor(Color.BLUE);
        g2d.fillRect(CHART_WIDTH - 150, 60, 10, 10);
        g2d.setColor(Color.BLACK);
        g2d.drawString("Divide & Conquer (ALG2)", CHART_WIDTH - 135, 70);
        
        // Draw y-axis label (rotated)
        AffineTransform original = g2d.getTransform();
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 12));
        g2d.rotate(-Math.PI/2, 15, CHART_HEIGHT / 2);
        g2d.drawString("Running Time (ms)", 15, CHART_HEIGHT / 2);
        g2d.setTransform(original);
        
        g2d.dispose();
        
        try {
            ImageIO.write(chart, "png", new File(outputPath));
            System.out.println("Chart saved to: " + outputPath);
        } catch (IOException e) {
            System.out.println("Error saving chart: " + e.getMessage());
        }
    }
    
    public static void generateEmpiricalVsPredictedChart(double[] empiricalTimes, double[] predictedTimes, 
    // Implementation #2: ClosestPairFinder.java (Continued)
    int[] nValues, String title, String outputPath) {
        BufferedImage chart = new BufferedImage(CHART_WIDTH, CHART_HEIGHT, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = chart.createGraphics();
        
        // Set up chart
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, CHART_WIDTH, CHART_HEIGHT);
        
        // Find maximum time for scaling
        double maxTime = 0;
        for (int i = 0; i < empiricalTimes.length; i++) {
            maxTime = Math.max(maxTime, Math.max(empiricalTimes[i], predictedTimes[i]));
        }
        
        // Add padding to max value
        maxTime *= 1.1;
        
        // Draw axes
        g2d.setColor(Color.BLACK);
        g2d.drawLine(MARGIN, CHART_HEIGHT - MARGIN, CHART_WIDTH - MARGIN, CHART_HEIGHT - MARGIN); // x-axis
        g2d.drawLine(MARGIN, MARGIN, MARGIN, CHART_HEIGHT - MARGIN); // y-axis
        
        // Draw grid lines and labels for y-axis
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 10));
        int numYGrids = 10;
        for (int i = 0; i <= numYGrids; i++) {
            int y = CHART_HEIGHT - MARGIN - (i * (CHART_HEIGHT - 2 * MARGIN) / numYGrids);
            g2d.setColor(new Color(220, 220, 220));
            g2d.drawLine(MARGIN, y, CHART_WIDTH - MARGIN, y);
            g2d.setColor(Color.BLACK);
            String label = String.format("%.2f", (i * maxTime / numYGrids));
            g2d.drawString(label, 5, y + 4);
        }
        
        // Plot data points and lines
        int dataWidth = CHART_WIDTH - 2 * MARGIN;
        int dataHeight = CHART_HEIGHT - 2 * MARGIN;
        
        // Draw X-axis labels and vertical grid lines
        for (int i = 0; i < nValues.length; i++) {
            int x = MARGIN + (i * dataWidth / (nValues.length - 1));
            g2d.setColor(new Color(220, 220, 220));
            g2d.drawLine(x, MARGIN, x, CHART_HEIGHT - MARGIN);
            g2d.setColor(Color.BLACK);
            g2d.drawString(String.valueOf(nValues[i]), x - 15, CHART_HEIGHT - MARGIN + 15);
        }
        
        // Plot Empirical data points and line
        g2d.setColor(Color.RED);
        int[] empiricalXPoints = new int[nValues.length];
        int[] empiricalYPoints = new int[nValues.length];
        
        for (int i = 0; i < nValues.length; i++) {
            empiricalXPoints[i] = MARGIN + (i * dataWidth / (nValues.length - 1));
            empiricalYPoints[i] = CHART_HEIGHT - MARGIN - (int)(empiricalTimes[i] * dataHeight / maxTime);
            g2d.fillOval(empiricalXPoints[i] - 3, empiricalYPoints[i] - 3, 6, 6);
        }
        
        for (int i = 0; i < nValues.length - 1; i++) {
            g2d.drawLine(empiricalXPoints[i], empiricalYPoints[i], empiricalXPoints[i+1], empiricalYPoints[i+1]);
        }
        
        // Plot Predicted data points and line
        g2d.setColor(Color.BLUE);
        int[] predictedXPoints = new int[nValues.length];
        int[] predictedYPoints = new int[nValues.length];
        
        for (int i = 0; i < nValues.length; i++) {
            predictedXPoints[i] = MARGIN + (i * dataWidth / (nValues.length - 1));
            predictedYPoints[i] = CHART_HEIGHT - MARGIN - (int)(predictedTimes[i] * dataHeight / maxTime);
            g2d.fillOval(predictedXPoints[i] - 3, predictedYPoints[i] - 3, 6, 6);
        }
        
        for (int i = 0; i < nValues.length - 1; i++) {
            g2d.drawLine(predictedXPoints[i], predictedYPoints[i], predictedXPoints[i+1], predictedYPoints[i+1]);
        }
        
        // Add title and labels
        g2d.setColor(Color.BLACK);
        g2d.setFont(new Font("SansSerif", Font.BOLD, 16));
        g2d.drawString(title, CHART_WIDTH / 2 - 170, 30);
        
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 12));
        g2d.drawString("Input Size (n)", CHART_WIDTH / 2 - 40, CHART_HEIGHT - 15);
        
        // Add legend
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 10));
        g2d.setColor(Color.RED);
        g2d.fillRect(CHART_WIDTH - 150, 40, 10, 10);
        g2d.setColor(Color.BLACK);
        g2d.drawString("Empirical", CHART_WIDTH - 135, 50);
        
        g2d.setColor(Color.BLUE);
        g2d.fillRect(CHART_WIDTH - 150, 60, 10, 10);
        g2d.setColor(Color.BLACK);
        g2d.drawString("Predicted", CHART_WIDTH - 135, 70);
        
        // Draw y-axis label (rotated)
        AffineTransform original = g2d.getTransform();
        g2d.setFont(new Font("SansSerif", Font.PLAIN, 12));
        g2d.rotate(-Math.PI/2, 15, CHART_HEIGHT / 2);
        g2d.drawString("Running Time (ms)", 15, CHART_HEIGHT / 2);
        g2d.setTransform(original);
        
        g2d.dispose();
        
        try {
            ImageIO.write(chart, "png", new File(outputPath));
            System.out.println("Chart saved to: " + outputPath);
        } catch (IOException e) {
            System.out.println("Error saving chart: " + e.getMessage());
        }
    }
}

class TableCreator {
    public static double[] generatePerformanceTable(double[] theoreticalRT, double[] empiricalRT, int[] sizes, String algorithm, String filename) {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filename))) {
            writer.println("Computing constant c for " + algorithm);
            writer.println("--------------------------------------------------");
            writer.println("n\tTheoreticalRT\tEmpiricalRT(ms)\tRatio\tPredictedRT");
            writer.println("--------------------------------------------------");
            
            // Calculate ratios
            double[] ratios = new double[sizes.length];
            for (int i = 0; i < sizes.length; i++) {
                ratios[i] = empiricalRT[i] / theoreticalRT[i];
                writer.printf("%d\t%.2e\t%.4f\t%.8e\t-\n", sizes[i], theoreticalRT[i], empiricalRT[i], ratios[i]);
            }
            
            // Find maximum ratio (avoiding outliers)
            Arrays.sort(ratios);
            // Use the average of middle values to avoid outliers
            double c = (ratios[ratios.length/2] + ratios[ratios.length/2 - 1]) / 2;
            if (ratios.length % 2 == 1) {
                c = ratios[ratios.length/2];
            }
            
            writer.println("--------------------------------------------------");
            writer.printf("Constant c = %.8e\n", c);
            writer.println("--------------------------------------------------");
            
            // Calculate predicted running times
            writer.println("\nn\tTheoreticalRT\tEmpiricalRT(ms)\tPredictedRT(ms)");
            writer.println("--------------------------------------------------");
            double[] predictedRT = new double[sizes.length];
            for (int i = 0; i < sizes.length; i++) {
                predictedRT[i] = c * theoreticalRT[i];
                writer.printf("%d\t%.2e\t%.4f\t%.4f\n", sizes[i], theoreticalRT[i], empiricalRT[i], predictedRT[i]);
            }
            
            System.out.println("Table saved to: " + filename);
            return predictedRT;
        } catch (IOException e) {
            System.out.println("Error creating table: " + e.getMessage());
            return new double[sizes.length]; // Return empty array on error
        }
    }
}

public class ClosestPairFinder {
    static final int ITERATIONS = 10; // Number of iterations (m)
    static final int MAX_RANGE = 100000; // Maximum coordinate range
    
    // ALG1: Brute Force algorithm to find closest pair
    public static PairResult findClosestPairBruteForce(Coordinate2D[] points) {
        int n = points.length;
        if (n < 2) return null;
        
        Coordinate2D minPoint1 = null, minPoint2 = null;
        double minDist = Double.MAX_VALUE;
        
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double dist = points[i].distanceTo(points[j]);
                if (dist < minDist) {
                    minDist = dist;
                    minPoint1 = points[i];
                    minPoint2 = points[j];
                }
            }
        }
        
        return new PairResult(minPoint1, minPoint2, minDist);
    }
    
    // ALG2: Divide and Conquer algorithm to find closest pair
    public static PairResult findClosestPairDivideConquer(Coordinate2D[] points) {
        // Sort points by x coordinate
        Coordinate2D[] sortedX = points.clone();
        Arrays.sort(sortedX, (p1, p2) -> Double.compare(p1.xPos, p2.xPos));
        
        return findClosestPairRecursive(sortedX, 0, sortedX.length - 1);
    }
    
    private static PairResult findClosestPairRecursive(Coordinate2D[] sortedX, int left, int right) {
        // Base cases
        if (right - left <= 2) {
            return findClosestPairBruteForce(Arrays.copyOfRange(sortedX, left, right + 1));
        }
        
        // Divide: find the middle point
        int mid = (left + right) / 2;
        double midX = sortedX[mid].xPos;
        
        // Conquer: find closest pairs in left and right halves
        PairResult leftPair = findClosestPairRecursive(sortedX, left, mid);
        PairResult rightPair = findClosestPairRecursive(sortedX, mid + 1, right);
        
        // Combine: determine minimum distance
        PairResult minPair = (leftPair.minDistance <= rightPair.minDistance) ? leftPair : rightPair;
        double delta = minPair.minDistance;
        
        // Collect points within delta of the middle line
        List<Coordinate2D> strip = new ArrayList<>();
        for (int i = left; i <= right; i++) {
            if (Math.abs(sortedX[i].xPos - midX) < delta) {
                strip.add(sortedX[i]);
            }
        }
        
        // Sort strip by y coordinate
        strip.sort((p1, p2) -> Double.compare(p1.yPos, p2.yPos));
        
        // Check pairs across the middle line
        for (int i = 0; i < strip.size(); i++) {
            for (int j = i + 1; j < strip.size() && (strip.get(j).yPos - strip.get(i).yPos) < delta; j++) {
                double dist = strip.get(i).distanceTo(strip.get(j));
                if (dist < delta) {
                    delta = dist;
                    minPair = new PairResult(strip.get(i), strip.get(j), dist);
                }
            }
        }
        
        return minPair;
    }
    
    // Generate array of unique random points
    public static Coordinate2D[] generateRandomPoints(int n) {
        Coordinate2D[] points = new Coordinate2D[n];
        Random random = new Random();
        Set<String> uniqueCoords = new HashSet<>();
        
        for (int i = 0; i < n; i++) {
            double x, y;
            String coordKey;
            
            do {
                x = random.nextDouble() * MAX_RANGE;
                y = random.nextDouble() * MAX_RANGE;
                coordKey = x + ":" + y;
            } while (uniqueCoords.contains(coordKey));
            
            uniqueCoords.add(coordKey);
            points[i] = new Coordinate2D(x, y, i);
        }
        
        return points;
    }
    
    public static void main(String[] args) {
        System.out.println("Starting Closest Pair Analysis...");
        
        // Define input sizes
        int[] sizes = new int[10];
        for (int i = 0; i < 10; i++) {
            sizes[i] = 10000 * (i + 1);
        }
        
        // Storage for timing results
        double[][] alg1RuntimesPerIteration = new double[ITERATIONS][sizes.length];
        double[][] alg2RuntimesPerIteration = new double[ITERATIONS][sizes.length];
        double[] alg1AvgRuntimes = new double[sizes.length];
        double[] alg2AvgRuntimes = new double[sizes.length];
        
        // Calculate theoretical running times
        double[] alg1TheoreticalRT = new double[sizes.length];
        double[] alg2TheoreticalRT = new double[sizes.length];
        
        for (int i = 0; i < sizes.length; i++) {
            alg1TheoreticalRT[i] = Math.pow(sizes[i], 2);
            alg2TheoreticalRT[i] = sizes[i] * Math.log(sizes[i]) / Math.log(2);
        }
        
        // Run experiments
        for (int sizeIndex = 0; sizeIndex < sizes.length; sizeIndex++) {
            int n = sizes[sizeIndex];
            System.out.println("\nProcessing input size n = " + n);
            
            for (int iter = 0; iter < ITERATIONS; iter++) {
                System.out.printf("  Iteration %d/%d... ", iter + 1, ITERATIONS);
                
                // Generate random points for this iteration
                Coordinate2D[] points = generateRandomPoints(n);
                
                // Run ALG1 (Brute Force)
                long startTime = System.nanoTime();
                PairResult alg1Result = findClosestPairBruteForce(points);
                long endTime = System.nanoTime();
                double alg1Time = (endTime - startTime) / 1_000_000.0; // Convert to milliseconds
                alg1RuntimesPerIteration[iter][sizeIndex] = alg1Time;
                
                // Run ALG2 (Divide and Conquer)
                startTime = System.nanoTime();
                PairResult alg2Result = findClosestPairDivideConquer(points);
                endTime = System.nanoTime();
                double alg2Time = (endTime - startTime) / 1_000_000.0; // Convert to milliseconds
                alg2RuntimesPerIteration[iter][sizeIndex] = alg2Time;
                
                // Verify results
                if (Math.abs(alg1Result.minDistance - alg2Result.minDistance) > 1e-10) {
                    System.out.println("Warning: Different results!");
                    System.out.println("  ALG1: " + alg1Result.minDistance);
                    System.out.println("  ALG2: " + alg2Result.minDistance);
                }
                
                System.out.printf("Done. ALG1: %.2f ms, ALG2: %.2f ms\n", alg1Time, alg2Time);
            }
            
            // Calculate average runtimes for this size
            double alg1Sum = 0, alg2Sum = 0;
            for (int iter = 0; iter < ITERATIONS; iter++) {
                alg1Sum += alg1RuntimesPerIteration[iter][sizeIndex];
                alg2Sum += alg2RuntimesPerIteration[iter][sizeIndex];
            }
            
            alg1AvgRuntimes[sizeIndex] = alg1Sum / ITERATIONS;
            alg2AvgRuntimes[sizeIndex] = alg2Sum / ITERATIONS;
            
            System.out.printf("Average results for n=%d:\n", n);
            System.out.printf("  ALG1 (Brute Force): %.2f ms\n", alg1AvgRuntimes[sizeIndex]);
            System.out.printf("  ALG2 (Divide & Conquer): %.2f ms\n", alg2AvgRuntimes[sizeIndex]);
        }
        
        // Generate tables and calculate predicted runtimes
        System.out.println("\nGenerating tables...");
        double[] alg1PredictedRT = TableCreator.generatePerformanceTable(
            alg1TheoreticalRT, alg1AvgRuntimes, sizes, "ALG1 (Brute Force)", "alg1_table.txt");
            
        double[] alg2PredictedRT = TableCreator.generatePerformanceTable(
            alg2TheoreticalRT, alg2AvgRuntimes, sizes, "ALG2 (Divide & Conquer)", "alg2_table.txt");
        
        // Generate charts
        System.out.println("\nGenerating charts...");
        ChartGenerator.generateRuntimeComparisonChart(
            alg1AvgRuntimes, alg2AvgRuntimes, sizes, "runtime_comparison_chart.png");
            
        ChartGenerator.generateEmpiricalVsPredictedChart(
            alg1AvgRuntimes, alg1PredictedRT, sizes, 
            "ALG1: Empirical vs. Predicted Running Time", "alg1_empirical_vs_predicted.png");
            
        ChartGenerator.generateEmpiricalVsPredictedChart(
            alg2AvgRuntimes, alg2PredictedRT, sizes,
            "ALG2: Empirical vs. Predicted Running Time", "alg2_empirical_vs_predicted.png");
        
        System.out.println("\nAnalysis complete. All tables and charts have been generated.");
    }
}