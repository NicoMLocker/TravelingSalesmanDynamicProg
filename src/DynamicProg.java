import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.List;


public class DynamicProg {

    static ThreadMXBean bean = ManagementFactory.getThreadMXBean();

    /* define constants */
    private final int N, start;
    private final int[][] distance;
    private List<Integer> tour = new ArrayList<>();
    private double minTourCost = Double.POSITIVE_INFINITY;
    private boolean ranSolver = false;
    static int numberOfTrials = 15;
    static int MAXINPUTSIZE = 10; //12 max for actual test
    static Random random = new Random();
    static String ResultsFolderPath = "/home/nicolocker/Results/"; // pathname to results folder
    static FileWriter resultsFile;
    static PrintWriter resultsWriter;

    public static void main(String[] args) {

        verifyWorks();
        System.out.println("\n");

        System.out.println("Running first full experiment...");
        runFullExperiment("DynamicProg1-CircularCost.txt");    //change all 3 to RandomCost, EuclideanCost or CircularCost depending on which is being used
        System.out.println("Running second full experiment...");
        runFullExperiment("DynamicProg2-CircularCost.txt");
        System.out.println("Running third full experiment...");
        runFullExperiment("DynamicProg3-CircularCost.txt");
    }

    public static void verifyWorks() {
        DynamicProg.GenerateRandomCostMatrix(10);
        DynamicProg.GenerateRandomEuclideanCostMatrix(10, 20);
        DynamicProg.GenerateRandomCircularGraphCostMatrix(10, 40);
    }


    public static void runFullExperiment(String resultsFileName) {

        try {
            resultsFile = new FileWriter(ResultsFolderPath + resultsFileName);
            resultsWriter = new PrintWriter(resultsFile);
        } catch (Exception e) {
            System.out.println("*****!!!!!  Had a problem opening the results file " + ResultsFolderPath + resultsFileName);
            return; // not very foolproof... but we do expect to be able to create/open the file...
        }

        ThreadCpuStopWatch TrialStopwatch = new ThreadCpuStopWatch(); // for timing an individual trial
        double prevTime = 0;

        resultsWriter.println("#NumOfVertices        AvgTime         Best Tour Cost"); // # marks a comment in gnuplot data
        resultsWriter.flush();

        for (int inputSize = 2; inputSize <= MAXINPUTSIZE; inputSize++) {
            System.out.println("Running test for input size " + inputSize + " ... ");
            System.gc();

            long batchElapsedTime = 0;
            int startNode = 0;

            int[][] matrix;
            for (long trial = 0; trial < numberOfTrials; trial++) {
                //matrix = DynamicProg.GenerateRandomCostMatrix(inputSize);      //** uncomment to test this one, comment the below
                //matrix = DynamicProg.GenerateRandomEuclideanCostMatrix(inputSize, 50);   //** uncomment to test this one, comment the above
                matrix = DynamicProg.GenerateRandomCircularGraphCostMatrix(inputSize, 50); //** uncomment to test this one, comment the above
                TrialStopwatch.start(); // *** uncomment this line if timing trials individually
                DynamicProg dynamicProg = new DynamicProg(startNode, matrix);
                dynamicProg.solve();
                batchElapsedTime = batchElapsedTime + TrialStopwatch.elapsedTime(); // *** uncomment this line if timing trials individually
            }

            double averageTimePerTrialInBatch = (double) batchElapsedTime / (double) numberOfTrials;

            double ratio = 0;
            if (prevTime > 0) {
                ratio = averageTimePerTrialInBatch / prevTime;
            }

            prevTime = averageTimePerTrialInBatch;

            /* print data for this size of input */
            resultsWriter.printf("%6d  %20.2f %13.2f\n", inputSize, averageTimePerTrialInBatch, ratio); // might as well make the columns look nice
            resultsWriter.flush();
            System.out.println(" ....done.");

        }
    }

    public static int[][] GenerateRandomCostMatrix ( int edgeCost){
        int[][] matrix = new int[edgeCost][edgeCost];

        for (int i = 0; i < edgeCost; i++) {
            //  System.out.println("++++");

            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                    //   System.out.println(matrix[i][j]);

                } else {
                    int temp = random.nextInt(100);
                    matrix[i][j] = temp;
                    //    System.out.println(matrix[i][j]);
                }
            }
        }
        System.out.println("Generate Random Cost Matrix:");
        printMatrix(matrix);

        return matrix;
    }

    public static int[][] GenerateRandomEuclideanCostMatrix ( int numOfVerticies, int coordinate){
        Vertex[] vertice = new Vertex[numOfVerticies];

        for (int i = 0; i < numOfVerticies; i++) {
            vertice[i] = new Vertex(random.nextInt(coordinate - 1), random.nextInt(coordinate - 1), i);
        }

        int[][] matrix = new int[numOfVerticies][numOfVerticies];

        for (int i = 0; i < numOfVerticies; i++) {
            //System.out.println("++++");
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                    //System.out.println(matrix[i][j]);
                } else {
                    matrix[i][j] = vertice[i].distance(vertice[j]);
                    matrix[j][i] = vertice[j].distance(vertice[i]);
                    //System.out.println(matrix[i][j]);
                }
            }
        }

        System.out.println("Generate Random Euclidean Cost Matrix:");
        printMatrix(matrix);

        return matrix;
    }

    public static int[][] GenerateRandomCircularGraphCostMatrix ( int numOfVerticies, int radius){
        Vertex[] vertice = new Vertex[numOfVerticies];

        double angle = 360 / numOfVerticies;
        double curAngle = 0;
        ArrayList<Vertex> sortVertices = new ArrayList<Vertex>();

        for (int i = 0; i < numOfVerticies; i++) {
            double radian = curAngle * Math.PI / 180;
            double x = Math.cos(radian) * radius;
            double y = Math.sin(radian) * radius;
            sortVertices.add(new Vertex(x, y, i));
            curAngle = curAngle + angle;
        }

        Collections.shuffle(sortVertices);
        vertice = sortVertices.toArray(vertice);

        int[][] matrix = new int[numOfVerticies][numOfVerticies];
        for (int i = 0; i < numOfVerticies; i++) {
            //System.out.println("++++");
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                    //System.out.println(matrix[i][j]);
                } else {
                    matrix[i][j] = vertice[i].distance(vertice[j]);
                    matrix[j][i] = vertice[j].distance(vertice[i]);
                    //System.out.println(matrix[i][j]);
                }
            }
        }

        System.out.println("Generate Circular Graph Cost Matrix:");
        printMatrix(matrix);

        return matrix;
    }

    public static void printMatrix ( int[][] matrix){
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%5d", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println("\n");
    }

    /*https://github.com/williamfiset/Algorithms/blob/master/com/williamfiset/algorithms/graphtheory/TspDynamicProgrammingIterative.java*/
    public DynamicProg(int start, int[][] distance) {
        N = distance.length;
        this.start = start;
        this.distance = distance;
    }

    public List<Integer> getTour() {
        if (!ranSolver) solve();
        return tour;
    }

    public double getTourCost() {
        if (!ranSolver) solve();
        return minTourCost;
    }

    public void solve() {

        if (ranSolver) return;

        final int END_STATE = (1 << N) - 1;
        Double[][] memo = new Double[N][1 << N];

        // Add all outgoing edges from the starting node to memo table.
        for (int end = 0; end < N; end++) {
            if (end == start) continue;
            memo[end][(1 << start) | (1 << end)] = Double.valueOf(distance[start][end]);
        }

        for (int r = 3; r <= N; r++) {
            for (int subset : combinations(r, N)) {
                if (notIn(start, subset)) continue;
                for (int next = 0; next < N; next++) {
                    if (next == start || notIn(next, subset)) continue;
                    int subsetWithoutNext = subset ^ (1 << next);
                    double minDist = Double.POSITIVE_INFINITY;
                    for (int end = 0; end < N; end++) {
                        if (end == start || end == next || notIn(end, subset)) continue;
                        double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
                        if (newDistance < minDist) {
                            minDist = newDistance;
                        }
                    }
                    memo[next][subset] = minDist;
                }
            }
        }

        // Connect tour back to starting node and minimize cost.
        for (int i = 0; i < N; i++) {
            if (i == start) continue;
            double tourCost = memo[i][END_STATE] + distance[i][start];
            if (tourCost < minTourCost) {
                minTourCost = tourCost;
            }
        }

        int lastIndex = start;
        int state = END_STATE;
        tour.add(start);

        for (int i = 1; i < N; i++) {
            int index = -1;
            for (int j = 0; j < N; j++) {
                if (j == start || notIn(j, state)) continue;
                if (index == -1) index = j;
                double prevDist = memo[index][state] + distance[index][lastIndex];
                double newDist = memo[j][state] + distance[j][lastIndex];
                if (newDist < prevDist) {
                    index = j;
                }
            }

            tour.add(index);
            state = state ^ (1 << index);
            lastIndex = index;
        }

        tour.add(start);
        Collections.reverse(tour);

        ranSolver = true;

        System.out.println("Tour: " + getTour());
        System.out.println("Tour Cost: " + getTourCost());
    }

    private static boolean notIn(int elem, int subset) {
        return ((1 << elem) & subset) == 0;
    }

    public static List<Integer> combinations(int r, int n) {
        List<Integer> subsets = new ArrayList<>();
        combinations(0, 0, r, n, subsets);
        return subsets;
    }

    private static void combinations(int set, int at, int r, int n, List<Integer> subsets) {

        // Return early if there are more elements left to select than what is available.
        int elementsLeftToPick = n - at;
        if (elementsLeftToPick < r) return;

        // We selected 'r' elements so we found a valid subset!
        if (r == 0) {
            subsets.add(set);
        } else {
            for (int i = at; i < n; i++) {
                // Try including this element
                set ^= (1 << i);

                combinations(set, i + 1, r - 1, n, subsets);

                // Backtrack and try the instance where we did not include this element
                set ^= (1 << i);
            }
        }
    }
}

