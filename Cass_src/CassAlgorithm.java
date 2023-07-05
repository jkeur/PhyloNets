//package cass;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.TimeoutException;
import java.util.stream.Collectors;

import static java.lang.Math.max;


/**
 * The Cass algorithm  for computing in minimum level network
 * Leo van Iersel, 2012
 * Modified by Hans Keur, 2023
 * -> Modified for creating the connected components (CC) in IG(C) and handling each CC separately + added timeout [min]
 * Compiled with Java 11
 *
 * Example usage: java org.cass.CassAlgorithm C:/documents/trees.txt --printretnum --timeout=5
 */
public class CassAlgorithm {
    public static final int FASE1 = 0;
    public static final int FASE2 = 1;
    private static final int MAX_CIRCLE_SIZE = 5;
    public static boolean PRINT_EXTRA_INFO = false;
    public static int DUMMY_NUMBER = 9999;
    public static int FOUND = 0; // number of networks found
    public static Vector<String> stringTaxa = new Vector<>();
    public static String leafShape = "circle";
    private static int K_REACHED = 0;

    /**
     * Stand-alone program
     *
     * @param args
     */
    public static void main(String[] args) throws Exception {
        boolean printClusters = false; // label each edge by the cluster its represents
        boolean colourEdges = false; // colour reticulation edges red
        boolean isOnlyOne = true; // compute only one network?
        boolean checkTrees = false; // only construct networks that display the trees?
        boolean printRetNum = false; // Print the level and reticulation number
        int timeout = 0; // Do not use a timeout
        ClusterSet clusterSet = new ClusterSet();

        BufferedReader br;

        if (args.length == 0)  // read from console
        {
            br = new BufferedReader(new InputStreamReader(System.in));
        } else   // read from file
        {
            String fileName = args[0];
            br = new BufferedReader(new FileReader(fileName));
        }

        String record;
        while ((record = br.readLine()) != null) {
            if (record.length() == 0 || record.startsWith("//") || record.startsWith("#")) {
                continue; //! ignore comments  and empty lines
            } else if (record.startsWith("."))  // end of input
                break;
            String[] clusterData = record.split(" ");
            Vector<Integer> cluster = addCluster(clusterData);
            clusterSet.addCluster(cluster);
        }

        for (int a = 1; a < args.length; a++) {
            String option1 = args[a];
            if (option1.equals("--printclusters")) {
                printClusters = true;
            } else if (option1.equals("--colouredges")) {
                colourEdges = true;
            } else if (option1.equals("--printretnum")) {
                printRetNum = true;
            } else if (option1.startsWith("--timeout=")) {
                timeout = 60 * Integer.parseInt(option1.split("=")[1]);
            } else {
                System.out.println("unknown option: " + option1);
            }
        }
        System.out.println(clusterSet.clusterVec);
        Vector<Vector<Vector<Integer>>> ccs = getConnectedComponents(clusterSet);

        // Sort the connected components from small to large
        ccs.sort(new Comparator<Vector<Vector<Integer>>>() {
            @Override
            public int compare(Vector<Vector<Integer>> cc1, Vector<Vector<Integer>> cc2) {
                // Sort based on the number of clusters per component
                if (cc1.size() < cc2.size()) {
                    return -1;
                } else if (cc1.size() > cc2.size()) {
                    return 1;
                }
                return 0;
            }
        });
//        printConnectedComponents(ccs);
//        System.exit(3);

        int k = 0;
        int r = 0;
        long start = System.currentTimeMillis();
        final long endTime = timeout > 0 ? start + 1000L * timeout : Long.MAX_VALUE;
        try {
            for (var cc : ccs) {  // For each connected component
                var CS = new ClusterSet();
                for (var cl : cc) {
                    CS.addCluster(cl);
                }

                // Run Cass for this connected component (defined in clusterSet)
                Vector<DiGraph> networks = minSL(CS, endTime, isOnlyOne, checkTrees);
                if (networks.size() > 0) {
                    var diGraph = (networks.elementAt(0));
                    diGraph.printDiGraph(printClusters, colourEdges);
                    if (printRetNum) {
                        System.out.println("cc: k=r=" + diGraph.level);
                    }
                    r += diGraph.level;
                    k = max(k, diGraph.level);
                } else {
                    System.err.println("No network found.");
                }
            }
        } catch (TimeoutException e) {
            System.err.println(e);
            k = max(k, K_REACHED);
            r += K_REACHED;
        }
        long tElapsed = System.currentTimeMillis() - start;
        System.out.println("t=" + (float) tElapsed / 1000);
        if (printRetNum) {
            System.out.println("k=" + k + ",r=" + r);
        }
    }

    private static Vector<Integer> addCluster(String[] cluster) {
        Vector<Integer> vecCluster = new Vector<>();
        for (int i = 0; i < cluster.length; i++) {
            Integer taxon;
            if (cluster[i].length() > CassAlgorithm.MAX_CIRCLE_SIZE) {
                CassAlgorithm.leafShape = "box";
            }
            if (stringTaxa.contains(cluster[i])) {
                taxon = Integer.valueOf(stringTaxa.indexOf(cluster[i]) + 1);
            } else {
                taxon = Integer.valueOf(stringTaxa.size() + 1);
                stringTaxa.add(cluster[i]);
            }
            vecCluster.add(taxon);
        }
        return vecCluster;
    }

//    public static Vector minSL(ClusterSet CS, ProgressListener progressListener, boolean isOnlyOne, boolean checkTrees) throws CanceledException {

    private static Vector<Vector<Vector<Integer>>> getConnectedComponents(ClusterSet CS) {
        Vector<Vector<Integer>> clusters = (Vector<Vector<Integer>>) CS.clusterVec.clone();
        Vector<Vector<Vector<Integer>>> components = new Vector<>();
        Vector<Vector<Integer>> comp = new Vector<>();
        // Remove the singleton clusters
        for (var c : CS.clusterVec) {
            if (c.size() == 1) {
                clusters.remove(c);
            }
        }
        Vector<Integer> cluster1 = clusters.firstElement();
        comp.add((Vector<Integer>) cluster1.clone());
        int clInd = 0;
        while (clusters.size() > 0) {
            clusters.remove(cluster1);
//            System.out.println(cluster1);
            for (var cluster2 : clusters) {
                int rel = ClusterSet.getRelation(cluster1, cluster2);
                if (rel == 2) {  // Incompatible
//                    System.out.println(cluster1 + " +-- " + cluster2);
                    if (!comp.contains(cluster2)) {
                        comp.add((Vector<Integer>) cluster2.clone());
                    }
                }
            }
            // Choose the next cluster
            clInd++;
            if (comp.size() >= clInd + 1) {
                cluster1 = comp.elementAt(clInd);
            } else {
                // A new connected component will be addressed. Note: clusters might be empty now
                if (comp.size() >= 2) {
//                    System.out.println("Adding comp. " + comp);
                    components.add((Vector<Vector<Integer>>) comp.clone());
//                } else {
//                    System.out.println("Neglect singleton " + comp);
                }
                clInd = 0;
                comp.clear();
                if (clusters.size() > 0) {
                    cluster1 = clusters.firstElement();
                    comp.add((Vector<Integer>) cluster1.clone());
                }
            }
        }

        return components;
    }

    private static void printConnectedComponents(Vector<Vector<Vector<Integer>>> components) {
        System.out.println(components.size() + " conn. comps. on the " + stringTaxa.size() + " taxa: " + stringTaxa);
        for (var comp : components) {
            Set taxa = new HashSet();
            for (Vector<Integer> cl : comp) {
                taxa.addAll(cl.stream().map(i -> stringTaxa.get(i - 1)).collect(Collectors.toList()));
            }
            System.out.println(comp.size() + " clusters on the " + taxa.size() + " taxa " + taxa.stream().sorted().collect(Collectors.toList()) + ":");
            for (var cl : comp) {
                System.out.println(cl.stream().map(i -> stringTaxa.get(i - 1)).collect(Collectors.toList()));
            }
            System.out.println();
        }
    }

    /**
     * Cass
     *
     * @param CS
     * @param endTime    End time in ms
     * @param isOnlyOne
     * @param checkTrees
     * @return
     * @throws TimeoutException
     */
    public static Vector<DiGraph> minSL(ClusterSet CS, long endTime, boolean isOnlyOne, boolean checkTrees) throws TimeoutException {
        int k = 0;
        FOUND = 0;
        Vector<DiGraph> output = new Vector<>();
        boolean found = false;

        // System.err.println("Collapsing");
        // collapse all maximal ST-sets
        Vector op = CS.collapse();
        ClusterSet collapsed = (ClusterSet) op.elementAt(0);
        Vector collapsedTaxa = (Vector) op.elementAt(1);
        Vector collapsedTaxaSets = (Vector) op.elementAt(2);

        while (!found) {
            System.err.println("Searching level " + k);
//            progressListener.incrementProgress();  // HANS
//            Vector networks = SL(collapsed, k, progressListener, isOnlyOne, checkTrees);
            Vector<DiGraph> networks = SL(collapsed, k, endTime, isOnlyOne, checkTrees);
            for (int j = 0; j < networks.size(); j++) {
                DiGraph colG = networks.elementAt(j);
                DiGraph G = uncollapse(colG, CS, collapsedTaxa, collapsedTaxaSets);
                G.cleanDiGraph();

                // check if network displays all clusters/trees, should not be necessary
                // this has to be done before postprocessing, since postprocessing does not update all attributes
                // boolean disp = G.displays(CS, k, checkTrees);

                // remove dummy leaves
                G.remDummyLeaves();
                G.cleanDiGraph();

                // contract edges if possible
                for (int i = 0; i < k; i++) {
                    // might have to repeat this a couple of times before everything has been cleaned up
                    G = G.postprocess();
                    G.cleanDiGraph();
                }
                G.level = k;

                //if(disp) {
                found = true;
                output.add(G);
                //} else {
                //    System.err.println("Incorrect network returned");
                //}
            }
            k++;
        }

        return output;
    }

    //    private static Vector SL(ClusterSet CS, int k, ProgressListener progressListener, boolean isOnlyOne, boolean checkTrees) throws CanceledException {
    private static Vector<DiGraph> SL(ClusterSet CS, int k, long endTime, boolean isOnlyOne, boolean checkTrees) throws TimeoutException {
        // int count = 0;

        // in this vector we will put all valid networks
        Vector<DiGraph> networks = new Vector<DiGraph>();

        if (k == 0) {
            // we try to create a tree
            // construct a tree containing precisely those clusters in CSp
            DiGraph T = buildTree(CS);
            if (T != null) {
                networks.add(T);
            }
            return networks;
        }

        Stack<StackObject> stack = new Stack<StackObject>();

        StackObject so = new StackObject(k);
        so.taxa = CS.taxa;
        so.type = FASE1;
        so.CS[0] = CS;
        so.collapsedTaxa[0] = null;
        so.collapsedTaxaSets[0] = null;
        so.removedLeaves[0] = null;

        stack.push(so);

        while (stack.size() != 0) {

            StackObject top = stack.peek();
            stack.pop();

            if (top.type == FASE1) {
                if (top.step == k) {
                    // build a tree and start the second phase
                    // construct a tree containing precisely those clusters in top.CS[k]
                    DiGraph T = buildTree(top.CS[k]);
                    if (T != null) {
                        // add the tree to the stackobject
                        // change the fase
                        top.type = FASE2;
                        top.network = T;
                        // put the stackobject back on the stack
                        stack.push(top);
                    }
                    continue;
                }

                // collapse all maximal ST-sets
                // IN THE FIRST STEP THERE IS NO NEED TO COLLAPSE
                Vector output = top.CS[top.step].collapse();
                ClusterSet collapsed = (ClusterSet) output.elementAt(0);
                Vector collapsedTaxa = (Vector) output.elementAt(1);
                Vector collapsedTaxaSets = (Vector) output.elementAt(2);

                // loop through all taxa
                for (int i = 0; i < collapsed.taxa.size() + 1; i++) {
                    ClusterSet CSp;
                    Integer x;
                    if (i == collapsed.taxa.size()) {
                        if (top.step == 0) {
                            // always remove a taxon in the first step
                            break;
                        }
                        // don't remove a leaf
                        // do not collapse
                        CSp = top.CS[top.step];
                        x = Integer.valueOf(CassAlgorithm.DUMMY_NUMBER);
                    } else {
                        // remove the taxon with index i
                        x = collapsed.taxa.elementAt(i);
                        CSp = collapsed.remLeaf(x);
                    }
                    // create a new stack object
                    StackObject virgin = new StackObject(k);
                    virgin.taxa = top.taxa;
                    virgin.type = FASE1;
                    virgin.step = top.step + 1;
                    for (int c = 0; c < virgin.step; c++) {
                        virgin.CS[c] = top.CS[c];
                        virgin.collapsedTaxa[c] = top.collapsedTaxa[c];
                        virgin.collapsedTaxaSets[c] = top.collapsedTaxaSets[c];
                        virgin.removedLeaves[c] = top.removedLeaves[c];
                    }
                    virgin.CS[virgin.step] = CSp;
                    virgin.collapsedTaxa[virgin.step] = collapsedTaxa;
                    virgin.collapsedTaxaSets[virgin.step] = collapsedTaxaSets;
                    virgin.removedLeaves[virgin.step] = x;
                    stack.push(virgin);
                }
            } else {
                // FASE 2
                if (top.step == 0) {
                    // a solution has been found

                    // we remove the root only if it has outdegree 1 and after removal the network still displays all clusters
                    // NOTE: if we keep an outdegree-1 root, the output network might not be simple, but I don't think that matters
                    if (top.network.outdeg == 1) {
                        if (((DiGraph) top.network.children.elementAt(0)).displays(top.CS[0], k, checkTrees)) {
                            top.network = ((DiGraph) top.network.children.elementAt(0));
                            top.network.indeg = 0;
                        }
                    }

                    // output this solution
                    networks.add(top.network);
                    FOUND++;
                    if (FOUND == 1) {
                        System.err.println("Network found");
                    }
                    if (isOnlyOne) {
                        return networks;
                    } else {
                        if (FOUND == 1) {
                            System.err.println("Starting search for alternative solutions");
                        }
                        continue;
                    }
                }
                // number all edges
                int[] num = new int[1];
                num[0] = 0;
                int e = top.network.numberEdges(num);
                top.network.cleanDiGraph();
                Integer x = top.removedLeaves[top.step];
                // if there are 2 cherries, we should only really sub cherries, unless the last leaf we removed was a dummy leaf, because then we haven't collapsed properly
                // if there are more cherries we can stop, unless the last leaf we removed was a dummy leaf, because then we haven't collapsed properly
                // if there's 1 cherry, we should sub that one and one arb other edge, unless the last leaf we removed was a dummy leaf, because then we haven't collapsed properly
                Vector cherries = top.network.findCherries();
                top.network.cleanDiGraph();
                if (cherries.size() > 2 && x.intValue() != CassAlgorithm.DUMMY_NUMBER) {
                    continue;
                }
                Vector cherry1 = new Vector();
                if (cherries.size() > 0) {
                    cherry1 = (Vector) cherries.elementAt(0);
                }
                Vector cherry2 = new Vector();
                if (cherries.size() > 1) {
                    cherry2 = (Vector) cherries.elementAt(1);
                }
                // loop through each pair of edges
                for (int e1 = 0; e1 < e; e1++) {
                    Integer edge1 = Integer.valueOf(e1);
                    // skip this edge if it is not a cherry-edge and there are two cherries
                    if (cherries.size() == 2 && !cherry1.contains(edge1) && !cherry2.contains(edge1) && x.intValue() != CassAlgorithm.DUMMY_NUMBER) {
                        continue;
                    }
                    for (int e2 = e1; e2 < e; e2++) {
                        Integer edge2 = Integer.valueOf(e2);
                        // if there is a cherry we have to hit it
                        if (cherries.size() > 0 && !cherry1.contains(edge1) && !cherry1.contains(edge2) && x.intValue() != CassAlgorithm.DUMMY_NUMBER) {
                            continue;
                        }
                        // if there are two cherries we have to hit both of them
                        if (cherries.size() == 2 && !cherry2.contains(edge1) && !cherry2.contains(edge2) && x.intValue() != CassAlgorithm.DUMMY_NUMBER) {
                            continue;
                        }
                        // if e1 = e2 we subdivide this edge twice
                        // copy top.network into a new DiGraph G
                        DiGraph G = top.network.cloneDiGraph();
                        top.network.cleanDiGraph();
                        // create a new reticulation
                        DiGraph r = new DiGraph();
                        r.indeg = 2;
                        // hang this reticulation below e1 and e2
                        G.hangBelowEdges(r, e1, e2);
                        G.cleanDiGraph();
                        // create a new leaf x below this reticulation
                        r.children.add(new DiGraph(x));
                        ((DiGraph) r.children.elementAt(0)).indeg = 1;
                        r.outdeg = 1;
                        // if G has cherries at this point we can forget about it
                        // but it shouldn't really have cherries if we're carefull
                        if (x.intValue() != CassAlgorithm.DUMMY_NUMBER) {
                            // ONLY UNCOLAPSE IF WE HAVEN'T ADDED A DUMMY LEAF
                            // OTHERWISE THE UNCOLLAPSE FUNCTION WON'T WORK PROPERLY
                            uncollapse(G, top.CS[top.step - 1], top.collapsedTaxa[top.step], top.collapsedTaxaSets[top.step]);
                        }
                        // check if the network displays all clusters
                        // (if we just added a dummy leaf, there is no need to check this)
                        if (x.intValue() == DUMMY_NUMBER || G.displays(top.CS[top.step - 1], k, checkTrees)) {
                            // create new stackobject
                            StackObject virgin = new StackObject(k);
                            virgin.taxa = top.taxa;
                            virgin.step = top.step - 1;
                            virgin.type = FASE2;
                            for (int c = 0; c < top.step; c++) {
                                virgin.CS[c] = top.CS[c];
                                virgin.collapsedTaxa[c] = top.collapsedTaxa[c];
                                virgin.collapsedTaxaSets[c] = top.collapsedTaxaSets[c];
                                virgin.removedLeaves[c] = top.removedLeaves[c];
                            }
                            virgin.network = G;
                            stack.push(virgin);
                        }
                    }
                    //if (progressListener != null)
                    //    progressListener.checkForCancel();
                }
            }
            /*
            if (count < 80) {
                System.err.print(".");
                count++;
            } else {
                System.err.println();
                count = 0;
            }
            */

//            if (progressListener != null)
//                progressListener.checkForCancel();
            if (System.currentTimeMillis() >= endTime) {
                System.out.println("Reached k=" + k);
                K_REACHED = k;
                throw new TimeoutException();
            }
        }
        return networks;
    }

    private static DiGraph buildTree(ClusterSet CS) {
        // find the maximal clusters
        Vector<Vector<Integer>> maxClusters = new Vector<>(0);
        for (int i = 0; i < CS.clusterVec.size(); i++) {
            Vector<Integer> cluster1 = CS.clusterVec.elementAt(i);
            boolean max = true;
            for (int j = 0; j < CS.clusterVec.size(); j++) {
                Vector<Integer> cluster2 = CS.clusterVec.elementAt(j);
                int rel = ClusterSet.getRelation(cluster1, cluster2);
                if (rel == 3) {
                    max = false;
                }
                if (rel == 2) {
                    return null;
                }
            }
            if (max) {
                maxClusters.add(cluster1);
            }
        }

        // create a vertex
        DiGraph G = new DiGraph(maxClusters.size());

        // the maximal clusters become children of the root
        for (int c = 0; c < maxClusters.size(); c++) {
            Vector<Integer> cluster = maxClusters.elementAt(c);
            if (cluster.size() == 1) {
                Integer label = cluster.elementAt(0);
                G.outdeg++;
                DiGraph leaf = new DiGraph(label);
                leaf.indeg = 1;
                G.children.add(leaf);
            } else {
                ClusterSet restricted = CS.restrict(cluster);
                G.outdeg++;
                DiGraph H = buildTree(restricted);
                H.indeg = 1;
                if (H != null) {
                    G.children.add(H);
                } else {
                    return null;
                }
            }
        }
        return G;
    }

    private static DiGraph uncollapse(DiGraph D, ClusterSet CS, Vector collapsedTaxa, Vector collapsedTaxaSets) {
        // we only uncollapse the last step!
        for (int i = 0; i < collapsedTaxa.size(); i++) {
            Integer x = (Integer) collapsedTaxa.elementAt(i);
            Vector cluster = (Vector) collapsedTaxaSets.elementAt(i);
            if (cluster.size() > 1) {
                // uncollapse this taxon
                ClusterSet restrictedCS = CS.restrict(cluster);
                DiGraph G = buildTree(restrictedCS);
                // replace the leaf by the subtree
                replaceLeafBySubtree(x, D, G);
                D.cleanDiGraph();
            }
        }
        return D;
    }

    private static DiGraph replaceLeafBySubtree(int x, DiGraph D, DiGraph G) {
        D.visited = true;
        for (int c = 0; c < D.outdeg; c++) {
            DiGraph child = ((DiGraph) D.children.elementAt(c));
            if (child.label.intValue() == x) {
                if (G.outdeg == 1) {
                    D.children.setElementAt(G.children.elementAt(0), c);
                } else {
                    D.children.setElementAt(G, c);
                }
            } else {
                if (!((DiGraph) D.children.elementAt(c)).visited) {
                    D.children.setElementAt(replaceLeafBySubtree(x, (DiGraph) D.children.elementAt(c), G), c);
                }
            }
        }
        return D;
    }

    static class StackObject {
        // the StackObject is a snapshot of the data in either fase 1 (collapse, remove leaf steps)
        //                                              or fase 2 (hang leaf back, decollapse steps)

        Vector<Integer> taxa; // all the taxa.... don't really need this
        int step; // the step we're at in the current fase
        int type; // either fase 1 or fase 2
        ClusterSet[] CS; // the cluster sets after 0,1,2,... "collapse, remove leaf"-steps
        Vector[] collapsedTaxaSets; // the sets of taxa we collapsed in step 0,1,2,...
        Vector[] collapsedTaxa; // the taxa we collapsed them into in step 0,1,2,...
        Integer[] removedLeaves; // the leaf we've removed in step 0,1,2,...

        // only for fase 2
        DiGraph network; // the solution in the current step of the second fase

        public StackObject(int k) {
            this.step = 0;
            CS = new ClusterSet[k + 1];
            collapsedTaxa = new Vector[k + 1];
            collapsedTaxaSets = new Vector[k + 1];
            removedLeaves = new Integer[k + 1];
        }
    }
}
