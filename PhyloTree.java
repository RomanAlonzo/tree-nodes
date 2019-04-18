/*
 * PhyloTree.java
 *
 * Defines a phylogenetic tree, which is a strictly binary tree
 * that represents inferred hierarchical relationships between species
 *
 * There are weights along each edge; the weight from parent to left child
 * is the same as parent to right child.
 *
 * Students may only use functionality provided in the packages
 *     java.lang
 *     java.util
 *     java.io
 *
 * Use of any additional Java Class Library components is not permitted
 *
 * Roman Alonzo
 *
 */

import java.lang.*;
import java.util.*;
import java.io.*;

public class PhyloTree {
    private PhyloTreeNode overallRoot;    // The actual root of the overall tree
    private int printingDepth;            // How many spaces to indent the deepest 
    // node when printing
    private int numSpecies = 0;
    private static ArrayList<Species> descendents = new ArrayList<Species>();


    // CONSTRUCTOR

    // PhyloTree
    // Pre-conditions:
    //        - speciesFile contains the path of a valid FASTA input file
    //        - printingDepth is a positive number
    // Post-conditions:
    //        - this.printingDepth has been set to printingDepth
    //        - A linked tree structure representing the inferred hierarchical
    //          species relationship has been created, and overallRoot points to
    //          the root of this tree
    // Notes:
    //        - A lot happens in this step!  See assignment description for details
    //          on the input format file and how to construct the tree
    //        - If you encounter a FileNotFoundException, print to standard error
    //          "Error: Unable to open file " + speciesFile
    //          and exit with status (return code) 1
    //    - Most of this should be accomplished by calls to loadSpeciesFile and buildTree
    public PhyloTree(String speciesFile, int printingDepth) {
        Scanner input;
        try {
            input = new Scanner(new File(speciesFile));
            this.numSpecies = 0;
            while (input.hasNext()) {
                if (input.next().contains(">")) {
                    numSpecies++;
                }
            }
            Species[] tree = loadSpeciesFile(speciesFile);
            buildTree(tree);
            descendents = new ArrayList<Species>(Arrays.asList(tree));
        } catch (FileNotFoundException e) {
            System.out.println("Error unable to open file " + speciesFile);
            System.exit(1);
        }
        this.printingDepth = printingDepth;
        return;
    }
    // ACCESSORS

    // getOverallRoot
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the overall root
    public PhyloTreeNode getOverallRoot() {
        return this.overallRoot;
    }

    // toString 
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns a string representation of the tree
    // Notes:
    //    - See assignment description for proper format
    //        (it will be a kind of reverse in-order [RNL] traversal)
    //    - Can be a simple wrapper around the following toString
    //    - Hint: StringBuilder is much faster than repeated concatenation
    public String toString() {
        return toString(this.getOverallRoot(), this.getWeightedHeight(), this.getWeightedHeight());
    }

    // toString 
    // Pre-conditions:
    //    - node points to the root of a tree you intend to print
    //    - weightedDepth is the sum of the edge weights from the
    //      overall root to the current root
    //    - maxDepth is the weighted depth of the overall tree
    // Post-conditions:
    //    - Returns a string representation of the tree
    // Notes:
    //    - See assignment description for proper format
    private String toString(PhyloTreeNode node, double weightedDepth, double maxDepth) {
        StringBuilder concat = new StringBuilder();
        String empty = "";
        if (node != null) {
            int k = (int) (this.printingDepth * (weightedDepth / maxDepth));
            if (!node.isLeaf()) {
                concat.append(toString(node.getRightChild(), weightedNodeDepth(node.getRightChild()), maxDepth));
                for (int i = 0; i < k; i++) {
                    concat.append(".");
                }
                empty = String.format("%.2f", node.getDistanceToChild());
                concat.append("[NONTERM " + empty + "]" + "\n");
                concat.append(toString(node.getLeftChild(), weightedNodeDepth(node.getLeftChild()), maxDepth));
            } else {
                for (int i = 0; i < k; i++) {
                    concat.append(".");
                }
                concat.append(node.getLabel() + "\n");
            }
        }
        return concat.toString();
    }

    //helper function for toString
    //handles negative infinity and recurses to find depth (if parent)
    public static double weightedNodeDepth(PhyloTreeNode node) {
        if (node == null) {
            return Double.NEGATIVE_INFINITY;
        }
        if (node.getParent() != null) {
            return node.getParent().getDistanceToChild() + weightedNodeDepth(node.getParent());
        } else {
            return 0.0;
        }
    }

    // toTreeString 
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns a string representation in tree format
    // Notes:
    //    - See assignment description for format details
    //    - Can be a simple wrapper around the following toTreeString
    public String toTreeString() {
        return toTreeString(this.getOverallRoot());
    }

    // toTreeString 
    // Pre-conditions:
    //    - node points to the root of a tree you intend to print
    // Post-conditions:
    //    - Returns a string representation in tree format
    // Notes:
    //    - See assignment description for proper format
    private String toTreeString(PhyloTreeNode node) {
        StringBuilder concat = new StringBuilder();
        String extra = "";
        if (node.isLeaf()) {
            extra = String.format("%.5f", node.getParent().getDistanceToChild());
            concat.append(node.getLabel() + ":" + extra);
        } else {
            concat.append("(");
            concat.append(toTreeString(node.getRightChild()));
            concat.append(",");
            concat.append(toTreeString(node.getLeftChild()));
            if (node == this.getOverallRoot()) {
                concat.append(")");
            } else {
                extra = String.format("%.5f", node.getDistanceToChild());
                concat.append("):" + extra);
            }
        }
        return concat.toString();
    }

    // getHeight
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the tree height as defined in class
    // Notes:
    //    - Can be a simple wrapper on nodeHeight
    public int getHeight() {
        return nodeHeight(this.getOverallRoot());
    }

    // getWeightedHeight
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the sum of the edge weights along the
    //      "longest" (highest weight) path from the root
    //      to any leaf node.
    // Notes:
    //   - Can be a simple wrapper for weightedNodeHeight
    public double getWeightedHeight() {
        return weightedNodeHeight(this.getOverallRoot());
    }

    // countAllSpecies
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the number of species in the tree
    // Notes:
    //    - Non-terminals do not represent species
    //    - This functionality is provided for you elsewhere
    //      just call the appropriate method
    public int countAllSpecies() {
        return this.getOverallRoot().getNumLeafs();
    }

    // getAllSpecies
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns an ArrayList containing all species in the tree
    // Notes:
    //    - Non-terminals do not represent species
    // Hint:
    //    - Call getAllDescendantSpecies
    public java.util.ArrayList<Species> getAllSpecies() {
        return this.descendents;
    }

    // findTreeNodeByLabel
    // Pre-conditions:
    //    - label is the label of a tree node you intend to find
    //    - Assumes labels are unique in the tree
    // Post-conditions:
    //    - If found: returns the PhyloTreeNode with the specified label
    //    - If not found: returns null
    public PhyloTreeNode findTreeNodeByLabel(String label) {
        return findTreeNodeByLabel(this.getOverallRoot(), label);
    }

    // findLeastCommonAncestor
    // Pre-conditions:
    //    - label1 and label2 are the labels of two species in the tree
    // Post-conditions:
    //    - If either node cannot be found: returns null
    //    - If both nodes can be found: returns the PhyloTreeNode of their
    //      common ancestor with the largest depth
    //      Put another way, the least common ancestor of nodes A and B
    //      is the only node in the tree where A is in the left tree
    //      and B is in the right tree (or vice-versa)
    // Notes:
    //    - Can be a wrapper around the static findLeastCommonAncestor
    public PhyloTreeNode findLeastCommonAncestor(String label1, String label2) {
        return findLeastCommonAncestor(findTreeNodeByLabel(label1), findTreeNodeByLabel(label2));
    }

    // findEvolutionaryDistance
    // Pre-conditions:
    //    - label1 and label2 are the labels of two species in the tree
    // Post-conditions:
    //    - If either node cannot be found: returns POSITIVE_INFINITY
    //    - If both nodes can be found: returns the sum of the weights 
    //      along the paths from their least common ancestor to each of
    //      the two nodes
    public double findEvolutionaryDistance(String label1, String label2) {
        PhyloTreeNode node1 = findTreeNodeByLabel(label1);
        PhyloTreeNode node2 = findTreeNodeByLabel(label2);
        if (node1 == null || node2 == null) {
            return Double.POSITIVE_INFINITY;
        } else {
            double sum = 0.0;
            while (nodeDepth(node1) != nodeDepth(node2)) {
                if (nodeDepth(node1) > nodeDepth(node2)) {
                    node1 = node1.getParent();
                    sum += node1.getDistanceToChild();
                } else {
                    node2 = node2.getParent();
                    sum += node2.getDistanceToChild();
                }
            }
            while (node2 != node1) {
                node1 = node1.getParent();
                sum += node1.getDistanceToChild();
                node2 = node2.getParent();
                sum += node2.getDistanceToChild();
            }
            return sum;
        }
    }

    // MODIFIER

    // buildTree
    // Pre-conditions:
    //    - species contains the set of species for which you want to infer
    //      a phylogenetic tree
    // Post-conditions:
    //    - A linked tree structure representing the inferred hierarchical
    //      species relationship has been created, and overallRoot points to
    //      the root of said tree
    // Notes:
    //    - A lot happens in this step!  See assignment description for details
    //      on how to construct the tree.
    //    - Be sure to use the tie-breaking conventions described in the pdf
    //    - Important hint: although the distances are defined recursively, you
    //      do NOT want to implement them recursively, as that would be very inefficient

    //make a forest of leaves
    //precompute and store all pairwise distances using species distance

    private void buildTree(Species[] species) {
        HashMap<String, PhyloTreeNode> forest = new HashMap<String, PhyloTreeNode>();
        MultiKeyMap distance = new MultiKeyMap();
        for (int i = 0; i < species.length; i++) {
            PhyloTreeNode node = new PhyloTreeNode(null, species[i]);
            forest.put(species[i].getName(), node);
        }
        for (int i = 0; i < species.length; i++) {
            for (int j = 0; j < species.length; j++) {
                distance.put(species[i].getName(), species[j].getName(), Species.distance(species[i], species[j]));
            }
        }
        double shortDistance;
        double distanceVar1 = 0.0;
        double distanceVar2 = 0.0;
        String space1 = "";
        String space2 = "";
        PhyloTreeNode node1 = null;
        PhyloTreeNode node2 = null;
        PhyloTreeNode parentNode = null;
        String parentLabel = "";
        while (forest.size() > 1) {
            shortDistance = Double.MAX_VALUE;
            for (Map.Entry<String, PhyloTreeNode> firstArg : forest.entrySet()) {
                for (Map.Entry<String, PhyloTreeNode> secondArg : forest.entrySet()) {
                    if (distance.get(firstArg.getKey(), secondArg.getKey()) != null) {
                        if (distance.get(firstArg.getKey(), secondArg.getKey()) < shortDistance && !firstArg.getKey().equals(secondArg.getKey())) {
                            shortDistance = distance.get(firstArg.getKey(), secondArg.getKey());
                            space1 = firstArg.getKey();
                            space2 = secondArg.getKey();
                            node1 = firstArg.getValue();
                            node2 = secondArg.getValue();
                        }
                    }
                }
            }
            forest.remove(space1);
            forest.remove(space2);
            if (space1.compareTo(space2) <= 0) {
                parentLabel = space1 + "+" + space2;
                parentNode = new PhyloTreeNode(parentLabel, null, node1, node2, shortDistance / 2);
            } else {
                parentLabel = space2 + "+" + space1;
                parentNode = new PhyloTreeNode(parentLabel, null, node2, node1, shortDistance / 2);
            }
            node1.setParent(parentNode);
            node2.setParent(parentNode);
            for (Map.Entry<String, PhyloTreeNode> thirdArg : forest.entrySet()) {
                if (!thirdArg.getKey().equals(parentLabel)) {
                    distanceVar1 = (((double) node1.getNumLeafs()) / ((double) (node2.getNumLeafs()) + ((double) node1.getNumLeafs()))) * distance.get(space1, thirdArg.getKey());
                    distanceVar2 = (((double) node2.getNumLeafs()) / ((double) (node2.getNumLeafs()) + ((double) node1.getNumLeafs()))) * distance.get(space2, thirdArg.getKey());
                    distance.put(parentLabel, thirdArg.getKey(), distanceVar1 + distanceVar2);
                }
            }
            forest.put(parentLabel, parentNode);
        }
        this.overallRoot = parentNode;
        return;
    }

    // STATIC

    // nodeDepth
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns -1
    //    - Else: returns the depth of the node within the overall tree
    public static int nodeDepth(PhyloTreeNode node) {
        if (node == null) {
            return -1;
        } else {
            if (node.getParent() == null) {
                return 0;
            } else {
                return 1 + nodeDepth(node.getParent());
            }
        }
    }


    // nodeHeight
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns -1
    //    - Else: returns the height subtree rooted at node
    public static int nodeHeight(PhyloTreeNode node) {
        if (node == null) {
            return -1;
        } else {
            if (node.isLeaf()) {
                return 0;
            } else {
                return 1 + Math.max(nodeHeight(node.getLeftChild()), nodeHeight(node.getRightChild()));
            }
        }
    }

    // weightedNodeHeight 
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns NEGATIVE_INFINITY
    //    - Else: returns the weighted height subtree rooted at node
    //     (i.e. the sum of the largest weight path from node
    //     to a leaf; this might NOT be the same as the sum of the weights
    //     along the longest path from the node to a leaf)
    public static double weightedNodeHeight(PhyloTreeNode node) {
        if (node != null) {
            if (!node.isLeaf()) {
                double rightSide = node.getDistanceToChild() + weightedNodeHeight(node.getRightChild());
                double leftSide = node.getDistanceToChild() + weightedNodeHeight(node.getLeftChild());
                if (rightSide >= leftSide) {
                    return rightSide;
                } else {
                    return leftSide;
                }
            } else {
                return 0.0;
            }
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

    // loadSpeciesFile
    // Pre-conditions:
    //    - filename contains the path of a valid FASTA input file
    // Post-conditions:
    //    - Creates and returns an array of species objects representing
    //      all valid species in the input file
    // Notes:
    //    - Species without names are skipped
    //    - See assignment description for details on the FASTA format
    // Hints:
    //    - Because the bar character ("|") denotes OR, you need to escape it
    //      if you want to use it to split a string, i.e. you can use "\\|" 
    public static Species[] loadSpeciesFile(String filename) {
        ArrayList<Species> array = new ArrayList<Species>();
        String[] nameArray = new String[7];
        String line = "";
        String name = "";
        String sequence = "";
        int inputLines = 0;

        try {
            Scanner input = new Scanner(new File(filename));
            line = input.next();
            inputLines++;
            while (input.hasNext()) {
                nameArray = line.split("\\|");
                name = nameArray[nameArray.length - 1];
                line = input.next();
                inputLines++;
                while (!line.contains(">") && input.hasNext()) {
                    sequence = sequence + line;
                    line = input.next();
                    inputLines++;
                }
                if (!input.hasNext()) {
                    sequence = sequence + line;
                }
                String[] finalSequence = new String[sequence.length()];
                finalSequence = sequence.split("");
                array.add(new Species(name, finalSequence));
            }
        } catch (FileNotFoundException e) {
            System.out.println("File not found: " + filename);
            System.exit(1);
        }
        Species[] specArray = array.toArray(new Species[array.size()]);
        descendents = array;
        String[] finalString = new String[specArray[5].getSequence().length];
        finalString = specArray[5].getSequence();
        return specArray;
    }

    // getAllDescendantSpecies
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - descendants is a non-null reference variable to an empty arraylist object
    // Post-conditions:
    //    - descendants is populated with all species in the subtree rooted at node
    //      in in-/pre-/post-order (they are equivalent here)
    private static void getAllDescendantSpecies(PhyloTreeNode node, java.util.ArrayList<Species> descendants) {
        for (int i = 0; i < descendents.size(); i++) {
            findTreeNodeByLabel(descendents.get(i).getName());
            if (!node.isLeaf()) {
                descendents.remove(descendents.get(i));
            }
        }
    }

    // findTreeNodeByLabel
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - label is the label of a tree node that you intend to locate
    // Post-conditions:
    //    - If no node with the label exists in the subtree, return null
    //    - Else: return the PhyloTreeNode with the specified label 
    // Notes:
    //    - Assumes labels are unique in the tree
    private static PhyloTreeNode findTreeNodeByLabel(PhyloTreeNode node, String label) {
        if (node.getLabel().equals(label)) {
            return node;
        }
        if (node.isLeaf()) {
            return null;
        } else {
            PhyloTreeNode right = findTreeNodeByLabel(node.getRightChild(), label);
            PhyloTreeNode left = findTreeNodeByLabel(node.getLeftChild(), label);
            if (left != null) {
                return left;
            } else {
                return right;
            }
        }
    }

    // findLeastCommonAncestor
    // Pre-conditions:
    //    - node1 and node2 point to nodes in the phylogenetic tree
    // Post-conditions:
    //    - If node1 or node2 are null, return null
    //    - Else: returns the PhyloTreeNode of their common ancestor 
    //      with the largest depth
    private static PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode node1, PhyloTreeNode node2) {
        if (node1 == null || node2 == null) {
            return null;
        } else {
            if (nodeDepth(node1) > nodeDepth(node2)) {
                return findLeastCommonAncestor(node1.getParent(), node2);
            } else if (nodeDepth(node1) < nodeDepth(node2)) {
                return findLeastCommonAncestor(node2.getParent(), node1);
            } else if (node1 == node2) {
                return node1;
            } else {
                return findLeastCommonAncestor(node1.getParent(), node2.getParent());
            }
        }
    }
}