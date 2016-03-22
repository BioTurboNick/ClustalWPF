using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.Tree
{
    class GuideTree<T>
    // Implements a rooted or unrooted tree, with methods for generating a tree from a distance matrix,
    // rooting and rerooting the tree, calcuating weights, and calculating steps for processing the objects
    // related by the tree.
    {
        public bool isRooted;
        TreeNode root; // Holds the root node for a rooted tree.
        List<TreeNode> nodesList = new List<TreeNode>();
        List<TreeNode> terminalNodesList = new List<TreeNode>();
        List<TreeLeaf> leavesList = new List<TreeLeaf>(); // stores the leaves in the same order as the matrices
        public double[,] similarityMatrix;
        List<Tuple<ReadOnlyCollection<T>, ReadOnlyCollection<T>, double>> steps = new List<Tuple<ReadOnlyCollection<T>, ReadOnlyCollection<T>, double>>(); // Stores the ordered arrangement of steps from terminal nodes up to the root, and the average similarity between them.

        public ReadOnlyCollection<Tuple<ReadOnlyCollection<T>, ReadOnlyCollection<T>, double>> Steps
        {
            get { return steps.AsReadOnly(); }
        }

        public ReadOnlyCollection<TreeLeaf> LeavesList
        {
            get { return leavesList.AsReadOnly(); }
        }

        public int leavesCount = 0;

        public TreeNode Root
        {
            get { return root; }
        }

        public static GuideTree<T> GetWeights(ref double[,] distanceMatrix, T[] leafDataList)
        {
            GuideTree<T> tree = GenerateTree(ref distanceMatrix, leafDataList);
            
            tree.GetWeights();

            return tree;
        }

        public static GuideTree<T> GetWeightsAndSteps(ref double[,] distanceMatrix, T[] leafDataList)
        {
            GuideTree<T> tree = GenerateTree(ref distanceMatrix, leafDataList);

            tree.GetWeightsAndSteps();

            return tree;
        }

        public void GetWeightsAndSteps()
        {
            GetWeights();

            CalculateSimilarities(); // Calculate a similarities matrix

            DetermineSteps();
        }

        public void GetWeights()
        // Calculates sequence weights from the guide tree and stores them in the tree.
        // For each leaf (sequence), Clustal sends to a secondary function to calculate the weight of a single sequence.
        // Then it normalizes the weights so that the sum of them is equal to the constant INT_SCALE_FACTOR, in order to convert
        // the floating-point values into integers. Unsure why this is necessary. I won't for the time being.
        // If there are only two sequences, it assigns equal weights.
        //
        // The single-sequence weight calculation is the distance to the sum of the branch lengths from the leaf to the root, each
        // time dividing by the 'order' of the node.

        {
            // Reset any previous weights
            foreach (TreeLeaf leaf in leavesList)
            {
                leaf.weight = 0;
            }


            // Step 1: Calculate sequence weights
            if (leavesCount > 2)
            {
                double sumWeights = 0;
                foreach (TreeLeaf leaf in leavesList)
                {
                    leaf.weight = leaf.parentBranch.BranchLength;
                    TreeNode parent = (TreeNode)leaf.parentBranch.Follow(leaf);
                    
                    while (parent.parentBranch != null)
                    {
                        leaf.weight += parent.parentBranch.BranchLength / parent.childCount;
                        parent = (TreeNode)parent.parentBranch.Follow(parent);
                    }
                    sumWeights += leaf.weight;
                }

                // Normalize the sequence weights
                if (sumWeights == 0)
                {
                    // For some reason all the sequence weights were 0, so we'll distribute evenly.
                    double weight = 1 / leavesCount;
                    
                    foreach (TreeLeaf leaf in leavesList)
                    {
                        leaf.weight = weight;
                    }
                }
                else
                {
                    // Note that Clustal adjusts very small or 0 weights up slightly. Unsure what effect not
                    // doing that will have.
                    double weightFactor = 1 / sumWeights;
                    foreach (TreeLeaf leaf in leavesList)
                    {
                        leaf.weight *= weightFactor;
                    }
                }
            }
            else
            {
                // If there are only two sequences, we can just set each to 0.5
                foreach (TreeLeaf leaf in leavesList)
                {
                    leaf.weight = 0.5;
                }
            }
                     
            // Step 2: Calculate similarities
            // Why is this necessary?


        }

        public void CalculateSimilarities()
        // Construct a similarities matrix
        // Clustal does this by storing the path to root for each leaf (storing the distance), then finds the common ancestor
        // for each pair, storing the distance from the other sequence to it.
        //
        // I think a better implementation would be to loop over the leaves to get their paths to root, then for each pair,
        // walk their paths until they diverge
        {
            similarityMatrix = new double[leavesCount, leavesCount];

            Dictionary<TreeLeaf, List<Tuple<TreeNode, double>>> pathsToRoot = new Dictionary<TreeLeaf, List<Tuple<TreeNode, double>>>();
            foreach (TreeLeaf leaf in leavesList)
            {
                List<Tuple<TreeNode, double>> pathToRoot = new List<Tuple<TreeNode, double>>();

                double distance = 0;
                TreeBranch parentBranch = leaf.parentBranch;
                TreeNode parentNode = (TreeNode)parentBranch.Follow(leaf);

                while (true)
                {
                    distance += parentBranch.BranchLength;
                    pathToRoot.Add(Tuple.Create(parentNode, distance));
                    parentBranch = parentNode.parentBranch;
                    if (parentBranch != null)
                    {
                        parentNode = (TreeNode)parentBranch.Follow(parentNode);
                    }
                    else
                    {
                        break;
                    }
                }
                pathsToRoot.Add(leaf, pathToRoot);
            }

            
            for (int i = 0; i < leavesCount - 1; i++)
            {
                TreeLeaf leafI = leavesList[i];
                List<Tuple<TreeNode, double>> pathToRootI = pathsToRoot[leafI];

                for (int j = i + 1; j < leavesCount; j++)
                {
                    TreeLeaf leafJ = leavesList[j];
                    List<Tuple<TreeNode, double>> pathToRootJ = pathsToRoot[leafJ];
                    bool found = false;

                    foreach (Tuple<TreeNode, double> pathMemberI in pathToRootI)
                    {
                        foreach (Tuple<TreeNode, double> pathmemberJ in pathToRootJ)
                        {
                            if (pathMemberI.Item1 == pathmemberJ.Item1)
                            {
                                found = true;
                                similarityMatrix[i, j] = 1 - (pathMemberI.Item2 + pathmemberJ.Item2);
                                break;
                            }
                        }
                        if (found) { break; }
                    }
                }
            }

            // Then Clustal forces any values less than 0.01 to be 0.01 (within the sequence ranges)
            // and sets small values above 1 to 1, but values above 1.1 go into an error handling method.
            // I won't implement this yet. Note that this is done while they are still in distances, not
            // similarities.

        }


        void DetermineSteps()
        // Determine the groups to be used in each step of the alignment.
        // Step through the tree to find the terminal leaf-only nodes and step back to join them together,
        // with each step stored in a tuple, and the leaf groups 1 and 2 stored in lists in members 1 and 2 of the
        // tuple, respectively.
        // Possible to parallelize, but it looks like the overhead is generally too much for the work being done.
        // may want to test that on larger sets, however.
        {
            if (!isRooted)
            {
                throw new InvalidOperationException("The tree must be rooted to determine alignment steps.");
            }

            Dictionary<TreeNode, List<TreeLeaf>> waitingList = new Dictionary<TreeNode, List<TreeLeaf>>(); // Store the members of one child of the node until the other child has been evaluated.

            steps = new List<Tuple<ReadOnlyCollection<T>, ReadOnlyCollection<T>, double>>();

            foreach (TreeNode terminalNode in terminalNodesList)
            {
                List<TreeLeaf> group1 = new List<TreeLeaf>();
                TreeLeaf leaf1 = (TreeLeaf)terminalNode.childBranch1.Follow(terminalNode);
                group1.Add(leaf1);
                List<TreeLeaf> group2 = new List<TreeLeaf>();
                TreeLeaf leaf2 = (TreeLeaf)terminalNode.childBranch2.Follow(terminalNode);
                group2.Add(leaf2);


                steps.Add(CreateStep(group1, group2));

                TreeNode currentNode = terminalNode;

                List<TreeLeaf> childLeaves = new List<TreeLeaf>();
                childLeaves.Add(leaf1);
                childLeaves.Add(leaf2);

                while (currentNode != root)
                {
                    TreeBranch parentBranch = currentNode.parentBranch;
                    TreeNode parentNode = (TreeNode)parentBranch.Follow(currentNode);
                    TreeBranch otherChildBranch;

                    group1 = new List<TreeLeaf>();
                    group2 = new List<TreeLeaf>();

                    if (parentNode.childBranch1 == parentBranch)
                    {
                        otherChildBranch = parentNode.childBranch2;
                    }
                    else // should be childBranch2
                    {
                        otherChildBranch = parentNode.childBranch1;
                    }

                    if (otherChildBranch.Follow(parentNode).GetType() == typeof(TreeLeaf))
                    // The other child is a leaf, so we're okay.
                    {
                        group1 = childLeaves;
                        TreeLeaf otherChildLeaf = (TreeLeaf)otherChildBranch.Follow(parentNode);
                        group2.Add(otherChildLeaf);
                    }
                    else if (waitingList.ContainsKey(parentNode))
                    // The other child is a node which was waiting for the other children to be evaluated.
                    {
                        group1 = childLeaves.ToList<TreeLeaf>();
                        group2 = waitingList[parentNode];
                    }
                    else
                    // The other child is a node. We'll need to save the current child and store the values
                    // for later. This also means we need to start over with another terminal node.
                    {
                        waitingList.Add(parentNode, childLeaves);
                        break;
                    }

                    steps.Add(CreateStep(group1, group2));

                    foreach (TreeLeaf leaf in group2)
                    {
                        childLeaves.Add(leaf);
                    }

                    currentNode = parentNode;
                }
            }
        }

        Tuple<ReadOnlyCollection<T>, ReadOnlyCollection<T>, double> CreateStep(List<TreeLeaf> group1, List<TreeLeaf> group2)
        {
            List<T> group1Data = new List<T>();
            List<T> group2Data = new List<T>();

            foreach (TreeLeaf leaf in group1)
            {
                group1Data.Add(leaf.Data);
            }
            foreach (TreeLeaf leaf in group2)
            {
                group2Data.Add(leaf.Data);
            }

            double similaritySum = 0;
            int totalComparisons = 0;
            double meanSimilarity = 0;
            foreach (Tree.GuideTree<T>.TreeLeaf leaf1 in group1)
            {
                foreach (Tree.GuideTree<T>.TreeLeaf leaf2 in group2)
                {
                    int i = leavesList.IndexOf(leaf1);
                    int j = leavesList.IndexOf(leaf2);

                    similaritySum += similarityMatrix[Math.Min(i, j), Math.Max(i, j)];
                    totalComparisons++;
                }
            }

            meanSimilarity = similaritySum / (double)totalComparisons;

            return Tuple.Create(group1Data.AsReadOnly(), group2Data.AsReadOnly(), meanSimilarity);
        }

        public static GuideTree<T> GenerateTree(ref double[,] distanceMatrix, T[] leafDataList)
        {
            GuideTree<T> tree;
            if (leafDataList.Count() >= 3)
            {
                bool useUPGMA = false; // Would be a setting from the interface.
                IClusteringAlgorithm<T> clusteringAlgorithm;

                if (useUPGMA)
                {
                    //cluster = new UPGMA();
                    throw new NotImplementedException("UPGMA algorithm is not yet implemented.");
                }
                else
                {
                    clusteringAlgorithm = new NJ<T>();
                    
                    tree = clusteringAlgorithm.GenerateTree(distanceMatrix, leafDataList);
                    tree.ReRoot(); // Root the unrooted tree
                }
                tree.CountNodeChildren(); // For each node, count the number of leaves downstream of the node.
            }
            else
            {
                tree = GenerateTwoLeafTree(distanceMatrix, leafDataList);
            }

            return tree;
        }

        static GuideTree<T> GenerateTwoLeafTree(double[,] distanceMatrix, T[] leafDataList)
        {
            GuideTree<T> tree = new GuideTree<T>();
            TreeNode rootNode = new TreeNode();
            
            TreeLeaf leaf1 = new TreeLeaf(leafDataList[0], tree);
            TreeLeaf leaf2 = new TreeLeaf(leafDataList[1], tree);

            double branchLength = 0.5 * distanceMatrix[0, 1];

            GuideTree<T>.TreeBranch.Build(rootNode, leaf1, branchLength);
            GuideTree<T>.TreeBranch.Build(rootNode, leaf2, branchLength);
            
            tree.root = rootNode;

            return tree;
        }

        void ReRoot()
        // Root an unrooted tree or adjust the root location after a rooted tree has beem modified
        {
            // Alternative algorithm:
            //   0. If the tree is already rooted, remove the root node.
            //   1. Find the average branch lengths extending from every node.
            //   2. Compute the differences using the largest branch
            //         Treat the branches as vectors leaving the node separated evenly by 120 degrees, and treat
            //         the largest branch length as the top vertical vector. Subtract the vertical components of the
            //         other two vectors from the largest to find the difference score.
            //         Since cos(60 deg) = 0.5, just multiply the other vectors' magnitudes by 0.5.
            //   3. Select the node with the smallest difference.
            //   4. Select the branch from that node with the largest average branch length.
            //      a. If two or three branches are tied, select the branch among the ties with the longest
            //         length to the next node.
            //      b. If they are still tied, arbitrarily pick the first.
            //   5. Insert the root between the identified nodes, and set the branch lengths
            //      satisfying x = (internode branch length - |difference|) / 2, where the side with the smaller average
            //      branch length is increased by the difference.
            //   6. Walk the tree from the new root to assert parent relationships as needed.

            if (isRooted)
            {
                double branchLength;
                TreeMember connection1, connection2;
                // The root should only connect to two members
                connection1 = root.Branches[0].Follow(root);
                connection2 = root.Branches[1].Follow(root);
                branchLength = root.Branches[0].BranchLength + root.Branches[1].BranchLength;

                root.Branches[0].Destroy();
                root.Branches[1].Destroy();

                TreeBranch.Build(connection1, connection2, branchLength);
            }

            TreeMember rootConnection1, rootConnection2;

            TreeNode minNode = new TreeNode();
            double minDifference = Double.MaxValue;
            TreeBranch minNodeMaxBranch = new TreeBranch();
            TreeBranch minNodeTiedMaxBranch = new TreeBranch();
            // Find the node with the smallest vector difference
            foreach (TreeNode node in nodesList) // parallelize?
            {
                Dictionary<TreeBranch, double> averageBranchLengths = new Dictionary<TreeBranch, double>();
                TreeBranch maxBranch = null;
                TreeBranch tiedMaxBranch = null;
                double maxAverageBranchLength = 0;
                foreach (TreeBranch branch in node.Branches)
                {
                    int numLeaves;
                    double averageBranchLength = GetAverageBranchLength(node, branch, true, out numLeaves);
                    averageBranchLengths.Add(branch, averageBranchLength);
                    if (averageBranchLength > maxAverageBranchLength)
                    {
                        maxAverageBranchLength = averageBranchLength;
                        maxBranch = branch;
                        tiedMaxBranch = null;
                    }
                    else if (averageBranchLength == maxAverageBranchLength)
                    {
                        tiedMaxBranch = branch;
                    }
                }
                
                double vectorDifference = maxAverageBranchLength;
                foreach (TreeBranch branch in node.Branches)
                {
                    if (branch != maxBranch)
                    {
                        vectorDifference -= 0.5 * averageBranchLengths[branch];
                    }
                }

                if (vectorDifference < minDifference)
                {
                    minDifference = vectorDifference;
                    minNode = node;
                    minNodeMaxBranch = maxBranch;
                    minNodeTiedMaxBranch = tiedMaxBranch;
                }
            }
            rootConnection1 = minNode;

            // Found the minimum node. Now we must decide which branch the root should replace.
            // In the case of ties (minDifference = 0 or minNodeTiedMaxBranch is not null), the root
            // will be placed along the longest adjacent branch. If tie isn't broken, it will be the
            // first branch used.
            TreeBranch rootBranch;
            if (minDifference == 0)
            {
                double maxBranchLength = 0;
                TreeBranch maxBranch = new TreeBranch();
                foreach (TreeBranch branch in minNode.Branches)
                {
                    double branchLength = branch.BranchLength;

                    if (branchLength > maxBranchLength)
                    {
                        maxBranchLength = branchLength;
                        maxBranch = branch;
                    }
                }
                rootBranch = maxBranch;
            }
            else if (minNodeTiedMaxBranch != null)
            {
                if (minNodeMaxBranch.BranchLength >= minNodeTiedMaxBranch.BranchLength)
                {
                    rootBranch = minNodeMaxBranch;
                }
                else
                {
                    rootBranch = minNodeTiedMaxBranch;
                }
            }
            else
            {
                rootBranch = minNodeMaxBranch;
            }
            rootConnection2 = rootBranch.Follow(rootConnection1);

            // Determine the lengths of the branches that will connect to the root.
            // This is determined by the equation:
            //     Li = (Aj - Ai + (Nj * L)) / (Nj + Ni)
            // Where "i" and "j" each represent one of the root connections, L indicates length,
            // A indicates average branch length, and N indicates number of leaves.
            int numLeaves1, numLeaves2;
            double averageBranchLength1 = GetAverageBranchLength(rootConnection1, rootBranch, false, out numLeaves1);
            double averageBranchLength2 = GetAverageBranchLength(rootConnection2, rootBranch, false, out numLeaves2);
            double branchDifference = averageBranchLength1 - averageBranchLength2;
            double rootBranchLength1, rootBranchLength2;
            if (branchDifference != 0)
            {
                rootBranchLength1 = (averageBranchLength2 - averageBranchLength1 + (numLeaves2 * rootBranch.BranchLength)) / (numLeaves1 + numLeaves2);
            }
            else
            {
                rootBranchLength1 = (numLeaves2 * rootBranch.BranchLength) / (numLeaves1 + numLeaves2);
            }
            rootBranchLength2 = rootBranch.BranchLength - rootBranchLength1;

            rootBranch.Destroy(); // Remove the branch that will be replaced by the root node.

            isRooted = true;
            root = new TreeNode(this);
            TreeBranch.Build(root, rootConnection1, rootBranchLength1);
            TreeBranch.Build(root, rootConnection2, rootBranchLength2);
            
            // Now we need to walk the tree from the root to ensure that all the downstream
            // nodes have the parent assigned correctly.
            AssertPaternity();
        }

        static double GetAverageBranchLength(TreeMember member, TreeBranch baseBranch, bool includeBranch, out int numLeaves)
        // From the given tree member, this will calculate the average branch length of the branches extending
        // from it. The switch "includeBranch" indicates whether the branch should be included in the numbers
        // (used for distance to the other end of the branch). If the TreeMember is a leaf, this function returns 0.
        {
            double averageBranchLength = 0;
            double sumBranchLength = GetSumBranchLength(member, baseBranch, includeBranch, out numLeaves);

            averageBranchLength = sumBranchLength / numLeaves;

            return averageBranchLength;
        }

        static double GetSumBranchLength(TreeMember member, TreeBranch baseBranch, bool includeBranch, out int numLeaves)
        // This is a recursive function
        //
        {
            // Algorithm:
            //   1. Store the branch length to the current member in "distance" (except the first if includeBranch is true).
            //   2. If the current member is a node, add the other branches to the stack and continue the loop.
            //      If the current member is a leaf, add 1 to the number of leaves, add the stemLength to the sumBranchLength, and subtract the length to this branch from the distance.
            //   3. 

            double sumBranchLength = 0;
            numLeaves = 0;
            Stack<Tuple<TreeMember, TreeBranch>> arguments = new Stack<Tuple<TreeMember, TreeBranch>>();
            double distance = 0;
            arguments.Push(Tuple.Create(baseBranch.Follow(member), baseBranch));
            TreeBranch prevBranch = new TreeBranch();
            bool stepBack = false;
            while (arguments.Count > 0)
            {
                Tuple<TreeMember, TreeBranch> nextArgument = arguments.Peek();
                TreeMember nextMember = nextArgument.Item1;
                TreeBranch nextBaseBranch = nextArgument.Item2;

                if (nextMember.GetType() == typeof(TreeNode))
                {
                    TreeNode node = (TreeNode)nextMember;
                    if (node.Branches.Contains(prevBranch))
                    // If the previously-used branch is a member of this node's branches, that means we've already used the
                    // branches of this node, and can discard it.
                    {
                        stepBack = true;
                    }
                    else
                    {
                        distance += nextBaseBranch.BranchLength;
                        foreach (TreeBranch branch in node.Branches)
                        {
                            if (branch != nextBaseBranch)
                            {
                                arguments.Push(Tuple.Create(branch.Follow(node), branch));
                            }
                        }
                    }
                }
                else
                // This is a leaf, so add the distance to the sum, increment the leaf count, and step backward
                {
                    distance += nextBaseBranch.BranchLength;
                    sumBranchLength += distance;
                    numLeaves++;
                    stepBack = true;
                }

                if (stepBack)
                {
                    prevBranch = nextBaseBranch;
                    distance -= nextBaseBranch.BranchLength;
                    arguments.Pop();
                    stepBack = false;
                }

            }

            if (!includeBranch)
            {   // If we didn't want to include the base branch length, correct for that here.
                sumBranchLength -= numLeaves * baseBranch.BranchLength;
            }

            return sumBranchLength;
        }

        void CountNodeChildren()
        // Walk the tree up from the terminal nodes and count the children of each node.
        {
            if (!isRooted)
            {
                throw new InvalidOperationException("The tree must be rooted to count children.");
            }

            foreach (TreeNode terminalNode in terminalNodesList)
            {
                terminalNode.childCount = 2;

                TreeNode currentNode = terminalNode;

                while (currentNode != root)
                {
                    TreeBranch parentBranch = currentNode.parentBranch;
                    TreeNode parentNode = (TreeNode)parentBranch.Follow(currentNode);
                    TreeBranch otherChildBranch;
                    
                    if (parentNode.childBranch1 == parentBranch)
                    {
                        otherChildBranch = parentNode.childBranch2;
                    }
                    else // should be childBranch2
                    {
                        otherChildBranch = parentNode.childBranch1;
                    }

                    if (otherChildBranch.Follow(parentNode).GetType() == typeof(TreeLeaf))
                    {
                        parentNode.childCount = currentNode.childCount + 1;
                    }
                    else if (parentNode.childCount > 0)
                    {
                        parentNode.childCount += currentNode.childCount;
                    }
                    else
                    {
                        // Need to wait for another terminal node to propogate up.
                        parentNode.childCount += currentNode.childCount;
                        break;
                    }

                    currentNode = parentNode;
                }
            }
        }

        void AssertPaternity()
        // For each branch out from the given root which is not a parent branch, assign that branch
        // as the child's parent branch.
        // As this could be deeply recursive and could result in a stack overflow, this must be coded as a loop.
        {
            Stack<Tuple<TreeMember, TreeBranch>> arguments = new Stack<Tuple<TreeMember, TreeBranch>>();
            TreeBranch prevBranch = new TreeBranch();
            bool stepBack = false;
            foreach (TreeBranch rootBranch in root.Branches)
            {
                arguments.Push(Tuple.Create(rootBranch.Follow(root), rootBranch));

                while (arguments.Count > 0)
                {
                    Tuple<TreeMember, TreeBranch> nextArgument = arguments.Peek();
                    TreeMember member = nextArgument.Item1;
                    TreeBranch parentBranch = nextArgument.Item2;
                    stepBack = false;

                    member.parentBranch = parentBranch;
                    ((TreeNode)parentBranch.Follow(member)).SetChild(parentBranch);
                    
                    if (member.GetType() == typeof(TreeNode))
                    {
                        TreeNode node = (TreeNode)member;
                        if (node.Branches.Contains(prevBranch))
                        // If the previously-used branch is a member of this node's branches, that means we've already used the
                        // branches of this node, and can discard it.
                        {
                            stepBack = true;
                        }
                        else
                        {
                            foreach (TreeBranch branch in node.Branches)
                            {
                                if (branch != parentBranch)
                                {
                                    arguments.Push(Tuple.Create(branch.Follow(node), branch));
                                }
                            }
                        }
                    }
                    else
                    {
                        stepBack = true;
                    }

                    if (stepBack)
                    {
                        prevBranch = parentBranch;
                        arguments.Pop();
                    }
                }
            }
        }
        
        public class TreeMember
        {
            public TreeBranch parentBranch;
            protected GuideTree<T> tree;

            public GuideTree<T> Tree
            {
                get { return tree; }
            }

            public TreeMember(GuideTree<T> newTree)
            {
                tree = newTree;
            }

            public TreeMember()
            // This shouldn't be used unless you don't intend to add the member to a tree.
            { }

            public TreeNode GetParent()
            {
                if (parentBranch != null)
                {
                    return (TreeNode)parentBranch.Follow(this);
                }
                else
                {
                    return null;
                }
            }

            public void RemoveParent()
            {
                parentBranch = null;
            }
        }

        public class TreeNode : TreeMember
        {
            List<TreeBranch> branches = new List<TreeBranch>(3);
            public int childCount = 0; // Holds the number of leaves downstream of the TreeNode in a rooted tree.
            public TreeBranch childBranch1;
            public TreeBranch childBranch2;

            public ReadOnlyCollection<TreeBranch> Branches
            {
                get { return branches.AsReadOnly(); }
            }
            
            public TreeNode(GuideTree<T> newTree)
                : base(newTree)
            {
                newTree.nodesList.Add(this);
            }

            public TreeNode()
            // This shouldn't be used unless you don't intend to add the member to a tree.
            { }

            public void AddBranch(TreeBranch branch)
            {
                // A node can only support three connections
                if (branches.Count == 3)
                {
                    throw new OverflowException("This node already has 3 connections.");
                }
                else
                {
                    // If the branch being added connects this node to a leaf, and one of the other branches also
                    // connects to a leaf, this is a terminal node. Add this node to the list of terminal nodes.
                    if (branch.Follow(this).GetType() == typeof(TreeLeaf))
                    {
                        foreach (TreeBranch otherBranch in branches)
                        {
                            if (otherBranch.Follow(this).GetType() == typeof(TreeLeaf))
                            {
                                tree.terminalNodesList.Add(this);
                            }
                        }
                    }

                    branches.Add(branch);
                }

            }

            public void RemoveBranch(TreeBranch branch)
            {
                // List's Remove method handles low-count and not-in-list conditions,
                // so I won't handle them here.
                // bool removeFlag = branches.Remove(branch); Not sure why I need the return code.
                branches.Remove(branch);

                if (childBranch1 == branch)
                {
                    childBranch1 = null;
                }
                else if (childBranch2 == branch)
                {
                    childBranch2 = null;
                }


                //if (!removeFlag)
                //{
                //    return false;
                //}

                if (parentBranch == branch)
                {
                    RemoveParent();
                }

                // If the branch removal results in no connections, the node can be removed from the tree.
                if (branches.Count == 0)
                {
                    tree.nodesList.Remove(this);
                }
                else if (true)
                {
                    // If the branch was a terminal node and the branch being removed was connected to a leaf,
                    // if only one other branch is connected to a leaf, this is no longer a terminal node.
                    // Remove this node from the list of terminal nodes.
                    if (branch.Follow(this).GetType() == typeof(TreeLeaf))
                    {
                        int leafCount = 0;
                        foreach (TreeBranch otherBranch in branches)
                        {
                            if (otherBranch.Follow(this).GetType() == typeof(TreeLeaf))
                            {
                                leafCount += 1;
                            }
                        }

                        if (leafCount <= 1)
                        {
                            tree.terminalNodesList.Remove(this);
                        }
                    }
                }
                //return true;

                
            }

            public TreeBranch GetBranch(TreeMember member)
            {
                foreach (TreeBranch branch in branches)
                {
                    if (branch.Follow(this) == member)
                    {
                        return branch;
                    }
                }
                return null;
            }

            public void SetChild(TreeBranch branch)
            {
                if (branches.Contains(branch))
                {
                    if (childBranch1 != branch && childBranch2 != branch) // If the branch isn't already a child
                    {
                        if (childBranch1 == null)
                        {
                            childBranch1 = branch;
                        }
                        else if (childBranch2 == null)
                        {
                            childBranch2 = branch;
                        }
                        else
                        {
                            throw new OverflowException("This node already has two children.");
                        }
                    }
                }
                else
                {
                    throw new ArgumentException("Cannot set a branch as a child of this node if it isn't connected to this node.");
                }
                
            }
        }

        public class TreeLeaf : TreeMember
        {
            T data;
            public double weight;

            public T Data
            {
                get { return data; }
            }

            public TreeLeaf(T leafData, GuideTree<T> newTree)
                : base(newTree)
            {
                data = leafData;
                newTree.leavesCount++;
                newTree.leavesList.Add(this);
            }
            
        }

        public class TreeBranch
        {
            double branchLength;
            TreeMember connection1;
            TreeMember connection2;

            public double BranchLength
            {
                get { return branchLength; }
            }

            TreeBranch(TreeMember member1, TreeMember member2, double length)
            {
                connection1 = member1;
                connection2 = member2;
                branchLength = length;
            }

            public TreeBranch()
            {
            }

            public static void Build(TreeMember member1, TreeMember member2, double length)
            {
                TreeBranch newBranch = new TreeBranch(member1, member2, length);
                
                if (member1.GetType() == typeof(TreeNode))
                {
                    ((TreeNode)member1).AddBranch(newBranch);
                }
                else
                {
                    member1.parentBranch = newBranch;
                }

                if (member2.GetType() == typeof(TreeNode))
                {
                    ((TreeNode)member2).AddBranch(newBranch);
                }
                else
                {
                    member2.parentBranch = newBranch;
                }
            }

            public void Destroy()
            {
                if (connection1.GetType() == typeof(TreeNode))
                {
                    ((TreeNode)connection1).RemoveBranch(this);
                }
                else
                {
                    connection1.RemoveParent();
                }

                if (connection2.GetType() == typeof(TreeNode))
                {
                    ((TreeNode)connection2).RemoveBranch(this);
                }
                else
                {
                    connection2.RemoveParent();
                }
            }

            public TreeMember Follow(TreeMember origin)
            // Returns the connection that the branch leads to, given the starting node
            {
                if (origin == connection1) { return connection2; }
                if (origin == connection2) { return connection1; }

                throw new ArgumentException("This branch does not connect to the member provided.");
            }
        }
    }

    // Clustal's RootedGuideTree (UPGMA) = Tree (NJ)
    // Clustal's ClusterTree does a lot of input/output stuff for some reason.
    // UPGMA produces a rooted tree. NJ produces an unrooted tree.

    // The clustering algorithm is set within the current UnRootedClusterTree or RootedClusterTree classes,
    // so does that mean that the rest of the code wouldn't care about the type of tree being produced/how it was produced?
    // If so, we'll just have one GuideTree class

    // Do I want to put the whole macromolecule into the TreeLeaf?
}
