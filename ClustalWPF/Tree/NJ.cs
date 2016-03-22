using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.Threading.Tasks;

namespace ClustalWPF.Tree
{
    class NJ<T> : IClusteringAlgorithm<T>
    // Implementation of the Neighbor-Joining (NJ) hierarchical clustering method.
    {
        public GuideTree<T> GenerateTree(double[,] distanceMatrix, T[] leafDataList)
        // Build and return an unrooted cluster tree from the pairwise distance matrix using the Neighbor Joining algorithm.
        // The number of sequences must be >= 3.
        {
            int numNodes = leafDataList.Count();
            double[] rowDistanceSums = new double[numNodes];
            double[] colDistanceSums = new double[numNodes];
            double[] distanceSums = new double[numNodes]; // Stores the sum of the distances of one sequence to every other sequence. Equal to the sum of the row and column of the distance matrix with the same index
            double rowTotal = 0;
            double[,] pairDistanceSums = new double[numNodes, numNodes]; // Stores the pre-computed values for every test node pair. The active on is given by the first dimension.
            double[,] notPairDistanceSums = new double[numNodes, numNodes]; // Stores the pre-computed values for the remainder of sequences for every test node pair.
            GuideTree<T>.TreeNode newNode;
            GuideTree<T> tree = new GuideTree<T>();

            LinkedList<Tuple<int, GuideTree<T>.TreeMember, double>> Nodes = new LinkedList<Tuple<int, GuideTree<T>.TreeMember, double>>(); // Stores the operational taxonomical units available and a mapping to the corresponding index in the distance matrix and the average branch length of the contained branches.

            // Initialize components necessary to direct the loop and pre-compute factors that are used
            // multiple times throughout the loop.

            // Populate the OTU array
            for (int i = 0; i < numNodes; i++)
            {
                GuideTree<T>.TreeLeaf leaf = new GuideTree<T>.TreeLeaf(leafDataList[i], tree);
                Nodes.AddLast(Tuple.Create<int, GuideTree<T>.TreeMember, double>(i, leaf, 0));
            }

            // Compute the row and column sums of the distance matrix and the total sums by sequence and in full.
            // We can skip the first column because it are always 0.
            for (int i = 0; i < numNodes - 1; i++)
            {
                for (int j = i + 1; j < numNodes; j++)
                {
                    double distanceIJ = distanceMatrix[i, j];
                    rowDistanceSums[i] += distanceIJ;
                    colDistanceSums[j] += distanceIJ;
                }
            }
            for (int i = 0; i < numNodes; i++)
            {
                distanceSums[i] = rowDistanceSums[i] + colDistanceSums[i];
                rowTotal += rowDistanceSums[i];
            }

            // Compute the pairDistanceSums. This is the sum of the distances from one OTU to every other OTU except its partner given by the second index.
            // Also compute the notPairDistanceSums. This is the sum of the distances from every other OTU to every other OTU, except the pair given by the indices.
            for (int i = 0; i < numNodes; i++)
            {
                for (int j = i + 1; j < numNodes; j++)
                {
                    double distanceIJ = distanceMatrix[i, j];
                    pairDistanceSums[i, j] = distanceSums[i] - distanceIJ;
                    pairDistanceSums[j, i] = distanceSums[j] - distanceIJ;

                    notPairDistanceSums[i, j] = rowTotal - rowDistanceSums[i] - rowDistanceSums[j] + distanceIJ - colDistanceSums[i] - colDistanceSums[j];
                }
            }

            // Main cycle. Loop until only 3 OTUs remain.
            for (; numNodes > 3; numNodes--)
            {
                // To save time in the main loop, pre-compute factors that we would otherwise use to divide, a slower operation
                double numOTULess2Factor = 1.0 / (double)(numNodes - 2);

                // Look for the neighbor pair which gives the smallest sum of branch lengths
                double minSumBranchLength = double.MaxValue; // Store the minimum sum branch length. Start high so the Min function works.
                double minPairDistance = 0; // Store the distance between the minimum OTU pair.
                LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> minI = new LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>>(Tuple.Create<int, GuideTree<T>.TreeMember, double>(0, null, 0)); // Store the sequence indexes at the minimum sum branch length and the OTU associated with it.
                LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> minJ = new LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>>(Tuple.Create<int, GuideTree<T>.TreeMember, double>(0, null, 0));

                //Could be parallelized
                for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuITuple = Nodes.First; otuITuple != Nodes.Last; otuITuple = otuITuple.Next)
                {
                    int i = otuITuple.Value.Item1;

                    for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuJTuple = otuITuple.Next; otuJTuple != null; otuJTuple = otuJTuple.Next)
                    {
                        int j = otuJTuple.Value.Item1;

                        // Sum of branch lengths formula:
                        //     X = Sum of the distances to every sequence from each member of the test pair excluding each other.
                        //     Y = Sum of the rows of the distance matrix excluding the test pair.
                        //     Z = Distance between the test pair.
                        //     N = Number of OTUs.
                        //
                        //     Sij = [(Xij + Xji) / (2 * (N - 2))]  +  [Z / 2]  +  [Y / (N - 2)]
                        // This is slightly simplified, using multiplication instead:
                        //     Sij = numOTULess2Factor * [0.5 * (Xij + Xji) + Y] + 0.5 * Z
                        double sumBranchLength = numOTULess2Factor * (0.5 * (pairDistanceSums[i, j] + pairDistanceSums[j, i]) + notPairDistanceSums[i, j]) + 0.5 * distanceMatrix[i, j];

                        // Decide whether to keep the branch pair.
                        // If the new pair results in a shorter total length, keep it.
                        // If it is a tie or greater, ignore it.
                        if (minSumBranchLength > sumBranchLength)
                        {
                            minI = otuITuple;
                            minJ = otuJTuple;
                            minSumBranchLength = sumBranchLength;
                        }
                    }
                }

                int minIIndex = minI.Value.Item1;
                int minJIndex = minJ.Value.Item1;
                GuideTree<T>.TreeMember minIOTU = minI.Value.Item2;
                GuideTree<T>.TreeMember minJOTU = minJ.Value.Item2;
                minPairDistance = distanceMatrix[minIIndex, minJIndex];
                double minIAverageBranchLength = minI.Value.Item3;
                double minJAverageBranchLength = minJ.Value.Item3;

                // Found the minimum pair.
                // Calculate the distance between each member of the minimum pair and their common node.
                // The distance formula for a single member of the pair is:
                //     X = Sum of the distances to every sequence from each member of the test pair excluding each other.
                //     N = Number of OTUs
                //
                //     Li = 0.5 * (Dij + numOTULess2Factor * [Xij - Xji])
                //
                // This length represents the length from this new node to the end of the branch, which is the average of the branches
                // emerging from it. Therefore if we are connecting a node, we must reduce these values by that as well
                //     Li = 0.5 * (Dij + numOTULess2Factor * [Xij - Xji]) - averageBranchLengthI
                double branchLengthI = 0.5 * (minPairDistance + numOTULess2Factor * (pairDistanceSums[minIIndex, minJIndex] - pairDistanceSums[minJIndex, minIIndex])) - minIAverageBranchLength;
                double branchLengthJ = 0.5 * (minPairDistance + numOTULess2Factor * (pairDistanceSums[minJIndex, minIIndex] - pairDistanceSums[minIIndex, minJIndex])) - minJAverageBranchLength;

                //
                // NOTE: Here, Clustal does a couple things that I haven't done
                // First, it takes very small positive and negative branch lengths (|x| < 0.0001) to 0
                //     Translation: "These sequences are so close to identical that we'll treat them as identical."
                //     May well be a correction for floating point rounding errors.
                // Then, after recording the branch lengths in the growing tree, it sets any negative
                //     or 0 minPairDistance to 0.000001 before adjusting the distance matrix and the others
                //     for the next steps. This would casue the average branch length used in the next round
                //     which uses the combined OTU will be non-0, reducing the new branch lengths by a small
                //     amount. I do not know why this is helpful in general.

                // The NJ method can sometimes produce negative branch lengths. The best practice
                // seems to be to set it to zero and adjust the sister branch accordingly, to keep the
                // distance between the two OTUs the same.
                //     Kuhner, M.K., and Felsenstein, J. (1994) A simulation comparison of phylogeny algorithms
                //     under equal and unequal evolutionary rates. /Mol. Biol. Evol./ 11(3):459-468.

                if (branchLengthI < 0)
                {
                    branchLengthJ += branchLengthI;
                    branchLengthI = 0;
                }
                else if (branchLengthJ < 0)
                {
                    branchLengthI += branchLengthJ;
                    branchLengthJ = 0;
                }

                double averageBranchLength = minPairDistance / 2;

                // Build the node for this minimum pair
                newNode = new GuideTree<T>.TreeNode(tree);
                GuideTree<T>.TreeBranch.Build(newNode, minIOTU, branchLengthI); // make sure that not passing these by ref is okay
                GuideTree<T>.TreeBranch.Build(newNode, minJOTU, branchLengthJ);
                
                // Remove these OTUs from the list and store the position in which to add the new one.
                // Since both are removed, we don't need to worry about them in the next loops through the list.
                Nodes.Remove(minJ);
                //LinkedListNode<Tuple<int, GuideTree.TreeMember>> afterMinINode = minI.Next;
                Nodes.AddAfter(minI, Tuple.Create<int, GuideTree<T>.TreeMember, double>(minIIndex, newNode, averageBranchLength));
                Nodes.Remove(minI);
              
                // Recalculate the distance matrix
                for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuI = Nodes.First; otuI != null; otuI = otuI.Next)
                {
                    double minIDistance, minJDistance;
                    int i = otuI.Value.Item1;

                    if (i == minIIndex) { continue; }

                    int minIRow = Math.Min(i, minIIndex);
                    int minICol = Math.Max(i, minIIndex);

                    minIDistance = distanceMatrix[minIRow, minICol];
                    minJDistance = distanceMatrix[Math.Min(i, minJIndex), Math.Max(i, minJIndex)];

                    // The formula for the distance between any other OTU z and the new combined OTU is the average of the two being combined:
                    //     D = 0.5 * (Diz + Djz)
                    distanceMatrix[minIRow, minICol] = 0.5 * (minIDistance + minJDistance);
                }
                //
                //
                // Attempt to combine the adjustments.
                // Can't parallelize because of using linked lists.
                if (numNodes > 4) // Don't run if we're on the last iteration.
                {
                    for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuI = Nodes.First; otuI != null; otuI = otuI.Next)
                    {
                        int i = otuI.Value.Item1;

                        if (i == minIIndex) { continue; }

                        pairDistanceSums[minIIndex, i] -= 0.5 * (pairDistanceSums[minIIndex, i] - pairDistanceSums[minJIndex, i]) + distanceMatrix[minIIndex, minJIndex];
                        pairDistanceSums[i, minIIndex] -= distanceMatrix[Math.Min(minJIndex, i), Math.Max(minJIndex, i)];

                        double minIDistance = distanceMatrix[Math.Min(minIIndex, i), Math.Max(minIIndex, i)];

                        for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuJ = Nodes.First; otuJ != null; otuJ = otuJ.Next)
                        {
                            int j = otuJ.Value.Item1;

                            if (i == j || j == minIIndex) { continue; }

                            if (j != minJIndex)
                            {
                                notPairDistanceSums[Math.Min(minIIndex, i), Math.Max(minIIndex, i)] -= distanceMatrix[Math.Min(minJIndex, j), Math.Max(minJIndex, j)];
                            }

                            pairDistanceSums[i, j] -= minIDistance;
                        }

                        for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuJ = otuI.Next; otuJ != null; otuJ = otuJ.Next)
                        {
                            int j = otuJ.Value.Item1;

                            if (i == j || j == minIIndex) { continue; }

                            notPairDistanceSums[i, j] -= distanceMatrix[minIIndex, minJIndex];

                            for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otuK = Nodes.First; otuK != null; otuK = otuK.Next)
                            {
                                int k = otuK.Value.Item1;

                                if (i == k || j == k || k == minIIndex) { continue; }

                                notPairDistanceSums[i, j] -= distanceMatrix[Math.Min(minIIndex, k), Math.Max(minIIndex, k)];
                            }
                        }
                    } 
                }
            }

            GuideTree<T>.TreeMember[] finalOTUs = new GuideTree<T>.TreeMember[3];
            int[] finalOTUIndices = new int[3];
            double[] finalOTUAvgBranchLength = new double[3];
            int index = 0;
            for (LinkedListNode<Tuple<int, GuideTree<T>.TreeMember, double>> otu = Nodes.First; otu != null; otu = otu.Next)
            {
                finalOTUs[index] = otu.Value.Item2;
                finalOTUIndices[index] = otu.Value.Item1;
                finalOTUAvgBranchLength[index] = otu.Value.Item3;
                index++;
            }

            // Connect the final three OTUs to a single node.
            double distance01 = distanceMatrix[Math.Min(finalOTUIndices[0], finalOTUIndices[1]), Math.Max(finalOTUIndices[0], finalOTUIndices[1])];
            double distance02 = distanceMatrix[Math.Min(finalOTUIndices[0], finalOTUIndices[2]), Math.Max(finalOTUIndices[0], finalOTUIndices[2])];
            double distance12 = distanceMatrix[Math.Min(finalOTUIndices[1], finalOTUIndices[2]), Math.Max(finalOTUIndices[1], finalOTUIndices[2])];

            double branchLength0 = 0.5 * (distance01 + distance02 - distance12) - finalOTUAvgBranchLength[0];
            double branchLength1 = 0.5 * (distance12 + distance01 - distance02) - finalOTUAvgBranchLength[1];
            double branchLength2 = 0.5 * (distance12 + distance02 - distance01) - finalOTUAvgBranchLength[2];
            // Clustal seemed to find a way to do this in fewer operations, but I've been having some trouble with that. But this is just being done once, so do we care?

            newNode = new GuideTree<T>.TreeNode(tree);
            GuideTree<T>.TreeBranch.Build(newNode, finalOTUs[0], branchLength0);
            GuideTree<T>.TreeBranch.Build(newNode, finalOTUs[1], branchLength1);
            GuideTree<T>.TreeBranch.Build(newNode, finalOTUs[2], branchLength2);

            return tree;
            
            // Note also that multiplication is faster than division, but there is virtually no speed difference between multiplication/addition/subtraction
        }
    }

}