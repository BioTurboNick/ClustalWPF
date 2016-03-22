using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.MultipleAlignment
{
    class MultipleAlignmentIteration
    {
        public void IterationOnTreeNode(List<AlignedMacromolecule> macromolecules)
        {
            // Clustal makes a new alignment object

            // then adds sequences from both profiles to it

            // runs "removeFirstIterate" and determine if it changed

            // if it has changed, recalculate the alignment lengths,
            // which are passed back out to the calling code

            

            // Combine the sequences

            DoRemoveFirstIteration(macromolecules);
        }

        public void DoRemoveFirstIteration(List<AlignedMacromolecule> macromolecules)
        {
            // Remove-first iteration strategy
            // Optimize the alignment score by progressively removing sequences
            // each time a sequence is removed, remove all-gap columns and do the profileAlignment again
            // if the fit is better, keep the new alignment. If not, discard it.


            // Rewritten:
            // 1. "Remove" a sequence from the profile by placing it in the first position
            //    We can do without this because we'll be using lists.
            // 2. Remove any gap-only columns from both profiles
            // 3. Calculate a simple distance matrix (find the percent identity and then convert to proportion nonidentical)
            // 4. Clustal here uses a temporary tree file, but unclear what it is doing with it.
            // 5. Call GetWeightsFromProfileAlign from the Tree Interface
            //    This function takes the alignment, the simple distance matrix, the profile tree names,
            //    the profile weights, and "false" for useTree1 and useTree2
            //    a. call GenerateTree(distanceMatrix) for each profile
            //    b. generate a new MultipleAlignment object
            //    c. Call MultipleAlignment.CalcPairwiseForProfileAlign
            //       i. uses MyersMillerProfileAlign
            // 

            // Oh wait, in this case, profile1 is the removed sequence, profile2 is the rest

            //*****************
            // Okay, here's what needs to happen.
            // I need to remove a macromolecule from the list and then remove all-gap columns from both.
            // Then I need to calculate a simple distance matrix for them, get their weights,
            // and then do an alignment on them and compare the new score to the old score.
            // ---> If the new score is better, save this state and score value.
            // Now we start fresh from the original alignment and remove the 2nd macromolecule...
            //
            // So I can loop through the main list of macromolecules
            // I will also maintain a second list of the same that I will remove and add
            // from each cycle.
            // Then within the loop I will make copies of everything and process them.
            // The copies are saved as "best so far" if the score improves.
            // If the score doesn't improve, they are left untouched.
            //
            // What I can't tell though is what is being saved between iterations that
            // makes it an iteration, rather than just a repeat.
            //     So this module receives a new Alignment object which is shared across iterations,
            //     Note that the alignments are "reset" prior to the MSA
            // It *seems* like the strategy should be:
            //     if removing one improves, then try again up to the number of iterations
            //     Alternatively, it might actaully mean removing each of the sequences in turn
            //     When no more improvement occurs, it stops and reverts to the alignment from
            //     the best-scoring round. The latter appears to make the most sense.
            //
            // So really, I shouldn't be starting fresh each round.
            // Which means its okay to add back the altered removed sequence, but how to
            // keep track of which are being removed?
            // Well as long as I'm not altering the main collection, the details will just
            // be the thing to change.
            // And I'll only need to make copies of the current alignment when I find a new best score.



            int iterations = 3; // will be a userparameter

            // First, make a copy of the alignment part of the macromolecules:
            List<AlignedMacromolecule> activeMacromolecules = new List<AlignedMacromolecule>();

            foreach (AlignedMacromolecule macromolecule in macromolecules)
            {
                activeMacromolecules.Add(macromolecule.Copy());
            }

            // Store the total macromolecules for the enumerator
            AlignedMacromolecule[] allMacromolecules = new AlignedMacromolecule[activeMacromolecules.Count];
            activeMacromolecules.CopyTo(allMacromolecules);

            bool improved = false; // Indicates whether any improvement has been made
            for (int i = 0; i < iterations; i++)
            {
                bool improvedThisIteration = false; // Indicates whether any improvement has been made in the current iteration
                foreach (AlignedMacromolecule removedMacromolecule in activeMacromolecules)
                {
                    activeMacromolecules.Remove(removedMacromolecule);

                    // Remove gaps from the removed sequence and gaps made unnecessary by the removal in the rest.
                    Alignment.RemoveRedundantGaps(activeMacromolecules);
                    removedMacromolecule.ClearGaps();

                    // Calculate simple distance matrix
                    int numSequences = macromolecules.Count;
                    double[,] distanceMatrix = new double[numSequences, numSequences];
                    for (int j = 0; j < macromolecules.Count; j++)
                    {
                        for (int k = 0; k < macromolecules.Count; k++)
                        {
                            double percentIdentity = CalculatePercentIdentity(allMacromolecules[j], allMacromolecules[k]);
                            distanceMatrix[j, k] = (100.0 - percentIdentity) / 100.0;
                        }
                    }
                    
                    // Calculate weights

                    // Here Clustal calls TreeInterface.GetWeightsForProfileAlign(Alignment, DistMatrix, treeName1, treeWeights1, treeName2, treeWeights2, numSeqs, profile1_numSeqs, useTree1, useTree2, success)
                    //     doesn't itself do anything but call next subroutine
                    //     --> GetWeightsForProfileAlignNJ(same as above)
                    //         if (!useTree1 && profile1_numSeqs >= 2) --> GenerateTreeFromDistMatNJ(DistMatrix, Alignment, Profile1_sequences, treeName1, success)
                    //         if (!useTree2 && profiel2_numSeqs >= 2) --> GenerateTreeFromDistMatNJ(DistMatrix, Alignment, Profile2_sequences, treeName2, success)
                    //         --> MSA.CalcPairwiseForProfileAlign(Alignment, DistMatrix)
                    //             Does a few things and then
                    //             --> MyersMillerProfileAlign.ProfileAlign -- this would be the normal one
                    //             And wraps up
                    //         if (profile1_numSeqs >= 2) --> Tree1.ReadTree
                    //         --> Tree1.CalcSeqWeights
                    //         if (profile2_numSeqs >= 2) --> Tree2.ReadTree
                    //         --> Tree2.CalcSeqWeights
                    //         convert distances to similarities
                    
                    // Tree.GuideTree<AlignedMacromolecule> tree = Tree.GuideTree<AlignedMacromolecule>.GetWeights(ref distanceMatrix, allMacromolecules);
                    // MultipleAlignment.PairwiseAlign(tree, distanceMatrix, allMacromolecules); // Right name? I think it is supposed to return aligned macromolecules




                    // "Reset" the profiles (does this do anything?)

                    
                    // Do multiple sequence alignment


                    // Check score
                    //     If better, save save current alignments



                    activeMacromolecules.Add(removedMacromolecule);
                }

                if (!improvedThisIteration)
                {
                    // If we haven't improved the score this past iteration, no point in continuing
                    break;
                }
            }

            if (improved)
            {
                // If we've found an improved alignment, commit the best one.

            }
        
        }
    }
}
