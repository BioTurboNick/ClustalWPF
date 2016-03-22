using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.PairwiseAlignment
{
    class FastPairwiseAlignment : IPairwiseAlignmentAlgorithm
    {
        static Dictionary<char, int> nucleicAcidCodeDict;
        static Dictionary<char, int> proteinCodeDict;

        static FastPairwiseAlignment()
        {
            // Initialize the nucleic acid code lookup dictionary
            nucleicAcidCodeDict = new Dictionary<char, int>(5);
            nucleicAcidCodeDict.Add('A', 0);
            nucleicAcidCodeDict.Add('C', 1);
            nucleicAcidCodeDict.Add('G', 2);
            nucleicAcidCodeDict.Add('T', 3);
            nucleicAcidCodeDict.Add('U', 3);

            // Initialize the protein code lookup dictionary
            proteinCodeDict = new Dictionary<char, int>(20);
            proteinCodeDict.Add('A', 0);
            proteinCodeDict.Add('C', 1);
            proteinCodeDict.Add('D', 2);
            proteinCodeDict.Add('E', 3);
            proteinCodeDict.Add('F', 4);
            proteinCodeDict.Add('G', 5);
            proteinCodeDict.Add('H', 6);
            proteinCodeDict.Add('I', 7);
            proteinCodeDict.Add('K', 8);
            proteinCodeDict.Add('L', 9);
            proteinCodeDict.Add('M', 10);
            proteinCodeDict.Add('N', 11);
            proteinCodeDict.Add('P', 12);
            proteinCodeDict.Add('Q', 13);
            proteinCodeDict.Add('R', 14);
            proteinCodeDict.Add('S', 15);
            proteinCodeDict.Add('T', 16);
            proteinCodeDict.Add('V', 17);
            proteinCodeDict.Add('W', 18);
            proteinCodeDict.Add('Y', 19);
        }

        public void PairwiseAlign(ref Alignment alignmentObject, ref double[,] distanceMatrix)
            // Performs a pairwise alignment on all sequences in the alignmentObject, and returns the distance scores in the upper-right
            // triangle of the distance matrix. These operations are performed in parallel.
        {
            int coupling = 2; // Will be a parameter. ktup in Clustal. Determines whether to combine the information from multiple adjacent bases during pairwise alignment.
            bool isNucleicAcid = alignmentObject.IsNucleicAcid;
            int numSequences = alignmentObject.NumberMacromolecules;

            // Loop over every sequence pair, encode them, and do a virtual pairwise alignment to find the highest scoring pairs.
            // Every iteration of both the outer loop and inner loop are independent, and so are done in parallel.
            double[,] parallelDistanceMatrix = distanceMatrix; // The anonymous method (parallel loop) cannot directly use a ref variable
            Alignment parallelAlignmentObject = alignmentObject;
            
            Parallel.For(0, numSequences, i =>
            {
                string sequenceI = parallelAlignmentObject.AlignedMacromolecules[i].Sequence;
                int sequenceILength = sequenceI.Length;
                int[] sequenceIEncoded;

                LinkedList<int>[] sequenceICodePositions = GetEncodedPositions(sequenceI, sequenceILength, isNucleicAcid, coupling, out sequenceIEncoded);

                Parallel.For(i + 1, numSequences, j =>
                {
                    //if (i != j)  // Aligning one sequence to itself would result in 0.
                    //{
                    string sequenceJ = parallelAlignmentObject.AlignedMacromolecules[j].Sequence;
                    int sequenceJLength = sequenceJ.Length;

                    LinkedList<int>[] sequenceJCodePositions = GetEncodedPositions(sequenceJ, sequenceJLength, isNucleicAcid, coupling);

                    int alignmentScore = AlignPair(sequenceIEncoded, sequenceICodePositions, sequenceJCodePositions, sequenceILength, sequenceJLength, coupling);

                    // Convert the alignment score to a percentage based on the smallest sequence length
                    //double normalizedScore = (double)alignmentScore / (double)(Math.Min(sequenceILength, sequenceJLength));// - coupling); //+ 1);
                    //temporarily doing it exactly as they do...
                    double normalizedScore = ((double)alignmentScore / (double)Math.Min(sequenceILength, sequenceJLength)) * 100.0;
                    normalizedScore = (100.0 - normalizedScore) / 100.0;

                    parallelDistanceMatrix[i, j] = normalizedScore; // Put the score into the distance matrix. (Only filling the top right triangle).
                    //}
                });
            });
        }

        static LinkedList<int>[] GetEncodedPositions(string sequence, int length, bool isNucleicAcid, int coupling, out int[] sequenceEncoded, bool storeEncoded = true)
        {
            int limit; // The maximum code value + 1
            int possibilities;
            int coupledLength = length - coupling + 1;
            Dictionary<char, int> sequenceCodeDict; // The dictionary allowing code lookup from a sequence charater.

            if (storeEncoded)
            {
                sequenceEncoded = new int[coupledLength]; // Check for efficiency: single reusable array
            }
            else
            {
                sequenceEncoded = new int[0];
            }
            

            if (isNucleicAcid)
            {
                possibilities = 4;
                limit = (int)Math.Pow(possibilities, coupling); // 4 possibilities for nucleic acids
                sequenceCodeDict = nucleicAcidCodeDict;
            }
            else
            {
                possibilities = 20;
                limit = (int)Math.Pow(possibilities, coupling); // 20 possibilities for proteins
                sequenceCodeDict = proteinCodeDict;
            }

            LinkedList<int>[] sequenceCodePositions = new LinkedList<int>[limit]; // Check for efficiency: single reusable array
            
            for (int i = 0; i < coupledLength; i++)
            {
                // Produce the encoded value based on the current sequence position and a number of following positions based on the coupling parameter
                int code = 0;
                for (int j = 0; j < coupling; j++)
                {
                    int seqCode = sequenceCodeDict[sequence[i + j]]; // Look up the code value for the sequence character.
                    code += seqCode * (int)Math.Pow(possibilities, j); // add the character into the coupled code value.
                }

                // If the code hasn't been found before, create a linked list for it at the array index equal to the code value and then append the sequence position to the end.
                // The linked list will allow quick traversal of the positions in the sequence matching the given code.
                if (sequenceCodePositions[code] == null)
                {
                    sequenceCodePositions[code] = new LinkedList<int>();
                }
                sequenceCodePositions[code].AddLast(i);

                // Build the encoded sequence for output if requested.
                if (storeEncoded)
                {
                    sequenceEncoded[i] = code;
                }
            }

            return sequenceCodePositions;
        }

        static LinkedList<int>[] GetEncodedPositions(string sequence, int length, bool isNucleicAcid, int coupling)
        {
            int[] dummy;
            return GetEncodedPositions(sequence, length, isNucleicAcid, coupling, out dummy, false);
        }

        static int AlignPair(int[] sequence1Encoded, LinkedList<int>[] sequence1CodePositions, LinkedList<int>[] sequence2CodePositions, int sequence1Length, int sequence2Length, int coupling)
        {
            // Essentially what this is doing is scoring a pair of positions as a match if the displacement
            // required to align them is approximately equal to one of the majority shifts.
            // Given the nucleic acid sequences, encoded with coupling = 2:
            //     ACTTCCGA:  4  13  15   7   5   9   2
            //     GACTTCCG:  2   4  13  15   7   5   9
            //
            // First, the algorithm figures out how much sequence 2 would have to shift to align every matching code.
            //               -1  -1  -1  -1  -1  -1   6
            //     These are counted up and stored in an array whose index is the displacement amount offset
            //     by the coupled length of sequence 2 - 1, = 6 here.
            //         0  0  0  0  0  6  0  0  0  0  0  0  1
            //
            // Next, the algorithm identifies the top non-zero displacements (determined by the parameter
            //     "numDiagonalsToUse", = 5 here) and flags all those positions within a certain tolerance
            //     (given by the parameter "halfWindowSize", = 2 here).
            //         F  F  F  T  T  T  T  T  F  F  T  T  T
            //     There may be a better way to do this.
            //
            // Now that we know this, we start looping through the encoded sequence 1 positions, and we look up
            //     the last match location in sequence 2 and ask whether the necessary displacement is okay.
            //            *
            //            4  13  15   7   5   9   2
            //         2  4  13  15   7   5   9
            //     This requires sequence 2 to shift -1 relative to sequence 1, which, offset as above, = 5,
            //         which is allowed.
            //     Because this first position matches, we assign a score = coupling = 2.
            //     We store this score and the positions of sequence 1 and 2 that were aligned as an "AlignedFragment" in a list.
            //         We also store this fragment at the index associated with its offset displacement.
            //     Then we look up the next match in sequence 2 for this code. Since there aren't any, we move
            //         to the next sequence 1 position.
            //                *
            //            4  13  15  7   5   9   2
            //         2  4  13  15  7   5   9
            //     We repeat the displacement-checking and initial score, then load the last fragment in the list,
            //         which is the maximum-scoring fragment yet found.
            //     We then check whether this previous fragment is located earlier in the sequence. Since it is,
            //         we want to sum the previous score with the current one. But because we are coupling
            //         multiple bases, if the new coupled position overlaps with the previous one, we reduce
            //         the score by the amount of overlap.
            //     Here, the overlap is 1, so the new score is 2 + 1 = 3, and this is now added to the end of the list.
            //     This continues until the last coupled position (here, = 9), where the score = 7.
            //     Now we move to the last base of sequence 1, which aligns with the first base of sequence 2.
            //                                *
            //         4  13  15  7   5   9   2
            //                                2  4  13  15   7   5   9
            //     This offset displacement = 12, which is allowed. However, this matched set (6, 0) is not located
            //         after the maximum scoring fragment (5, 6), so we walk backward in the fragments collection to find
            //         one that is less than it. There won't be any with this sequence. So we keep the score of 2
            //         and insert this fragment in the list just before the first score greater than it.
            //     Two situations not covered by this example:
            //         If the test fragment is located prior to the current matched pair, but they aren't
            //         displaced by the same amount, we look back to see whether there is a prior fragment registered
            //         with the same displacement. We try to both build on this one if one is found or start a new
            //         fragment if not, and also build on the test fragment and include a gap penalty. We keep the
            //         maximum of these two.
            //
            // Finally, we return the maximum score found.

            // Trying this for another sequence pair, with halfWindowSize = 5:
            //     ACTTTCCG -> 4 13 15 15  7  5  9
            //     ACTCTCCG -> 4 13  7 13  7  5  9
            //     Displ: 4 = |  6 = |||||
            //     Slopes:     T  T  T  T  T  T  T  T  T  T  T  T  F
            //
            //     Round 1: 4 13 15 15  7  5  9
            //              4 13  7 13  7  5  9
            //              *
            //         Fragment added: 2, 0, 1 (score, seq1, seq2)
            //         Displacment[6] added: 2, 0, 1
            //         Max Fragment: 2, 0, 1
            //
            //     Round 2:       4 13 15 15  7  5  9
            //              4 13  7 13  7  5  9
            //                       *
            //         This position is (1, 3); the previous fragment was (0, 0), so the max fragment is not
            //         prior to it.
            //         Fragment added: 2, 1, 3
            //         Displacement[4] added: 2, 1, 3
            //         Max Fragment: 2, 1, 3
            //
            //     Round 3: 4 13 15 15  7  5  9
            //              4 13  7 13  7  5  9
            //                 *
            //         This position is (1, 1); the previous fragment was (1, 3), but this is not prior to it.
            //         The next previous fragment was (0, 0), and this one is prior; it is also the same displacement.
            //         Fragment added: 3, 1, 1
            //         Displacment[6] replaced: 3, 1, 1
            //         Max Fragment: 3, 1, 1
            //
            //     Round 4: There is no 15 in sequence 2
            //
            //     Round 5: There is no 15 in sequence 2
            //
            //     Round 6: 4 13 15 15  7  5  9
            //              4 13  7 13  7  5  9
            //                          *
            //         This position is (4, 4); the previous fragment was (1, 1), is prior, and is the same displacement.
            //         This position is also greater than 1 position away from the last, so we add a score of 2.
            //         Fragment added: 5, 4, 4
            //         Displacement[6] replaced: 5, 4, 4
            //         Max Fragment: 5, 4, 4
            //
            //     Round 7: 4 13 15 15  7  5  9
            //                    4 13  7 13  7  5  9
            //                          *
            //         This position is (4, 2); the previous fragment was (5, 4) but is not prior.
            //         The next prior is (1, 1), but it also is not prior.
            //         The next prior is (0, 0), and this is prior, but not the same displacement.
            //         Displacement = 8. No prior displacement is there.
            //         Subtotal1 = 2  Subtotal2 = 2 - 3 + 2 = 1
            //         Fragment inserted: 2, 4, 2
            //         Displacement[8] added: 2, 4, 2
            //         Max Fragment: 5, 4, 4
            //
            //     Round 8: 4 13 15 15  7  5  9
            //              4 13  7 13  7  5  9
            //                             *
            //         This position is (5, 5); the previous fragment was (4, 4), is prior and the same displacement.
            //         Fragment added: 6, 5, 5
            //         Displacement[6] replaced: 6, 5, 5
            //         Max Fragment: 6, 5, 5
            //
            //     Round 9: 4 13 15 15  7  5  9
            //              4 13  7 13  7  5  9
            //                                *
            //         This position is (6, 6); the previous fragment was (5, 5), is prior and the same displacement.
            //         Fragment added: 7, 6, 6
            //         Displacement[6] replaced: 7, 6, 6
            //         Max Fragment: 7, 6, 6
            //
            //     *************Return Score = 7
            //
            // Trying this for another sequence pair, with halfWindowSize = 5:
            //     ACTGGGTAGA -> 4 13 11 10 10 14  3  8  2
            //     ACTTCCGA   -> 4 13 15  7  5  9  2
            //     Displ: 6 = ||  8 = |
            //     Slopes:       F  T  T  T  T  T  T  T  T  T  T  T  T
            //
            //     Round 1: 4 13 11 10 10 14  3  8  2
            //              4 13 15  7  5  9  2
            //              *
            //         Fragment added: 2, 0, 0 (score, seq1, seq2)
            //         Displacment[6] added: 2, 0, 0
            //         Max Fragment: 2, 0, 0
            //
            //     Round 2: 4 13 11 10 10 14  3  8  2
            //              4 13 15  7  5  9  2
            //                 *
            //         This position is (1, 1); the previous fragment was (0, 0), so the max fragment is
            //         prior to it.
            //         Fragment added: 3, 1, 1 (score, seq1, seq2)
            //         Displacment[6] replaced: 3, 0, 0
            //         Max Fragment: 3, 1, 1
            //
            //     Round 3: 4 13 11 10 10 14  3  8  2
            //                    4 13 15  7  5  9  2
            //                                      *
            //         This position is (8, 6); the prevoius fragment was (1,1), so the max fragment is
            //         prior to it, but is not the same displacement.
            //         Subtotal 1: 2
            //         Subtotal 2: 3 - 3 + 2
            //         Fragment added: 2, 8, 6
            //         Displacement[8] added: 2, 8, 6
            //         Max Fragment: 3, 1, 1
            //
            //     Return 3

            //         This should become 0.75 in the distance matrix. = 1-0.75 = 0.25 score ratio = 1/4 = 2/8 = 3/16
            //         This highly suggests that Clustal is returning "2" for this sequence pair

            int limit = sequence1CodePositions.Count();
            int sequence1CoupledLength = sequence1Length - coupling + 1;
            int sequence2CoupledLength = sequence2Length - coupling + 1;
            int displacementIndexOffset = sequence2CoupledLength - 1;
            int displacementsArrayLength = displacementIndexOffset + sequence1CoupledLength;
            int[] displacements = new int[displacementsArrayLength]; // consider array memory usage. Though it is faster to create a new array than to loop through a current array to reset to 0.
            int numDiagonalsToUse = 5; // This would be the "signif" parameter from Clustal
            int halfWindowSize = 5; // This would be the "window" parameter from Clustal
            int gapPenalty = 3; // "windowGap" from Clustal; reduces score when a gap is introduced.
            bool[] slopes = new bool[displacementsArrayLength];
            LinkedList<AlignmentFragment> fragmentsList = new LinkedList<AlignmentFragment>();
            
            // For all the possible alignments of one code in one sequence and one code in another, count the
            // displacement that would be required to align them and store the counts in an array.
            for (int code = 0; code < limit; code++)
            {
                LinkedList<int> sequence1CodePosition = sequence1CodePositions[code];
                if (sequence1CodePosition == null) { continue; } // If it is null, none of this code is in the sequence.

                // Loop over all the positions in sequence 1 with this code value
                for (LinkedListNode<int> seq1PosNode = sequence1CodePosition.First; seq1PosNode != null; seq1PosNode = seq1PosNode.Next)
                {
                    LinkedList<int> sequence2CodePosition = sequence2CodePositions[code];
                    if (sequence2CodePosition == null) { break; } // If it is null, none of this code is in the sequence.

                    // Loop over all the positions in sequence 2 with this code value
                    for (LinkedListNode<int> seq2PosNode = sequence2CodePosition.First; seq2PosNode != null; seq2PosNode = seq2PosNode.Next)
                    {
                        // Calculate how much sequence 2 would need to shift to align these matched codes and store it in the displacements array, which is structured as follows:
                        // ----- 
                        // -------     Sequence 2 is not shifted at all. (Diff = Seq1Pos(0) - Seq2Pos(0) + DisplacementIndexOffset(6) = 6)
                        //
                        //       -----
                        // -------     Sequence 2 is shifted all the way to the left. (Diff = Seq1Pos(0) - Seq2Pos(6) + DisplacementIndexOffset(6) = 0)
                        //
                        // -----
                        //     ------- Sequence 2 is shifted all the way to the right. (Diff = Seq1Pos(4) - Seq2Pos(0) + DisplacementIndexOffset(6) = 10)
                        int displacementIndex = seq1PosNode.Value - seq2PosNode.Value + displacementIndexOffset;
                        displacements[displacementIndex]++;
                    }
                }
            }

            // We want to sort the displacements array so we can look up the top displacements,
            // But we also want to remember what the original indices were, so we make an array of the same size
            // populated with the indices and co-sort it based on the displacements matrix
            int[] displacementIndices = new int[displacementsArrayLength];
            for (int i = 0; i < displacementsArrayLength; i++)
            {
                displacementIndices[i] = i;
            }
            Array.Sort(displacements, displacementIndices);
            
            int displacementsArrayUpperBound = displacementsArrayLength - 1;
            for (int i = displacementsArrayUpperBound; i > displacementsArrayUpperBound - numDiagonalsToUse; i--)
            {
                if (displacements[i] == 0) { break; } // If we've hit a 0, we can stop.

                // Get the location of this displacement hit and find the bounds of the window around it.
                int displacementIndex = displacementIndices[i];
                int leftBound = Math.Max(0, displacementIndex - halfWindowSize); // Lock this to 0 if it would otherwise be negative.
                int rightBound = Math.Min(displacementsArrayUpperBound, displacementIndex + halfWindowSize); // Lock this to the length if it would otherwise be too large.

                // For each position within the bounds, set that position in the slopes array to 1. Also, I need a better name for the slopes array.
                for (int j = leftBound; j <= rightBound; j++)
                {
                    slopes[j] = true;
                }
            }

            LinkedListNode<AlignmentFragment>[] displacementFragmentNodes = new LinkedListNode<AlignmentFragment>[displacementsArrayLength];
            int couplingLess1 = coupling - 1; // Pre-calculate this value to save time in the loop.

            // Loop through the encoded sequence 1, and then look up the last position in sequence 2 with the same code and step through that set of positions.
            for (int i = 0; i < sequence1CoupledLength; i++)
            {
                int code = sequence1Encoded[i];

                LinkedList<int> sequence2CodePosition = sequence2CodePositions[code];
                if (sequence2CodePosition == null) { continue; } // If it is null, none of this code is in the sequence.

                // Loop over all the positions in sequence 2 with this code value
                //for (LinkedListNode<int> seq2PosNode = sequence2CodePosition.First; seq2PosNode != null; seq2PosNode = seq2PosNode.Next)
                for (LinkedListNode<int> seq2PosNode = sequence2CodePosition.Last; seq2PosNode != null; seq2PosNode = seq2PosNode.Previous) // Trying the reverse order, which is actually what Clustal does.
                {
                    int displacementIndex = i - seq2PosNode.Value + displacementIndexOffset;

                    if (slopes[displacementIndex]) // We only want to do something if this particular displacement was previously flagged.
                    {
                        int fragmentScore = coupling; // Stores the current score of the fragment.
                        LinkedListNode<AlignmentFragment> fragmentNode = fragmentsList.Last;
                        int seq2Pos = seq2PosNode.Value;

                        // A-Loop
                        // This won't execute the first time through.
                        while (fragmentNode != null)
                        {
                            int newSeq1Pos = fragmentNode.Value.Seq1Pos; //get the matched-pair sequence positions of the next fragment
                            int newSeq2Pos = fragmentNode.Value.Seq2Pos;
                            if (FragmentIsPrior(i, seq2Pos, newSeq1Pos, newSeq2Pos, coupling))
                            {
                                if (i - seq2Pos == newSeq1Pos - newSeq2Pos)
                                {
                                    fragmentScore = fragmentNode.Value.FragmentScore;
                                    if (i > newSeq1Pos + coupling - 1)
                                    {
                                        fragmentScore += coupling;
                                    }
                                    else
                                    {
                                        int remainder = i - newSeq1Pos;
                                        fragmentScore += remainder;
                                    }
                                }
                                else
                                {
                                    int subtotal1;
                                    if (displacementFragmentNodes[displacementIndex] == null)
                                    {
                                        subtotal1 = coupling;
                                    }
                                    else
                                    {
                                        LinkedListNode<AlignmentFragment> testFragmentNode = displacementFragmentNodes[displacementIndex];
                                        AlignmentFragment testFragment = testFragmentNode.Value;
                                        int testFragmentScore = testFragment.FragmentScore;
                                        int testSeq1Pos = testFragment.Seq1Pos;
                                        if (i > testSeq1Pos + couplingLess1)
                                        {
                                            subtotal1 = testFragmentScore + coupling;
                                        }
                                        else
                                        {
                                            int remainder = i - testSeq1Pos;
                                            subtotal1 = testFragmentScore + remainder;
                                        }
                                    }
                                    int subtotal2 = fragmentNode.Value.FragmentScore - gapPenalty + coupling;
                                    fragmentScore = Math.Max(subtotal1, subtotal2);
                                }
                                //fragmentNode = null; // I need to check how this is working... it is possible that there is no need to set this to null.
                                //Actually, I think I can just break here.
                                break;
                            }
                            else
                            {
                                fragmentNode = fragmentNode.Previous; // I originally had this as .Next, but the in Clustal, next was always finding the next lowest score, which would be previous here.
                            }
                        }

                        // Add the fragment to the list

                        AlignmentFragment currFragment = new AlignmentFragment(fragmentScore, i, seq2Pos);
                        LinkedListNode<AlignmentFragment> currFragmentNode;

                        if (fragmentsList.Last == null || fragmentsList.Last.Value.FragmentScore <= fragmentScore)
                        {
                            // If the list is empty or the last entry in the list has a smaller or equal score
                            fragmentsList.AddLast(currFragment);
                            currFragmentNode = fragmentsList.Last;
                        }
                        else
                        {
                            // From the highest score (end of the list), find the fragment in the list which
                            // has a lower or equal score and insert it after.
                            LinkedListNode<AlignmentFragment> testFragmentNode = fragmentsList.Last;
                            for (; testFragmentNode.Value.FragmentScore > fragmentScore; testFragmentNode = testFragmentNode.Previous);

                            fragmentsList.AddAfter(testFragmentNode, currFragment);
                            currFragmentNode = testFragmentNode.Next;
                        }
                        displacementFragmentNodes[displacementIndex] = currFragmentNode;
                    }
                }
            }

            if (fragmentsList.Last != null)
            {
                return fragmentsList.Last.Value.FragmentScore;
            }
            else
            {
                return 0;
            }           
        }

        //static Array FillArray(Array inputArray, object value)
        //{
        //    Parallel.For(0, inputArray.Length, i =>
        //    {
        //        inputArray.SetValue(value, i);
        //    });
        //    //for (int i = 0; i < inputArray.Length; i++)
        //    //{
        //    //    inputArray.SetValue(value, i);
        //    //}

        //    return inputArray;
        //}

        struct AlignmentFragment
        {
            public int FragmentScore;
            public int Seq1Pos;
            public int Seq2Pos;

            public AlignmentFragment(int fragmentScore, int seq1Pos, int seq2Pos)
            {
                FragmentScore = fragmentScore;
                Seq1Pos = seq1Pos;
                Seq2Pos = seq2Pos;
            }
        }
        //                               8       6       1       1
        static bool FragmentIsPrior(int a1, int b1, int a2, int b2, int coupling)
        // Determine whether the second set of positions are located prior to the first set of positions
        {
            if ((a1 - b1 == a2 - b2 && a2 < a1) ||
               (a2 + coupling - 1 < a1 && b2 + coupling - 1 < b1))
               // (8 - 6 == 1 - 1 && 1 < 8) || (1 + 2 - 1 < 8 && 1 + 2 - 1 < 6)
               // (    2 == 0     &&   T  ) || (        2 < 8 &&         2 < 6)
               // (       F       &&   T  ) || (          T   &&           T  )
               //                  F        ||                 T
               //                            T
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    
    }
}
