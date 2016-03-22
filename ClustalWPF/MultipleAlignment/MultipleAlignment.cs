using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;

namespace ClustalWPF.MultipleAlignment
{
    class MultipleAlignment
    {
        public static void AlignSequences(ref Alignment alignmentObject, SubstitutionMatrix.SubstitutionMatrixSeries subMatrix, ref double[,] similarityMatrix, ReadOnlyCollection<Tree.GuideTree<AlignedMacromolecule>.TreeLeaf> leavesList, ReadOnlyCollection<Tuple<ReadOnlyCollection<AlignedMacromolecule>, ReadOnlyCollection<AlignedMacromolecule>, double>> steps)
        {
            IMultipleAlignmentAlgorithm alignmentAlgorithm = new MyersMillerProfileAlign();
            int numSteps = steps.Count();
                        
            // Here, Clustal checks the similarity matrix and looks up the most closely-related sequence
            // for each sequence.
            // But it could have used the distance matrix and found the minimum distance. Consider switching to avoid unnecessary calculations

            foreach (Tuple<ReadOnlyCollection<AlignedMacromolecule>, ReadOnlyCollection<AlignedMacromolecule>, double> step in steps)
            {
                // Here, Clustal would check to make sure that at least one pair of sequences isn't too divergent.
                // It also records which ones are too divergent so that the alignment doesn't use them.
                
                double score = alignmentAlgorithm.Align(ref alignmentObject, subMatrix, step);
            }
        }
    }
}
