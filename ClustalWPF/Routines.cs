using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.IO;
using System.Windows;
using System.Windows.Media;

namespace ClustalWPF
{
    internal class Routines
    {
        public const char oldGap = '-';
        public const char newGap = '_';
        public static readonly char[] nucleotideCodes = new char[] { 'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y' };
        public static readonly char[] nucleotideCodesPlusGaps = new char[] { 'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', oldGap, newGap };
        public static readonly char[] nondegenerateNucleotideCodes = new char[] { 'A', 'C', 'G', 'T', 'U' };
        public static readonly char[] degenerateNucleotideCodes = new char[] { 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'X', 'Y' };
        public static readonly char[] aminoAcidCodes = new char[] { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z' };
        public static readonly char[] aminoAcidCodesPlusGaps = new char[] { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', oldGap, newGap };
        public static readonly char[] nondegenerateAminoAcidCodes = new char[] { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
        public static readonly char[] degenerateAminoAcidCodes = new char[] { 'B', 'X', 'Z' };


        internal static bool ValidateConverterInput(object values)
        // Call this at the start of every IValueConverter or IMultiValueConverter subroutine.
        // When the data binding's sources are set at runtime, the designer will be sending
        // DependencyProperty.UnsetValue instead. This will check all values passed to it and
        // return true if all values are set, or false if at least one value is unset.
        {
            try // Try to treat the input object as an array
            {
                object[] valuesArray = (object[])values;

                foreach (object value in valuesArray)
                {
                    if (value == System.Windows.DependencyProperty.UnsetValue)
                    {
                        return false;
                    }
                }
            }
            catch // Input object was not an array, so treat as a single object.
            {
                if (values == System.Windows.DependencyProperty.UnsetValue)
                {
                    return false;
                }
            }
            return true;

        }

        // TODO Currently I am dynamically allocating arrays. I should just figure out how big to make them and then work with that.

        internal static double GetTextDisplayLength(string text, FontFamily fontFamily, FontStyle fontStyle, FontWeight fontWeight, FontStretch fontStretch, int fontSize)
        {
            Typeface typeface = new Typeface(fontFamily, fontStyle, fontWeight, fontStretch);
            FormattedText ft = new FormattedText(text, System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, typeface, fontSize, Brushes.White);
            return ft.WidthIncludingTrailingWhitespace;
        }

        internal static bool IsNucleicAcid(String sequence)
            // redo this
        // When it isn't otherwise specified what macromolecule type is being used, we can guess based on the code used.
        // We count the number of sequence positions which are a standard nucleic acid (A, C, G, T). If a code unique
        // to nucleic acid (B, U) or unique to protein (E, F, L, P, Q, X) is encountered, return true or false, respectively.
        // If no unique codes are observed, if >85% are standard nucleic acid codes, return true; else return false.
        {
            int NACount = 0;
            foreach (char seqChar in sequence.ToUpperInvariant())
            {
                if (seqChar == 'A' || seqChar == 'C' || seqChar == 'G' || seqChar == 'T')
                {
                    NACount++;
                }
                else if (seqChar == 'U' || seqChar == 'B')
                {
                    return true; // If even a single sequence code is unique to DNA, treat it like DNA.
                }
                else if (seqChar == 'E' || seqChar == 'F' || seqChar == 'L' || seqChar == 'P' || seqChar == 'Q' || seqChar == 'X')
                {
                    return false; // If even a single sequence code is unique to protein, treat it like a protein.
                }
            }

            if (NACount / sequence.Length >= 0.85) // Treat the sequence as nucleic acid with some degenerate nucleic acid codes.
            {
                return true;
            }
            else // Treat sequence as protein.
            {
                return false;
            }
        }

        internal static void Align(ref Alignment alignmentObject, ref BackgroundWorker myBackgroundWorker)
        // Runs complete Clustal multiple sequence alignment on all loaded molecules
        {
            // Currently not set up to save files during alignment.
            // Also not set to use the "stats" option

            if (alignmentObject.IsNucleicAcid)
            {
                // Load alignment parameters for nucleic acids
            }
            else
            {
                // Load alignment parameters for protein
            }

            // Step 1: Do pairwise alignments and generate the distance matrix
            myBackgroundWorker.ReportProgress(0, "Starting pairwise alignments...");

            int numSequences = alignmentObject.NumberMacromolecules;
            double[,] distanceMatrix = new double[numSequences, numSequences];

            // This will select between fast and full pairwise alignment. Just implementing fast alignment to start.
            PairwiseAlignment.IPairwiseAlignmentAlgorithm pairwiseAlignmentAlgorithm;
            if (true)
            {
                // This sub basically calculates the score of the best alignment between each pair of sequences
                pairwiseAlignmentAlgorithm = new PairwiseAlignment.FastPairwiseAlignment();
            }
            else
            {
                throw new NotImplementedException("Full Pairwise Alignment is not yet implemented");
                // pairwiseAlignmentAlgorithm = new PairwiseAlignment.FullPairwiseAlignment();
            }
            pairwiseAlignmentAlgorithm.PairwiseAlign(ref alignmentObject, ref distanceMatrix);

            // Step 2: Generate a phylogenetic tree of the sequences from the distance matrix
            myBackgroundWorker.ReportProgress(0, "Generating phylogenetic tree...");

            AlignedMacromolecule[] macromolecules = new AlignedMacromolecule[alignmentObject.NumberMacromolecules];
            ReadOnlyCollection<AlignedMacromolecule> alignedMacromolecules = alignmentObject.AlignedMacromolecules;
            for (int i = 0; i < alignmentObject.NumberMacromolecules; i++)
            {
                macromolecules[i] = alignedMacromolecules[i];
            }

            Tree.GuideTree<AlignedMacromolecule> tree = Tree.GuideTree<AlignedMacromolecule>.GetWeightsAndSteps(ref distanceMatrix, macromolecules);

            // Transfer the weights from the tree to the macromolecules and fetch the similarities and steps
            foreach (Tree.GuideTree<AlignedMacromolecule>.TreeLeaf leaf in tree.LeavesList)
            {
                AlignedMacromolecule macromolecule = (AlignedMacromolecule)leaf.Data;
                macromolecule.SequenceWeight = leaf.weight;
            }

            ReadOnlyCollection<Tuple<ReadOnlyCollection<AlignedMacromolecule>, ReadOnlyCollection<AlignedMacromolecule>, double>> alignmentSteps = tree.Steps; 
            
            // Step 3: Multiple sequence alignment
            myBackgroundWorker.ReportProgress(0, "Starting multiple sequence alignment...");

            SubstitutionMatrix.SubstitutionMatrixSeries subMatrixType;

            if (alignmentObject.IsNucleicAcid) // Replaced by more advanced logic later
            {
                subMatrixType = new SubstitutionMatrix.NucleotideIdentity();
            }
            else
            {
                subMatrixType = new SubstitutionMatrix.BLOSUM();
            }
            
            MultipleAlignment.MultipleAlignment.AlignSequences(ref alignmentObject, subMatrixType, ref tree.similarityMatrix, tree.LeavesList, alignmentSteps);
            
            
        }
    }

    public struct FontSettings
    // Fully describes a font
    {
        public Typeface Typeface;
        public int FontSize;        

        public FontSettings(Typeface newTypeface, int fontSize)
        {
            Typeface = newTypeface;
            FontSize = fontSize;
        }
    }

    //public struct AlignedMacromolecule
    //// Associates a macromolecule with an array indicating the aligned position of each member of its sequence.
    //{
    //    Macromolecule macromolecule;
    //    int[] alignedPositions;
    //    double sequenceWeight;

    //    public Macromolecule Macromolecule
    //    {
    //        get { return macromolecule; }
    //        set { macromolecule = value; }
    //    }

    //    public int[] AlignedPositions
    //    {
    //        get { return alignedPositions; }
    //        set { alignedPositions = value; }
    //    }

    //    public double SequenceWeight
    //    {
    //        get { return sequenceWeight; }
    //        set { sequenceWeight = value; }
    //    }

    //    public AlignedMacromolecule(Macromolecule newMacromolecule, int[] newAlignedPositions)
    //    {
    //        macromolecule = newMacromolecule;
    //        alignedPositions = newAlignedPositions;
    //        sequenceWeight = 0;
    //    }
    //}

    public class AlignedMacromolecule : Macromolecule
    {
        int[] alignedPositions;
        List<int> newGaps; // Store the locations of gaps which have been newly inserted
        double sequenceWeight;

        public int[] AlignedPositions
        {
            get { return alignedPositions; }
            set { alignedPositions = value; }
        }

        public List<int> NewGaps
        {
            get { return newGaps; }
            set { newGaps = value; }
        }

        public double SequenceWeight
        {
            get { return sequenceWeight; }
            set { sequenceWeight = value; }
        }

        public AlignedMacromolecule(Macromolecule newMacromolecule, int[] newAlignedPositions)
            : base(newMacromolecule)
        {
            alignedPositions = newAlignedPositions;
            sequenceWeight = 0;
            newGaps = new List<int>();
        }


        public AlignedMacromolecule(Macromolecule newMacromolecule)
            : base(newMacromolecule)
        {
            alignedPositions = new int[newMacromolecule.Sequence.Length];

            for (int i = 0; i < newMacromolecule.Sequence.Length; i++)
            {
                alignedPositions[i] = i;
            }

            sequenceWeight = 0;
        }

        public AlignedMacromolecule Copy()
        {
            //Produces a copy of the alignment information only, and sets sequemceWeight to 0
            AlignedMacromolecule copy = new AlignedMacromolecule(this, (int[])alignedPositions.Clone());

            foreach (int newGap in newGaps)
            {
                copy.newGaps.Add(newGap);
            }

            sequenceWeight = 0;

            return copy;
        }

        public void ClearGaps()
        {
            alignedPositions = new int[Sequence.Length];
            newGaps = new List<int>();
        }
    }

    enum ReturnCodes
    {
        OK,
        NoFileName,
        FileNotFound,
        NoSequencesInFile,
        EmptySequencesRemoved,
        NoNonEmptySequencesInFile,
        FileErrors,
        UnrecognizedFileType,
        IOException,
        MixedToNA,
        MixedToProtein
    }

    
}
