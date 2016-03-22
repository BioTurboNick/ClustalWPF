using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using ClustalWPF;

namespace ClustalWPF
{
    internal class Alignment: INotifyPropertyChanged //, IList<AlignedMacromolecule>
    // This class defines an Alignment object which holds the currently-loaded sequences and methods to manipulate them.
    {
        // A sequence alignment is really a collection of characters which have been assigned to positions relative to each other.
        // The alignment information is necessary for display and saving.
        //List<Tuple<Macromolecule, int[]>> alignedMacromolecules = new List<Tuple<Macromolecule, int[]>>(); // stores the sequences for the original unaligned sequences coupled to an array specifying the aligned locations of the sequence.
        List<AlignedMacromolecule> alignedMacromolecules = new List<AlignedMacromolecule>(); // stores the sequences for the original unaligned sequences and an array specifying the aligned locations of the sequence.
        int alignmentLength = 0;
        int maxAlignmentLength = 0; // Stores the maximum possible size of the alignment, which is 2 * longest sequence - 2;
        bool isNucleicAcid = false; // As currently coded, the most recently added macromolecule (or first of a range being added) sets this.
        bool isChanged = false; // Specifies whether the alignment has changed since it was last saved.
        //int count = 0;
        
        public ReadOnlyCollection<AlignedMacromolecule> AlignedMacromolecules // This property creates a read-only list. Anything that wants to add or remove a member of this list must go through AddAlignedMacromolecule or RemoveAlignedMacromolecule.
        {
            get { return alignedMacromolecules.AsReadOnly(); }
        }

        public int AlignmentLength // Setting this property affects the display and alignment mode, not the sequences loaded.
        {
            get { return alignmentLength; }
            set { alignmentLength = value; NotifyPropertyChanged("AlignmentLength"); }
        }

        public bool IsNucleicAcid
        {
            get { return isNucleicAcid; }
            set { isNucleicAcid = value; NotifyPropertyChanged("IsNucleicAcid"); }
        }

        public bool IsChanged
        {
            get { return isChanged; }
            set { isChanged = value; NotifyPropertyChanged("IsChanged"); }
        }

        public int NumberMacromolecules
        {
            get { return AlignedMacromolecules.Count; }
        }

        public int MaxAlignmentLength
        {
            get { return maxAlignmentLength; }
        }

        //internal void NudgeSequences(List<int> sequenceIndices, int positionStart, int positionEnd, int nudgeAmount)
        //// Moves the sequences given by the sequenceIndices at the given absolute position in relation to the others by the amount given in nudgeAmount.
        //{
        //    NudgeSequences(alignedMacromolecules, sequenceIndices, positionStart, positionEnd, nudgeAmount);
            
        //    // Ensure that the minimum offset is 0.
        //    NormalizeAlignedPositions();
        //    // Recalculate the length of the alignment.
        //    CalculateAlignmentLength();
        //}

        static void NudgeSequences(List<AlignedMacromolecule> macromolecules, int positionStart, int positionEnd, int nudgeAmount)
        // Moves the sequences given by the sequenceIndices at the given absolute position in relation to the others by the amount given in nudgeAmount.
        {
            foreach (AlignedMacromolecule macromolecule in macromolecules)
            {
                //int[] seqAlignedPos = alignedMacromolecules[index].Item2;
                int[] seqAlignedPos = macromolecule.AlignedPositions;
                int sequenceLength = seqAlignedPos.Length;

                // Find the first offset within or after the nudged region.
                // If none of the offsets are within or beyond the start position,
                // first will default to a value beyond the offset list.
                int first = sequenceLength;
                for (int i = 0; i < sequenceLength; i++)
                {
                    if (seqAlignedPos[i] >= positionStart)
                    {
                        first = i;
                        break;
                    }
                }

                int last = sequenceLength - 1;
                if (first < sequenceLength) // Only do if first is before the end of the sequence.
                {
                    // Find the first offset before the end of the nudged region.
                    // If none of the offsets are before the end position,
                    // last will equal the index of the last member of the string.
                    // if last is also prior to the first sequence, this loop will assign to it -1.
                    for (int i = first; i < sequenceLength; i++)
                    {
                        if (seqAlignedPos[i] > positionEnd)
                        {
                            last = i - 1;
                            break;
                        }
                    }
                }
                // Positions have been found.
                // Special cases:
                // first == 0, last == -1:
                //     the shifting positions are before the start of the sequence.
                // first == sequenceOffset.Count, last == sequenceOffset.Count - 1:
                //     the shifting positions are after the end of the sequence.

                if (last < first) // occurs when the positions are beyond an end of the sequence and thus all should be moved
                {
                    NudgeWholeSequence(macromolecule, nudgeAmount);
                }
                else
                {
                    // Adjust the offsets of the bases within the region
                    for (int i = first; i <= last; i++)
                    {
                        seqAlignedPos[i] += nudgeAmount;
                    }

                    // Fix collisions at the edges. Causes gaps in the direction of movement to be closed up.
                    for (int i = first; i > 0; i--)
                    {
                        if (seqAlignedPos[i] <= seqAlignedPos[i - 1])
                        {
                            seqAlignedPos[i - 1] = seqAlignedPos[i] - 1;
                        }
                        else
                        {
                            break;
                        }
                    }

                    for (int i = last; i < sequenceLength - 1; i++)
                    {
                        if (seqAlignedPos[i] >= seqAlignedPos[i + 1])
                        {
                            seqAlignedPos[i + 1] = seqAlignedPos[i] + 1;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }

        static void NudgeSequences(List<AlignedMacromolecule> macromolecules, int positionStart, int nudgeAmount)
        // nudge entire sequence to the right of the start position.
        {
            NudgeSequences(macromolecules, positionStart, int.MaxValue, nudgeAmount);
        }

        static void NudgeWholeSequence(AlignedMacromolecule macromolecule, int nudgeAmount)
        // Moves the entire sequences given by the sequenceIndices by the nudgeAmount.
        {
            //int[] seqAlignedPos = alignedMacromolecules[index].Item2;
            int[] seqAlignedPos = macromolecule.AlignedPositions;

            for (int i = 0; i < seqAlignedPos.Length; i++)
            {
                seqAlignedPos[i] += nudgeAmount;
            }
        }

        void NormalizeAlignedPositions()
        // Sets the minimum offset to 0 and recalculates the maximum alignment length.
        {
            throw new NotImplementedException();
            //// Determine the minimum offset value
            //int minAlignedPos = 0;
            ////foreach (Tuple<Macromolecule, int[]> seqTuple in alignedMacromolecules)
            //foreach (AlignedMacromolecule alignedMacromolecule in alignedMacromolecules)
            //{
            //    minAlignedPos = Math.Min(minAlignedPos, alignedMacromolecule.AlignedPositions[0]);
            //}

            //// If the minimum offset value is not zero, adjust the offsets so that it is.
            //if (minAlignedPos != 0)
            //{
            //    for (int i = 0; i < alignedMacromolecules.Count; i++)
            //    {
            //        NudgeWholeSequence(i, -minAlignedPos);
            //    }

            //    // Recalculate the length of the alignment
            //    CalculateAlignmentLength();
            //}
        }

        void CalculateAlignmentLength()
        // Calculate the length of the full alignment
        // Should make this automatic somehow...
        {
            int maxLength = 0;
            //foreach (Tuple<Macromolecule, int[]> alignedMacromolecule in alignedMacromolecules)
            foreach (AlignedMacromolecule alignedMacromolecule in alignedMacromolecules)
            {
                maxLength = Math.Max(maxLength, alignedMacromolecule.AlignedPositions.Max() + 1);
            }
            alignmentLength = maxLength;
            NotifyPropertyChanged("AlignmentLength");
        }

        internal void Clear()
        // Reset the alignment
        {
            alignedMacromolecules.Clear();
            alignmentLength = 0;
            isChanged = false;
        }

        void RecalculateMaximumLength()
        {
            int maxLength = 0;

            foreach (AlignedMacromolecule macromolecule in alignedMacromolecules)
            {
                maxLength = Math.Max(maxLength, macromolecule.Sequence.Length);
            }

            maxAlignmentLength = 2 * maxLength - 2;
        }

        //internal void AddAlignedMacromolecule(Tuple<Macromolecule, int[]> newAlignedMacromolecule, bool doAdjustments = true)
        internal void AddAlignedMacromolecule(AlignedMacromolecule newAlignedMacromolecule, bool doAdjustments = true)
        // Provides an interface to add an aligned macromolecule to the list.
        // Adjust the name if it is not unique, adjust the alignmentLength, and add it to the list.
        {
            //Macromolecule newMacromolecule = newAlignedMacromolecule.Item1;
            Macromolecule newMacromolecule = newAlignedMacromolecule;

            // Rename the macromolecule if another with the same name is already loaded.
            // Appends a "1" if already exists; then if this new name matches another, add "2", etc.
            // Could lead to some weirdness if the loaded sequences were "seq1" "seq4" "seq4", which would
            // become "seq1", "seq4", "seq4 1", but this is tolerable.
            int nameInteger = 1;
            for (int i = 0; i < alignedMacromolecules.Count; i++)
            {
                //Macromolecule testMacromolecule = alignedMacromolecules[i].Item1;
                Macromolecule testMacromolecule = alignedMacromolecules[i];
                if (newMacromolecule.Name == testMacromolecule.Name) // same name
                {
                    newMacromolecule.Name += " " + nameInteger.ToString();
                    nameInteger++;
                }
            }

            alignedMacromolecules.Add(newAlignedMacromolecule);

            if (doAdjustments)
            { // Make into another sub so the Add routines can be synced?
                NormalizeAlignedPositions();
                CalculateAlignmentLength();
                NotifyPropertyChanged("AlignedMacromolecules");
                //isNucleicAcid = newAlignedMacromolecule.Item1.IsNucleicAcid;
                isNucleicAcid = newAlignedMacromolecule.IsNucleicAcid;
                NotifyPropertyChanged("IsNucleicAcid");
                RecalculateMaximumLength();
            }
        }

        //internal void AddAlignedMacromolecules(List<Tuple<Macromolecule, int[]>> newAlignedMacromolecules)
        internal void AddAlignedMacromolecules(List<AlignedMacromolecule> newAlignedMacromolecules)
        // Provides an interface to add multiple aligned macromolecule to the list.
        {
            //foreach (Tuple<Macromolecule, int[]> newAlignedMacromolecule in newAlignedMacromolecules)
            foreach (AlignedMacromolecule newAlignedMacromolecule in newAlignedMacromolecules)
            {
                AddAlignedMacromolecule(newAlignedMacromolecule, false);
            }

            NormalizeAlignedPositions();
            CalculateAlignmentLength();
            NotifyPropertyChanged("AlignedMacromolecules");
            //IsNucleicAcid = newAlignedMacromolecules[0].Item1.IsNucleicAcid;
            IsNucleicAcid = newAlignedMacromolecules[0].IsNucleicAcid;
            NotifyPropertyChanged("IsNucleicAcid");
            RecalculateMaximumLength();
        }

        //internal void RemoveAlignedMacromolecule(Tuple<Macromolecule, int[]> alignedMacromoleculeToRemove, bool doAdjustments = true)
        internal void RemoveAlignedMacromolecule(AlignedMacromolecule alignedMacromoleculeToRemove, bool doAdjustments = true)
        // Provides an interface to remove an aligned macromolecule from the list.
        {
            alignedMacromolecules.Remove(alignedMacromoleculeToRemove);

            if (doAdjustments)
            {
                NormalizeAlignedPositions();
                CalculateAlignmentLength();
                NotifyPropertyChanged("AlignedMacromolecules");
                RecalculateMaximumLength();
            }
        }

        //internal void RemoveAlignedMacromolecules(List<Tuple<Macromolecule, int[]>> alignedMacromoleculesToRemove)
        internal void RemoveAlignedMacromolecules(List<AlignedMacromolecule> alignedMacromoleculesToRemove)
        // Provides an interface to remove an aligned macromolecule from the list.
        {
            //foreach (Tuple<Macromolecule, int[]> alignedMacromoleculeToRemove in alignedMacromolecules)
            foreach (AlignedMacromolecule alignedMacromoleculeToRemove in alignedMacromolecules)
            {
                RemoveAlignedMacromolecule(alignedMacromoleculeToRemove, false);
            }

            NormalizeAlignedPositions();
            CalculateAlignmentLength();
            NotifyPropertyChanged("AlignedMacromolecules");
            RecalculateMaximumLength();
        }

        internal static string ToAlignedSequence(string sequence, int[] alignedPositions, int alignmentLength)
        // Returns a string representing the fully aligned sequence given an original sequence, a matching array
        // representing the aligned positions, and the length of the fully aligned sequence.
        {
            StringBuilder alignedSequenceStringBuilder = new StringBuilder();

            for (int seqIndex = 0, i = 0; seqIndex < alignmentLength; seqIndex++)
            {
                if (i < alignedPositions.Length && alignedPositions[i] == seqIndex)
                {
                    alignedSequenceStringBuilder.Append(sequence[i]);
                    i++;
                }
                else
                {
                    alignedSequenceStringBuilder.Append('-');
                }
            }

            return alignedSequenceStringBuilder.ToString();
        }

        internal static string ToSequence(string alignedSequence, out int[] alignedPositions)
        // Returns a string containing the no-gap sequence and as an output parameter an array of
        // integers representing the aligned positions of the sequence.
        {
            StringBuilder sequenceStringBuilder = new StringBuilder();
            List<int> alignedPosList = new List<int>();

            for (int seqIndex = 0; seqIndex < alignedSequence.Length; seqIndex++)
            {
                if (alignedSequence[seqIndex] != '-')
                {
                    sequenceStringBuilder.Append(alignedSequence[seqIndex]);
                    alignedPosList.Add(seqIndex);
                }
            }

            alignedPositions = alignedPosList.ToArray();

            return sequenceStringBuilder.ToString();
        }

        public event PropertyChangedEventHandler PropertyChanged;

        private void NotifyPropertyChanged(string info)
        {
            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs(info));
            }
        }

        internal static void RemoveRedundantGaps(List<AlignedMacromolecule> macromolecules)
        {
            // Remove any gaps which appear in all of the provided sequences


            // Find the maximum length
            int maxLength = 0;
            foreach (AlignedMacromolecule macromolecule in macromolecules)
            {
                maxLength = Math.Max(maxLength, macromolecule.AlignedPositions.Last());
            }
            
            bool[] residuePresent = new bool[maxLength];

            // Find all positions for which residues exist
            foreach (AlignedMacromolecule macromolecule in macromolecules)
            {
                for (int i = 0; i < macromolecule.AlignedPositions.Length; i++)
                {
                    residuePresent[macromolecule.AlignedPositions[i]] = true;
                }
            }

            // Nudge the sequences to remove unnecessary gaps
            for (int i = 0; i < maxLength; i++)
            {
                int offset = 0; // keep track of the number of shifts done
                if (!residuePresent[i])
                {
                    NudgeSequences(macromolecules, i + 1 - offset, -1);
                }
            }
        }
    }
}
