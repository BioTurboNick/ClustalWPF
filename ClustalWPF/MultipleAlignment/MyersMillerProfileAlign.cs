using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Runtime.CompilerServices;
using System.Linq;
using System.Text;

namespace ClustalWPF.MultipleAlignment
{
    class MyersMillerProfileAlign : IMultipleAlignmentAlgorithm
    {
        int maxLength = 0; // Maximum length of the alignment;
        int lastPrint = 0;
        int printPtr = 1;
        int[] displ;
        int[] alnPath1;
        int[] alnPath2;
        int alignmentLength;
        Profile group1Profile = new Profile();
        Profile group2Profile = new Profile();
        char[] residueCodes = Routines.nucleotideCodesPlusGaps; // should be one with gaps, need to do programatically

        public double Align(ref Alignment alignmentObject, SubstitutionMatrix.SubstitutionMatrixSeries subMatrixClass, Tuple<ReadOnlyCollection<AlignedMacromolecule>, ReadOnlyCollection<AlignedMacromolecule>, double> step)
        {
            bool switchGroups;
            ReadOnlyCollection<AlignedMacromolecule> group1;
            ReadOnlyCollection<AlignedMacromolecule> group2;
            double meanSimilarity = step.Item3;

            int numSeqs1 = step.Item1.Count();
            int numSeqs2 = step.Item2.Count();
            int numSeqs = alignmentObject.NumberMacromolecules;

            if (numSeqs1 == 0 || numSeqs2 == 0)
                // If one of the groups has no alignable sequences, can't do anything!
            {
                return 0;
            }

                        
            // Make the group with the most sequences group 1
            // Figure out the structure penalties here
            if (switchGroups = numSeqs2 > numSeqs1)
            {
                group1 = step.Item2;
                group2 = step.Item1;
            }
            else
            {
                group1 = step.Item1;
                group2 = step.Item2;
            }           

            // Make the first group
            int group1MaxLength = 0;
            foreach (AlignedMacromolecule macromolecule in group1)
            {
                group1MaxLength = Math.Max(group1MaxLength, macromolecule.AlignedPositions.Last() + 1);
            }


            // Make the second group
            int group2MaxLength = 0;
            foreach (AlignedMacromolecule macromolecule in group2)
            {
                group2MaxLength = Math.Max(group2MaxLength, macromolecule.AlignedPositions.Last() + 1);
            }

            maxLength = group1MaxLength + group2MaxLength + 2; // Clustal sets the Alignment parameter here.
            displ = new int[maxLength + 1];
            alnPath1 = new int[maxLength + 1];
            alnPath2 = new int[maxLength + 1];

            
            // Calculate the real length of profiles, removing gaps:

            int group1NoGapsLength = 0;
            foreach (AlignedMacromolecule macromolecule in group1)
            {
                group1NoGapsLength += macromolecule.Sequence.Length;
            }
            group1NoGapsLength = (int)(group1NoGapsLength / (double)numSeqs1);

            int group2NoGapsLength = 0;
            foreach (AlignedMacromolecule macromolecule in group2)
            {
                group2NoGapsLength += macromolecule.Sequence.Length;
            }
            group2NoGapsLength = (int)(group2NoGapsLength / (double)numSeqs2);

            int minNoGapsLength = Math.Min(group1NoGapsLength, group2NoGapsLength);
            int maxNoGapsLength = Math.Max(group1NoGapsLength, group2NoGapsLength);
           
            // I'm just going to reproduce the Clustal code directly. May recode later once I fully understand what it does.
            
            // There is quite literally no point to the scaleVals.scale value, and no point to setting this here.
            // Tuple<double, double> scaleVals = Tuple.Create(1.0, 100.0); // scaleVals.scale and scaleVals.intScale, respectively.
                            
            double gapOpeningCoefficient, gapExtensionCoefficient;
            double gapOpenPenalty = 10.0; // Will be user selectable parameter
            double gapExtendPenalty = 0.2; // Will be user selectable parameter

            //Will be user option
            bool useNegative = false;

            SubstitutionMatrix.SubstitutionMatrix subMatrix = subMatrixClass.GetMatrix(meanSimilarity, minNoGapsLength, useNegative);

            double gapPenaltyScaleFactor = subMatrixClass.GetScaleFactor(meanSimilarity, useNegative);

            // Then it sets DNA vs. Protein parameters
            if (alignmentObject.IsNucleicAcid)
            {
                gapOpeningCoefficient = 100.0 * gapOpenPenalty * gapPenaltyScaleFactor;
                gapExtensionCoefficient = 100.0 * gapExtendPenalty * gapPenaltyScaleFactor;
            }
            else
            {
                gapOpeningCoefficient = CalculateProteinGapOpeningCoefficient(gapOpenPenalty, group1NoGapsLength, group2NoGapsLength, subMatrix.AverageScore, gapPenaltyScaleFactor);
                gapExtensionCoefficient = 100.0 * gapExtendPenalty;
            }            

            // We need one profile with substitution matrix information and one without!
            // But this will change when we have the LE scoring function. (Whatever that is.)

            // Calculate the profile arrays. The first group, but not the second, incorporates substitution information.
            group1Profile = Profile.Calculate(group1, group1MaxLength, gapOpeningCoefficient, gapExtensionCoefficient, subMatrix);
            group2Profile = Profile.Calculate(group2, group2MaxLength, gapOpeningCoefficient, gapExtensionCoefficient, null); //new SubstitutionMatrix());

            // User Myers and Miller to align the two sequences
            double score = ProgDiff(0, 0, group1MaxLength, group2MaxLength, group1Profile.GapOpeningPenalties[0], group1Profile.GapExtensionPenalties[0]);

            alignmentLength = ProgTracepath();

            addGGaps();

            group1MaxLength = alignmentLength;

            if (true) // Clustal: DoRemoveFirstIteration() == TREE ???
            {
                // Combine the sequences and submit for iterations 
            }



            return 0;
        }

        private void addGGaps()
        {
            foreach (AlignedMacromolecule macromolecule in group1Profile.Macromolecules)
            {
                for (int i = 0; i < alignmentLength; i++)
                {
                    if (alnPath1[i] == 1)
                    {
                        // move logic to Alignment? Or to AlignedMacromolecule?
                        for (int j = 0; j < macromolecule.Sequence.Length; j++)
                        {
                            if (macromolecule.AlignedPositions[j] > i - 1)
                            {
                                macromolecule.AlignedPositions[j]++;
                                macromolecule.NewGaps.Add(i - 1);
                            }
                        }
                    }
                }
            }

            foreach (AlignedMacromolecule macromolecule in group2Profile.Macromolecules)
            {
                for (int i = 0; i < alignmentLength; i++)
                {
                    if (alnPath2[i] == 1)
                    {
                        // move logic to Alignment? Or to AlignedMacromolecule?
                        for (int j = 0; j < macromolecule.Sequence.Length; j++)
                        {
                            if (macromolecule.AlignedPositions[j] > i - 1)
                            {
                                macromolecule.AlignedPositions[j]++;
                                macromolecule.NewGaps.Add(i - 1);
                            }
                        }
                    }
                }
            }
        }

        private int ProgTracepath()
        {
            int k, pos, toDo;
            int alignLen;
            pos = 0;

            toDo = printPtr - 1;

            for (int i = 1; i <= toDo; i++)
            {
                if (displ[i] == 0)
                {
                    alnPath1[pos] = 2; // indicates that the position should be recorded as-is
                    alnPath2[pos] = 2;
                    pos++;
                }
                else
                {
                    if ((k = displ[i]) > 0)
                    {
                        for (int j = 0; j <= k - 1; j++)
                        {
                            alnPath2[pos + j] = 2;
                            alnPath1[pos + j] = 1; // indicates that the position should have a gap
                        }
                        pos += k;
                    }
                    else
                    {
                        k = (displ[i] < 0) ? displ[i] * -1 : displ[i];
                        for (int j = 0; j <= k - 1; j++)
                        {
                            alnPath1[pos + j] = 2;
                            alnPath2[pos + j] = 1;
                        }
                        pos += k;
                    }
                }
            }

            alignLen = pos;
            return alignLen;
        }

        private double CalculateProteinGapOpeningCoefficient(double gapOpenPenalty, int group1NoGapsLength, int group2NoGapsLength, double subMatrixAverageScore, double scaleFactor)
        {
            double logmin = 1.0;
            double logdiff = 1.0;
            double gapOpeningCoefficient;

            if (group1NoGapsLength == 0 || group2NoGapsLength == 0)
            {
                int minNoGapsLength = Math.Min(group1NoGapsLength, group2NoGapsLength);
                logmin = 1.0 / Math.Log10((double)minNoGapsLength);
                if (group1NoGapsLength != group2NoGapsLength)
                {
                    int maxNoGapsLength = Math.Max(group1NoGapsLength, group2NoGapsLength);
                    logdiff = 1.0 + 0.5 * Math.Log10((double)minNoGapsLength / (double)maxNoGapsLength);
                }
                else
                {
                    logdiff = 1.0;
                }
                if (logdiff < 0.9)
                {
                    logdiff = 0.9;
                }
            }

            bool isNegative = false; // Will be user parameter
            if (isNegative)
            {
                gapOpeningCoefficient = 100.0 * gapOpenPenalty;
            }
            else
            {
                if (subMatrixAverageScore <= 0)
                {
                    gapOpeningCoefficient = 100.0 * gapOpenPenalty + logmin;
                }
                else
                {
                    gapOpeningCoefficient = scaleFactor * subMatrixAverageScore * (gapOpenPenalty / (logdiff * logmin));
                }
            }

            return gapOpeningCoefficient;
        }

        int pd_t, pd_tl, pd_g, pd_h;
        int pd_i, pd_j;
        int pd_f, pd_e, pd_s;

        // sb1 = sb2 = 0    se1 = length1    se2 = length2
       //  sb1, sb2, se1 - sb1, se2 - sb2, start gap open penalty, end gap open penalty

        // whew, might actually need this to be fully recursive with bits accessible...
        double t, tl, g, h;
        double hh, f, e, s;
        double[] HH;
        double[] DD;
        double[] RR;
        double[] SS;
        double[] gS;
        double ProgDiff(int A, int B, int M, int N, double startGapOpenPenalty, double endGapOpenPenalty)
        {
            // for now pretty much just reproducing the Clustal code because I
            // don't entirely understand it.
            // I'll need to check the indices to make sure they are looking up
            // the right things.

            int midi, midj, gapType;
            double midh;
            //double hh;

            if (N <= 0) // if profile B is empty...
            {
                if (M > 0) // ...and if profile A is not empty
                {
                    // Delete the residues between A[1] and A[M]
                    ProgDelete(M); //?
                }
                return -GetGapPenalty1(A, B, M);
            }

            if (M <= 1) // if profile 1 is small...
            {
                if (M <= 0) // if profile 1 is empty...
                {
                    // Insert residues B[1] to B[N]
                    ProgAdd(N); //?
                    return -GetGapPenalty2(A, B, N);
                }

                if (startGapOpenPenalty == 0)
                {
                    midh = -GetGapPenalty1(A + 1, B + 1, N);
                }
                else
                {
                    midh = -GetGapPenalty2(A + 1, B, 1) - GetGapPenalty1(A + 1, B + 1, N);
                }

                midj = 0;

                for (int j = 1; j <= N; j++)
                {
                    hh = -GetGapPenalty1(A, B + 1, j - 1) + GetProfileScore(A + 1, B + j) - GetGapPenalty1(A + 1, B + j + 1, N - j);
                    if (hh > midh)
                    {
                        midh = hh;
                        midj = j;
                    }
                }

                if (midj == 0)
                {
                    ProgAdd(N);
                    ProgDelete(1);
                }
                else
                {
                    if (midj > 1)
                    {
                        ProgAdd(midj - 1);
                    }
                    ProgAlign();
                    if (midj < N)
                    {
                        ProgAdd(N - midj);
                    }
                }
                return midh;
            }

            midi = M / 2; // Divide profile 1 in half

            /* In a forward phase, calculate all HH[j] and DD[j] */
           
            HH = new double[maxLength + 1];
            DD = new double[maxLength + 1];
            //double t, tl;

            HH[0] = 0;
            t = -GetGapOpeningPenalty1(A, B + 1);
            tl = -GetGapExtensionPenalty1(A, B + 1);

            for (int j = 1; j <= N; j++)
            {
                HH[j] = t = t + tl;
                DD[j] = t - GetGapOpeningPenalty2(A + 1, B + j);
            }

            if (startGapOpenPenalty == 0)
            {
                t = 0;
            }
            else
            {
                t = -GetGapOpeningPenalty2(A + 1, B);
            }
            tl = -GetGapExtensionPenalty2(A + 1, B);

            //double f, s, e, g, h;

            for (int i = 1; i <= midi; i++)
            {
                s = HH[0];
                HH[0] = hh = t = t + tl;
                f = t - GetGapOpeningPenalty1(A + i, B + 1);

                for (int j = 1; j <= N; j++)
                {
                    g = GetGapOpeningPenalty1(A + i, B + j);
                    h = GetGapExtensionPenalty1(A + i, B + j);

                    if ((hh = hh - g - h) > (f = f - h))
                    {
                        f = hh;
                    }

                    g = GetGapOpeningPenalty2(A + i, B + j);
                    h = GetGapExtensionPenalty2(A + i, B + j);

                    if ((hh = HH[j] - g - h) > (e = DD[j] - h))
                    {
                        e = hh;
                    }

                    hh = s + GetProfileScore(A + i, B + j);

                    if (f > hh)
                    {
                        hh = f;
                    }

                    if (e > hh)
                    {
                        hh = e;
                    }

                    s = HH[j];
                    HH[j] = hh;
                    DD[j] = e;
                }
            }

            DD[0] = HH[0];

            // In a reverse phase, calculate all RR[j] and SS[j]

            RR = new double[maxLength + 1];
            SS = new double[maxLength + 1];
            gS = new double[maxLength + 1];

            RR[N] = 0;
            tl = 0;

            for (int j = N - 1; j >= 0; j--)
            {
                g = -GetGapOpeningPenalty1(A + M, B + j + 1);
                tl -= GetGapExtensionPenalty1(A + M, B + j + 1);
                RR[j] = g + tl;
                SS[j] = RR[j] - GetGapOpeningPenalty2(A + M, B + j);
                gS[j] = GetGapOpeningPenalty2(A + M, B + j);
            }

            tl = 0;
            for (int i = M - 1; i >= midi; i--)
            {
                s = RR[N];
                if (endGapOpenPenalty == 0)
                {
                    g = 0;
                }
                else
                {
                    g = -GetGapOpeningPenalty2(A + i + 1, B + N);
                }
                tl -= GetGapExtensionPenalty2(A + i + 1, B + N);
                RR[N] = hh = g + tl;
                t = GetGapOpeningPenalty1(A + i, B + N);
                f = RR[N] - t;

                for (int j = N - 1; j >= 0; j--)
                {
                    g = GetGapOpeningPenalty1(A + i, B + j + 1);
                    h = GetGapExtensionPenalty1(A + i, B + j + 1);

                    if ((hh = hh - g - h) > (f = f - h - g + t))
                    {
                        f = hh;
                    }

                    t = g;
                    g = GetGapOpeningPenalty2(A + i + 1, B + j);
                    h = GetGapExtensionPenalty2(A + i + 1, B + j);
                    hh = RR[j] - g - h;

                    if (i == M - 1)
                    {
                        e = SS[j] - h;
                    }
                    else
                    {
                        e = SS[j] - h - g + GetGapOpeningPenalty2(A + i + 2, B + j);
                        gS[j] = g;
                    }

                    if (hh > e)
                    {
                        e = hh;
                    }
                    hh = s + GetProfileScore(A + i + 1, B + j + 1);

                    if (f > hh)
                    {
                        hh = f;
                    }

                    if (e > hh)
                    {
                        hh = e;
                    }

                    s = RR[j];
                    RR[j] = hh;
                    SS[j] = e;
                }
            }

            SS[N] = RR[N];
            gS[N] = GetGapOpeningPenalty2(A + midi + 1, B + N);

            // Find midj, such that HH[j] + RR[j] or DD[j] + SS[j] + gap is the maximum
            // My note: Okay, what I see here is that HH/RR is tracking the situation without a gap,
            // and DD/SS/gap is tracking the alignment if a gap has been opened.

            midh = HH[0] + RR[0];
            midj = 0;
            gapType = 1; // 1 would presumably be new gaps, 2 old? 

            for (int j = 0; j <= N; j++)
            {
                hh = HH[j] + RR[j];
                if (hh >= midh)
                {
                    if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j]))
                    {
                        midh = hh;
                        midj = j;
                    }
                }
            }

            for (int j = N; j >= 0; j--)
            {
                hh = DD[j] + SS[j] + gS[j];
                if (hh > midh)
                {
                    midh = hh;
                    midj = j;
                    gapType = 2;
                }
            }

            // Conquer recursively around midpoint

            if (gapType == 1)
            {
                ProgDiff(A, B, midi, midj, startGapOpenPenalty, 1);
                ProgDiff(A + midi, B + midj, M - midi, N - midj, 1, endGapOpenPenalty);
            }
            else // gapType == 2
            {
                ProgDiff(A, B, midi - 1, midj, startGapOpenPenalty, 0);
                ProgDelete(2);
                ProgDiff(A + midi + 1, B + midj, M - midi - 1, N - midj, 0, endGapOpenPenalty);
            }

            return midh; // return the score of the best alignment
            
        }

        

        double GetGapPenalty1(int i, int j, int k)
        {
            double gapPenalty;

            //User parameter
            bool endGapPenalties = false;

            if (k <= 0)
            {
                return 0;
            }

            if (!endGapPenalties && (i == 0 || i == group1Profile.ProfileLength))
            {
                return 0;
            }

            gapPenalty = group2Profile.GapOpeningPenalties[j] + group1Profile.GapOpeningPenalties[i];
            for (int l = 0; l < k && l + j < group2Profile.ProfileLength; l++)
            {
                gapPenalty += group2Profile.GapExtensionPenalties[l + j];
            }

            return gapPenalty;
        }

        double GetGapPenalty2(int i, int j, int k)
        {
            double gapPenalty;

            //User parameter
            bool endGapPenalties = false;

            if (k <= 0)
            {
                return 0;
            }

            if (!endGapPenalties && (j == 0 || j == group2Profile.ProfileLength))
            {
                return 0;
            }

            gapPenalty = group1Profile.GapOpeningPenalties[i] + group2Profile.GapOpeningPenalties[j];
            for (int l = 0; l < k && l + i < group1Profile.ProfileLength; l++)
            {
                gapPenalty += group1Profile.GapExtensionPenalties[l + i];
            }

            return gapPenalty;
        }

        // Clustal makes this an inline function. I may want to check
        // if using this .NET 4.5-specific option helps or not.
        [MethodImplAttribute(MethodImplOptions.AggressiveInlining)]
        double GetGapOpeningPenalty1(int i, int j)
        {
            double penalty;

            //User parameter
            bool endGapPenalties = false;

            if (!endGapPenalties && (i == 0 || i == group1Profile.ProfileLength))
            {
                return 0;
            }

            penalty = group2Profile.GapOpeningPenalties[j] + group1Profile.GapOpeningPenalties[i];

            return penalty;
        }

        // Clustal makes this an inline function. I may want to check
        // if using this .NET 4.5-specific option helps or not.
        [MethodImplAttribute(MethodImplOptions.AggressiveInlining)]
        double GetGapOpeningPenalty2(int i, int j)
        {
            double penalty;

            //User parameter
            bool endGapPenalties = false;

            if (!endGapPenalties && (j == 0 || j == group2Profile.ProfileLength))
            {
                return 0;
            }

            penalty = group1Profile.GapOpeningPenalties[i] + group2Profile.GapOpeningPenalties[j];

            return penalty;
        }

        // Clustal makes this an inline function. I may want to check
        // if using this .NET 4.5-specific option helps or not.
        [MethodImplAttribute(MethodImplOptions.AggressiveInlining)]
        double GetGapExtensionPenalty1(int i, int j)
        {
            double penalty;

            //User parameter
            bool endGapPenalties = false;

            if (!endGapPenalties && (i == 0 || i == group1Profile.ProfileLength))
            {
                return 0;
            }

            penalty = group2Profile.GapExtensionPenalties[j];

            return penalty;
        }

        // Clustal makes this an inline function. I may want to check
        // if using this .NET 4.5-specific option helps or not.
        [MethodImplAttribute(MethodImplOptions.AggressiveInlining)]
        double GetGapExtensionPenalty2(int i, int j)
        {
            double penalty;

            //User parameter
            bool endGapPenalties = false;

            if (!endGapPenalties && (j == 0 || j == group2Profile.ProfileLength))
            {
                return 0;
            }

            penalty = group1Profile.GapExtensionPenalties[i];

            return penalty;
        }

        // Clustal makes this an inline function. I may want to check
        // if using this .NET 4.5-specific option helps or not.
        [MethodImplAttribute(MethodImplOptions.AggressiveInlining)]
        double GetProfileScore(int n, int m)
        {
            double score = 0;

            // Since this is called so often, I should store the bound externally
            foreach (char residue in residueCodes) // make sure this is one with gaps
            {
                score += group1Profile.ResidueWeights[residue][n] * group2Profile.ResidueWeights[residue][m];
            }
            
            return score / 10;
        }

        void ProgDelete(int k)
        {
            if (lastPrint < 0)
            {
                lastPrint = displ[printPtr - 1] -= k;
            }
            else
            {
                lastPrint = displ[printPtr++] = -k;
            }
        }

        void ProgAdd(int k)
        {
            if (lastPrint < 0)
            {
                displ[printPtr - 1] = k;
                displ[printPtr++] = k;
            }
            else
            {
                lastPrint = displ[printPtr++] = k;
            }
        }

        void ProgAlign()
        {
            displ[printPtr++] = lastPrint = 0;
        }
    }
}
