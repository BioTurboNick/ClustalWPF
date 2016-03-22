using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;

namespace ClustalWPF.MultipleAlignment
{
    class Profile
    {
        // The profile is how informataion about the sequences are stored for use in alignement.
        // It mainly consists of a two-dimensional matrix providing information about the sequences in
        // the group forming the profile:
        //   - Each column corresponds to a position in the alignment
        //   - Each row corresponds to a possible residue type (ordered alphabetically by single-letter code),
        //     followed by a row for old gaps and a row for new gaps
        // Gap opening and extension penalties are stored in separate arrays. These contain extra entries
        // which store the value for the start and end of the profile.
        // Another array stores the number of sequences with a gap at each position, and if it is a protein,
        // it also indicates the distance of non-gap residues to the gap within the gapDistance parameter.

        ReadOnlyCollection<AlignedMacromolecule> profileMacromolecules;
        int profileLength;
        Dictionary<char, double[]> residueWeights; // temporarily reconfiguring to have the residues be 1-profileLength, with 0 being initial infor
        double[] newGapWeights;
        double[] oldGapWeights;
        double[] gapOpeningPenalties; // temporarily reconfiguring to have the residues be 1-profileLength, with 0 being initial infor
        double[] gapExtensionPenalties; // temporarily reconfiguring to have the residues be 1-profileLength, with 0 being initial infor
        int[] gapPositions;

        //User parameters
        bool useResiduePenalties;
        bool useHydrophilicPenalties;
        bool useVariabilityPenalties;

        // Pascarella and Argos residue specific gap modification factors from
        // http://www.ncbi.nlm.nih.gov/pubmed/7984417
        // But note that these values are transformed by (2 - x) * 100
        static int[] pascarellaFactors = new int[] { 87, 87, 104, 69, 80, 139, 100, 104, 68, 79, 71, 137, 126, 93, 128, 124, 111, 75, 100, 77 };
        static Dictionary<char, int> pascarellaFactorsDict = new Dictionary<char, int>();

        static char[] hydrophilicResidues = new char[] { 'G', 'P', 'S', 'N', 'D', 'Q', 'E', 'K', 'R' }; // Should be obtained from central location

        public int ProfileLength
        {
            get { return profileLength; }
        }

        public Dictionary<char, double[]> ResidueWeights
        {
            get { return residueWeights; }
        }

        public double[] GapOpeningPenalties
        {
            get { return gapOpeningPenalties; }
        }

        public double[] GapExtensionPenalties
        {
            get { return gapExtensionPenalties; }
        }

        public ReadOnlyCollection<AlignedMacromolecule> Macromolecules
        {
            get { return profileMacromolecules; }
        }

        Profile(int length, char[] residues, ReadOnlyCollection<AlignedMacromolecule> macromolecules) // need to figure out how to better specify these
        {
            profileLength = length;
            profileMacromolecules = macromolecules;
            residueWeights = NewProfileDict(residues, length + 1);
            newGapWeights = new double[length];
            oldGapWeights = new double[length];
            gapOpeningPenalties = new double[length + 1];
            gapExtensionPenalties = new double[length + 1];
            gapPositions = new int[length];
        }

        public Profile()
        { }

        static Profile()
        {
            for (int i = 0; i < Routines.nondegenerateAminoAcidCodes.Length; i++)
            {
                pascarellaFactorsDict.Add(Routines.nondegenerateAminoAcidCodes[i], pascarellaFactors[i]);
            }
        }

        public static Profile Calculate(ReadOnlyCollection<AlignedMacromolecule> macromolecules, int length, double gapOpeningCoefficient, double gapExtensionCoefficient, SubstitutionMatrix.SubstitutionMatrix subMatrix)
        {
            // For the moment, ignoring struct penalties and gap penalty masks

            char[] residues;

            // Configure parameters
            bool isNucleicAcid = macromolecules[0].IsNucleicAcid; // More global parameter to access?
            if (isNucleicAcid) 
            {
                residues = Routines.nucleotideCodesPlusGaps;
            }
            else
            {
                residues = Routines.aminoAcidCodesPlusGaps;
            }

            Profile profile = new Profile(length, residues, macromolecules);
            
            profile.CountGaps();

            profile.CalculateGapPenalties(gapOpeningCoefficient, gapExtensionCoefficient);

            if (subMatrix == null) // The substitution matrix is empty
            {
                profile.CalculateResidueWeights(residues);
            }
            else
            {
                profile.CalculateResidueWeights(residues, ref subMatrix);
            }
            
            return profile;
        }

        void CalculateGapPenalties(double gapOpeningCoefficient, double gapExtensionCoefficient)
        {
            bool isNucleicAcid = profileMacromolecules[0].IsNucleicAcid;

            // User parameters
            int gapDistance = 4;
            bool useStructurePenalties = false;
            useVariabilityPenalties = true;
            useResiduePenalties = true;
            useHydrophilicPenalties = true;
            bool endGapPenalties = true;

            // user parameter?
            double reducedGap = 1.0;
            
            double[] variabilityPenalties = new double[0];
            double[] residuePenalties = new double[0];
            double[] hydrophilicPenalties = new double[0];

            if (!isNucleicAcid)
            {
                if (useStructurePenalties)
                {
                    useVariabilityPenalties = useHydrophilicPenalties = useResiduePenalties = false;
                }
                else if (useVariabilityPenalties && profileMacromolecules.Count == 2)
                {
                    if (CalculatePercentIdentity(profileMacromolecules[0], profileMacromolecules[1]) > 0.6)
                    {
                        useHydrophilicPenalties = useResiduePenalties = false;
                    }
                }
                else
                {
                    useVariabilityPenalties = false;
                }

                if (useVariabilityPenalties)
                {
                    variabilityPenalties = new double[profileLength];
                    variabilityPenalties = CalculateVariabilityPenalties();
                }

                if (useResiduePenalties)
                {
                    residuePenalties = new double[profileLength];
                    residuePenalties = CalculateResiduePenalties();
                }

                if (useHydrophilicPenalties)
                {
                    hydrophilicPenalties = new double[profileLength];
                    hydrophilicPenalties = CalculateHydrophilicPenalties();
                }
            }

            for (int i = 0; i < profileLength; i++)
            {
                if (gapPositions[i] <= 0) // For all non-gap positions
                {
                    if (isNucleicAcid)
                    {
                        gapOpeningPenalties[i + 1] = gapOpeningCoefficient;
                    }
                    else
                    {
                        gapOpeningPenalties[i + 1] = CalculateLocalPenalty(gapOpeningCoefficient, i, ref residuePenalties, ref hydrophilicPenalties, ref variabilityPenalties);
                    }
                    gapExtensionPenalties[i + 1] = gapExtensionCoefficient;

                    // Increase the gap penalties near existing gaps

                    if (gapPositions[i] < 0 && !useStructurePenalties)
                    {
                        gapOpeningPenalties[i + 1] = gapOpeningPenalties[i + 1] * (2.0 + 2.0 * (double)(gapDistance + gapPositions[i]) / (double)gapDistance);
                    }
                }
                else
                {
                    double scale = ((double)(profileMacromolecules.Count - gapPositions[i]) / (double)(profileMacromolecules.Count)) * reducedGap;
                    gapOpeningPenalties[i + 1] = scale * gapOpeningCoefficient;
                    gapExtensionPenalties[i + 1] = 0.5 * gapExtensionCoefficient;
                }

                // Apply the gap penalty mask

                if (useStructurePenalties)
                {
                    throw new NotImplementedException();
                    // I'm not sure what to do with these yet.
                }

                // Make sure no penalty is zero, even for all-gap positions

                if (gapOpeningPenalties[i + 1] < 1)
                {
                    gapOpeningPenalties[i + 1] = 1;
                }
                if (gapExtensionPenalties[i + 1] < 1)
                {
                    gapExtensionPenalties[i + 1] = 1;
                }
            }

            // Set the penalties at the beginning and end of the profile
            if (endGapPenalties) // Why not also check useEndGap?
            {
                gapOpeningPenalties[0] = gapOpeningCoefficient;
                gapExtensionPenalties[0] = gapExtensionCoefficient;
            }
            else
            {
                gapOpeningPenalties[0] = gapOpeningPenalties[profileLength] = 0;
                gapExtensionPenalties[0] = gapExtensionPenalties[profileLength] = 0;
            }
        }

        void CountGaps()
        {
            // Will be user parameters:
            bool useEndGaps = true;
            bool endGapPenalties = true;
            int gapDistance = 4;

            bool isNucleicAcid = profileMacromolecules[0].IsNucleicAcid;

            foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
            {
                // count the gaps within the profile (Clustal looks for base < 0 or > maxAA)
                for (int i = 1; i < macromolecule.Sequence.Length; i++)
                {
                    if (macromolecule.AlignedPositions[i] - macromolecule.AlignedPositions[i - 1] > 1)
                    {
                        for (int j = macromolecule.AlignedPositions[i - 1] + 1; j < macromolecule.AlignedPositions[i]; j++)
                        {
                            gapPositions[j]++;
                        }
                    }
                }

                if (useEndGaps || endGapPenalties)
                {
                    // add in end gaps
                    for (int i = 0; i < macromolecule.AlignedPositions[0]; i++)
                    {
                        gapPositions[i]++;
                    }
                    for (int i = profileLength - 1; i > macromolecule.AlignedPositions[macromolecule.Sequence.Length - 1]; i--)
                    {
                        gapPositions[i]++;
                    }
                }
            }

            if (!isNucleicAcid && gapDistance > 0)
            {
                // As near as I can tell, this ends up marking those non-gap residues within gapDistance - 1 of the start of a gap
                // So the gap sequence        0  1  1  1  0  0  0  0  0  0  0  2  2  2  2  0  0  1  1  1  0  0  with gapDistance = 4 would give
                // would give gapPositions = -1  1  1  1 -1 -2 -3  0 -3 -2 -1  2  2  2  2 -1 -1  1  1  1 -1 -2

                for (int i = 0; i < profileLength; i++)
                {
                    if (gapPositions[i] > 0)
                    {
                        for (int j = 1; j < gapDistance; j++) // Mark positions prior to the gap
                        {
                            int position = i - j;
                            if (position >= 0 && gapPositions[position] == 0) // We only want to mark non-gap positions
                            {
                                if (gapPositions[position] >= 0)
                                {
                                    gapPositions[position] = -j;
                                }
                                else
                                {
                                    gapPositions[position] = Math.Max(gapPositions[position], -j); // use the smallest absolute value
                                }
                            }
                        }

                        while (i < profileLength && gapPositions[i] > 0) // Scan until the end of the gap, copying gap information over
                        {
                            i++;
                        }

                        int iAdjust = 0; // Keep track of how far i can be advanced after this loop
                        for (int j = 1; j < gapDistance; j++) // Mark positions following the gap, same strategy as above
                        {
                            int position = i + j;
                            if (position < profileLength)
                            {
                                if (gapPositions[position] == 0)
                                {
                                    if (gapPositions[position] >= 0)
                                    {
                                        gapPositions[position] = -j;
                                    }
                                    else
                                    {
                                        gapPositions[position] = Math.Max(gapPositions[position], -j);
                                    }
                                    iAdjust++;
                                }
                                else // found another gap, can stop prematurely
                                {
                                    break;
                                }
                            }
                        }
                        i += iAdjust;
                    }
                }
            }
        }

        double[] CalculateVariabilityPenalties() // Only used if there are only two sequences in the profile.
        {
            int windowExtent = 5; // window size is +/- 5 at each position.
            int variabilityScoreLowerLimit = 50;

            double[] penalties = new double[profileLength];
            double[] maxPossiblePenalties = new double[profileLength];

            for (int i = 0, j = 0; i < profileMacromolecules[0].Sequence.Length && j < profileMacromolecules[1].Sequence.Length; i++, j++)
            {
                bool endFound = false;

                if (profileMacromolecules[0].AlignedPositions[i] - profileMacromolecules[1].AlignedPositions[j] != 0)
                {
                    endFound = FindNextAlignedPositions(ref i, ref j, profileMacromolecules[0].AlignedPositions, profileMacromolecules[1].AlignedPositions);
                }

                if (!endFound)
                {
                    if (profileMacromolecules[0].Sequence[i] == profileMacromolecules[1].Sequence[j])
                    {
                        for (int k = Math.Max(0, profileMacromolecules[0].AlignedPositions[i] - windowExtent); k < Math.Min(profileLength, profileMacromolecules[0].AlignedPositions[i] + windowExtent); k++)
                        {
                            // An amino acid match counts for every position within the window.
                            penalties[k]++;
                            maxPossiblePenalties[k]++;
                        }
                    }
                    else
                    {
                        for (int k = Math.Max(0, profileMacromolecules[0].AlignedPositions[i] - windowExtent); k < Math.Min(profileLength, profileMacromolecules[0].AlignedPositions[i] + windowExtent); k++)
                        {
                            maxPossiblePenalties[k]++;
                        }
                    }
                }
            }

            // Clustal says that at this point, penalties are between -maxPossibleScore and +maxPossibleScore
            // However, I don't know where the negatives would come in.

            for (int i = 0; i < profileLength; i++)
            {
                penalties[i] += maxPossiblePenalties[i];
                if (maxPossiblePenalties[i] > 0)
                {
                    penalties[i] = (penalties[i] * 100) / (2 * maxPossiblePenalties[i]); // Note: integer division! Change?
                }
                else
                {
                    penalties[i] = 100;
                }

                if (penalties[i] < variabilityScoreLowerLimit)
                {
                    penalties[i] = variabilityScoreLowerLimit;
                }
            }

            // Clustal says that now the scores are between the lower limit and 100.

            return penalties;
        }

        double[] CalculateResiduePenalties()
        {
            double[] penalties = new double[profileLength];
            foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
            {
                for (int i = 0; i < macromolecule.Sequence.Length; i++)
                {
                    penalties[macromolecule.AlignedPositions[i]] += (180 - pascarellaFactorsDict[macromolecule.Sequence[i]]);
                }
            }
            return penalties;
        }

        double[] CalculateHydrophilicPenalties()
        {
            bool[] isHydrophilic = new bool[profileLength];
            double[] penalties = new double[profileLength];
            foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
            {
                for (int i = 0; i < macromolecule.Sequence.Length; i++)
                {
                    if (hydrophilicResidues.Contains(macromolecule.Sequence[i]))
                    {
                        isHydrophilic[macromolecule.AlignedPositions[i]] = true;
                    }
                }

                for (int i = 0; i < profileLength; i++)
                {
                    if (isHydrophilic[i])
                    {
                        int startStretch = i;
                        while (i < profileLength && isHydrophilic[i])
                        {
                            i++;
                        }
                        int afterStretch = i;

                        if (afterStretch - startStretch > 3)
                        {
                            for (int j = startStretch; j < afterStretch; j++)
                            {
                                penalties[j] += 100;
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < profileLength; i++)
            {
                penalties[i] = penalties[i] / (double)(profileMacromolecules.Count);
            }

            return penalties;
        }

        public double CalculateLocalPenalty(double penalty, int n, ref double[] residuePenalties, ref double[] hydrophilicPenalties, ref double[] varibilityPenalties)
        {
            bool usedHydrophilic = false;
            double localPenalty;

            localPenalty = 1.0;
            if (useVariabilityPenalties)
            {
                localPenalty *= varibilityPenalties[n] / 100.0;
            }

            if (useHydrophilicPenalties && (usedHydrophilic = hydrophilicPenalties[n] > 0))
            {
                localPenalty *= 0.5;
            }

            if (useResiduePenalties && !usedHydrophilic)
            {
                localPenalty *= residuePenalties[n] / 100.0;
            }

            localPenalty *= penalty;
            return localPenalty;
        }

        static double CalculatePercentIdentity(AlignedMacromolecule macromolecule1, AlignedMacromolecule macromolecule2)
        {
            double percentIdentity = 0;
            int matched = 0;

            for (int i = 0, j = 0; i < macromolecule1.Sequence.Length && j < macromolecule2.Sequence.Length; i++, j++)
            {
                bool endFound = false;

                if (macromolecule1.AlignedPositions[i] - macromolecule2.AlignedPositions[j] != 0)
                {
                    endFound = FindNextAlignedPositions(ref i, ref j, macromolecule1.AlignedPositions, macromolecule2.AlignedPositions);
                }

                if (!endFound && macromolecule1.Sequence[i] == macromolecule2.Sequence[j])
                {
                    matched++;
                }
            }

            percentIdentity = 100 * (double)matched / (double)macromolecule1.Sequence.Length;

            return percentIdentity;
        }

        static bool FindNextAlignedPositions(ref int i, ref int j, int[] alignedPositions1, int[] alignedPositions2)
            // returns true if it reaches the end of the sequences before finding another aligned position
        {
            bool endFound = false;
            int difference = 0;

            do
            {
                while (!endFound && (difference = alignedPositions1[i] - alignedPositions2[j]) < 0)
                {
                    // There were gaps in sequence 1 but not in sequence 2. Move sequence 2 position up until they align or pass.
                    j++;
                    endFound = j >= alignedPositions2.Length;
                }

                while (!endFound && (difference = alignedPositions1[i] - alignedPositions2[j]) > 0)
                {
                    // There were gaps in sequence 2 but not in sequence 1. Move sequence 1 position up until they align or pass.
                    i++;
                    endFound = i >= alignedPositions1.Length;
                }
            } while (difference != 0 && !endFound);

            return endFound;
        }

        protected void CalculateResidueWeights(char[] residues)
        {
            // Calculate the maximum possible score value
            double maxScore = 0;
            foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
            {
                maxScore += macromolecule.SequenceWeight;
            }

            if (maxScore > 0)
            {
                Dictionary<char, double[]> residueSums = NewProfileDict(residues, profileLength + 1);

                int previousPos = 0;

                foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
                {
                    for (int i = 0; i < macromolecule.Sequence.Length; i++)
                    {
                        int currentPos = macromolecule.AlignedPositions[i];

                        residueSums[macromolecule.Sequence[i]][currentPos] += macromolecule.SequenceWeight;
                        
                        if (currentPos > previousPos + 1)
                        {
                            // There was a gap
                            for (int j = previousPos + 1; j < currentPos; j++)
                            {
                                if (macromolecule.NewGaps.Contains(j))
                                {
                                    residueSums[Routines.newGap][j] += macromolecule.SequenceWeight;
                                }
                                else
                                {
                                    residueSums[Routines.oldGap][j] += macromolecule.SequenceWeight;
                                }
                            }
                        }
                    }
                }

                foreach (char residue in residues)
                {
                    double[] singleResidueWeights = residueWeights[residue];
                    double[] singleResidueSums = residueSums[residue];
                    for (int i = 0; i < profileLength; i++)
                    {
                        singleResidueWeights[i + 1] = 10 * residueSums[residue][i] / maxScore; // This should modify the base dictionary array, but check
                    }
                }
            }
        }

        protected void CalculateResidueWeights(char[] residues, ref SubstitutionMatrix.SubstitutionMatrix subMatrix)
        {
            double maxScore;
            Dictionary<char, double[]> weighting;

            CalculateResidueWeightsCommon(residues, out maxScore, out weighting);

            // Clustal does some things by looking up the gaps in the substitution matrix, but it appears
            // that it will always be equal to 0.
            for (int i = 0; i < profileLength; i++)
            {
                if (gapPositions[i] < profileMacromolecules.Count) // If at least one sequence has a residue here.
                {
                    double scale = (double)(profileMacromolecules.Count - gapPositions[i]) / (double)(profileMacromolecules.Count);

                    foreach (char residue1 in residues)
                    {
                        double score = 0;
                        foreach (char residue2 in residues)
                        {
                            score += weighting[residue2][i] * subMatrix[residue2][residue1];
                        }
                        residueWeights[residue1][i + 1] = (score / maxScore) * scale;
                    }

                    double oldGapScore = 0; // though if these really are always 0, I should just get rid of them or they're pointless operations.
                    double newGapScore = 0;
                    foreach (char residue in residues)
                    {
                        oldGapScore += weighting[residue][i] * subMatrix[residue][Routines.oldGap];
                        newGapScore += weighting[residue][i] * subMatrix[residue][Routines.newGap];
                    }
                    residueWeights[Routines.oldGap][i + 1] = (oldGapScore / maxScore) * scale;
                    residueWeights[Routines.newGap][i + 1] = (newGapScore / maxScore) * scale;
                }
            }
        }

        static Dictionary<char, double[]> NewProfileDict(char[] residues, int profileLength)
        {
            Dictionary<char, double[]> newDict = new Dictionary<char, double[]>(residues.Length);
            foreach (char residue in residues)
            {
                newDict.Add(residue, new double[profileLength]);
            }

            return newDict;
        }

        void CalculateResidueWeightsCommon(char[] residues, out double maxScore, out Dictionary<char, double[]> weights)
        {
            // Calculate the maximum possible score value
            maxScore = 0;
            foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
            {
                maxScore += macromolecule.SequenceWeight;
            }

            weights = NewProfileDict(residues, profileLength);

            if (maxScore > 0)
            {
                foreach (AlignedMacromolecule macromolecule in profileMacromolecules)
                {
                    int previousPos = 0;
                    for (int i = 0; i < macromolecule.Sequence.Length; i++)
                    {
                        int currentPos = macromolecule.AlignedPositions[i];
                        weights[macromolecule.Sequence[i]][currentPos] += macromolecule.SequenceWeight;

                        if (currentPos > previousPos + 1)
                        {
                            // There was a gap
                            for (int j = previousPos + 1; j < currentPos; j++)
                            {
                                if (macromolecule.NewGaps.Contains(j))
                                {
                                    weights[Routines.newGap][j] += macromolecule.SequenceWeight;
                                }
                                else
                                {
                                    weights[Routines.oldGap][j] += macromolecule.SequenceWeight;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
