using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    public class SubstitutionMatrix
        // Defines a substitution matrix that allows value lookup by residue character
        // Clustal expands the matrix to include both types of gaps as columns too. However, they set all of these to '0'.
        // So there's really no need to program them into the matrix.
    {
        Dictionary<char, Dictionary<char, double>> matrix;
        string name;
        double averageScore;
        char[] residueCodes;
        
        const int gapResidueValue = 0;
        const int gapGapValue = 0;

        public Dictionary<char, double> this[char key]
        {
            get
            {
                return matrix[key];
            }
        }

        public string Name
        {
            get { return name; }
        }

        public double AverageScore
        {
            get { return averageScore; }
        }

        public SubstitutionMatrix()
        { }

        SubstitutionMatrix(string newName, Dictionary<char, Dictionary<char, double>> newMatrix, char[] newResidueCodes)
        {
            name = newName;
            matrix = newMatrix;
            residueCodes = newResidueCodes;
        }

        //public static SubstitutionMatrix Create(string name, List<List<double>> matrixInfo, char[] residueCodes)
        //// matrixInfo is a list of lists, with each sub-list containing a line from the substituion matrix
        //// file, which is stored as the bottom-left triangle of the symmetrical matrix.
        //{
        //    int maxResidues = residueCodes.Length;
        //    bool useNegativeMatrix = false;
            
        //    Dictionary<char, Dictionary<char, double>> matrix = new Dictionary<char, Dictionary<char, double>>(residueCodes.Length + 2);
        //    Dictionary<char, double> oldGapRow = new Dictionary<char, double>(residueCodes.Length + 2);
        //    Dictionary<char, double> newGapRow = new Dictionary<char, double>(residueCodes.Length + 2);

        //    for (int i = 0; i < maxResidues; i++)
        //    {
        //        int limit = matrixInfo[i].Count;
        //        Dictionary<char, double> row = new Dictionary<char, double>(residueCodes.Length + 2);
        //        for (int j = 0; j < limit; j++)
        //        {
        //            row.Add(residueCodes[j], matrixInfo[i][j] * 100);
        //        }
        //        for (int j = limit; j < maxResidues; j++)
        //        {
        //            row.Add(residueCodes[j], matrixInfo[j][i] * 100);
        //        }
        //        row.Add(Routines.oldGap, gapResidueValue);
        //        row.Add(Routines.newGap, gapResidueValue);

        //        matrix.Add(residueCodes[i], row);

        //        oldGapRow.Add(residueCodes[i], gapResidueValue);
        //        newGapRow.Add(residueCodes[i], gapResidueValue);
        //    }

        //    oldGapRow.Add(Routines.oldGap, gapGapValue);
        //    oldGapRow.Add(Routines.newGap, gapGapValue);
        //    newGapRow.Add(Routines.oldGap, gapGapValue);
        //    newGapRow.Add(Routines.newGap, gapGapValue);

        //    matrix.Add(Routines.oldGap, oldGapRow);
        //    matrix.Add(Routines.newGap, newGapRow);
            
        //    SubstitutionMatrix subMatrix = new SubstitutionMatrix(name, matrix, residueCodes);

        //    subMatrix.CalculateMatrixAverageScore();

        //    if (!useNegativeMatrix)
        //    {
        //        subMatrix.MakePositive();
        //    }

            




        //    return subMatrix;
        //}

        public static SubstitutionMatrix Create(string name, double[] matrixInfo, char[] residueCodes)
        // matrixInfo is a single array containing the bottom-left triangle of the symmetrical matrix.
        {
            int maxResidues = residueCodes.Length;
            bool useNegativeMatrix = false;

            Dictionary<char, Dictionary<char, double>> matrix = new Dictionary<char, Dictionary<char, double>>(residueCodes.Length + 2);
            Dictionary<char, double> oldGapRow = new Dictionary<char, double>(residueCodes.Length + 2);
            Dictionary<char, double> newGapRow = new Dictionary<char, double>(residueCodes.Length + 2);

            int rowStart = 0;
            for (int i = 0; i < maxResidues; i++) // rows
            {
                Dictionary<char, double> row = new Dictionary<char, double>(residueCodes.Length + 2);
                for (int j = 0; j <= i; j++) // Populate the row
                {
                    row.Add(residueCodes[j], matrixInfo[rowStart + j] * 100);
                }
                matrix.Add(residueCodes[i], row);

                for (int j = 0; j < i; j++) // Populate the column too
                {
                    matrix[residueCodes[j]].Add(residueCodes[i], matrixInfo[rowStart + i]);
                }

                rowStart += i + 1;

                oldGapRow.Add(residueCodes[i], gapResidueValue);
                newGapRow.Add(residueCodes[i], gapResidueValue);
            }

            foreach (KeyValuePair<char, Dictionary<char, double>> row in matrix) // Add gap columns
            {
                row.Value.Add(Routines.oldGap, gapResidueValue);
                row.Value.Add(Routines.newGap, gapResidueValue);
            }

            oldGapRow.Add(Routines.oldGap, gapGapValue);
            oldGapRow.Add(Routines.newGap, gapGapValue);
            newGapRow.Add(Routines.oldGap, gapGapValue);
            newGapRow.Add(Routines.newGap, gapGapValue);

            matrix.Add(Routines.oldGap, oldGapRow);
            matrix.Add(Routines.newGap, newGapRow);

            SubstitutionMatrix subMatrix = new SubstitutionMatrix(name, matrix, residueCodes);

            subMatrix.CalculateMatrixAverageScore();

            if (!useNegativeMatrix)
            {
                subMatrix.MakePositive();
            }

            return subMatrix;
        }

        void CalculateMatrixAverageScore()
        {
            int totalResidues = residueCodes.Length;

            double sumDifferent = 0;

            for (int i = 0; i < totalResidues; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    sumDifferent += matrix[residueCodes[i]][residueCodes[j]];
                }
            }

            averageScore = -sumDifferent / ((totalResidues * totalResidues - totalResidues) / 2);
        }

        void MakePositive()
        {
            // Make the matrix positive by adding -(minValue) to every position (excluding gap columns/rows)
            double minValue = GetMatrixMinimumValue();

            for (int i = 0; i < residueCodes.Length; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    matrix[residueCodes[i]][residueCodes[j]] -= minValue;
                }
            }
        }

        double GetMatrixMinimumValue()
        {
            double minValue = 0;

            for (int i = 0; i < residueCodes.Length; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    minValue = Math.Min(minValue, matrix[residueCodes[i]][residueCodes[j]]);
                }
            }

            return minValue;
        }

    }
}
