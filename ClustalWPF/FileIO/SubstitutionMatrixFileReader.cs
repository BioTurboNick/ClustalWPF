using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using ClustalWPF.SubstitutionMatrix;

namespace ClustalWPF.FileIO
{
    class SubstitutionMatrixFileReader
    // This class contains methods pertaining to reading files.

    // I'm not sure I like the structure of this.
    {



        //public static List<SubstitutionMatrixSeries> LoadSubstituionMatrices(bool isNucleicAcid)
        //{
        //    int line = 0;
        //    string fileName = "";
        //    List<SubstitutionMatrixSeries> matricies = new List<SubstitutionMatrixSeries>();
            
        //    char[] residueCodes;

        //    if (isNucleicAcid)
        //    {
        //        residueCodes = Routines.nucleotideCodes;
        //        fileName = System.Environment.CurrentDirectory + Path.DirectorySeparatorChar + "Nucleotide Substitution Matrices";
        //    }
        //    else
        //    {
        //        residueCodes = Routines.aminoAcidCodes;
        //        fileName = System.Environment.CurrentDirectory + Path.DirectorySeparatorChar + "Amino Acid Substitution Matrices";
        //    }

        //    List<string> fileLines = ReadFile(fileName);

        //    if (fileLines.Count == 0)
        //    {
        //        throw new Exception(); // Not best way to handle not finding the file
        //    }

        //    while (line < fileLines.Count)
        //    {
        //        string currentLine = fileLines[line];

        //        if (currentLine.StartsWith("[")) // Found the start of a new matrix
        //        {
        //            string name = currentLine.Substring(1, currentLine.IndexOf(']') - 1);

        //            List<List<double>> matrixInfo = new List<List<double>>(residueCodes.Length);

        //            // Read out the matrix
        //            line++;
                    

        //            while (line < fileLines.Count)
        //            {
        //                currentLine = fileLines[line];
        //                if (currentLine.Trim().Length == 0 || currentLine.StartsWith("["))
        //                {
        //                    break;
        //                }

        //                string[] splitLine = currentLine.Split(new char[] { ',', ' ' }, StringSplitOptions.RemoveEmptyEntries);
        //                List<double> matrixLine = new List<double>(splitLine.Length);

        //                foreach (string value in splitLine)
        //                {
        //                    matrixLine.Add(double.Parse(value));
        //                }

        //                matrixInfo.Add(matrixLine);

        //                line++;
        //            }

        //            matrices.Add(SubstitutionMatrix.SubstitutionMatrix.Create(name, matrixInfo, residueCodes));
        //        }

        //        line++;
        //    }

        //    return matrices;
        //}

    }

}
