using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class BLOSUM : SubstitutionMatrixSeries
    {
        static Dictionary<string, SubstitutionMatrix> matrices = new Dictionary<string, SubstitutionMatrix>();

        public override SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            // intScale = 100 -- This is used to produce values on the same scale
            // as each other.
            if (useNegative) // Clustal: || !getDistanceTree
            {
                // scale 0.75 -- these values are passed back to the profile alignment
                // Seems like this should be handled there instead of here
                return matrices["BLOSUM40"];
            }
            else if (percentIdentity > 80) 
            {
                // scale 0.75
                return matrices["BLOSUM80"];
            }
            else if (percentIdentity > 60)
            {
                // scale 0.75
                return matrices["BLOSUM62x2"];
            }
            else if (percentIdentity > 40)
            {
                // scale 0.75
                return matrices["BLOSUM45"];
            }
            else if (percentIdentity > 30)
            {
                // scale 0.5
                return matrices["BLOSUM45"];
            }
            else if (percentIdentity > 20)
            {
                // scale 0.6
                return matrices["BLOSUM45"];
            }
            else
            {
                // scale 0.6
                return matrices["BLOSUM30"];
            }
        }

        public override double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            if (useNegative || percentIdentity > 40) // || !getDistanceTree
            {
                return 0.75;
            }
            else if (percentIdentity > 30)
            {
                return 0.5;
            }
            else
            {
                return 0.6;
            }
        }

        static BLOSUM()
        {
            matrices.Add("BLOSUM30", SubstitutionMatrix.Create("BLOSUM30", MatrixData.BLOSUM30, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM35", SubstitutionMatrix.Create("BLOSUM35", MatrixData.BLOSUM35, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM40", SubstitutionMatrix.Create("BLOSUM40", MatrixData.BLOSUM40, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM45", SubstitutionMatrix.Create("BLOSUM45", MatrixData.BLOSUM45, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM50", SubstitutionMatrix.Create("BLOSUM50", MatrixData.BLOSUM50, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM55", SubstitutionMatrix.Create("BLOSUM55", MatrixData.BLOSUM55, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM62", SubstitutionMatrix.Create("BLOSUM62", MatrixData.BLOSUM62, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM62x2", SubstitutionMatrix.Create("BLOSUM62x2", MatrixData.BLOSUM62x2, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM65", SubstitutionMatrix.Create("BLOSUM65", MatrixData.BLOSUM65, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM70", SubstitutionMatrix.Create("BLOSUM70", MatrixData.BLOSUM70, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM75", SubstitutionMatrix.Create("BLOSUM75", MatrixData.BLOSUM75, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM80", SubstitutionMatrix.Create("BLOSUM80", MatrixData.BLOSUM80, Routines.aminoAcidCodes));
            matrices.Add("BLOSUM90", SubstitutionMatrix.Create("BLOSUM90", MatrixData.BLOSUM90, Routines.aminoAcidCodes));
        }
    }
}
