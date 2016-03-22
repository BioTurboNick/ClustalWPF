using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class Gonnet : SubstitutionMatrixSeries
    {
        static Dictionary<string, SubstitutionMatrix> matrices = new Dictionary<string, SubstitutionMatrix>();

        public override SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            if (useNegative) // Clustal: || !getDistanceTree
            {
                return matrices["Gonnet250"];
            }
            else if (percentIdentity > 35)
            {
                return matrices["Gonnet80"];
            }
            else if (percentIdentity > 25)
            {
                if (minLength < 100)
                {
                    return matrices["Gonnet250"];
                }
                else
                {
                    return matrices["Gonnet120"];
                }
            }
            else
            {
                if (minLength < 100)
                {
                    return matrices["Gonnet350"];
                }
                else
                {
                    return matrices["Gonnet160"];
                }
            }
        }

        public override double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            if (!useNegative && percentIdentity > 35) // (!useNegative && getDistanceTree)
            {
                return 0.25;
            }
            else
            {
                return 0.5;
            }
        }

        static Gonnet()
        {
            matrices.Add("Gonnet40", SubstitutionMatrix.Create("Gonnet40", MatrixData.Gonnet40, Routines.aminoAcidCodes));
            matrices.Add("Gonnet80", SubstitutionMatrix.Create("Gonnet80", MatrixData.Gonnet80, Routines.aminoAcidCodes));
            matrices.Add("Gonnet120", SubstitutionMatrix.Create("Gonnet120", MatrixData.Gonnet120, Routines.aminoAcidCodes));
            matrices.Add("Gonnet160", SubstitutionMatrix.Create("Gonnet160", MatrixData.Gonnet160, Routines.aminoAcidCodes));
            matrices.Add("Gonnet250", SubstitutionMatrix.Create("Gonnet250", MatrixData.Gonnet250, Routines.aminoAcidCodes));
            matrices.Add("Gonnet300", SubstitutionMatrix.Create("Gonnet300", MatrixData.Gonnet300, Routines.aminoAcidCodes));
            matrices.Add("Gonnet350", SubstitutionMatrix.Create("Gonnet350", MatrixData.Gonnet350, Routines.aminoAcidCodes));
        }
    }
}
