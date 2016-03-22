using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class PAM : SubstitutionMatrix
    {
        static Dictionary<string, SubstitutionMatrix> matrices = new Dictionary<string, SubstitutionMatrix>();

        public static SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            // intScale = 100
            if (useNegative) // Clustal: || !getDistanceTree
            {
                // scale 0.75
                return matrices["PAM120"];
            }
            else if (percentIdentity > 80)
            {
                // scale 0.75
                return matrices["PAM20"];
            }
            else if (percentIdentity > 60)
            {
                // scale 0.75
                return matrices["PAM60"];
            }
            else if (percentIdentity > 40)
            {
                // scale 0.75
                return matrices["PAM120"];
            }
            else
            {
                // scale 0.75
                return matrices["PAM350"];
            }
        }

        public double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            return 0.75;
        }

        static PAM()
        {
            matrices.Add("PAM20", SubstitutionMatrix.Create("PAM20", MatrixData.PAM20, Routines.aminoAcidCodes));
            matrices.Add("PAM60", SubstitutionMatrix.Create("PAM60", MatrixData.PAM60, Routines.aminoAcidCodes));
            matrices.Add("PAM120", SubstitutionMatrix.Create("PAM120", MatrixData.PAM120, Routines.aminoAcidCodes));
            matrices.Add("PAM160", SubstitutionMatrix.Create("PAM160", MatrixData.PAM160, Routines.aminoAcidCodes));
            matrices.Add("PAM250", SubstitutionMatrix.Create("PAM250", MatrixData.PAM250, Routines.aminoAcidCodes));
            matrices.Add("PAM350", SubstitutionMatrix.Create("PAM350", MatrixData.PAM350, Routines.aminoAcidCodes));
        }
    }
}
