using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class MD : SubstitutionMatrixSeries
    {
        static Dictionary<string, SubstitutionMatrix> matrices = new Dictionary<string, SubstitutionMatrix>();

        public override SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            throw new NotImplementedException();
            // Clustal doesn't use this series, so there's no logic available.
        }

        public override double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            throw new NotImplementedException();
        }

        static MD()
        {
            matrices.Add("MD40", SubstitutionMatrix.Create("MD40", MatrixData.MD40, Routines.aminoAcidCodes));
            matrices.Add("MD120", SubstitutionMatrix.Create("MD120", MatrixData.MD120, Routines.aminoAcidCodes));
            matrices.Add("MD250", SubstitutionMatrix.Create("MD250", MatrixData.MD250, Routines.aminoAcidCodes));
            matrices.Add("MD350", SubstitutionMatrix.Create("MD350", MatrixData.MD350, Routines.aminoAcidCodes));
        }
    }
}
