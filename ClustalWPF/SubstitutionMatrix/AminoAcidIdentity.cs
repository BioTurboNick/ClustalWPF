using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class AminoAcidIdentity : SubstitutionMatrixSeries
    {
        static SubstitutionMatrix matrix;

        public override SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            return matrix;
        }

        public override double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            return 1;
        }

        static AminoAcidIdentity()
        {
            matrix = SubstitutionMatrix.Create("Identity", MatrixData.AminoAcidIdentity, Routines.aminoAcidCodes);
        }
    }
}
