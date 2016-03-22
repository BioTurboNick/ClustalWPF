using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class IUB : SubstitutionMatrixSeries
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

        static IUB() // Also referred to as SW_GAP in Clustal code
        {
            matrix = SubstitutionMatrix.Create("IUB", MatrixData.IUB, Routines.nucleotideCodes);
        }
    }
}
