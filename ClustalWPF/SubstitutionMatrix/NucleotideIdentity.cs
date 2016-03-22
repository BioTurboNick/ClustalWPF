using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class NucleotideIdentity : SubstitutionMatrixSeries
    {
        static SubstitutionMatrix matrix;

        public override SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            return matrix;
        }

        public override double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            return 0.66;
        }

        static NucleotideIdentity()  // Also referred to as ClustalVDNA in Clustal code, or ClustalW (1.6)
        {
            matrix = SubstitutionMatrix.Create("Identity", MatrixData.NucleotideIdentity, Routines.nucleotideCodes);
        }
    }
}
