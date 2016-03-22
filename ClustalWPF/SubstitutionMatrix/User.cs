using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    class User : SubstitutionMatrixSeries
    {
        Dictionary<string, SubstitutionMatrix> matrices;

        public override SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative)
        {
            // Clustal uses an algorithm to select a user matrix and scale value
            // It also loads in a special format.
            throw new NotImplementedException();
        }

        public override double GetScaleFactor(double percentIdentity, bool useNegative)
        {
            throw new NotImplementedException();
        }
    }
}
