using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClustalWPF.SubstitutionMatrix
{
    public abstract class SubstitutionMatrixSeries
    {
        public abstract SubstitutionMatrix GetMatrix(double percentIdentity, double minLength, bool useNegative);

        public abstract double GetScaleFactor(double percentIdentity, bool useNegative);
    }
}
