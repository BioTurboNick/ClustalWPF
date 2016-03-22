using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;

namespace ClustalWPF.PairwiseAlignment
{
    interface IPairwiseAlignmentAlgorithm
    {
        void PairwiseAlign(ref Alignment alignmentObject, ref double[,] distanceMatrix);
    }
}
