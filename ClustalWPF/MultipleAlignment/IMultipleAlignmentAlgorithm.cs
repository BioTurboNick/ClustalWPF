using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;

namespace ClustalWPF.MultipleAlignment
{
    interface IMultipleAlignmentAlgorithm
    {
        double Align(ref Alignment alignmentObject, SubstitutionMatrix.SubstitutionMatrixSeries subMatrixClass, Tuple<ReadOnlyCollection<AlignedMacromolecule>, ReadOnlyCollection<AlignedMacromolecule>, double> step);
    }
}
