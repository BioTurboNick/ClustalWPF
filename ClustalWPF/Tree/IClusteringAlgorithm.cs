using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;

namespace ClustalWPF.Tree
{
    interface IClusteringAlgorithm<T>
    {
        GuideTree<T> GenerateTree(double[,] distanceMatrix, T[] leafDataList);
    }
}
