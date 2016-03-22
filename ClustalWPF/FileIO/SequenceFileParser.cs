using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace ClustalWPF.FileIO
{
    
    internal abstract class SequenceFileParser
    {
        protected StreamReader fileReader;
        protected InputFileTypes fileType;

        //internal abstract ReturnCodes ReadFile(out List<Tuple<Macromolecule, int[]>> loadedMacromolecules);
        internal abstract ReturnCodes ReadFile(out List<AlignedMacromolecule> loadedMacromolecules);
        
    }
}
