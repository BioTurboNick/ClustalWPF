using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace ClustalWPF.FileIO
{
    internal class ClustalFileParser : SequenceFileParser
    // This class defines a parser object which reads a Clustal file and outputs the contained information.
    // 
    // Format:
    //// CLUSTAL 2.1 multiple sequence alignment
    //// 
    //// Sequence1     <aligned sequence 1-60> 60
    //// Sequence2     <aligned sequence 1-60> 60
    //// Sequence3     <aligned sequence 1-60> 60
    ////               <optional degree of conservation line>
    //// 
    //// Sequence1     <sequence bases 61-120> 120
    //// Sequence2     <sequence bases 61-120> 120
    //// Sequence3     <sequence bases 61-120> 120
    ////               <optional degree of conservation line>
    //
    // Format notes:
    // The number at the end of the line is optional.
    // The degree of conservation line uses:
    //     * = identical
    //     : = conserved substitutions
    //     . = semi-conserved substitutions
    // These elements are only included to enhance human-readability of the file. We don't care about them.
    {
        internal ClustalFileParser(ref StreamReader newFileReader, InputFileTypes newFileType)
        {
            fileReader = newFileReader;
            fileType = newFileType;
        }

        //internal override ReturnCodes ReadFile(out List<Tuple<Macromolecule, int[]>> loadedAlignedMacromolecules)
        internal override ReturnCodes ReadFile(out List<AlignedMacromolecule> loadedAlignedMacromolecules)
        {
            //loadedAlignedMacromolecules = new List<Tuple<Macromolecule, int[]>>();
            loadedAlignedMacromolecules = new List<AlignedMacromolecule>();
            ReturnCodes returnCode = ReturnCodes.OK;
            string lineIn = "";
            Dictionary<string, StringBuilder> alignedSequenceBuilderDict = new Dictionary<string, StringBuilder>();

            // Find the first non-whitespace line. This is the header, and can be ignored.
            while (!fileReader.EndOfStream)
            {
                lineIn = fileReader.ReadLine();
                if (!String.IsNullOrWhiteSpace(lineIn))
                {
                    break;
                }
            }

            // Scan through the remainder of the file to read the sequences.
            while (!fileReader.EndOfStream)
            {
                lineIn = fileReader.ReadLine();

                // Discard uninformative lines
                if (IsClustalBlankLine(lineIn))
                {
                    continue;
                }

                // Get the first two space-separated elements, name and up to 60 bases of sequence
                string[] lineElements = lineIn.Split(new char[]{ ' ' }, 2, StringSplitOptions.RemoveEmptyEntries);
                
                string name = lineElements[0];
                string sequenceLine = lineElements[1];

                if (!alignedSequenceBuilderDict.ContainsKey(name))
                {
                    alignedSequenceBuilderDict.Add(name, new StringBuilder(sequenceLine));
                }
                else
                {
                    alignedSequenceBuilderDict[name].Append(sequenceLine);
                }
            }

            // Now that the sequences have been read, load them into the Macromolecule object,
            // create the aligned positions array, and determine if it is a protein or DNA.
            foreach (KeyValuePair<string, StringBuilder> alignedSeqBuilderItem in alignedSequenceBuilderDict)
            {
                Macromolecule macromolecule = new Macromolecule();
                int[] alignedPos;

                macromolecule.Name = alignedSeqBuilderItem.Key;
                macromolecule.Sequence = Alignment.ToSequence(alignedSeqBuilderItem.Value.ToString(), out alignedPos);
                //alignedMacromolecule.Macromolecule.Sequence = Alignment.ToSequence(alignedSeqBuilderItem.Value.ToString(), out alignedPos);
                macromolecule.IsNucleicAcid = Routines.IsNucleicAcid(macromolecule.Sequence);
                
                //loadedAlignedMacromolecules.Add(Tuple.Create(macromolecule, alignedPos));
                loadedAlignedMacromolecules.Add(new AlignedMacromolecule(macromolecule, alignedPos));
            }
            
            return returnCode;
        }

        bool IsClustalBlankLine(string line)
        // This function checks whether a line should be discarded.
        // Only lines that contain letters or '-' are meaningful.
        {
            if (line.Length == 0 || line.StartsWith("!", StringComparison.Ordinal))
            {
                return true;
            }

            foreach (char lineChar in line)
            {
                if (Char.IsLetter(lineChar) || lineChar == '-')
                {
                // Line is informative
                    return false;
                }
            }

            return true;
        }
    }
}
