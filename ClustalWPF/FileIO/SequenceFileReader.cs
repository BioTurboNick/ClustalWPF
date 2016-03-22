using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace ClustalWPF.FileIO
{
    class SequenceFileReader
    // This class contains methods pertaining to reading files.
    {

        internal static ReturnCodes LoadSequencesFromFiles(ref Alignment alignmentObject, bool append, string[] fileNames, out List<Tuple<string, int, ReturnCodes>> fileErrors)
        // This method provides the main interface for loading sequences. It takes an array of file names and calls LoadSequencesFromFile.
        // It then checks the return codes and the number of sequences loaded, then clears the alignment if not appending, then adds the
        // new sequences to the alignment. In addition to the main return code, it returns as an output parameter a list containing the error codes.
        {
            ReturnCodes returnCode = ReturnCodes.OK;
            fileErrors = new List<Tuple<string, int, ReturnCodes>>();
            //List<Tuple<Macromolecule, int[]>> loadedAlignedMacromolecules = new List<Tuple<Macromolecule, int[]>>();
            List<AlignedMacromolecule> loadedAlignedMacromolecules = new List<AlignedMacromolecule>();

            foreach (string fileName in fileNames)
            {
                //List<Tuple<Macromolecule, int[]>> fileAlignedMacromolecules;
                List<AlignedMacromolecule> fileAlignedMacromolecules;
                ReturnCodes fileReturnCode = ReturnCodes.OK;
                fileReturnCode = LoadSequencesFromFile(fileName, out fileAlignedMacromolecules);

                if (fileReturnCode != ReturnCodes.OK)
                {
                    // Add any errors to the array
                    fileErrors.Add(Tuple.Create(fileName, fileAlignedMacromolecules.Count, fileReturnCode));
                    returnCode = ReturnCodes.FileErrors;
                }
                
                if (fileReturnCode == ReturnCodes.OK || fileReturnCode == ReturnCodes.EmptySequencesRemoved)
                {
                    // These codes mean it is okay to add their macromolecules to the main list.
                    //foreach (Tuple<Macromolecule, int[]> alignedMacromolecule in fileAlignedMacromolecules)
                    foreach (AlignedMacromolecule alignedMacromolecule in fileAlignedMacromolecules)
                    {
                        loadedAlignedMacromolecules.Add(alignedMacromolecule);
                    }
                }
            }
            
            // Ensure that some sequences were loaded.
            if (loadedAlignedMacromolecules.Count > 0)
            {   
                // Check that all are the same type, force the minority to conform to the majority, and append to the error list.
                int numChanged;
                ReturnCodes typeReturnCode = EnsureAllSameType(loadedAlignedMacromolecules, out numChanged);
                if (typeReturnCode != ReturnCodes.OK) // 
                {
                    fileErrors.Add(Tuple.Create("DNA to Protein", numChanged, typeReturnCode));
                    returnCode = ReturnCodes.FileErrors;
                }

                if (!append)
                {
                    alignmentObject.Clear();
                }

                // It's now okay to add macromolecules to the alignment object
                alignmentObject.AddAlignedMacromolecules(loadedAlignedMacromolecules);               
            }

            return returnCode;
        }
        
        //static ReturnCodes LoadSequencesFromFile(string fileName, out List<Tuple<Macromolecule, int[]>> loadedAlignedMacromolecules)
        static ReturnCodes LoadSequencesFromFile(string fileName, out List<AlignedMacromolecule> loadedAlignedMacromolecules)
        // Opens the file given by fileName, reads it, and passes the loaded macromolecules back as an output parameter.
        // In addition, it replaces blank names with "Unnamed Sequence" and removes 
        {
            ReturnCodes returnCode = ReturnCodes.OK;
            SequenceFileParser fileParser;
            StreamReader fileReader;
            //loadedAlignedMacromolecules = new List<Tuple<Macromolecule, int[]>>();
            loadedAlignedMacromolecules = new List<AlignedMacromolecule>();

            // Check that fileName is not empty
            if (fileName.Length == 0) // Note: checking the length of a string as 0 is more efficient than comparing it to an empty string.
            {
                return ReturnCodes.NoFileName;
            }

            // IO operations are vulnerable to file system errors.
            try
            {
                // Open the file
                fileReader = new StreamReader(fileName);

                try
                {
                    // Determine the file type to choose the correct parser.
                    fileParser = DetermineFileType(ref fileReader);

                    // Read the file
                    returnCode = fileParser.ReadFile(out loadedAlignedMacromolecules);
                }
                catch (System.IO.FileNotFoundException ex)
                {
                    return ReturnCodes.FileNotFound;
                }
                catch (System.IO.IOException ex)
                {
                    return ReturnCodes.IOException;
                }
                finally
                {
                    // The file is open, so we need to close it
                    fileReader.Close();
                }
            }
            catch (FileNotFoundException ex)
            {
                // File wasn't opened, so just return the error code
                return ReturnCodes.FileNotFound;
            }            

            if (returnCode != ReturnCodes.OK)
            {
                return returnCode;
            }

            // Check that at least 1 sequence was loaded
            if (loadedAlignedMacromolecules.Count == 0)
            {
                return ReturnCodes.NoSequencesInFile;
            }

            // Check each sequence to ensure it is not empty.
            // Set the return code to reflect this.
            for (int i = 0; i < loadedAlignedMacromolecules.Count; i++)
            {
                //Tuple<Macromolecule, int[]> alignedMacromolecule = loadedAlignedMacromolecules[i];
                AlignedMacromolecule alignedMacromolecule = loadedAlignedMacromolecules[i];

                //if (alignedMacromolecule.Item1.Sequence.Length == 0) // empty sequence
                if (alignedMacromolecule.Sequence.Length == 0) // empty sequence
                {
                    returnCode = ReturnCodes.EmptySequencesRemoved;
                    loadedAlignedMacromolecules.Remove(alignedMacromolecule);
                }
            }

            // Check that we still have at least 1 sequence
            if (loadedAlignedMacromolecules.Count == 0)
            {
                return ReturnCodes.NoNonEmptySequencesInFile;
            }

            // Set any no-name molecules to "Unnamed Sequence"
            //foreach (Tuple<Macromolecule, int[]> alignedMacromolecule in loadedAlignedMacromolecules)
            foreach (AlignedMacromolecule alignedMacromolecule in loadedAlignedMacromolecules)
            {
                //Macromolecule macromolecule = alignedMacromolecule.Item1;
               if (alignedMacromolecule.Name.Length == 0) // no name
                {
                    alignedMacromolecule.Name = "Unnamed Sequence";
                }  
            }
            
            return returnCode;
        }

        public static SequenceFileParser DetermineFileType(ref StreamReader fileReader)
        // Read the first lines of a file to determine what type it is so that the appropriate parser can be assigned.
        {
            string lineIn = "";
            SequenceFileParser fileParser;
            InputFileTypes fileType = InputFileTypes.None; // To be removed once Clustal method is replicated.

            // Find the first non-whitespace line
            while (!fileReader.EndOfStream)
            {
                lineIn = fileReader.ReadLine();
                if (!String.IsNullOrWhiteSpace(lineIn))
                {
                    break;
                }
            }

            // Trim trailing whitespace.
            lineIn.TrimEnd();

            // Evaluate the first line to determine the file type.
            if (lineIn.StartsWith("CLUSTAL", StringComparison.OrdinalIgnoreCase))
            {
                // Should go first because it will be one of the most frequent files opened.
                fileParser = new ClustalFileParser(ref fileReader, InputFileTypes.Clustal);
            }
            else if (lineIn.StartsWith(">", StringComparison.Ordinal))
            {
                // FASTA/Pearson is probably the next-most common type. Also indicates PIR.
                if (lineIn.Length >= 3 && lineIn[3] == ';')
                {
                    fileType = InputFileTypes.PIR;
                }
                else
                {
                    fileType = InputFileTypes.FASTAPearson;
                }
                throw new NotImplementedException();
            }
            else if (lineIn.StartsWith("ID", StringComparison.OrdinalIgnoreCase))
            {
                fileType = InputFileTypes.EMBLSwissProt;
                throw new NotImplementedException();
            }
            else if (lineIn.StartsWith("LOCUS", StringComparison.OrdinalIgnoreCase))
            {
                fileType = InputFileTypes.GenBank;
                throw new NotImplementedException();
            }
            else if (lineIn.StartsWith("!!AA_MULTIPLE_ALIGNMENT", StringComparison.OrdinalIgnoreCase) ||
                lineIn.StartsWith("!!NA_MULTIPLE_ALIGNMENT", StringComparison.OrdinalIgnoreCase) ||
                lineIn.StartsWith("MSF", StringComparison.OrdinalIgnoreCase) ||
                lineIn.StartsWith("PILEUP", StringComparison.OrdinalIgnoreCase))
            {
                fileType = InputFileTypes.MSF;
                throw new NotImplementedException();
            }
            else if (lineIn.StartsWith("!!RICH_SEQUENCE", StringComparison.OrdinalIgnoreCase))
            {
                fileType = InputFileTypes.RSF;
                throw new NotImplementedException();
            }
            else if (lineIn.StartsWith("\"", StringComparison.Ordinal) ||
                lineIn.StartsWith("%", StringComparison.Ordinal) ||
                lineIn.StartsWith("#", StringComparison.Ordinal))
            {
                fileType = InputFileTypes.GDE;
                throw new NotImplementedException();

                // Note: % = Protein, # = DNA
            }
            else if (lineIn.StartsWith("#NEXUS", StringComparison.OrdinalIgnoreCase))
            {
                fileType = InputFileTypes.NEXUS;
                throw new NotImplementedException();
            }
            else
            {
                fileType = InputFileTypes.None;
                throw new NotImplementedException();
            }

            // Rewind the FileReader
            fileReader.DiscardBufferedData();
            fileReader.BaseStream.Position = 0;

            return fileParser;
        }

        //static ReturnCodes EnsureAllSameType(List<Tuple<Macromolecule, int[]>> alignedMacromolecules, out int numChanged)
        static ReturnCodes EnsureAllSameType(List<AlignedMacromolecule> alignedMacromolecules, out int numChanged)
        // Counts the number of loaded molecules which are NA and Protein and forces the minority to become the majority.
        // Also returns a return code if it does this, and outputs the number changed in the
        {   
            int proteinCount = 0;
            int NACount = 0;
            numChanged = 0;
            bool forceDNA;
            //foreach (Tuple<Macromolecule, int[]> alignedMacromolecule in alignedMacromolecules)
            foreach (AlignedMacromolecule alignedMacromolecule in alignedMacromolecules)
            {
                //if (alignedMacromolecule.Item1.IsNucleicAcid)
                if (alignedMacromolecule.IsNucleicAcid)
                {
                    NACount++;
                }
                else
                {
                    proteinCount++;
                }
            }

            if (NACount > 0 && proteinCount >= NACount) // Though I doubt it'll ever happen, we'll break the tie to favor protein.
            {
                numChanged = NACount;
                //foreach (Tuple<Macromolecule, int[]> alignedMacromolecule in alignedMacromolecules)
                foreach (AlignedMacromolecule alignedMacromolecule in alignedMacromolecules)
                {
                    //alignedMacromolecule.Item1.IsNucleicAcid = false;
                    alignedMacromolecule.IsNucleicAcid = false;
                }
                return ReturnCodes.MixedToProtein;
            }
            else if (proteinCount > 0 && NACount > proteinCount)
            {
                numChanged = proteinCount;
                //foreach (Tuple<Macromolecule, int[]> alignedMacromolecule in alignedMacromolecules)
                foreach (AlignedMacromolecule alignedMacromolecule in alignedMacromolecules)
                {
                    //alignedMacromolecule.Item1.IsNucleicAcid = true;
                    alignedMacromolecule.IsNucleicAcid = true;
                }
                return ReturnCodes.MixedToNA;
            }

            return ReturnCodes.OK;
        }
        
    }

    enum InputFileTypes
    {
        None,
        FASTAPearson,
        EMBLSwissProt,
        GenBank,
        Clustal,
        MSF,
        RSF,
        NEXUS,
        PIR,
        GDE
    }
}
