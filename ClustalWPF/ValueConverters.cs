using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Data;
using System.Windows.Media;

namespace ClustalWPF
{
    internal sealed class SequenceAlignmentConverter : IMultiValueConverter
    // This class provides for conversion from a sequence string and aligned position pair to an aligned string.
    // Convert and ConvertBack provide the interfaces for XAML data binding converters.
    {
        public object Convert(object[] values, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            if (!Routines.ValidateConverterInput(values))
            {
                return null;
            }

            string sequence = (string)values[0];
            int[] alignedPositions = (int[])values[1];
            int alignmentLength = (int)values[2];

            string alignedSequence = Alignment.ToAlignedSequence(sequence, alignedPositions, alignmentLength);

            return alignedSequence;
        }

        public object[] ConvertBack(object value, Type[] targetTypes, object parameter, System.Globalization.CultureInfo culture)
        {
            if (!Routines.ValidateConverterInput(value))
            {
                return null;
            }

            string alignedSequence = (string)value;
            int[] alignedPositions;
            
            string sequence = Alignment.ToSequence(alignedSequence, out alignedPositions);

            object[] returnArray = new object[3];

            returnArray[0] = sequence;
            returnArray[1] = alignedPositions;
            returnArray[2] = alignedSequence.Length;

            return returnArray;
        }
    }

    internal sealed class NamesLengthConverter : IMultiValueConverter
    {
     public object Convert(object[] values, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            if (!Routines.ValidateConverterInput(values))
            {
                return null;
            }
            
            IList<AlignedMacromolecule> alignedMacromolecules = (IList<AlignedMacromolecule>)values[0];
            FontSettings fontSettings = (FontSettings)values[1];
            double adjustment = System.Convert.ToDouble(parameter);
                        
            double maxDisplayLength = 0;

            foreach (AlignedMacromolecule alignedMacromolecule in alignedMacromolecules)
            {
                string testMacromoleculeName = alignedMacromolecule.Name;

                double testDisplayLength = Routines.GetTextDisplayLength(testMacromoleculeName, fontSettings.Typeface.FontFamily, fontSettings.Typeface.Style, fontSettings.Typeface.Weight, fontSettings.Typeface.Stretch, fontSettings.FontSize);
                
                if (testDisplayLength > maxDisplayLength)
                {
                    maxDisplayLength = testDisplayLength;
                }
            }

            return maxDisplayLength + adjustment;
        }

        public object[] ConvertBack(object value, Type[] targetTypes, object parameter, System.Globalization.CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
