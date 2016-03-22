using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Diagnostics;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using Microsoft.Win32;
using ClustalWPF.FileIO;

namespace ClustalWPF
{
    // Name idea: mAlign, MAlign

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public static RoutedCommand Align = new RoutedCommand();
        BackgroundWorker alignBackgroundWorker = new BackgroundWorker();
        //Stopwatch x = new Stopwatch();

        public FontSettings NamesDisplayFont
        {
            get { return (FontSettings)GetValue(NamesDisplayFontProperty); }
            set { SetValue(NamesDisplayFontProperty, value); }
        }

        public static readonly DependencyProperty NamesDisplayFontProperty =
            DependencyProperty.Register("NamesDisplayFont", typeof(FontSettings), typeof(MainWindow), new UIPropertyMetadata(new FontSettings(new Typeface("Arial"), 12)));

        public string StatusMessage
        {
            get { return (string)GetValue(StatusMessageProperty); }
        }

        public static readonly DependencyProperty StatusMessageProperty =
            DependencyProperty.Register("StatusMessage", typeof(string), typeof(MainWindow), new UIPropertyMetadata(""));



        public MainWindow()
        {
            InitializeComponent();

            TestAlign();
            //ApplicationCommands.Open.Execute(null, this);

            //Initialize the alignBackgroundWorker
            alignBackgroundWorker.DoWork += AlignBackgroundWorkerDoWork;
            alignBackgroundWorker.RunWorkerCompleted += AlignBackgroundWorkerCompleted;
            alignBackgroundWorker.ProgressChanged += AlignBackgroundWorkerProgressChanged;
            alignBackgroundWorker.WorkerReportsProgress = true;

            
        }

        void TestAlign()
        {
            Alignment x = new Alignment();
            Macromolecule y1 = new Macromolecule();
            y1.Sequence = "ACTGGGTAGA";
            y1.Name = "TestSeq1";
            y1.identifier = 1;
            int[] y11 = new int[10] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            //Macromolecule y2 = new Macromolecule();
            //y2.Sequence = "CTGGTAG";
            //y2.Name = "TestSeq1";
            //y2.identifier = 1;
            //int[] y22 = new int[7] { 1, 2, 4, 5, 6, 7, 8 };

            x.AlignmentLength = 10; // Need to set up Alignment to do this automatically when a sequence is added.

            //x.AddAlignedMacromolecule(Tuple.Create(y1, y11));
            x.AddAlignedMacromolecule(new AlignedMacromolecule(y1, y11));
            //x.AddAlignedMacromolecule(Tuple.Create(y2, y22));
            
            //.NudgeSequences(new List<int>{1}, 6, 6, 4);

            DataContext = x;
        }

        private void OpenCommandExecuted(object sender, ExecutedRoutedEventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Title = "Load sequence files";
            ofd.Multiselect = true;
            bool? success = ofd.ShowDialog();            
            
            if (success.Value)
            {
                //x.Start();
                Alignment currentAlignment = (Alignment)DataContext;
                ReturnCodes returnCode = ReturnCodes.OK;
                List<Tuple<string, int, ReturnCodes>> fileErrors;
                bool append = false;
                //x.Stop(); MessageBox.Show("Time to initialize varaibles was " + x.ElapsedMilliseconds + " ms."); Environment.Exit(1);
                
                if (currentAlignment.NumberMacromolecules > 0)
                {
                    MessageBoxResult appendChoice = MessageBox.Show("Do you want to append these sequences to the current alignment?", "", MessageBoxButton.YesNoCancel, MessageBoxImage.Question);
                    
                    switch (appendChoice)
                    {
                        case MessageBoxResult.Yes:
                            // User wants to append.
                            append = true;
                            break;
                        case MessageBoxResult.No:
                            // User doesn't want to append.
                            append = false;

                            bool cancelAction;
                            CheckToSaveChanges(out cancelAction);
                            if (cancelAction)
                            {
                                // User wants to cancel opening the file
                                return;
                            }
                            break;
                        case MessageBoxResult.Cancel:
                            // User wants to cancel opening the file.
                            return;
                    }
                }

                returnCode = SequenceFileReader.LoadSequencesFromFiles(ref currentAlignment, append, ofd.FileNames, out fileErrors);

                if (returnCode == ReturnCodes.FileErrors)
                {
                    StringBuilder errorMessage = new StringBuilder();

                    errorMessage.AppendLine("The program encountered an error loading your sequence(s):");

                    foreach (Tuple<string, int, ReturnCodes> fileError in fileErrors)
                    {
                        string fileName = fileError.Item1;
                        int numSeq = fileError.Item2;
                        ReturnCodes errorCode = fileError.Item3;
                        errorMessage.Append(' ', 4);
                        errorMessage.Append(fileName);
                        
                        switch (errorCode)
                        {
                            case ReturnCodes.EmptySequencesRemoved:
                                errorMessage.Append(" had empty sequences removed; ");
                                errorMessage.AppendFormat("{0} loaded.", numSeq);
                                break;
                            case ReturnCodes.FileNotFound:
                                errorMessage.Append(" could not be found or was inaccessible.");
                                break;
                            case ReturnCodes.NoNonEmptySequencesInFile:
                                errorMessage.Append(" contained only empty sequences.");
                                break;
                            case ReturnCodes.NoSequencesInFile:
                                errorMessage.Append(" was a recognized type, but contained no sequences.");
                                break;
                            case ReturnCodes.UnrecognizedFileType:
                                errorMessage.Append(" was not a recognized file type.");
                                break;
                            case ReturnCodes.IOException:
                                errorMessage.Append(" threw an unspecified IO error.");
                                break;
                            case ReturnCodes.MixedToNA:
                                errorMessage.Append("Some sequences were detected as protein, but since most loaded were detected as NA, they were forced to be NA; ");
                                errorMessage.AppendFormat("{0} sequences were affected.", numSeq);
                                break;
                            case ReturnCodes.MixedToProtein:
                                errorMessage.Append("Some sequences were detected as NA, but since most loaded were detected as protein, they were forced to be protein; ");
                                errorMessage.AppendFormat("{0} sequences were affected.", numSeq);
                                break;
                        }
                        errorMessage.AppendLine();
                    }
                    ErrorWindow error = new ErrorWindow();
                    error.DataContext = errorMessage.ToString();
                    error.Owner = this;
                    error.ShowDialog();
                }
            }
        }

        private void OpenCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
        // Determine whether the command can be executed or not.
        {
            e.CanExecute = true;
        }

        private void AlignCommandExecuted(object sender, ExecutedRoutedEventArgs e)
        {
            // Adjust anything necessary to show that work is being done, prevent the user form mucking it up.
            alignBackgroundWorker.RunWorkerAsync(DataContext);

            //////use old methods
            ////Alignment aligntemp = (Alignment)DataContext;
            ////ClustalW.Clustal.align(ref aligntemp);
        }

        private void AlignBackgroundWorkerDoWork(object sender, DoWorkEventArgs e)
        // Do the alignment in the background
        {
            Alignment currentAlignment = (Alignment)e.Argument;
            
            Routines.Align(ref currentAlignment, ref alignBackgroundWorker);
        }

        private void AlignBackgroundWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        // Reset things that prevented the user from mucking things up done in AlignCommandExecuted
        {

        }

        private void AlignBackgroundWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            SetValue(StatusMessageProperty, (string)e.UserState);
        }

        private void AlignCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
        // More than one sequence must be loaded for this command to be executed.
        {
            if (((Alignment)DataContext).NumberMacromolecules > 1)
            {
                e.CanExecute = true;
            }
            else
            {
                e.CanExecute = false;
            }
        }

        private void CheckToSaveChanges(out bool cancelAction)
        // Check whether the current alignment has changed, and if so, ask if the user
        // would like to save changes. If the user cancels or closes the box, return a value
        // to tell the calling code to abort.
        {
            cancelAction = false;
            if (((Alignment)DataContext).IsChanged)
            {
                MessageBoxResult save = MessageBox.Show("The current sequence alignment has unsaved changes. Would you like to save them?", "", MessageBoxButton.YesNoCancel, MessageBoxImage.Question);

                switch (save)
                {
                    case MessageBoxResult.Yes:
                        // User wants to save the current alignment.
                        ApplicationCommands.Save.Execute(null, this);
                        break;
                    case MessageBoxResult.No:
                        // User wants to discard the current alignment.
                        break;
                    case MessageBoxResult.Cancel:
                        // User wants to abort.
                        cancelAction = true;
                        break;
                }
            }
        }
    }

}
