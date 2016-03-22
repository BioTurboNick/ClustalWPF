using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ClustalWPF
{
    public class Macromolecule
    // This class defines a Macromolecule object which holds all the information and methods related to a macromolecule.
    {
        string sequence = "";
        string name = "";
        internal string title = "";
        internal ulong identifier = 0;
        bool isNucleicAcid = false;

        public string Sequence
        {
            get { return sequence; }
            set { sequence = value; }
        }

        public string Name
        {
            get { return name; }
            set { name = value; }
        }

        public bool IsNucleicAcid
        {
            get { return isNucleicAcid; }
            set { isNucleicAcid = value; }
        }

        protected Macromolecule(Macromolecule newMacromolecule)
        {
            sequence = newMacromolecule.Sequence;
            name = newMacromolecule.Name;
            isNucleicAcid = newMacromolecule.IsNucleicAcid;
        }

        public Macromolecule()
        {}
    }
}
