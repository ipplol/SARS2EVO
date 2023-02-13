using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SampleMutlist转氨基酸挑Spike
{
    class MutationAnno //突变注释
    {
        public int Position;
        public string Titv;
        public string Gene;
        public string NSS;
        public string AminoChange;
    }
    class Program
    {
        static void Main(string[] args)
        {
            string workfold = "***";// e.g. "//NAS8500/g/VariationMutation/VirusAntibodyEscape"
            Dictionary<string, MutationAnno> DicMutAnno = new Dictionary<string, MutationAnno>();
            StreamReader read = new StreamReader(workfold + "/Data/AllMutationAnnomutresult.tsv");
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                MutationAnno newa = new MutationAnno();
                newa.Position = Convert.ToInt32(line1[2]);
                newa.AminoChange = line1[6];
                newa.Gene = line1[5];
                newa.NSS = line1[7];
                newa.Titv = line1[1];
                DicMutAnno.Add(line1[2] + line1[4], newa);
                line = read.ReadLine();
            }
            read.Close();

            int i, j, k;
            StreamWriter write = new StreamWriter(workfold + "/Data/sampled_sample_mutlist.transfer.tsv");
            write.Write("Seq\tOriginal\tCount\tSpikeNuc\tCount\tRBDNuc\tCount\tNTDNuc\tCount\tGenomeAA\tCount\tSpikeAA\tCount\tRBDAA\tCount\tNTDAA\tCount\tORF1abAA\tcount\tNGeneAA\tCount\n");
            //write.Write("Seq\tOriginal\tCount\tRBDNuc\tCount\tRBDAA\tCount\n");

            read = new StreamReader(workfold + "/Data/sampled_sample_mutlist.tsv");
            line = read.ReadLine();
            while(line!=null)
            {
                line = line.Replace(' ', ',').Replace('\t', ':');
                if (!line.Contains(":")) line += ":";
                string[] line1 = line.Split(':');//:
                string[] mut = line1[1].Split(',');//,
                string output = line1[0] + "\t" + line1[1] +"\t" + (line1[1].Length-line1[1].Replace("/","").Length) + "\t";
                int count = 0;
                string transfered = "";

                //SpikeNuc
                count = 0;
                transfered = "";
                for (i=0;i<mut.Length;i++)
                    if(mut[i]!="")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 21563 && pos <= 25384)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            transfered += mut[i] + " ";
                            count++;
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //RBDNuc
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 22517 && pos <= 23185)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            transfered += mut[i] + " ";
                            count++;
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //NTDNuc
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 21599 && pos <= 22477)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            transfered += mut[i] + " ";
                            count++;
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //GenomeAA
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                        if (RefVar[0].Length > 1 || RefVar[1].Length > 1)
                        {
                            transfered += Convert.ToString(pos) + "inDel ";
                            count++;
                        }
                        else
                            if (DicMutAnno.ContainsKey(mut[i]) && DicMutAnno[mut[i]].NSS == "NS")
                        {
                            transfered += DicMutAnno[mut[i]].AminoChange + " ";
                            count++;
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //SpikeAA
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 21563 && pos <= 25384)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            if (RefVar[0].Length > 1 || RefVar[1].Length > 1)
                            {
                                transfered += Convert.ToString(pos) + "inDel ";
                                count++;
                            }
                            else
                                if (DicMutAnno.ContainsKey(mut[i]) && DicMutAnno[mut[i]].NSS == "NS")
                            {
                                transfered += DicMutAnno[mut[i]].AminoChange + " ";
                                count++;
                            }
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //RBDAA
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 22517 && pos <= 23185)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            if (RefVar[0].Length > 1 || RefVar[1].Length > 1)
                            {
                                transfered += Convert.ToString(pos) + "inDel ";
                                count++;
                            }
                            else
                                if (DicMutAnno.ContainsKey(mut[i]) && DicMutAnno[mut[i]].NSS == "NS")
                            {
                                transfered += DicMutAnno[mut[i]].AminoChange + " ";
                                count++;
                            }
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" +Convert.ToString(count)+ "\t";

                //NTDAA
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 21599 && pos <= 22477)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            if (RefVar[0].Length > 1 || RefVar[1].Length > 1)
                            {
                                transfered += Convert.ToString(pos) + "inDel ";
                                count++;
                            }
                            else
                                if (DicMutAnno.ContainsKey(mut[i]) && DicMutAnno[mut[i]].NSS == "NS")
                            {
                                transfered += DicMutAnno[mut[i]].AminoChange + " ";
                                count++;
                            }
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //Orf1abAA
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 266 && pos <= 21552)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            if (RefVar[0].Length > 1 || RefVar[1].Length > 1)
                            {
                                transfered += Convert.ToString(pos) + "inDel ";
                                count++;
                            }
                            else
                                if (DicMutAnno.ContainsKey(mut[i]) && DicMutAnno[mut[i]].NSS == "NS")
                            {
                                transfered += DicMutAnno[mut[i]].AminoChange + " ";
                                count++;
                            }
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count) + "\t";

                //NGeneAA
                count = 0;
                transfered = "";
                for (i = 0; i < mut.Length; i++)
                    if (mut[i] != "")
                    {
                        k = 0;
                        while (mut[i][k] <= '9' && mut[i][k] >= '0') k++;
                        int pos = Convert.ToInt32(mut[i].Substring(0, k));
                        if (pos >= 28274 && pos <= 29533)
                        {
                            string[] RefVar = mut[i].Substring(k, mut[i].Length - k).Split('/');
                            if (RefVar[0].Length > 1 || RefVar[1].Length > 1)
                            {
                                transfered += Convert.ToString(pos) + "inDel ";
                                count++;
                            }
                            else
                                if (DicMutAnno.ContainsKey(mut[i]) && DicMutAnno[mut[i]].NSS == "NS")
                            {
                                transfered += DicMutAnno[mut[i]].AminoChange + " ";
                                count++;
                            }
                        }
                    }
                if (transfered != "") transfered = transfered.Substring(0, transfered.Length - 1);
                output += transfered + "\t" + Convert.ToString(count);

                write.Write(output + "\n");
                line = read.ReadLine();
            }
            read.Close();
            write.Close();
        }
    }
}
