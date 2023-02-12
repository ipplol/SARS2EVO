using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 计算指定突变在进化枝的RoHo
{
    public class Mutation
    {
        public string MutationName;
        public List<double> RoHoScores = new List<double>();
        public List<string> MutDates = new List<string>();
    }
    public class Clade
    {
        public string CladeName;
        public Dictionary<string,Mutation> MutationList = new Dictionary<string, Mutation>();
    }
    public class Node
    {
        public string Lineage;
        public string EarliestDate;
    }
    class Program
    {
        static List<Clade> CladeList = new List<Clade>();
        static Dictionary<string, string> MutRefDic = new Dictionary<string, string>();//Spike nuc mutation 2 aa mutation
        static Dictionary<string, Node> NodeDic = new Dictionary<string, Node>();
        static string Workfold = "***";
        static void ReadinMutAndClade()
        {
            int i, j, k;
            StreamReader read = new StreamReader(Workfold + "/Usher1128/RoHo/RoHoMutationList.tsv");
            List<string> tmpm = new List<string>();
            string line = read.ReadLine();
            while(line!=null)
            {
                tmpm.Add(line);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(Workfold + "/Usher1128/RoHo/RoHoCladeList.tsv");
            line = read.ReadLine();
            while(line!=null)
            {
                Clade newc = new Clade();
                newc.CladeName = line;
                for(i=0;i<tmpm.Count;i++)
                {
                    Mutation newm = new Mutation();
                    newm.MutationName = tmpm[i];
                    if(!newc.MutationList.ContainsKey(tmpm[i]))
                        newc.MutationList.Add(tmpm[i], newm);
                }
                CladeList.Add(newc);
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinMutationMap()
        {
            StreamReader read = new StreamReader(Workfold + "/Data/AllMutationAnnomutresult.tsv");
            int i, j, k;
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                string nuc = line1[4][0] + line1[2] + line1[4][2];
                string AA = "OutOfRange";
                if (line1[5] == "S")
                    AA = line1[6].Substring(1, line1[6].Length - 1);
                MutRefDic.Add(nuc, AA);
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinNodeinfo()
        {
            StreamReader read = new StreamReader(Workfold + "/Usher1128/Data/MAT1128.json.nodeinfo");
            int i, j, k;
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if(line1[0].Contains("node_"))
                {
                    Node newd = new Node();
                    newd.EarliestDate = line1[1];
                    newd.Lineage = line1[3];
                    NodeDic.Add(line1[0], newd);
                }
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinRoHo()
        {
            StreamReader read = new StreamReader(Workfold + "/Usher1128/RoHo/RoHo.tsv");
            int i, j, k;
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                string AAmut = MutRefDic[line1[0]];
                string clade = "NA";
                if(NodeDic.ContainsKey(line1[1]))clade = NodeDic[line1[1]].Lineage;
                if (clade.Contains("BM") || clade.Contains("BL"))
                    clade = "BA.2";
                if (clade.Contains("BA.5")|| clade.Contains("BF")|| clade.Contains("BE"))
                    clade = "BA.4";

                //定位突变枝
                for (k = 0; k < CladeList.Count; k++)
                {
                    //先找能不能对上的
                    if (clade == CladeList[k].CladeName)
                        break;
                }
                if (k >= CladeList.Count)
                {
                    string[] tmp1 = clade.Split('.');//再找Lineage的上一级
                    string tmps = "";
                    for (j = 0; j < tmp1.Length - 1; j++) tmps += tmp1[j] + ".";
                    if (tmps != "") tmps = tmps.Substring(0, tmps.Length - 1);
                    for (k = 0; k < CladeList.Count; k++)
                        if (tmps == CladeList[k].CladeName)
                            break;
                }

                if(k<CladeList.Count()&&CladeList[k].MutationList.ContainsKey(AAmut))
                {
                    CladeList[k].MutationList[AAmut].RoHoScores.Add(Convert.ToDouble(line1[6]));
                    CladeList[k].MutationList[AAmut].MutDates.Add(NodeDic[line1[1]].EarliestDate);
                }
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void WriteResult()
        {
            StreamWriter write = new StreamWriter(Workfold + "/Usher1128/RoHo/MutationRoHo.txt");
            write.WriteLine("Clade\tMutation\tDate\tRoHo\tGroup");
            int i, j, k;
            for(i=0;i<CladeList.Count;i++)
            {
                foreach(string val in CladeList[i].MutationList.Keys)
                {
                    CladeList[i].MutationList[val].RoHoScores.Sort();
                    for (j = 0; j < CladeList[i].MutationList[val].RoHoScores.Count; j++)
                    {
                        string group = "NA";
                        if (CladeList[i].MutationList[val].RoHoScores[j] > 0) group = "Positive";
                        if (CladeList[i].MutationList[val].RoHoScores[j] < 0) group = "Negative";
                        if(group!="NA")
                            write.WriteLine(CladeList[i].CladeName + "\t" + CladeList[i].MutationList[val].MutationName + "\t" + CladeList[i].MutationList[val].MutDates[j] + "\t" + CladeList[i].MutationList[val].RoHoScores[j] + "\t" + group);
                    }
                }
            }
            write.Close();

            write = new StreamWriter(Workfold + "/Usher1128/RoHo/MutationRoHo_average.txt");
            for (i = 0; i < CladeList.Count; i++)
            {
                foreach (string val in CladeList[i].MutationList.Keys)
                {
                    if (CladeList[i].MutationList[val].RoHoScores.Count >= 3)
                    {
                        double meanRoHo = 0;
                        for (j = 0; j < CladeList[i].MutationList[val].RoHoScores.Count; j++)
                            meanRoHo += CladeList[i].MutationList[val].RoHoScores[j];
                        meanRoHo /= CladeList[i].MutationList[val].RoHoScores.Count;
                        write.WriteLine(CladeList[i].CladeName + "\t" + CladeList[i].MutationList[val].MutationName + "\t" + meanRoHo);
                    }
                    else
                        write.WriteLine(CladeList[i].CladeName + "\t" + CladeList[i].MutationList[val].MutationName + "\tNA");
                }
            }
            write.Close();
            return;
        }
        static void Main(string[] args)
        {
            ReadinMutAndClade();
            ReadinMutationMap();
            ReadinNodeinfo();
            ReadinRoHo();
            WriteResult();
        }
    }
}
