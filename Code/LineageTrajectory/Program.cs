using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 汇总Lineage数据_准备画图用的数据
{
    public class Lineage
    {
        public int NodeInd;//节点编号
        public string Name;//Lineage名称
        public int FatherInd;//父节点编号
        public List<int> ChildList = new List<int>();//子节点列表

        public double DMS_Expr_Bind=0;//DMS RBD表达与ACE结合力
        public double Escape_normalized=0;//逃逸得分
        public double Delta_ACEbind = 0;//相比于父节点的变化
        public double Delta_ESCscore = 0;

        public string CollectionDate5P;//前5%序列的采样时间
        public int CountryNumber;//在多少国家检出
        public int SeqNumber;//序列数量
        public List<string> MutationList = new List<string>();//突变列表
        public List<string> NewMutList = new List<string>();//相比于父节点的新增突变
        public string NewMutString = "";
        public List<string> BackMutList = new List<string>();//相比于父节点的回复突变
        public string BackMutString = "";
        public string VOC = "Other";
        //public string DMSGroup = "Wuhan-Hu-1_v1";//基于什么数据算DMS得分
        public string DMSGroup = "Wuhan-Hu-1_v1";//基于什么数据算DMS得分
        public bool pick4plot = false;//是否用于画图
        public int MCAwithExtraMut = 0;//拥有额外突变的最近祖先节点
    }
    public class AAchangeDMS
    {
        public double deltaExpr = 0; // single mutation score
        public double deltaBind = 0;
    }
    public class DMS_bindexpr
    {
        public Dictionary<string, AAchangeDMS> MutScore = new Dictionary<string, AAchangeDMS>(); //single variant score
    }
    public class AAchangeEscape
    {
        public Dictionary<string, double> MutEscape = new Dictionary<string, double>();
    }
    class Program
    {
        static List<Lineage> LineageList = new List<Lineage>();
        static Dictionary<string, int> LineageMap = new Dictionary<string, int>();
        static string workfold = "***";
        static Dictionary<string, DMS_bindexpr> DMSScore = new Dictionary<string, DMS_bindexpr>();//final variants score from bloom
        static Dictionary<string, AAchangeEscape> EscapeScore = new Dictionary<string, AAchangeEscape>();//variants的逃逸得分
        static string ImmuneBackground = "WT";//用什么突变背景算逃逸 WT BA.2
        static void ReadinRelations()//读入Lineage节点间的关系
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/Data/lineageRelations.tsv");
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                Lineage newl = new Lineage();
                newl.NodeInd = Convert.ToInt32(line1[0]);
                newl.Name = line1[1];
                newl.FatherInd = Convert.ToInt32(line1[2]);
                string[] child = line1[4].Split(',');
                for (i = 0; i < child.Count(); i++)
                    if(child[i]!="")
                        newl.ChildList.Add(Convert.ToInt32(child[i]));
                LineageMap.Add(newl.Name, newl.NodeInd);
                LineageList.Add(newl);
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinDateAAmut()//读入节点的采样时间和突变，计算与父节点的差异
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/Data/Lineage.Date.AAmut");
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (LineageMap.ContainsKey(line1[0]))
                {
                    int ind = LineageMap[line1[0]];
                    LineageList[ind].CollectionDate5P = line1[1];
                    LineageList[ind].CountryNumber = Convert.ToInt32(line1[4]);
                    LineageList[ind].SeqNumber = Convert.ToInt32(line1[5]);
                    string[] mut = line1[3].Split(',');
                    for (i = 0; i < mut.Count(); i++)
                        if (mut[i] != "")
                            LineageList[ind].MutationList.Add(mut[i].Substring(1, mut[i].Length - 1));
                    LineageList[ind].MutationList.Sort();
                }
                line = read.ReadLine();
            }

            //assign每个节点所属的voc，计算差异
            LineageList[LineageMap["B.1.1.7"]].VOC = "Alpha";
            LineageList[LineageMap["B.1.351"]].VOC = "Beta";
            LineageList[LineageMap["P.1"]].VOC = "Gamma";
            LineageList[LineageMap["B.1.617.2"]].VOC = "Delta";
            LineageList[LineageMap["BA.1"]].VOC = "Omicron";
            LineageList[LineageMap["BA.2"]].VOC = "Omicron";
            LineageList[LineageMap["BA.3"]].VOC = "Omicron";
            LineageList[LineageMap["BA.4"]].VOC = "Omicron";
            LineageList[LineageMap["BA.5"]].VOC = "Omicron";
            //assign每个节点所属的DMS数据
            LineageList[LineageMap["B.1.1.7"]].DMSGroup = "Alpha";
            LineageList[LineageMap["B.1.351"]].DMSGroup = "Beta";
            LineageList[LineageMap["P.1"]].DMSGroup = "Beta";
            LineageList[LineageMap["B.1.525"]].DMSGroup = "Eta";
            LineageList[LineageMap["B.1.617.2"]].DMSGroup = "Delta";
            LineageList[LineageMap["BA.1"]].DMSGroup = "Omicron_BA1";
            LineageList[LineageMap["BA.2"]].DMSGroup = "Omicron_BA2";
            LineageList[LineageMap["BA.3"]].DMSGroup = "Omicron_BA2";
            LineageList[LineageMap["BA.4"]].DMSGroup = "Omicron_BA2";
            LineageList[LineageMap["BA.5"]].DMSGroup = "Omicron_BA2";
            for (i = 1; i < LineageList.Count; i++)
            {
                if (LineageList[i].FatherInd != -1)
                {
                    if(LineageList[i].VOC=="Other")
                        LineageList[i].VOC = LineageList[LineageList[i].FatherInd].VOC;
                    if (LineageList[i].DMSGroup == "Wuhan-Hu-1_v1")
                        LineageList[i].DMSGroup = LineageList[LineageList[i].FatherInd].DMSGroup;
                    for (j = 0; j < LineageList[i].MutationList.Count(); j++)
                        if (!LineageList[LineageList[i].FatherInd].MutationList.Contains(LineageList[i].MutationList[j]))
                        { 
                            LineageList[i].NewMutList.Add(LineageList[i].MutationList[j]);
                            LineageList[i].NewMutString += LineageList[i].MutationList[j];
                        }
                    for (j = 0; j < LineageList[LineageList[i].FatherInd].MutationList.Count(); j++)
                        if (!LineageList[i].MutationList.Contains(LineageList[LineageList[i].FatherInd].MutationList[j]))
                        { 
                            LineageList[i].BackMutList.Add(LineageList[LineageList[i].FatherInd].MutationList[j]);
                            LineageList[i].BackMutString += LineageList[LineageList[i].FatherInd].MutationList[j];
                        }
                }
            }
            return;
        }
        static void ReadinDMS()//读入DMS 结合力与表达与逃逸
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/Data/final_variant_scores.csv");
            string AA20 = "ACDEFGHIKLMNPQRSTVWY";
            DMSScore.Add("Omicron_BA1", new DMS_bindexpr());
            DMSScore.Add("Omicron_BA2", new DMS_bindexpr());
            DMSScore.Add("Wuhan-Hu-1_v1", new DMS_bindexpr());
            DMSScore.Add("Wuhan-Hu-1_v2", new DMS_bindexpr());
            DMSScore.Add("Beta", new DMS_bindexpr());
            DMSScore.Add("Eta", new DMS_bindexpr());
            DMSScore.Add("Alpha", new DMS_bindexpr());
            DMSScore.Add("Delta", new DMS_bindexpr());
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split(',');
                AAchangeDMS newdms = new AAchangeDMS();
                if (line1[6] != "NA") newdms.deltaBind = Convert.ToDouble(line1[6]);
                if (line1[13] != "NA") newdms.deltaExpr = Convert.ToDouble(line1[13]);
                DMSScore[line1[0]].MutScore.Add(line1[4].Substring(1, 4), newdms);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(workfold + "/Data/EscapeScore_PKU_WT.single.txt");
            EscapeScore.Add("WT", new AAchangeEscape());
            line = read.ReadLine();line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                EscapeScore["WT"].MutEscape.Add(line1[0], Convert.ToDouble(line1[1]));
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(workfold + "/Data/EscapeScore_PKU_BA.1.single.txt");
            EscapeScore.Add("BA.1", new AAchangeEscape());
            line = read.ReadLine(); line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                EscapeScore["BA.1"].MutEscape.Add(line1[0], Convert.ToDouble(line1[1]));
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void CalculateLineageScore()//计算每个节点的表达结合力，和免疫逃逸力
        {
            int i, j, k;
            LineageList[0].DMS_Expr_Bind = 0;
            LineageList[0].Escape_normalized = 0;
            for(i=1;i<LineageList.Count();i++)
            {
                if (LineageList[i].FatherInd == -1) continue;
                LineageList[i].DMS_Expr_Bind = LineageList[LineageList[i].FatherInd].DMS_Expr_Bind;
                LineageList[i].Escape_normalized = LineageList[LineageList[i].FatherInd].Escape_normalized;
                if(LineageList[i].NewMutList.Count()!=0)//有新突变
                {
                    for(j=0;j< LineageList[i].NewMutList.Count();j++)
                    {
                        LineageList[i].DMS_Expr_Bind += DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].NewMutList[j]].deltaBind;
                        LineageList[i].DMS_Expr_Bind += DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].NewMutList[j]].deltaExpr;
                        LineageList[i].Delta_ACEbind += DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].NewMutList[j]].deltaBind;
                        LineageList[i].Delta_ACEbind += DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].NewMutList[j]].deltaExpr;
                        
                        LineageList[i].Escape_normalized += EscapeScore[ImmuneBackground].MutEscape[LineageList[i].NewMutList[j]];
                        LineageList[i].Delta_ESCscore += EscapeScore[ImmuneBackground].MutEscape[LineageList[i].NewMutList[j]];
                    }
                }
                if (LineageList[i].BackMutList.Count() != 0)//有回复突变
                {
                    for (j = 0; j < LineageList[i].BackMutList.Count(); j++)
                    {
                        LineageList[i].DMS_Expr_Bind -= DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].BackMutList[j]].deltaBind;
                        LineageList[i].DMS_Expr_Bind -= DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].BackMutList[j]].deltaExpr;
                        LineageList[i].Delta_ACEbind -= DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].BackMutList[j]].deltaBind;
                        LineageList[i].Delta_ACEbind -= DMSScore[LineageList[LineageList[i].FatherInd].DMSGroup].MutScore[LineageList[i].BackMutList[j]].deltaExpr;

                        LineageList[i].Escape_normalized -= EscapeScore[ImmuneBackground].MutEscape[LineageList[i].BackMutList[j]];
                        LineageList[i].Delta_ESCscore -= EscapeScore[ImmuneBackground].MutEscape[LineageList[i].BackMutList[j]];
                    }
                }
            }
        }
        static void WriteLineages()//输出所有的节点信息
        {
            StreamWriter write = new StreamWriter(workfold + "/LineageInfo.tsv");
            int i, j, k;
            write.WriteLine("Lineage\tCollectionDate5P\tNewMutCount\tNewMut70P\tMut70PCount\tMut70P\tDMS_Infect\tDMS_Escape\tCountryNumber\tSeqNumber\tVOC\tNodeInd\tFatherInd\tChildNum");
            for(i=0;i<LineageList.Count;i++)
            {
                if (LineageList[i].SeqNumber > 0)
                {
                    string output = LineageList[i].Name + "\t";
                    output += LineageList[i].CollectionDate5P + "\t";
                    output += LineageList[i].NewMutList.Count() + "\t";
                    string newmut = "";
                    for (j = 0; j < LineageList[i].NewMutList.Count(); j++)
                        newmut += LineageList[i].NewMutList[j] + ",";
                    if (newmut != "") newmut = newmut.Substring(0, newmut.Length - 1);
                    output += newmut + "\t";
                    output += LineageList[i].MutationList.Count() + "\t";
                    string mut = "";
                    for (j = 0; j < LineageList[i].MutationList.Count(); j++)
                        mut += LineageList[i].MutationList[j] + ",";
                    if (mut != "") mut = mut.Substring(0, mut.Length - 1);
                    output += mut + "\t";
                    output += Convert.ToString(LineageList[i].DMS_Expr_Bind) + "\t";
                    output += Convert.ToString(LineageList[i].Escape_normalized) + "\t";
                    output += Convert.ToString(LineageList[i].CountryNumber) + "\t";
                    output += Convert.ToString(LineageList[i].SeqNumber) + "\t";
                    output += LineageList[i].VOC + "\t";
                    output += Convert.ToString(LineageList[i].NodeInd) + "\t";
                    output += Convert.ToString(LineageList[i].FatherInd) + "\t";
                    output += Convert.ToString(LineageList[i].ChildList.Count);
                    string child = "";
                    for (j = 0; j < LineageList[i].ChildList.Count(); j++)
                        child += LineageList[i].ChildList[j] + ",";
                    if (child != "") child = child.Substring(0, child.Length - 1);
                    //output += child + "\t";
                    
                    write.WriteLine(output);
                }
            }
            write.Close();
            return;
        }
        static void PickLineages()//挑选画图用的Lineage
        {
            int i, j, k;
            for(i=LineageList.Count-1;i>=0;i--)
            {
                if (LineageList[i].NewMutList.Count > 0 && LineageList[i].FatherInd != -1)
                    //if (LineageList[i].SeqNumber >= 1000 || LineageList[i].CountryNumber >= 20 || LineageList[i].VOC == "Omicron")
                        LineageList[i].pick4plot = true;
                if (LineageList[i].pick4plot == true && LineageList[LineageList[i].FatherInd].FatherInd != -1)
                    LineageList[LineageList[i].FatherInd].pick4plot = true;
            }
            StreamWriter write = new StreamWriter(workfold + "/LineagePlotData_" + ImmuneBackground + ".txt");
            write.WriteLine("Lineage\tVOC\tDate5P\tCountryNumber\tSeqNumber\tDMS_Infect\tDMS_Escape\tBranchGroup\tMidpoint\tNewMut\tBackMut\tΔESC\tΔACE");
            for(i=1;i<LineageList.Count;i++)
            {
                if (LineageList[i].pick4plot == true)
                {
                    string output = LineageList[i].Name + "\t";
                    output += LineageList[i].VOC + "\t";
                    output += Convert.ToString(LineageList[i].CollectionDate5P) + "\t";
                    output += Convert.ToString(LineageList[i].CountryNumber) + "\t";
                    output += Convert.ToString(LineageList[i].SeqNumber) + "\t";
                    output += Convert.ToString(LineageList[i].DMS_Expr_Bind) + "\t";
                    output += Convert.ToString(LineageList[i].Escape_normalized) + "\t";
                    output += Convert.ToString(LineageList[i].Name) + "\t";
                    if (LineageList[i].NewMutList.Count > 0) output += "Lineage\t";
                    else output += "Midpoint\t";
                    output += LineageList[i].NewMutString + "\t";
                    output += LineageList[i].BackMutString + "\t";
                    output += Convert.ToString(LineageList[i].Delta_ESCscore) + "\t" + Convert.ToString(LineageList[i].Delta_ACEbind);
                    write.WriteLine(output);

                    int fatherInd = LineageList[i].FatherInd;
                    output = LineageList[fatherInd].Name + "\t";
                    output += LineageList[fatherInd].VOC + "\t";
                    output += Convert.ToString(LineageList[fatherInd].CollectionDate5P) + "\t";
                    output += Convert.ToString(LineageList[fatherInd].CountryNumber) + "\t";
                    output += Convert.ToString(LineageList[fatherInd].SeqNumber) + "\t";
                    output += Convert.ToString(LineageList[fatherInd].DMS_Expr_Bind) + "\t";
                    output += Convert.ToString(LineageList[fatherInd].Escape_normalized) + "\t";
                    output += Convert.ToString(LineageList[i].Name) + "\t";
                    if (LineageList[fatherInd].NewMutList.Count > 0) output += "Lineage\t";
                    else output += "Midpoint\t";
                    output += LineageList[fatherInd].NewMutString + "\t";
                    output += LineageList[fatherInd].BackMutString + "\t";
                    output += Convert.ToString(LineageList[fatherInd].Delta_ESCscore) + "\t" + Convert.ToString(LineageList[fatherInd].Delta_ACEbind);
                    write.WriteLine(output);
                }
            }
            write.Close();
        }
        static void Main(string[] args)
        {
            ReadinRelations();//读入Lineage节点间的关系
            ReadinDateAAmut();//读入节点的采样时间和突变，计算与父节点的差异
            ReadinDMS();//读入DMS 结合力与表达与逃逸
            CalculateLineageScore();//计算每个节点的表达结合力，和免疫逃逸力
            WriteLineages();//输出所有的节点信息
            PickLineages();//挑选画图用的Lineage
        }
    }
}
