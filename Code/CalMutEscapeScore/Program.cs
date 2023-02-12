using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Diagnostics;

namespace 适应性突变估计_北大算法_单位点
{
    public class Antibody
    {
        public string Name;
        public Dictionary<string, double> NeutralDic = new Dictionary<string, double>();//variant, IC50
        public Dictionary<string, double> Escape = new Dictionary<string, double>();//pos mut, escape score
        public double MaxEscape = 0;
        public string source;//wt BA.1
        public string group;//12 group
    }
    class Program
    {
        static Dictionary<string, double> Dic_DMS_Exp = new Dictionary<string, double>();//e^(bind+expr)
        static Dictionary<string, int> Dic_Mut_Codon = new Dictionary<string, int>();//突变是否可以通过单个碱基突变到达
        static Dictionary<string, Antibody> Dic_Antibody = new Dictionary<string, Antibody>();//name, antibody
        static List<string> variantList = new List<string>();//测了中和能力的变异株
        static List<string> groupList = new List<string>();//12类抗体
        static string AA20 = "ACDEFGHIKLMNPQRSTVWY";
        static string RBD331_531_AA;
        static string RBD331_531_Nuc;
        List<string> AAMutationList = new List<string>();
        static string Workfold = "***/Usher1128";
        static void Readin()//读入数据
        {
            int i, j, k;

            //读入DMS
            StreamReader read = new StreamReader(Workfold + "/Data/single_mut_effects.csv");
            string line;
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split(',');
                string AAmut = line1[1] + line1[3].Substring(1, 1);
                double bindExpr = 0;
                if (line1[8] != "NA") bindExpr += Convert.ToDouble(line1[8]);
                if (line1[11] != "NA") bindExpr += Convert.ToDouble(line1[11]);
                if (bindExpr >= 0)
                    Dic_DMS_Exp.Add(AAmut, Math.Exp(bindExpr));//【取这俩加和的e幂】
                else
                    Dic_DMS_Exp.Add(AAmut, 1);
                line = read.ReadLine();
            }
            read.Close();


            //读入抗体中和
            read = new StreamReader(Workfold + "/Data/NeutralWTBA125_Cross.txt");
            line = read.ReadLine();
            string[] linev = line.Split('\t');
            for (i = 4; i < 9; i++)
                variantList.Add(linev[i]);
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Antibody newa = new Antibody();
                newa.Name = line1[0];
                double IC50 = 1;
                for (i = 4; i < 9; i++)
                    if (line1[i] != "NA" && line1[i] != "--")
                    {
                        IC50 = Convert.ToDouble(line1[i]);
                        if (IC50 > 1) IC50 = 1;
                        if (IC50 < 0.0005) IC50 = 0.0005;
                        newa.NeutralDic.Add(variantList[i - 4], IC50);
                    }
                    else
                        newa.NeutralDic.Add(variantList[i - 4], 1);
                newa.source = line1[3];
                newa.group = line1[2];
                if (!groupList.Contains(newa.group)) groupList.Add(newa.group);
                Dic_Antibody.Add(line1[0], newa);
                line = read.ReadLine();
            }
            read.Close();

            //读入逃逸得分
            read = new StreamReader(Workfold + "/Data/use_abs_res.csv");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split(',');
                string AAmut = line1[1] + line1[2];
                double Escape = Convert.ToDouble(line1[3]);
                if (Dic_Antibody.ContainsKey(line1[0]))
                {
                    Dic_Antibody[line1[0]].Escape.Add(AAmut, Escape);
                    if (Escape > Dic_Antibody[line1[0]].MaxEscape) Dic_Antibody[line1[0]].MaxEscape = Escape;
                }
                line = read.ReadLine();
            }
            read.Close();

            //读入RBD序列，计算密码子
            read = new StreamReader(Workfold + "/Data/RBD_331_531.txt");
            RBD331_531_AA = read.ReadLine();
            RBD331_531_Nuc = read.ReadLine();
            read.Close();
            Dictionary<string, string> CodonMap = new Dictionary<string, string>();
            read = new StreamReader(Workfold + "/Data/mimazi.txt");
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                CodonMap.Add(line1[0], line1[1]);
                line = read.ReadLine();
            }
            read.Close();
            string ATCG = "ATCG";
            for (i = 0; i < 201; i++)
            {
                string refCodon = RBD331_531_Nuc.Substring(i * 3, 3);
                for (j = 0; j < 3; j++)
                {
                    for (k = 0; k < 4; k++)
                    {
                        char[] mutCodon = refCodon.ToCharArray();
                        mutCodon[j] = ATCG[k];
                        string tmps = Convert.ToString(i + 331) + CodonMap[new string(mutCodon)];
                        if (!Dic_Mut_Codon.ContainsKey(tmps))
                            Dic_Mut_Codon.Add(tmps, 1);
                    }
                }
            }

            StreamWriter writeDMS = new StreamWriter(Workfold + "/Data/bindExpr.txt");
            for (i = 331; i <= 531; i++)
            {
                double sum = 0;
                int count = 0;
                for (j = 0; j < AA20.Length; j++)
                {
                    string AAmut = Convert.ToString(i) + Convert.ToString(AA20[j]);
                    double S_dms = 1;
                    if (Dic_DMS_Exp.ContainsKey(AAmut))
                        S_dms = Dic_DMS_Exp[AAmut];
                    if (Dic_Mut_Codon.ContainsKey(AAmut))
                    {
                        sum += S_dms;
                        count++;
                    }
                }
                writeDMS.WriteLine(Convert.ToString(i) + "\t" + Convert.ToString(sum / count));
            }
            writeDMS.Close();
        }
        static List<double> EscapeScoreCalculator(string Variant, string Antibody)
        {
            int i, j, k;
            List<double> siteScoreList = new List<double>();
            for (i = 331; i <= 531; i++)
            {
                for (j = 0; j < AA20.Length; j++)
                {
                    string AAmut = Convert.ToString(i) + Convert.ToString(AA20[j]);

                    double S_neutral = 0;
                    if (Dic_Antibody[Antibody].NeutralDic[Variant] != 0)
                        S_neutral = -1 * Math.Log10(Dic_Antibody[Antibody].NeutralDic[Variant]);

                    double S_escape = 0;
                    if (Dic_Antibody[Antibody].MaxEscape != 0 && Dic_Antibody[Antibody].Escape.ContainsKey(AAmut))
                        S_escape = Dic_Antibody[Antibody].Escape[AAmut] / Dic_Antibody[Antibody].MaxEscape;

                    double S = S_neutral * S_escape;
                    //double S = S_escape;

                    siteScoreList.Add(S);
                }
            }
            return siteScoreList;
        }
        static void AllAntibodyEscapeScoreCalculator()
        {
            string variantTEsc = "BA.1";
            StreamWriter write = new StreamWriter(Workfold + "/Data/EscapeScore_PKU_" + variantTEsc + ".single.txt");
            int i, j, k, tmpi,tmpj;
            List<double> TotalEscapeScore = new List<double>();
            List<string> AAmut = new List<string>();
            //先算所有抗体的得分
            for (i = 331; i <= 531; i++)
                for (j = 0; j < AA20.Length; j++)
                {
                    TotalEscapeScore.Add(0);
                    AAmut.Add(Convert.ToString(i) + AA20[j]);
                }
            foreach (string val in Dic_Antibody.Keys)
            {
                if (Dic_Antibody[val].source == "WT convalescents") //variantTEsc
                {
                    List<double> tmpl = new List<double>(EscapeScoreCalculator(variantTEsc, val));
                    for (i = 0; i < tmpl.Count; i++)
                        TotalEscapeScore[i] += tmpl[i];
                    if (Dic_Antibody[val].source == "BA.1 convalescents")
                        Console.WriteLine("?");
                }
            }
            write.WriteLine("Pos\t" + variantTEsc);
            for (i = 0; i < TotalEscapeScore.Count; i++)
            {
                string output = AAmut[i] + "\t";
                output += Convert.ToString(TotalEscapeScore[i]);
                write.WriteLine(output);
            }
            write.Close();

            //再算每一类的
            groupList.Sort();
            write = new StreamWriter(Workfold + "/Data/EscapeScore_PKU_" + variantTEsc + ".12.txt");
            string outputline = "Mut";
            for (i = 0; i < groupList.Count; i++) outputline += "\t" + groupList[i];
            write.WriteLine(outputline);
            List<string> outgroup = new List<string>();
            for (i = 331; i <= 531; i++)
                for (j = 0; j < AA20.Length; j++)
                {
                    outgroup.Add(Convert.ToString(i) + AA20[j]);
                }
            for (k = 0; k < groupList.Count; k++)
            {
                TotalEscapeScore = new List<double>();
                for (tmpi = 331; tmpi <= 531; tmpi++)
                    for (tmpj = 0; tmpj < AA20.Length; tmpj++)
                        TotalEscapeScore.Add(0);
                foreach (string val in Dic_Antibody.Keys)
                {
                    if (Dic_Antibody[val].group == groupList[k] && Dic_Antibody[val].source == variantTEsc)
                    {
                        List<double> tmpl = new List<double>(EscapeScoreCalculator(variantTEsc, val));
                        for (tmpi = 0; tmpi < tmpl.Count; tmpi++)
                            if(tmpl[tmpi]>0)
                            TotalEscapeScore[tmpi] += tmpl[tmpi];
                    }
                }
                for (i = 0; i < TotalEscapeScore.Count; i++)
                    outgroup[i] += "\t" + Convert.ToString(TotalEscapeScore[i]);
            }
            for (i = 0; i < outgroup.Count; i++)
                write.WriteLine(outgroup[i]);
                write.Close();
        }
        static void Main(string[] args)
        {
            Readin();//读入数据
            AllAntibodyEscapeScoreCalculator();
            return;
        }
    }
}
