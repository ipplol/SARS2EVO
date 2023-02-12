using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 不同进化枝突变针对12类抗体的逃逸分布_只算逃不逃逸
{
    public class Mutation
    {
        public string MutName;
        public List<double> EscScores = new List<double>();
        public List<int> Esc12TF = new List<int>();//是否12类抗体的逃逸突变，1是，0否
        public bool NotEscMut = true;//是否为非逃逸突变
    }
    public class Clade
    {
        public string CladeName;
        public int CladeRBDTotalEvent = 0;
        public int TotalEscapeMutation = 0;
        public int TotalPositiveDMSMutation = 0;
        public int TotalESC_DMSMutation = 0;
        public List<int> Esc12Group = new List<int>();
    }
    class Program
    {
        static string workfold = "***/Usher1128";
        static Dictionary<string, Mutation> MutationList = new Dictionary<string, Mutation>();
        static Dictionary<string, double> DMS_WT = new Dictionary<string, double>();
        static Dictionary<string, double> DMS_BA1 = new Dictionary<string, double>();
        static Dictionary<string, double> DMS_BA2 = new Dictionary<string, double>();
        static List<Clade> CladeList = new List<Clade>();
        static void WriteResult()
        {
            StreamWriter write = new StreamWriter(workfold + "/Group12/CladeRBDMut12Group.txt");
            int i, j, k;
            write.WriteLine("Clade\tA\tB\tC\tD1\tD2\tE1\tE2.1\tE2.2\tE3\tF1\tF2\tF3");
            for(i=0;i<CladeList.Count;i++)
            {
                string output = CladeList[i].CladeName;
                for (j = 0; j < 12; j++)
                    output += "\t" + Convert.ToString((double)CladeList[i].Esc12Group[j] / CladeList[i].CladeRBDTotalEvent);
                write.WriteLine(output);
            }
            write.Close();

            write = new StreamWriter(workfold + "/Group12/CladeRBDEscMutProp.txt");
            write.WriteLine("Clade\tTotalEvent\tEscMut\tDMSPositiveMut\tEsc&DMSMut");
            for (i = 0; i < CladeList.Count; i++)
            {

                write.WriteLine(CladeList[i].CladeName + "\t" + Convert.ToString(CladeList[i].CladeRBDTotalEvent) + "\t" + Convert.ToString(CladeList[i].TotalEscapeMutation) + "\t" + Convert.ToString(CladeList[i].TotalPositiveDMSMutation) + "\t" + Convert.ToString(CladeList[i].TotalESC_DMSMutation));
            }
            write.Close();

            return;
        }
        static void EscapeDistribution()
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/CladeRBDMutRate.txt");
            //StreamReader read = new StreamReader(workfold + "/Group12/random.txt");
            //Clade\tCladeMutNum\tMut\tESC\tEventCount
            //需要去重复
            string line = read.ReadLine();
            line = read.ReadLine();
            string tmp = "#TBD#";
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (line1[0] != tmp)
                {
                    Clade newc = new Clade();
                    newc.CladeName = line1[0];
                    newc.CladeRBDTotalEvent = Convert.ToInt32(line1[1]);
                    for (i = 0; i < 12; i++) newc.Esc12Group.Add(0);
                    CladeList.Add(newc);
                    tmp = line1[0];
                }
                int indE = 0, indD = 0;
                int n = 0;
                for(i=0;i<12;i++)
                {
                    CladeList[CladeList.Count - 1].Esc12Group[i] += MutationList[line1[2]].Esc12TF[i] * Convert.ToInt32(line1[4]);
                    n += MutationList[line1[2]].Esc12TF[i];
                }
                if (n != 0)
                { CladeList[CladeList.Count - 1].TotalEscapeMutation += Convert.ToInt32(line1[4]); indE = 1; }

                if (line1[0].Contains("BA.1"))
                {
                    if (DMS_BA1[line1[2]] > 0)
                    { CladeList[CladeList.Count - 1].TotalPositiveDMSMutation += Convert.ToInt32(line1[4]); indD = 1; }
                }
                else if (line1[0].Contains("BA"))
                {
                    if (DMS_BA2[line1[2]] > 0)
                    { CladeList[CladeList.Count - 1].TotalPositiveDMSMutation += Convert.ToInt32(line1[4]); indD = 1; }
                }
                else if (DMS_WT[line1[2]] > 0)
                { CladeList[CladeList.Count - 1].TotalPositiveDMSMutation += Convert.ToInt32(line1[4]); indD = 1; }

                if(indE+indD==2)
                {
                    CladeList[CladeList.Count - 1].TotalESC_DMSMutation += Convert.ToInt32(line1[4]);
                }
                line = read.ReadLine();
            }
            return;
        }
        static void ReadinAndCalAverage3()
        {
            int mutationCount = 0;
            List<double> ScoreSum = new List<double>();
            StreamReader read = new StreamReader(workfold + "/Data/EscapeScore_PKU_WT.12EscOnly.txt");
            StreamWriter write = new StreamWriter(workfold + "/Data/MutationEsc12.txt");
            write.WriteLine("Mutation\tA\tB\tC\tD1\tD2\tE1\tE2.1\tE2.2\tE3\tF1\tF2\tF3");
            int i, j, k;
            for (i = 1; i <= 12; i++) ScoreSum.Add(0);
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                Mutation newa = new Mutation();
                newa.MutName = line1[0];
                for(i=1;i<line1.Length;i++)
                {
                    newa.EscScores.Add(Convert.ToDouble(line1[i]));
                    ScoreSum[i - 1] += Convert.ToDouble(line1[i]);
                }
                MutationList.Add(newa.MutName,newa);
                mutationCount++;
                line = read.ReadLine();
            }
            for (i = 0; i < 12; i++)
                ScoreSum[i] = ScoreSum[i] / mutationCount * 3;
            //阈值为3倍平均数

            foreach(string val in MutationList.Keys)
            {
                string output = val;
                for(j=0;j<12;j++)
                {
                    if (MutationList[val].EscScores[j] >= ScoreSum[j])
                    {
                        MutationList[val].Esc12TF.Add(1);
                        MutationList[val].NotEscMut = false;
                        output += "\t1";
                    }
                    else
                    {
                        MutationList[val].Esc12TF.Add(0);
                        output += "\t0";
                    }
                }
                write.WriteLine(output);
            }
            read.Close();
            write.Close();
            return;
        }
        static void ReadinDMS()
        {
            StreamReader read = new StreamReader(workfold + "/Group12/DMS.txt");
            int i, j, k;
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (line1[3] != "NA") DMS_WT.Add(line1[0], Convert.ToDouble(line1[3]));
                else DMS_WT.Add(line1[0], 0);
                if (line1[4] != "NA") DMS_BA1.Add(line1[0], Convert.ToDouble(line1[4]));
                else DMS_BA1.Add(line1[0], 0);
                if (line1[5] != "NA") DMS_BA2.Add(line1[0], Convert.ToDouble(line1[5]));
                else DMS_BA2.Add(line1[0], 0);
                line = read.ReadLine();
            }
            read.Close();
        }
        static void Main(string[] args)
        {
            ReadinDMS();
            ReadinAndCalAverage3();//读入突变逃逸得分，计算突变是否逃逸12类抗体
            EscapeDistribution();//统计分布
            WriteResult();
            return;
        }
    }
}
