using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 不同进化枝在出现后的突变率差异
{
    class SpikeSite //RBD NTD Other
    {
        public int SitePos;//氨基酸位置
        public int EventCount = 0;//该突变发生次数
        public int EventWithThisMut = 0;//以该突变为背景的突变事件数
        public double MutRate = 0;
        public double MutProp = 0;
        public double Prop = 0;
    }
    class MutationAnno //突变注释
    {
        public int Position;
        public string Titv;
        public string Gene;
        public string NSS;
        public string AminoChange;
        public string Context3;
        public string Context5;
        public string CodonPos;
    }
    class GeneWindow
    {
        public string GeneName;
        public Dictionary<string, SpikeSite> Gene_SiteList = new Dictionary<string, SpikeSite>();
        public int CladeMutEventNum = 0;//该时间段发生在该进化枝下的事件数
    }
    class TimeWindow
    {
        public string WindowName;//clade
        public int DateStart;
        public int DateEnd;
        public Dictionary<string, SpikeSite> SpikeSiteList = new Dictionary<string, SpikeSite>();
        public int TotalMutEventNum = 0;//该时间段总共突变事件数
        public int CladeMutEventNum = 0;//该时间段发生在该进化枝下的事件数
        public int CladeRBDMutEventNum = 0;
        public int NS_MutEventNum = 0;
        public int S_MutEventNum = 0;
        public int backmutNum = 0;

        public Dictionary<string, GeneWindow> GeneMutDic = new Dictionary<string, GeneWindow>();//每一个基因的突变情况
        public List<string> topN = new List<string>();
    }
    class Program
    {
        static string location = "global";
        //static string location = "United Kingdom END England";
        static string workfold = "***/Usher1128";
        static Dictionary<string, MutationAnno> MutAnnoMap = new Dictionary<string, MutationAnno>();
        static List<TimeWindow> timeWindowList = new List<TimeWindow>();
        static string AA20 = "ACDEFGHIKLMNPQRSTVWY";
        static string RBD331_531;
        static string RBDnuc;
        static Dictionary<string, string> CodonMap = new Dictionary<string, string>();//密码子
        static Dictionary<string, double> EscapeScore_WT_12 = new Dictionary<string, double>();
        static Dictionary<string, double> EscapeScore_WT_Total = new Dictionary<string, double>();
        static Dictionary<string, double> EscapeScore_BA1_12 = new Dictionary<string, double>();
        static Dictionary<string, double> EscapeScore_BA1_Total = new Dictionary<string, double>();
        static Dictionary<string, bool> IfMutationEscape = new Dictionary<string, bool>();//这个突变是否是逃逸突变
        static List<double> EnrichmentScoreCurve = new List<double>();
        static string[] Group12 = "A,B,C,D1,D2,E1,E2.1,E2.2,E3,F1,F2,F3".Split(',');
        
        static void ReadMutAnno()//读入突变注释
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/Data/AllMutationAnnomutresult.tsv");
            string line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                string[] nuc = line1[4].Split('/');
                string key = nuc[0] + line1[2] + nuc[1];
                MutationAnno ma = new MutationAnno();
                ma.Position = Convert.ToInt32(line1[2]);
                ma.Titv = line1[1];
                ma.Gene = line1[5];
                ma.NSS = line1[7];
                ma.AminoChange = line1[6];
                ma.Context3 = line1[8];
                ma.Context5 = line1[9];
                ma.CodonPos = line1[10];
                MutAnnoMap.Add(key, ma);
                line = read.ReadLine();
            }
            read.Close();
            read = new StreamReader(workfold + "/Data/RBD_331_531.txt");
            RBD331_531 = read.ReadLine();
            RBDnuc = read.ReadLine();
            read.Close();

            //读入密码子
            read = new StreamReader(workfold + "/Data/mimazi.txt");
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                CodonMap.Add(line1[0], line1[1]);
                line = read.ReadLine();
            }
            read.Close();

            return;
        }
        static string CombineMutation(List<string> backgroundmut, string newmut, int AApos)
        {
            char[] refnuc = RBDnuc.ToCharArray();
            int i, j, k;
            for(i=0;i<backgroundmut.Count;i++)
            {
                int pos = Convert.ToInt32(backgroundmut[i].Substring(1, backgroundmut[i].Length - 2));
                if (pos >= 22553 && pos <= 23155)
                    refnuc[pos - 22553] = backgroundmut[i][backgroundmut[i].Length - 1];
            }
            j = Convert.ToInt32(newmut.Substring(1, newmut.Length - 2));
            if (j >= 22553 && j <= 23155)
                refnuc[j - 22553] = newmut[newmut.Length - 1];
            return Convert.ToString(AApos) + CodonMap["" + refnuc[(AApos - 331) * 3]+refnuc[(AApos - 331) * 3+1]+refnuc[(AApos - 331) * 3+2]];
        }
        static void ReadMutEvent()//读入突变事件
        {
            StreamWriter writetmp = new StreamWriter(workfold + "/BackMut.txt");
            int i, j, k, tmpi =0;

            //读入突变事件
            StreamReader read = new StreamReader(workfold + "/Data/MAT1128.json.mutevent");
            read.ReadLine();
            string line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                if (location == "global" || location.ToUpper().Contains(line1[6].ToUpper()))
                {
                    if (MutAnnoMap[line1[3]].NSS == "S")
                    {
                        //定位突变枝
                        string[] lineage1 = line1[4].Split('.');
                        for (k = 0; k < timeWindowList.Count; k++)
                        {
                            //if (line1[8].Contains(timeWindowList[k].WindowName))//先找能不能对上的
                            //    break;
                            if (line1[4] == timeWindowList[k].WindowName)
                                break;
                        }
                        if (k >= timeWindowList.Count)
                        {
                            string[] tmp1 = line1[4].Split('.');//再找Lineage的上一级
                            string tmps = "";
                            for (j = 0; j < tmp1.Length - 1; j++) tmps += tmp1[j] + ".";
                            if (tmps != "") tmps = tmps.Substring(0, tmps.Length - 1);
                            for (k = 0; k < timeWindowList.Count; k++)
                                if (tmps == timeWindowList[k].WindowName)
                                    break;
                        }
                        int AApos = -1;
                        int startp = 0, endp = MutAnnoMap[line1[3]].AminoChange.Length - 1;
                        if (MutAnnoMap[line1[3]].AminoChange != "intergenic" && MutAnnoMap[line1[3]].AminoChange != "Overlap")
                        {
                            while (!(MutAnnoMap[line1[3]].AminoChange[startp] <= '9' && MutAnnoMap[line1[3]].AminoChange[startp] >= '0')) startp++;
                            while (!(MutAnnoMap[line1[3]].AminoChange[endp] <= '9' && MutAnnoMap[line1[3]].AminoChange[endp] >= '0')) endp--;
                            AApos = Convert.ToInt32(MutAnnoMap[line1[3]].AminoChange.Substring(startp, endp - startp + 1));
                        }
                        if (k < timeWindowList.Count) timeWindowList[k].S_MutEventNum++;
                    }


                    if (MutAnnoMap[line1[3]].NSS == "NS")
                    {
                        if (line1[5].Length < 7)
                        {
                            line = read.ReadLine();
                            continue;
                        }

                        //定位采集时间
                        int month = Convert.ToInt32(line1[5].Substring(0, 4)) * 100 + Convert.ToInt32(line1[5].Substring(5, 2));

                        //定位突变枝
                        /*if ("19A19B20A20B20C20D20E (EU1)20F20G".Contains(line1[8]))
                            line1[8] = "WT";*/
                        string[] lineage1 = line1[4].Split('.');
                        for (k = 0; k < timeWindowList.Count; k++)
                        {
                            //if (line1[8].Contains(timeWindowList[k].WindowName))//先找能不能对上的
                            //    break;
                            if (line1[4]  == timeWindowList[k].WindowName)
                                break;
                        }
                        if (k >= timeWindowList.Count)
                        {
                            string[] tmp1 = line1[4].Split('.');//再找Lineage的上一级
                            string tmps = "";
                            for (j = 0; j < tmp1.Length - 1; j++) tmps += tmp1[j] + ".";
                            if (tmps != "") tmps = tmps.Substring(0, tmps.Length - 1);
                            for (k = 0; k < timeWindowList.Count; k++)
                                if (tmps == timeWindowList[k].WindowName)
                                    break;
                        }

                        if (k < timeWindowList.Count && month >= timeWindowList[k].DateStart && month <= timeWindowList[k].DateEnd )
                        {
                            timeWindowList[k].NS_MutEventNum++;
                            int AApos = -1;
                            int startp = 0, endp = MutAnnoMap[line1[3]].AminoChange.Length - 1;
                            if (MutAnnoMap[line1[3]].AminoChange != "intergenic" && MutAnnoMap[line1[3]].AminoChange != "Overlap")
                            {
                                while (!(MutAnnoMap[line1[3]].AminoChange[startp] <= '9' && MutAnnoMap[line1[3]].AminoChange[startp] >= '0')) startp++;
                                while (!(MutAnnoMap[line1[3]].AminoChange[endp] <= '9' && MutAnnoMap[line1[3]].AminoChange[endp] >= '0')) endp--;
                                AApos = Convert.ToInt32(MutAnnoMap[line1[3]].AminoChange.Substring(startp, endp - startp + 1));
                            }

                            string[] nucmut = line1[2].Split(',');

                            //测试这个突变是不是回复突变
                            string backmuttest = "" + line1[3][line1[3].Length - 1] + line1[3].Substring(1, line1[3].Length - 2) + line1[3][0];
                            if (nucmut.Contains(backmuttest))
                            {
                                writetmp.WriteLine(timeWindowList[k].WindowName + "\t" + MutAnnoMap[line1[3]].AminoChange);
                                timeWindowList[k].backmutNum++;
                                line = read.ReadLine(); continue;
                            }

                            if (MutAnnoMap[line1[3]].Gene == "S")//只关注Spike上的突变
                            {
                                timeWindowList[k].CladeMutEventNum++;//该突变分支下发生的突变+1

                                //转换背景突变为AA突变
                                List<string> aamut = new List<string>();
                                List<string> aamut_pos = new List<string>();
                                for (i = 0; i < nucmut.Count() - 1; i++)
                                    if (nucmut[i] != "" && MutAnnoMap[nucmut[i]].NSS == "NS" && MutAnnoMap[nucmut[i]].Gene == "S")
                                    {
                                        if (!aamut.Contains(MutAnnoMap[nucmut[i]].AminoChange.Substring(1, MutAnnoMap[nucmut[i]].AminoChange.Length - 1)))
                                        {
                                            aamut.Add(MutAnnoMap[nucmut[i]].AminoChange.Substring(1, MutAnnoMap[nucmut[i]].AminoChange.Length - 1));
                                            aamut_pos.Add(MutAnnoMap[nucmut[i]].AminoChange.Substring(1, MutAnnoMap[nucmut[i]].AminoChange.Length - 2));
                                        }
                                    }
                                for (i = 0; i < aamut.Count; i++)
                                    if (timeWindowList[k].SpikeSiteList.ContainsKey(aamut[i]))
                                        timeWindowList[k].SpikeSiteList[aamut[i]].EventWithThisMut++;
                                //ref氨基酸
                                for (i = 331; i <= 531; i++)
                                    if (!aamut_pos.Contains(Convert.ToString(i)))
                                        if (timeWindowList[k].SpikeSiteList.ContainsKey(Convert.ToString(i) + RBD331_531[i - 331]))
                                            timeWindowList[k].SpikeSiteList[Convert.ToString(i) + RBD331_531[i - 331]].EventWithThisMut++;

                                string key = MutAnnoMap[line1[3]].AminoChange.Substring(1, MutAnnoMap[line1[3]].AminoChange.Length - 1);
                                //RBD NS突变
                                if (MutAnnoMap[line1[3]].NSS == "NS" && MutAnnoMap[line1[3]].Gene == "S" && !aamut.Contains(key) && AApos >= 331 && AApos <= 531)
                                {
                                    key = CombineMutation(nucmut.ToList(), line1[3], AApos);
                                    timeWindowList[k].SpikeSiteList[key].EventCount++;
                                    timeWindowList[k].CladeRBDMutEventNum++;
                                }
                            }

                            foreach(string val in timeWindowList[k].GeneMutDic.Keys)
                            {
                                if (val.Contains(MutAnnoMap[line1[3]].Gene))
                                {
                                    //转换背景突变为AA突变
                                    List<string> aamut = new List<string>();
                                    List<string> aamut_pos = new List<string>();
                                    for (i = 0; i < nucmut.Count() - 1; i++)
                                        if (nucmut[i] != "" && MutAnnoMap[nucmut[i]].NSS == "NS" && val.Contains(MutAnnoMap[nucmut[i]].Gene))
                                        {
                                            if (!aamut.Contains(MutAnnoMap[nucmut[i]].AminoChange.Substring(1, MutAnnoMap[nucmut[i]].AminoChange.Length - 1)))
                                            {
                                                aamut.Add(MutAnnoMap[nucmut[i]].AminoChange.Substring(1, MutAnnoMap[nucmut[i]].AminoChange.Length - 1));
                                                aamut_pos.Add(MutAnnoMap[nucmut[i]].AminoChange.Substring(1, MutAnnoMap[nucmut[i]].AminoChange.Length - 2));
                                            }
                                        }
                                    for (i = 0; i < aamut.Count; i++)
                                        if (timeWindowList[k].GeneMutDic[val].Gene_SiteList.ContainsKey(aamut[i]))
                                            timeWindowList[k].GeneMutDic[val].Gene_SiteList[aamut[i]].EventWithThisMut++;
                                    timeWindowList[k].GeneMutDic[val].CladeMutEventNum++;//该突变分支下发生的突变+1
                                    string key = MutAnnoMap[line1[3]].AminoChange.Substring(1, MutAnnoMap[line1[3]].AminoChange.Length - 1);
                                    timeWindowList[k].GeneMutDic[val].Gene_SiteList[key].EventCount++;
                                }
                            }

                        }
                    }
                }
                line = read.ReadLine();
            }

            writetmp.Close();
            read.Close();
            return;
        }

        static void CalMutRateOutputForGenes(string Gene)//计算突变率，输出
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(workfold + "/GeneTopMut/CladeMutRate_" + Gene + ".txt");
            StreamWriter writep = new StreamWriter(workfold + "/GeneTopMut/CladeHeatmap_" + Gene + ".txt");
            StreamWriter writeBarplot = new StreamWriter(workfold + "/GeneTopMut/CladeBarplot_" + Gene + ".txt");
            string output = "Clade\tCladeMutNum\tMut\tEventCount";
            write.WriteLine(output);
            for (j = 0; j < timeWindowList.Count; j++)
            {
                foreach (string valk in timeWindowList[j].GeneMutDic[Gene].Gene_SiteList.Keys)
                {
                    timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].Prop = (double)timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].EventWithThisMut / timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum;
                    double mutRate = (double)timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].EventCount / timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum;
                    if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].Prop < 0.9 && timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum > 100 && mutRate >= 0.00001)
                    {
                        timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].MutRate = mutRate;
                    }
                }
            }
            string title = "Mut";
            //输出每个窗口区域前
            List<string> RBD5 = new List<string>();
            for (j = 0; j < timeWindowList.Count; j++)
            {
                title += "\t" + timeWindowList[j].WindowName;
                List<string> siteID = new List<string>();
                for (i = 1; i <= timeWindowList[j].GeneMutDic[Gene].Gene_SiteList.Count/20; i++)
                    for (k = 0; k < 20; k++)
                        siteID.Add(Convert.ToString(i) + AA20[k]);
                int v1, v2, tmpi;
                string tmps;
                for (v1 = 0; v1 < siteID.Count && v1 < 30; v1++)
                    for (v2 = v1 + 1; v2 < siteID.Count; v2++)
                        if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[siteID[v1]].MutRate < timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[siteID[v2]].MutRate)
                        { tmps = siteID[v1]; siteID[v1] = siteID[v2]; siteID[v2] = tmps; }

                timeWindowList[j].topN.Clear();
                for (i = 0; i < 5; i++)//【each clade top N】
                {
                    timeWindowList[j].topN.Add(siteID[i]);
                    if (!RBD5.Contains(siteID[i]) && timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[siteID[i]].MutRate > 0) RBD5.Add(siteID[i]);
                }
            }
            writep.WriteLine(title);
            RBD5.Sort();
            //输出
            for (i = 0; i < RBD5.Count(); i++)
            {
                title = RBD5[i];
                bool ind = true;
                for (j = 0; j < timeWindowList.Count; j++)
                {
                    if ((double)timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[RBD5[i]].Prop > 0.9)
                    { title += "\tNA"; ind = false; }
                    else if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[RBD5[i]].MutRate != -1)
                    {
                        title += "\t" + Convert.ToString(timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[RBD5[i]].MutRate);
                    }
                    else
                    {
                        title += "\t-6";
                    }
                }
                if(ind)
                    writep.WriteLine(title);
            }

            //输出所有突变
            for (j = 0; j < timeWindowList.Count; j++)
            {
                foreach (string vali in timeWindowList[j].GeneMutDic[Gene].Gene_SiteList.Keys)
                    if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[vali].EventCount != 0)
                        write.WriteLine(timeWindowList[j].WindowName + "\t" + Convert.ToString(timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum) + "\t" + vali + "\t" + Convert.ToString(timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[vali].EventCount));
            }
            write.Close();

            //输出topN的物种组成图
            writeBarplot.WriteLine("Clade\tMut\tPercent");
            for (k = 0; k < timeWindowList.Count; k++)
            {
                double otherMutNum = 0;
                for (i = 0; i < timeWindowList[k].topN.Count(); i++)
                {
                    otherMutNum += timeWindowList[k].GeneMutDic[Gene].Gene_SiteList[timeWindowList[k].topN[i]].EventCount;
                }
                writeBarplot.WriteLine(timeWindowList[k].WindowName + "\tZOther\t" + Convert.ToString(otherMutNum / timeWindowList[k].GeneMutDic[Gene].CladeMutEventNum));
            }
            writeBarplot.Close();
            write.Close();
            writep.Close();
            return;
        }

        static void CalMutRateOutputForNTD_Others(string Gene)//计算突变率，输出
        {
            int i, j, k;
            StreamWriter writep = new StreamWriter(workfold + "/GeneTopMut/CladeHeatmap_" + "OtherSpike" + ".txt");
            string output = "Clade\tCladeMutNum\tMut\tEventCount";
            for (j = 0; j < timeWindowList.Count; j++)
            {
                foreach (string valk in timeWindowList[j].GeneMutDic[Gene].Gene_SiteList.Keys)
                {
                    timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].Prop = (double)timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].EventWithThisMut / timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum;
                    double mutRate = (double)timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].EventCount / timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum;
                    if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].Prop < 0.9 && timeWindowList[j].GeneMutDic[Gene].CladeMutEventNum > 100 && mutRate >= 0.00001)
                    {
                        timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[valk].MutRate = mutRate;
                    }
                }
            }
            string title = "Mut";
            //输出每个窗口区域前
            List<string> RBD5 = new List<string>();
            for (j = 0; j < timeWindowList.Count; j++)
            {
                title += "\t" + timeWindowList[j].WindowName;
                List<string> siteID = new List<string>();
                for (i = 1; i <= 13; i++)
                    for (k = 0; k < 20; k++)
                        siteID.Add(Convert.ToString(i) + AA20[k]);
                for (i = 305; i <= 331; i++)
                    for (k = 0; k < 20; k++)
                        siteID.Add(Convert.ToString(i) + AA20[k]);
                for (i = 531; i <= 669; i++)
                    for (k = 0; k < 20; k++)
                        siteID.Add(Convert.ToString(i) + AA20[k]);
                int v1, v2, tmpi;
                string tmps;
                for (v1 = 0; v1 < siteID.Count && v1 < 30; v1++)
                    for (v2 = v1 + 1; v2 < siteID.Count; v2++)
                        if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[siteID[v1]].MutRate < timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[siteID[v2]].MutRate)
                        { tmps = siteID[v1]; siteID[v1] = siteID[v2]; siteID[v2] = tmps; }

                timeWindowList[j].topN.Clear();
                for (i = 0; i < 5; i++)//【每个clade取前几】
                {
                    timeWindowList[j].topN.Add(siteID[i]);
                    if (!RBD5.Contains(siteID[i]) && timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[siteID[i]].MutRate > 0) RBD5.Add(siteID[i]);
                }
            }
            writep.WriteLine(title);
            RBD5.Sort();
            //输出
            for (i = 0; i < RBD5.Count(); i++)
            {
                title = RBD5[i];
                bool ind = true;
                for (j = 0; j < timeWindowList.Count; j++)
                {
                    if ((double)timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[RBD5[i]].Prop > 0.9)
                    { title += "\tNA"; ind = false; }
                    else if (timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[RBD5[i]].MutRate != -1)
                    {
                        title += "\t" + Convert.ToString(timeWindowList[j].GeneMutDic[Gene].Gene_SiteList[RBD5[i]].MutRate);
                    }
                    else
                    {
                        title += "\t-6";
                    }
                }
                if (ind)
                    writep.WriteLine(title);
            }
            writep.Close();
            return;
        }

        static void CalMutRateOutput()//计算突变率，输出
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(workfold + "/CladeRBDMutRate.txt");
            StreamWriter write12 = new StreamWriter(workfold + "/CladeRBDMut_ESC12.txt");
            StreamWriter writep = new StreamWriter(workfold + "/CladeHeatmap.txt");
            StreamWriter writeBarplot = new StreamWriter(workfold + "/CladeBarplot.txt");
            string output = "Clade\tCladeMutNum\tMut\tESC\tEventCount";
            write.WriteLine(output);
            write12.WriteLine("Clade\tA\tB\tC\tD1\tD2\tE1\tE2.1\tE2.2\tE3\tF1\tF2\tF3");
            for (i = 331; i <= 531; i++)
            {
                if(i==352)
                    output = Convert.ToString(i);
                for (j = 0; j < timeWindowList.Count; j++)
                {
                    for (k = 0; k < 20; k++)
                    {
                        string key = Convert.ToString(i) + AA20[k];
                        //double mutRate = (double)timeWindowList[j].SpikeSiteList[key].EventCount / timeWindowList[j].CladeMutEventNum;
                        double mutRate = (double)timeWindowList[j].SpikeSiteList[key].EventCount / timeWindowList[j].CladeRBDMutEventNum;

                        if (timeWindowList[j].SpikeSiteList[key].Prop < 0.9 && timeWindowList[j].CladeMutEventNum > 100 && mutRate >= 0.00001)
                        {
                            timeWindowList[j].SpikeSiteList[key].MutRate = mutRate;
                        }
                    }
                }

            }
            
            string title = "Mut";
            //输出每个窗口RBD区域前10
            List<string> RBD5 = new List<string>();
            for (j = 0; j < timeWindowList.Count; j++)
            {
                title += "\t" + timeWindowList[j].WindowName;
                List<string> siteID = new List<string>();
                for (i = 331; i <= 531; i++)
                    for (k = 0; k < 20; k++)
                        siteID.Add(Convert.ToString(i) + AA20[k]);
                int v1, v2, tmpi;
                string tmps;
                for (v1 = 0; v1 < siteID.Count; v1++)
                    for (v2 = v1 + 1; v2 < siteID.Count; v2++)
                        if (timeWindowList[j].SpikeSiteList[siteID[v1]].MutRate < timeWindowList[j].SpikeSiteList[siteID[v2]].MutRate)
                        { tmps = siteID[v1]; siteID[v1] = siteID[v2]; siteID[v2] = tmps; }

                timeWindowList[j].topN.Clear();
                for (i = 0; i < 10; i++)//【每个clade取前几】
                {
                    timeWindowList[j].topN.Add(siteID[i]);
                    if (!RBD5.Contains(siteID[i]) && timeWindowList[j].SpikeSiteList[siteID[i]].MutRate > 0) RBD5.Add(siteID[i]);
                }

                //12类抗体的逃逸程度中位数
                write12.Write(timeWindowList[j].WindowName);
                for (tmpi = 0; tmpi < Group12.Length; tmpi++)
                {
                    List<double> scape = new List<double>();
                    double average = 0;
                    for (i = 0; i < siteID.Count(); i++)//【输出所有突变】siteID.Count()
                    {
                        for (k = 0; k < timeWindowList[j].SpikeSiteList[siteID[i]].EventCount; k++)
                        {
                            scape.Add(EscapeScore_WT_12[siteID[i] + "_" + Group12[tmpi]]);
                            average += EscapeScore_WT_12[siteID[i] + "_" + Group12[tmpi]];
                        }
                    }
                    scape.Sort();
                    average /= siteID.Count();
                    //write.Write("\t" + Convert.ToString(scape[scape.Count / 2]));
                    write12.Write("\t" + Convert.ToString(average));
                }
                write12.Write("\n");
            }
            
            writep.WriteLine(title);
            //手动读取给定的位点
            //RBD5.Clear();
            /*StreamReader read = new StreamReader(workfold + "/Mutations.txt");
            string line = read.ReadLine();
            while (line != null) 
            {
                //if(RBD5.Contains(line))RBD5.Remove(line); 
                if (!RBD5.Contains(line)) RBD5.Add(line);
                line = read.ReadLine(); 
            }*/

            RBD5.Sort();
            //输出
            for (i = 0; i < RBD5.Count(); i++)
            {
                title = RBD5[i];
                for (j = 0; j < timeWindowList.Count; j++)
                {
                    output = timeWindowList[j].WindowName + "\t";
                    output += RBD5[i] + "\t";
                    output += Convert.ToString((double)timeWindowList[j].SpikeSiteList[RBD5[i]].EventWithThisMut / timeWindowList[j].CladeMutEventNum);
                    if (timeWindowList[j].SpikeSiteList[RBD5[i]].MutRate != -1)
                        output += "\t" + Convert.ToString(timeWindowList[j].SpikeSiteList[RBD5[i]].MutRate);
                    else
                        output += "\t0";
                    //if((double)timeWindowList[j].SpikeSiteList[RBD5[i]].EventWithThisMut / timeWindowList[j].CladeMutEventNum > 0.01 || timeWindowList[j].SpikeSiteList[RBD5[i]].MutRate != -1)
                    //    write.WriteLine(output);

                    if ((double)timeWindowList[j].SpikeSiteList[RBD5[i]].EventWithThisMut / timeWindowList[j].CladeMutEventNum > 0.9)
                    { title += "\tNA"; }
                    else if (timeWindowList[j].SpikeSiteList[RBD5[i]].MutRate != -1)
                    {
                        title += "\t" + Convert.ToString(timeWindowList[j].SpikeSiteList[RBD5[i]].MutRate);
                    }
                    else
                    {
                        title += "\t-6";
                    }
                }
                writep.WriteLine(title);
            }

            //输出所有突变
            for (j = 0; j < timeWindowList.Count; j++)
            {
                foreach (string vali in timeWindowList[j].SpikeSiteList.Keys)
                    //for(k=0;k< timeWindowList[j].SpikeSiteList[vali].EventCount;k++)
                    if(timeWindowList[j].SpikeSiteList[vali].EventCount!=0)
                        write.WriteLine(timeWindowList[j].WindowName + "\t" + Convert.ToString(timeWindowList[j].CladeRBDMutEventNum) + "\t" + vali + "\t" + EscapeScore_WT_Total[vali] + "\t" + Convert.ToString(timeWindowList[j].SpikeSiteList[vali].EventCount));
            }
            write.Close();

            //输出topN的物种组成图
            writeBarplot.WriteLine("Clade\tMut\tPercent\tNSS\tBackMutProp");
            for (k = 0; k < timeWindowList.Count; k++)
            {
                double otherMutNum = timeWindowList[k].CladeRBDMutEventNum;
                for (i = 0; i < timeWindowList[k].topN.Count(); i++)
                {
                    otherMutNum -= timeWindowList[k].SpikeSiteList[timeWindowList[k].topN[i]].EventCount;
                }
                writeBarplot.WriteLine(timeWindowList[k].WindowName + "\tZOther\t" + Convert.ToString( 1 - otherMutNum / timeWindowList[k].CladeRBDMutEventNum) + "\t" + Convert.ToString((double)timeWindowList[k].NS_MutEventNum / timeWindowList[k].S_MutEventNum) + "\t" + Convert.ToString((double)timeWindowList[k].backmutNum / timeWindowList[k].NS_MutEventNum));
            }

            writeBarplot.Close();
            write.Close();
            writep.Close();
            write12.Close();
            return;
        }
        static void ReadClade()
        {
            StreamReader readClade = new StreamReader(workfold + "/CladeList.tsv");
            int i, j, k;
            string line = readClade.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                TimeWindow newt = new TimeWindow();
                newt.WindowName = line1[0];
                newt.DateStart = Convert.ToInt32(line1[1]);
                newt.DateEnd = Convert.ToInt32(line1[2]);
                for (i = 331; i <= 531; i++)
                    for (k = 0; k < 20; k++)
                    {
                        SpikeSite news = new SpikeSite();
                        news.SitePos = i;
                        newt.SpikeSiteList.Add(Convert.ToString(i) + AA20[k], news);
                    }

                //初始化基因列表
                string[] gene = "orf1aorf1b,S,ORF3a,E,M,ORF6,ORF7a,ORF7b,ORF8,N,ORF10".Split(',');
                string[] aalength = "7096,1274,276,76,223,62,122,44,122,420,39".Split(',');
                for (j = 0; j < gene.Length; j++)
                {
                    GeneWindow newg = new GeneWindow();
                    newg.GeneName = gene[j];
                    for (i = 1; i <= Convert.ToInt32(aalength[j]); i++)
                        for (k = 0; k < 20; k++)
                        {
                            SpikeSite news = new SpikeSite();
                            news.SitePos = i;
                            newg.Gene_SiteList.Add(Convert.ToString(i) + AA20[k], news);
                        }
                    newt.GeneMutDic.Add(gene[j], newg);
                }
                
                timeWindowList.Add(newt);
                line = readClade.ReadLine();
            }
            readClade.Close();
        }
        static void ReadEscapeScore()
        {
            StreamReader read = new StreamReader(workfold + "/Data/EscapeScore_PKU_WT.12EscOnly.txt");
            int i, j, k;
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                for (i = 0; i < Group12.Length; i++)
                    EscapeScore_WT_12.Add(line1[0] + "_" + Group12[i], Convert.ToDouble(line1[i+1]));
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(workfold + "/Data/EscapeScore_PKU_WT.single.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                EscapeScore_WT_Total.Add(line1[0], Convert.ToDouble(line1[1]));
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(workfold + "/Data/MutationEsc12.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                bool esc = false;
                for (i = 1; i < line1.Count(); i++)
                    if (line1[i] == "1") esc = true;
                IfMutationEscape.Add(line1[0], esc);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(workfold + "/Data/EscapeScore_PKU_BA.1.single.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                EscapeScore_BA1_Total.Add(line1[0], Convert.ToDouble(line1[1]));
                line = read.ReadLine();
            }
            read.Close();
        }
        static void ESC_GSEA()
        {
            StreamWriter write_curve = new StreamWriter(workfold + "/GSEA/Mutation_ES_Curve.txt");
            StreamWriter write_es_pvalue = new StreamWriter(workfold + "/GSEA/Mutation_ES_Pvalue.txt");
            StreamWriter write_escape_score = new StreamWriter(workfold + "/GSEA/Mutation_escape_score.txt");
            write_curve.WriteLine("Clade\tMutNum\tES");
            write_escape_score.WriteLine("Clade\tMutNum\tESC");
            write_es_pvalue.WriteLine("Clade\tES\tRank\tPvalue");
            int i, j, k;
            for(i=0;i<timeWindowList.Count;i++)
            {
                List<string> mutationlist = new List<string>();
                List<int> count = new List<int>();//发生次数
                List<bool> escape = new List<bool>();//是否逃逸突变
                int totalMutNum = 0;
                int escapeMutNum = 0;
                double sumEscScore = 0;
                foreach(string val in timeWindowList[i].SpikeSiteList.Keys)
                {
                    if(timeWindowList[i].SpikeSiteList[val].EventCount!=0)
                    {
                        totalMutNum++;
                        mutationlist.Add(val);
                        count.Add(timeWindowList[i].SpikeSiteList[val].EventCount);
                        escape.Add(IfMutationEscape[val]);
                        if (IfMutationEscape[val])
                        {
                            escapeMutNum++;
                            sumEscScore += EscapeScore_WT_Total[val];//这里切换WT或者BA.1
                        }
                    }
                }
                int tmpi; string tmps; bool tmpb;
                
                //开始排序 从大到小
                for(j=0;j<mutationlist.Count;j++)
                    for(k=j+1;k<mutationlist.Count;k++)
                        if(count[j]<count[k])
                        {
                            tmps = mutationlist[j];mutationlist[j] = mutationlist[k];mutationlist[k] = tmps;
                            tmpi = count[j];count[j] = count[k];count[k] = tmpi;
                            tmpb = escape[j];escape[j] = escape[k];escape[k] = tmpb;
                        }

                //计算es
                double ES = CalculateEnrichmentScore(mutationlist, escape, sumEscScore, (double)1 / (totalMutNum - escapeMutNum));
                int rank = 0;
                for (k = 0; k < EnrichmentScoreCurve.Count; k++)
                {
                    if (EnrichmentScoreCurve[k] == ES)
                        rank = k;
                    write_curve.WriteLine(timeWindowList[i].WindowName + "\t" + Convert.ToString((double)k / EnrichmentScoreCurve.Count) + "\t" + EnrichmentScoreCurve[k]);
                    if(escape[k])
                        write_escape_score.WriteLine(timeWindowList[i].WindowName + "\t" + Convert.ToString((double)k / EnrichmentScoreCurve.Count) + "\t" + EscapeScore_WT_Total[mutationlist[k]]);
                }

                //随机1000次，计算pvalue
                double pvalue = 0;
                for(k=0;k<50000;k++)
                {
                    List<string> tmp_mutationlist = new List<string>(mutationlist);
                    List<bool> tmp_escape = new List<bool>(escape);
                    List<string> rdm_mutationlist = new List<string>();
                    List<bool> rdm_escape = new List<bool>();
                    Random rdm = new Random((int)(k*Math.E));
                    while(tmp_mutationlist.Count>0)
                    {
                        int rdnnumber = rdm.Next(0, tmp_mutationlist.Count);
                        rdm_mutationlist.Add(tmp_mutationlist[rdnnumber]);
                        rdm_escape.Add(tmp_escape[rdnnumber]);
                        tmp_mutationlist.RemoveRange(rdnnumber, 1);
                        tmp_escape.RemoveRange(rdnnumber, 1);
                    }
                    if (CalculateEnrichmentScore(rdm_mutationlist, rdm_escape, sumEscScore, (double)1 / (totalMutNum - escapeMutNum)) >= ES)
                        pvalue++;
                }
                pvalue /= 50000;
                write_es_pvalue.WriteLine(timeWindowList[i].WindowName + "\t" + Convert.ToString(ES) + "\t" + Convert.ToString((double)rank/mutationlist.Count) + "\t" + Convert.ToString(pvalue));
            }
            write_curve.Close();
            write_es_pvalue.Close();
            write_escape_score.Close();
        }
        static double CalculateEnrichmentScore(List<string> mutationlist, List<bool> escape, double sumEscScore, double penaltyPoints)
        {
            int i, j, k;
            double es = 0,esmax = 0;
            EnrichmentScoreCurve = new List<double>();
            for(i=0;i<mutationlist.Count;i++)
            {
                if (escape[i])
                {
                    es += (double)EscapeScore_WT_Total[mutationlist[i]] / sumEscScore;
                }
                else
                    es -= penaltyPoints;
                if (es > esmax) esmax = es;
                EnrichmentScoreCurve.Add(es);
            }
            return esmax;
        }
        static void Main(string[] args)
        {
            ReadClade();//读入要计算的进化枝
            ReadMutAnno();//读入突变注释
            ReadMutEvent();//读入突变事件
            ReadEscapeScore();//读入逃逸得分
            CalMutRateOutput();//计算突变率，输出

            string[] gene = "orf1aorf1b,S,ORF3a,E,M,ORF6,ORF7a,ORF7b,ORF8,N,ORF10".Split(',');
            for(int i=0;i<gene.Length;i++)
                CalMutRateOutputForGenes(gene[i]);
            //CalMutRateOutputForNTD_Others("S");

            ESC_GSEA();//计算突变逃逸得分的GSEA
        }
    }
}
