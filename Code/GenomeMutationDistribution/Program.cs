using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 计算不同时间_地点的突变率
{
    class SpikeRegion //RBD NTD Other
    {
        public int RBDTotal = 0;
        public int NTDTotal = 0;
        public int OtherTotal = 0;
        public int RBDNonZeroSite = 0;
        public int NTDNonZeroSite = 0;
        public int OtherNonZeroSite = 0;
        public Dictionary<int, int> RBDCount = new Dictionary<int, int>();//计数每个氨基酸位置的突变发生次数
        public Dictionary<int, int> NTDCount = new Dictionary<int, int>();//计数每个氨基酸位置的突变发生次数
        public Dictionary<int, int> OtherCount = new Dictionary<int, int>();//计数每个氨基酸位置的突变发生次数
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
    class TimeMutationRate //时间间隔内的突变率
    {
        public int CollectionMonth;
        public int TotalMutEventNum = 0;
        public int TotalNSMutEventNum = 0;
        public Dictionary<string, int> GeneMutEventNum = new Dictionary<string, int>();
        public Dictionary<string, int> GeneNSMutEventNum = new Dictionary<string, int>();
    }
    class Program
    {
        static string location = "global";
        //static string location = "United Kingdom END England";
        static string workfold = "***/Usher1128";
        static Dictionary<string, MutationAnno> MutAnnoMap = new Dictionary<string, MutationAnno>();
        static Dictionary<string, TimeMutationRate> TimeMutRateList = new Dictionary<string, TimeMutationRate>();
        static Dictionary<string, SpikeRegion> TimeSpikeList = new Dictionary<string, SpikeRegion>();
        static List<string> GeneList = new List<string>();
        static List<double> GeneLengthPercent = new List<double>();//按长度归一化
        static Dictionary<string, int> MutationCount = new Dictionary<string, int>();//突变发生次数
        static void ReadMutAnno()//读入突变注释
        {
            StreamReader read = new StreamReader(workfold + "/Data/AllMutationAnnomutresult.tsv");
            string line = read.ReadLine();
            while(line!=null)
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
                if (!GeneList.Contains(line1[5]))
                {
                    GeneList.Add(line1[5]);
                    GeneLengthPercent.Add(1);
                }
                else
                    GeneLengthPercent[GeneList.IndexOf(line1[5])]++;
                line = read.ReadLine();
            }
            read.Close();

            int i, j;
            double k = 0;
            for (i = 0; i < GeneLengthPercent.Count; i++)
                k += GeneLengthPercent[i];
            for (i = 0; i < GeneLengthPercent.Count; i++)
                GeneLengthPercent[i] /= k;
            GeneList.Add("RBD");
            GeneLengthPercent.Add(0.020165);
            GeneList.Add("NTD");
            GeneLengthPercent.Add(0.029395);
            GeneList.Add("SpikeOther");
            GeneLengthPercent.Add(0.078253);
            return;
        }
        static void ReadMutEvent()//读入突变事件
        {
            int i, j, k;
            //初始化
            Dictionary<string, int> ModGeneMutEventNum = new Dictionary<string, int>();
            for(i=0;i<GeneList.Count();i++)
            {
                ModGeneMutEventNum.Add(GeneList[i], 0);
            }
            for (i=2020;i<=2022;i++)
            {
                for(j=1;j<=12;j++)
                {
                    TimeMutationRate newt = new TimeMutationRate();
                    string key;
                    if (j < 10)
                        key = Convert.ToString(i) + "-0" + Convert.ToString(j);
                    else
                        key = Convert.ToString(i) + "-" + Convert.ToString(j);
                    newt.CollectionMonth = i * 100 + j;
                    newt.TotalMutEventNum=0;
                    newt.TotalNSMutEventNum = 0;
                    newt.GeneMutEventNum = new Dictionary<string, int>(ModGeneMutEventNum);
                    newt.GeneNSMutEventNum = new Dictionary<string, int>(ModGeneMutEventNum);
                    TimeMutRateList.Add(key, newt);

                    SpikeRegion news = new SpikeRegion();
                    for (k = 1; k <= 12; k++) news.OtherCount.Add(k, 0);
                    for (k = 13; k <= 305; k++) news.NTDCount.Add(k, 0);
                    for (k = 306; k <= 330; k++) news.OtherCount.Add(k, 0);
                    for (k = 331; k <= 531; k++) news.RBDCount.Add(k, 0);
                    for (k = 532; k <= 1274; k++) news.OtherCount.Add(k, 0);
                    TimeSpikeList.Add(key, news);
                }
            }

            //读入突变事件
            StreamReader read = new StreamReader(workfold + "/Data/MAT1128.json.mutevent");
            read.ReadLine();
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if(location == "global" || location.ToUpper().Contains(line1[6].ToUpper()))
                {
                    if (line1[5].Length < 7)
                    {
                        line = read.ReadLine();
                        continue;
                    }

                    string[] nucmut = line1[2].Split(',');
                    //测试这个突变是不是回复突变
                    /*string backmuttest = "" + line1[3][line1[3].Length - 1] + line1[3].Substring(1, line1[3].Length - 2) + line1[3][0];
                    if (nucmut.Contains(backmuttest))
                    { line = read.ReadLine(); continue; }*/

                    string month = line1[5].Substring(0, 7);
                    if (TimeMutRateList.ContainsKey(month))
                    {
                        int AApos = -1;
                        int startp = 0, endp = MutAnnoMap[line1[3]].AminoChange.Length - 1;
                        if (MutAnnoMap[line1[3]].AminoChange != "intergenic" && MutAnnoMap[line1[3]].AminoChange != "Overlap")
                        {
                            while (!(MutAnnoMap[line1[3]].AminoChange[startp] <= '9' && MutAnnoMap[line1[3]].AminoChange[startp] >= '0')) startp++;
                            while (!(MutAnnoMap[line1[3]].AminoChange[endp] <= '9' && MutAnnoMap[line1[3]].AminoChange[endp] >= '0')) endp--;
                            AApos = Convert.ToInt32(MutAnnoMap[line1[3]].AminoChange.Substring(startp, endp - startp + 1));
                        }

                        TimeMutRateList[month].TotalMutEventNum++;
                        TimeMutRateList[month].GeneMutEventNum[MutAnnoMap[line1[3]].Gene]++;
                        if (MutAnnoMap[line1[3]].Gene == "S")
                        {
                            if (AApos >= 331 && AApos <= 531)
                                TimeMutRateList[month].GeneMutEventNum["RBD"]++;
                            else if (AApos >= 13 && AApos <= 305)
                                TimeMutRateList[month].GeneMutEventNum["NTD"]++;
                            else
                                TimeMutRateList[month].GeneMutEventNum["SpikeOther"]++;
                        }
                        //NS突变
                        if (MutAnnoMap[line1[3]].NSS == "NS")
                        {
                            TimeMutRateList[month].TotalNSMutEventNum++;
                            TimeMutRateList[month].GeneNSMutEventNum[MutAnnoMap[line1[3]].Gene]++;
                            if (MutAnnoMap[line1[3]].Gene == "S")
                            {
                                if (AApos >= 331 && AApos <= 531)
                                {
                                    TimeMutRateList[month].GeneNSMutEventNum["RBD"]++;
                                    TimeSpikeList[month].RBDCount[AApos]++;
                                    TimeSpikeList[month].RBDTotal++;
                                }
                                else if (AApos >= 13 && AApos <= 305)
                                {
                                    TimeMutRateList[month].GeneNSMutEventNum["NTD"]++;
                                    TimeSpikeList[month].NTDCount[AApos]++;
                                    TimeSpikeList[month].NTDTotal++;
                                }
                                else
                                {
                                    TimeMutRateList[month].GeneNSMutEventNum["SpikeOther"]++;
                                    TimeSpikeList[month].OtherCount[AApos]++;
                                    TimeSpikeList[month].OtherTotal++;
                                }
                            }
                        }
                        //计数突变发生次数
                        if (MutationCount.ContainsKey(line1[3]))
                            MutationCount[line1[3]]++;
                        else
                            MutationCount.Add(line1[3], 1);
                    }
                }
                line = read.ReadLine();
            }

            //叠加 调整为三个月滑窗 345 678 91011 1212
            for (i = 2020; i <= 2022; i++)
            {
                for (j = 1; j <= 12; j++)
                {
                    string key;
                    if (j < 10)
                        key = Convert.ToString(i) + "-0" + Convert.ToString(j);
                    else
                        key = Convert.ToString(i) + "-" + Convert.ToString(j);
                    if(key.Substring(5,2)=="03" || key.Substring(5, 2) == "06" || key.Substring(5, 2) == "09")
                    {
                        //后俩月的加过来
                        string key1, key2;
                        k = j + 1;
                        if (k < 10)
                            key1 = Convert.ToString(i) + "-0" + Convert.ToString(k);
                        else
                            key1 = Convert.ToString(i) + "-" + Convert.ToString(k);
                        k++;
                        if (k < 10)
                            key2 = Convert.ToString(i) + "-0" + Convert.ToString(k);
                        else
                            key2 = Convert.ToString(i) + "-" + Convert.ToString(k);
                        TimeMutRateList[key].TotalMutEventNum += TimeMutRateList[key1].TotalMutEventNum + TimeMutRateList[key2].TotalMutEventNum;
                        TimeMutRateList[key].TotalNSMutEventNum += TimeMutRateList[key1].TotalNSMutEventNum + TimeMutRateList[key2].TotalNSMutEventNum;
                        for (k = 0; k < GeneList.Count; k++)
                        { 
                            TimeMutRateList[key].GeneMutEventNum[GeneList[k]] += TimeMutRateList[key1].GeneMutEventNum[GeneList[k]] + TimeMutRateList[key2].GeneMutEventNum[GeneList[k]];
                            TimeMutRateList[key].GeneNSMutEventNum[GeneList[k]] += TimeMutRateList[key1].GeneNSMutEventNum[GeneList[k]] + TimeMutRateList[key2].GeneNSMutEventNum[GeneList[k]];
                        }
                        for (k = 1; k <= 12; k++) TimeSpikeList[key].OtherCount[k] += TimeSpikeList[key1].OtherCount[k] + TimeSpikeList[key2].OtherCount[k];
                        for (k = 13; k <= 305; k++) TimeSpikeList[key].NTDCount[k] += TimeSpikeList[key1].NTDCount[k] + TimeSpikeList[key2].NTDCount[k];
                        for (k = 306; k <= 330; k++) TimeSpikeList[key].OtherCount[k] += TimeSpikeList[key1].OtherCount[k] + TimeSpikeList[key2].OtherCount[k];
                        for (k = 331; k <= 531; k++) TimeSpikeList[key].RBDCount[k] += TimeSpikeList[key1].RBDCount[k] + TimeSpikeList[key2].RBDCount[k];
                        for (k = 532; k <= 1274; k++) TimeSpikeList[key].OtherCount[k] += TimeSpikeList[key1].OtherCount[k] + TimeSpikeList[key2].OtherCount[k];
                        TimeSpikeList[key].RBDTotal += TimeSpikeList[key1].RBDTotal + TimeSpikeList[key2].RBDTotal;
                        TimeSpikeList[key].NTDTotal += TimeSpikeList[key1].NTDTotal + TimeSpikeList[key2].NTDTotal;
                        TimeSpikeList[key].OtherTotal += TimeSpikeList[key1].OtherTotal + TimeSpikeList[key2].OtherTotal;
                    }
                    else
                    {
                        if (key.Substring(5, 2) == "12" && i!=2022)
                        {
                            string key1, key2;
                            key1 = Convert.ToString(i+1) + "-01";
                            key2 = Convert.ToString(i+1) + "-02";
                            TimeMutRateList[key].TotalMutEventNum += TimeMutRateList[key1].TotalMutEventNum + TimeMutRateList[key2].TotalMutEventNum;
                            TimeMutRateList[key].TotalNSMutEventNum += TimeMutRateList[key1].TotalNSMutEventNum + TimeMutRateList[key2].TotalNSMutEventNum;
                            for (k = 0; k < GeneList.Count; k++)
                            {
                                TimeMutRateList[key].GeneMutEventNum[GeneList[k]] += TimeMutRateList[key1].GeneMutEventNum[GeneList[k]] + TimeMutRateList[key2].GeneMutEventNum[GeneList[k]];
                                TimeMutRateList[key].GeneNSMutEventNum[GeneList[k]] += TimeMutRateList[key1].GeneNSMutEventNum[GeneList[k]] + TimeMutRateList[key2].GeneNSMutEventNum[GeneList[k]];
                            }
                            for (k = 1; k <= 12; k++) TimeSpikeList[key].OtherCount[k] += TimeSpikeList[key1].OtherCount[k] + TimeSpikeList[key2].OtherCount[k];
                            for (k = 13; k <= 305; k++) TimeSpikeList[key].NTDCount[k] += TimeSpikeList[key1].NTDCount[k] + TimeSpikeList[key2].NTDCount[k];
                            for (k = 306; k <= 330; k++) TimeSpikeList[key].OtherCount[k] += TimeSpikeList[key1].OtherCount[k] + TimeSpikeList[key2].OtherCount[k];
                            for (k = 331; k <= 531; k++) TimeSpikeList[key].RBDCount[k] += TimeSpikeList[key1].RBDCount[k] + TimeSpikeList[key2].RBDCount[k];
                            for (k = 532; k <= 1274; k++) TimeSpikeList[key].OtherCount[k] += TimeSpikeList[key1].OtherCount[k] + TimeSpikeList[key2].OtherCount[k];
                            TimeSpikeList[key].RBDTotal += TimeSpikeList[key1].RBDTotal + TimeSpikeList[key2].RBDTotal;
                            TimeSpikeList[key].NTDTotal += TimeSpikeList[key1].NTDTotal + TimeSpikeList[key2].NTDTotal;
                            TimeSpikeList[key].OtherTotal += TimeSpikeList[key1].OtherTotal + TimeSpikeList[key2].OtherTotal;
                        }
                    }
                }
            }

            //统计Spike三个位置有多少位点发生了突变
            foreach(string key in TimeSpikeList.Keys)
            {
                foreach (int key1 in TimeSpikeList[key].OtherCount.Keys)
                    if (TimeSpikeList[key].OtherCount[key1] != 0) TimeSpikeList[key].OtherNonZeroSite++;
                foreach (int key1 in TimeSpikeList[key].NTDCount.Keys)
                    if (TimeSpikeList[key].NTDCount[key1] != 0) TimeSpikeList[key].NTDNonZeroSite++;
                foreach (int key1 in TimeSpikeList[key].RBDCount.Keys)
                    if (TimeSpikeList[key].RBDCount[key1] != 0) TimeSpikeList[key].RBDNonZeroSite++;
            }

            read.Close();
            return;
        }
        static void CalMutRateOutput()//计算突变率，输出
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(workfold + "/MutationRate.txt");
            StreamWriter writeNSp = new StreamWriter(workfold + "/NSPercent.txt");
            StreamWriter writeNSMutRate = new StreamWriter(workfold + "/NSMutRate.txt");
            StreamWriter writeSMutRate = new StreamWriter(workfold + "/SMutRate.txt");
            StreamWriter writeMutCount = new StreamWriter(workfold + "/MutationCount.txt");
            StreamWriter writeSpikeMutCount = new StreamWriter(workfold + "/SpikeMutationCount.txt");
            StreamWriter writeSpike50 = new StreamWriter(workfold + "/Spike50.txt");

            string output = "Gene";
            foreach (string val in TimeMutRateList.Keys)
                if(val!="2020-01"&& val != "2020-02"  && val != "2022-12")
                    if (val.Contains("03") || val.Contains("06") || val.Contains("09") || val.Contains("12"))
                        output += "\t" + Convert.ToInt32(TimeMutRateList[val].CollectionMonth);
            write.WriteLine(output);
            writeNSp.WriteLine(output);
            writeNSMutRate.WriteLine(output);
            writeSMutRate.WriteLine(output);
            writeSpike50.WriteLine("Month\tRegion\tM50");

            for (i=0;i<GeneList.Count;i++)
            {
                if (GeneList[i] == "Overlap")
                    continue;

                //各个基因突变的百分比
                output = GeneList[i];
                foreach (string val in TimeMutRateList.Keys)
                    if (val != "2020-01" && val != "2020-02"  && val != "2022-12")
                        if (val.Contains("03") || val.Contains("06") || val.Contains("09") || val.Contains("12"))
                            output += "\t" + Convert.ToString((double)TimeMutRateList[val].GeneMutEventNum[GeneList[i]] / Math.Max(1,TimeMutRateList[val].TotalMutEventNum) / GeneLengthPercent[i]);
                write.WriteLine(output);

                //各个基因NS突变的百分比
                output = GeneList[i];
                foreach (string val in TimeMutRateList.Keys)
                    if (val != "2020-01" && val != "2020-02" && val != "2022-12")
                        if (val.Contains("03") || val.Contains("06") || val.Contains("09") || val.Contains("12"))
                            output += "\t" + Convert.ToString((double)TimeMutRateList[val].GeneNSMutEventNum[GeneList[i]] / Math.Max(1, TimeMutRateList[val].TotalNSMutEventNum) / GeneLengthPercent[i]);
                writeNSMutRate.WriteLine(output);

                //各个基因S突变的百分比
                output = GeneList[i];
                foreach (string val in TimeMutRateList.Keys)
                    if (val != "2020-01" && val != "2020-02" && val != "2022-12")
                        if (val.Contains("03") || val.Contains("06") || val.Contains("09") || val.Contains("12"))
                            output += "\t" + Convert.ToString(((double)TimeMutRateList[val].GeneMutEventNum[GeneList[i]] - (double)TimeMutRateList[val].GeneNSMutEventNum[GeneList[i]]) / Math.Max(1, TimeMutRateList[val].TotalMutEventNum - TimeMutRateList[val].TotalNSMutEventNum) / GeneLengthPercent[i]);
                writeSMutRate.WriteLine(output);

                //各个基因突变的NS占比 %
                output = GeneList[i];
                foreach (string val in TimeMutRateList.Keys)
                    if (val != "2020-01" && val != "2020-02" && val != "2022-12")
                        if (val.Contains("03") || val.Contains("06") || val.Contains("09") || val.Contains("12"))
                            output += "\t" + Convert.ToString((double)TimeMutRateList[val].GeneNSMutEventNum[GeneList[i]] / Math.Max(1, TimeMutRateList[val].GeneMutEventNum[GeneList[i]]));
                writeNSp.WriteLine(output);

            }

            foreach (string val1 in MutationCount.Keys)
                writeMutCount.WriteLine(val1 + "\t" + Convert.ToString(MutationCount[val1]));

            int TotalRBDN = 0,N;
            List<int> TotalRBDList = new List<int>();
            for (k = 331; k <= 531; k++)
                TotalRBDList.Add(0);

            foreach (string val in TimeMutRateList.Keys)
                if (val != "2020-01" && val != "2020-02" && val != "2022-12")
                    if (val.Contains("03") || val.Contains("06") || val.Contains("09") || val.Contains("12"))
                    {
                        List<int> tmpOtherList = new List<int>();
                        List<int> tmpRBDList = new List<int>();
                        List<int> tmpNTDList = new List<int>();
                        for (k = 1; k <= 12; k++) { tmpOtherList.Add(TimeSpikeList[val].OtherCount[k]); writeSpikeMutCount.WriteLine(val + "\tOther\t" + Convert.ToString(k) + "\t" + TimeSpikeList[val].OtherCount[k]); }
                        for (k = 13; k <= 305; k++) { tmpNTDList.Add(TimeSpikeList[val].NTDCount[k]); writeSpikeMutCount.WriteLine(val + "\tNTD\t" + Convert.ToString(k) + "\t" + TimeSpikeList[val].NTDCount[k]); }
                        for (k = 306; k <= 330; k++) { tmpOtherList.Add(TimeSpikeList[val].OtherCount[k]); writeSpikeMutCount.WriteLine(val + "\tOther\t" + Convert.ToString(k) + "\t" + TimeSpikeList[val].OtherCount[k]); }
                        for (k = 331; k <= 531; k++) 
                        { 
                            tmpRBDList.Add(TimeSpikeList[val].RBDCount[k]); writeSpikeMutCount.WriteLine(val + "\tRBD\t" + Convert.ToString(k) + "\t" + TimeSpikeList[val].RBDCount[k]); 
                            TotalRBDList[k-331]+= TimeSpikeList[val].RBDCount[k];
                            TotalRBDN += TimeSpikeList[val].RBDCount[k];
                        }
                        for (k = 532; k <= 1274; k++) { tmpOtherList.Add(TimeSpikeList[val].OtherCount[k]); writeSpikeMutCount.WriteLine(val + "\tOther\t" + Convert.ToString(k) + "\t" + TimeSpikeList[val].OtherCount[k]); }

                        int tmpi;
                        tmpOtherList.Sort();
                        tmpRBDList.Sort();
                        tmpNTDList.Sort();
                        N = 0;
                        for (i = tmpOtherList.Count() - 1; i >= 0; i--)
                        {
                            if (tmpOtherList[i] + N >= TimeSpikeList[val].OtherTotal / 2)
                            { //writeSpike50.WriteLine(val + "\tOther\t" + Convert.ToString((double)(tmpOtherList.Count() - i) / TimeSpikeList[val].OtherNonZeroSite)); break;
                              }
                            else
                            
                                N += tmpOtherList[i];
                                //if (i == tmpOtherList.Count() - 4)
                                //    writeSpike50.WriteLine(val + "\tOther\t" + Convert.ToString((double)N / TimeSpikeList[val].OtherTotal));
                            
                        }
                        N = 0;
                            for (i = tmpNTDList.Count() - 1; i >= 0; i--)
                            if (tmpNTDList[i] + N >= TimeSpikeList[val].NTDTotal / 2)
                            {// writeSpike50.WriteLine(val + "\tNTD\t" + Convert.ToString((double)(tmpNTDList.Count() - i) / TimeSpikeList[val].NTDNonZeroSite)); break;
                             }
                            else
                            {
                                N += tmpNTDList[i];
                                //if (i == tmpNTDList.Count() - 4)
                                //    writeSpike50.WriteLine(val + "\tNTD\t" + Convert.ToString((double)N / TimeSpikeList[val].NTDTotal)); 
                            }
                        N = 0;
                        for (i = tmpRBDList.Count() - 1; i >= 0; i--)
                        {
                                        if (tmpRBDList[i] + N >= TimeSpikeList[val].RBDTotal / 2)
                                        { writeSpike50.WriteLine(val + "\tRBD\t" + Convert.ToString((double)(tmpRBDList.Count() - i) / TimeSpikeList[val].RBDNonZeroSite)) ; 
                                            //writeSpike50.WriteLine(val + "\tRBD\t" + Convert.ToString(tmpRBDList.Count() - i) + "\t" + Convert.ToString(TimeSpikeList[val].RBDNonZeroSite)); 
                                            break;
                                        }
                                        else
                                            N += tmpRBDList[i];
                                //if (i == tmpRBDList.Count() - 4)
                                //    writeSpike50.WriteLine(val + "\tRBD\t" + Convert.ToString((double)N / TimeSpikeList[val].RBDTotal)); 
                        }
                    }
            int RBDNonZeroSite = 0;
            for (k = 331; k <= 531; k++)
                if (TotalRBDList[k-331] != 0) RBDNonZeroSite++;
            TotalRBDList.Sort();
            N = 0;
            for (i = TotalRBDList.Count() - 1; i >= 0; i--)
            {
                if (TotalRBDList[i] + N >=  TotalRBDN / 2)
                {
                    writeSpike50.WriteLine("Total\tRBD\t" + Convert.ToString((double)(TotalRBDList.Count() - i) / RBDNonZeroSite));
                    break;
                }
                else
                    N += TotalRBDList[i];
            }


                write.Close();
            writeNSp.Close();
            writeNSMutRate.Close();
            writeSMutRate.Close();
            writeMutCount.Close();
            writeSpike50.Close();
            return;
        }
        static void Main(string[] args)
        {
            ReadMutAnno();//读入突变注释
            ReadMutEvent();//读入突变事件
            CalMutRateOutput();//计算突变率，输出
        }
    }
}
