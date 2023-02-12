using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 统计每一个lineage的突变和meta分布
{
    public class Lineage
    {
        public string lineageName;
        public List<int> collectionDate = new List<int>();
        public int totalSeqCount = 0;
        public Dictionary<string, double> Dic_mutProp = new Dictionary<string, double>();
        public List<string> countryList = new List<string>();
    }
    class Program
    {
        static Dictionary<string, Lineage> lineageDic = new Dictionary<string, Lineage>();
        static void Main(string[] args)
        {
            int i, j, k;
            StreamReader read = new StreamReader("***/sample_mut_loc_time.tsv.AAmut");
            StreamWriter write = new StreamWriter("***/Data/Lineage.Date.AAmut");
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if(line1.Length > 4 && lineageDic.ContainsKey(line1[4]))
                {
                    string[] date = line1[3].Split('-');
                    if (date.Count() == 3)
                    {
                        int collectionDate = Convert.ToInt32(date[0]) * 10000 + Convert.ToInt32(date[1]) * 100 + Convert.ToInt32(date[2]);
                        lineageDic[line1[4]].collectionDate.Add(collectionDate);
                        string[] mut = line1[1].Split(' ');
                        for (i = 0; i < mut.Length; i++)
                        {
                            if (lineageDic[line1[4]].Dic_mutProp.ContainsKey(mut[i]))
                                lineageDic[line1[4]].Dic_mutProp[mut[i]]++;
                            else
                                lineageDic[line1[4]].Dic_mutProp.Add(mut[i], 1);
                        }
                        if (!lineageDic[line1[4]].countryList.Contains(line1[2]))
                            lineageDic[line1[4]].countryList.Add(line1[2]);
                        lineageDic[line1[4]].totalSeqCount++;
                    }
                }
                else
                {
                    Lineage newl = new Lineage();
                    string[] date = line1[3].Split('-');
                    if (date.Count() == 3)
                    {
                        int collectionDate = Convert.ToInt32(date[0]) * 10000 + Convert.ToInt32(date[1]) * 100 + Convert.ToInt32(date[2]);
                        newl.collectionDate.Add(collectionDate);
                        string[] mut = line1[1].Split(' ');
                        for (i = 0; i < mut.Length; i++)
                        {
                            if (newl.Dic_mutProp.ContainsKey(mut[i]))
                                newl.Dic_mutProp[mut[i]]++;
                            else
                                newl.Dic_mutProp.Add(mut[i], 1);
                        }
                        newl.totalSeqCount++;
                        newl.countryList.Add(line1[2]);
                        lineageDic.Add(line1[4], newl);
                    }
                }
                line = read.ReadLine();
            }
            read.Close();

            //每个lineage 采样时间排序，算每个突变的占比 输出
            write.WriteLine("Lineage\tCollectionDate5P\tCollectionDate25P\tMut70P\tCountryNumber\tSeqNumber");
            foreach (string val in lineageDic.Keys)
            {
                string output = val + "\t";
                lineageDic[val].collectionDate.Sort();
                output += Convert.ToString(lineageDic[val].collectionDate[lineageDic[val].collectionDate.Count / 20]) + "\t";
                output += Convert.ToString(lineageDic[val].collectionDate[lineageDic[val].collectionDate.Count / 4]) + "\t";
                string mut = "";
                foreach (string mal in lineageDic[val].Dic_mutProp.Keys)
                    if (lineageDic[val].Dic_mutProp[mal] / lineageDic[val].totalSeqCount >= 0.7)
                        mut += mal + ",";
                if (mut != "") mut = mut.Substring(0, mut.Length - 1);
                output += mut + "\t";
                output += Convert.ToString(lineageDic[val].countryList.Count) + "\t";
                output += Convert.ToString(lineageDic[val].totalSeqCount);
                write.WriteLine(output);
            }
            write.Close();
        }
    }
}
