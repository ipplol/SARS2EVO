using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

namespace SampleMutlist过滤高质量序列并抽样
{
    class Metadata
    {
        public string CollectionDate;
        public string CollectionMonth;
        public string Lineage;
        public string VOC;
        public string Location;
    }
    class Sequence
    {
        public string Accession;
        public Metadata Meta;
        public string Mutation;
    }
    class Program
    {
        static Dictionary<string, Metadata> Metadata_Dic = new Dictionary<string, Metadata>();
        static List<Sequence> Sequence_List = new List<Sequence>();
        static string workfold = "***/Data";
        static void ReadinMutlist()
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/sample_mutlist.tsv");
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                string[] line2 = line1[0].Split('|');
                if(Metadata_Dic.ContainsKey(line2[1]))
                {
                    Sequence news = new Sequence();
                    news.Accession = line2[1];
                    news.Mutation = line1[1];
                    news.Meta = Metadata_Dic[line2[1]];
                    Sequence_List.Add(news);
                }
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinMetadata()
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/public-latest.metadata.tsv");
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                //if(line1[5]=="Complete" && line1[7]=="High")
                //{
                    Metadata newa = new Metadata();
                    newa.CollectionDate = line1[2];
                    if (line1[2].Split('-').Length == 3)
                        newa.CollectionMonth = line1[2].Substring(0, 7);
                    else
                        newa.CollectionMonth = "NA";
                    newa.Lineage = line1[8];
                    newa.Location = line1[3];
                    if(line1[1]!="null" && line1[1] != "NULL" && line1[1] != "" && line1[1] != " " && !Metadata_Dic.ContainsKey(line1[1]))
                        Metadata_Dic.Add(line1[1], newa);
                    if (line1[0] != "null" && line1[0] != "NULL" && line1[0] != "" && line1[0] != " " && !Metadata_Dic.ContainsKey(line1[0]))
                        Metadata_Dic.Add(line1[0], newa);
                //}
                line = read.ReadLine();
            }
            return;
            read.Close();
        }
        static void SamplingByMonth()//根据采样月份抽样
        {
            StreamWriter writeseq = new StreamWriter(workfold + "/sampled_sample_mutlist.tsv");
            StreamWriter writemeta = new StreamWriter(workfold + "/sampled_metadata.tsv");
            int i, j, k;
            List<string> monthList = new List<string>();
            List<List<int>> sampledInd = new List<List<int>>();
            monthList.Add("2019-12");
            sampledInd.Add(new List<int>());
            for(i=2020;i<=2022;i++)
                for(j=1;j<=12;j++)
                {
                    if (j < 10)
                        monthList.Add(Convert.ToString(i) + "-0" + Convert.ToString(j));
                    else
                        monthList.Add(Convert.ToString(i) + "-" + Convert.ToString(j));
                    sampledInd.Add(new List<int>());
                }
            for(i=0;i<Sequence_List.Count;i++)
            {
                if (monthList.Contains(Sequence_List[i].Meta.CollectionMonth))
                    sampledInd[monthList.IndexOf(Sequence_List[i].Meta.CollectionMonth)].Add(i);
            }
            //开始抽样 每个滑窗抽200
            for (i = 0; i < sampledInd.Count; i++)
            {
                if(sampledInd[i].Count<200)
                {
                    for(j=0;j<sampledInd[i].Count;j++)
                    {
                        writeseq.WriteLine(Sequence_List[sampledInd[i][j]].Accession + "\t" + Sequence_List[sampledInd[i][j]].Mutation);
                        writemeta.WriteLine(Sequence_List[sampledInd[i][j]].Accession + "\t" + Sequence_List[sampledInd[i][j]].Meta.Lineage + "\t" + Sequence_List[sampledInd[i][j]].Accession + "\t" + Sequence_List[sampledInd[i][j]].Meta.CollectionDate);
                    }
                }
                else
                {
                    k = sampledInd[i].Count / 200;
                    for (j = 0; j < sampledInd[i].Count; j+=k)
                    {
                        writeseq.WriteLine(Sequence_List[sampledInd[i][j]].Accession + "\t" + Sequence_List[sampledInd[i][j]].Mutation);
                        writemeta.WriteLine(Sequence_List[sampledInd[i][j]].Accession + "\t" + Sequence_List[sampledInd[i][j]].Meta.Lineage + "\t" + Sequence_List[sampledInd[i][j]].Accession + "\t" + Sequence_List[sampledInd[i][j]].Meta.CollectionDate);
                    }
                }
            }
            writemeta.Close();
            writeseq.Close();
            return;
        }
        static void Main(string[] args)
        {
            ReadinMetadata();
            ReadinMutlist();//and find metadata
            SamplingByMonth();//根据采样月份抽样
        }
    }
}
