using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 还原lineage上下关系
{
    public class Lineage
    {
        public int NodeInd;
        public string LineageName;
        public string ParentName = "NA";
        public int ParentInd = -1;
        public List<int> ChildListInd = new List<int>();
    }
    class Program
    {
        static List<Lineage> LineageList = new List<Lineage>();
        static Dictionary<string, int> LineageMap = new Dictionary<string, int>();
        static void Main(string[] args)
        {
            StreamReader read = new StreamReader("***/lineages-website-master/data/lineages.yml");
            StreamWriter write = new StreamWriter("***/Data/lineageRelations.tsv");
            write.WriteLine("NodeInd\tName\tFatherInd\tChildNum\tChildInd");
            int i, j, k;
            string line = read.ReadLine();
            while(line!=null)
            {
                if(line[0]=='-')
                {
                    string[] line1 = line.Split(' ');
                    Lineage newl = new Lineage();
                    newl.LineageName = line1[2];
                    newl.NodeInd = LineageList.Count();
                    line = read.ReadLine();
                    while(line != null && line[0] != '-')
                    {
                        if(line.Contains("parent"))
                        {
                            string[] line2 = line.Split(' ');
                            newl.ParentName = line2[3];
                        }
                        line = read.ReadLine();
                    }
                    LineageList.Add(newl);
                    LineageMap.Add(newl.LineageName, LineageList.Count - 1);
                }
            }
            read.Close();
            for(i=0;i<LineageList.Count;i++)
            {
                if(LineageList[i].ParentName!="NA" && LineageMap.ContainsKey(LineageList[i].ParentName))
                {
                    LineageList[i].ParentInd = LineageMap[LineageList[i].ParentName];
                    LineageList[LineageList[i].ParentInd].ChildListInd.Add(LineageList[i].NodeInd);
                }
            }
            for(i = 0; i < LineageList.Count; i++)
            {
                string output = Convert.ToString(i) + "\t" + LineageList[i].LineageName + "\t";
                output += Convert.ToString(LineageList[i].ParentInd) + "\t" + Convert.ToString(LineageList[i].ChildListInd.Count()) + "\t";
                string childind = "";
                for (j = 0; j < LineageList[i].ChildListInd.Count(); j++)
                    childind += LineageList[i].ChildListInd[j] + ",";
                if (childind != "")
                    childind = childind.Substring(0, childind.Length - 1);
                output += childind;
                write.WriteLine(output);
            }
            write.Close();
        }
    }
}
