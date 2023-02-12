using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Json2Mutation_Csharp
{
    public class Node
    {
        public string name;
        public bool IsLeaf = true;
        public List<string> Nuc_mutation = new List<string>();
        public List<string> Total_Nuc_mutation = new List<string>();
        public string CollectionDate = null;
        public string Location;
        public string Lineage;
        public string NextstrainClade;
        public string accessionID = "NA";
        public List<int> ChildernIndList = new List<int>();//子节点编号
        public int SameSequenceNumber = 0;//和这个节点一样的叶子有多少
    }
    public class Metadata
    {
        public string name;
        public string CollectionDate;
        public string Location;
        public string Lineage;
        public string NextstrainClade = "Other";
    }
    class Program
    {
        static List<string> JsonLine = new List<string>(234217729);
        static List<Node> NodesList = new List<Node>();
        static Dictionary<string, Metadata> MetadataDic = new Dictionary<string, Metadata>();
        static string Workfold = "***";
        static bool DateCompare(string A, string B)//A>B true else false
        {
            int i;
            for (i = 0; i < A.Length && i < B.Length; i++)
            {
                if (A[i] > B[i])
                    return true;
            }
            return false;
        }
        static List<string> FindEarliestDate(List<string> DateList, List<string> LocationList)//从list中找出最早采样时间
        {
            string tmps;
            int i, j;
            for (i = 0; i < DateList.Count; i++)
            {
                for (j = i + 1; j < DateList.Count; j++)
                {
                    if (DateCompare(DateList[i], DateList[j]))
                    {
                        tmps = DateList[i]; DateList[i] = DateList[j]; DateList[j] = tmps;
                        tmps = LocationList[i]; LocationList[i] = LocationList[j]; LocationList[j] = tmps;
                    }
                }
            }
            List<string> RT = new List<string>();
            RT.Add(DateList[0]);
            RT.Add(LocationList[0]);
            return RT;
        }
        static void AssignDateToNodes()//给树中间的节点分配采样时间
        {
            /*
             * 树上存在完全一致的序列，这些序列难以判断发生突变，因为他们到前一个节点的距离都是0
             * 前一个节点，也是树中间节点的突变，可以代表这些序列的突变
             * 用这些序列的最早采样时间，作为前一个节点的采样时间
             */
            int i, j, k;
            for (i = 0; i < NodesList.Count; i++)
            {
                if (!NodesList[i].IsLeaf)
                {
                    List<string> DateList = new List<string>();
                    List<string> LocationList = new List<string>();
                    for (j = 0; j < NodesList[i].ChildernIndList.Count; j++)
                    {
                        //该中间节点有叶子节点的子代且子代没有额外突变
                        if (NodesList[NodesList[i].ChildernIndList[j]].IsLeaf && NodesList[NodesList[i].ChildernIndList[j]].Nuc_mutation.Count == 0)
                            if (NodesList[NodesList[i].ChildernIndList[j]].CollectionDate != null)
                            {
                                DateList.Add(NodesList[NodesList[i].ChildernIndList[j]].CollectionDate);
                                LocationList.Add(NodesList[NodesList[i].ChildernIndList[j]].Location);
                            }
                        NodesList[i].SameSequenceNumber = DateList.Count;
                    }
                    if (DateList.Count > 0)
                    {
                        List<string> DateLocation = FindEarliestDate(DateList, LocationList);
                        NodesList[i].CollectionDate = DateLocation[0];
                        NodesList[i].Location = DateLocation[1];
                        NodesList[i].Lineage = NodesList[NodesList[i].ChildernIndList[0]].Lineage;
                        NodesList[i].NextstrainClade = NodesList[NodesList[i].ChildernIndList[0]].NextstrainClade;
                    }

                }
            }

            return;
        }
        static void LeafFindMetadata()//给叶子节点找metadata
        {
            int i, j, k;
            for (i = 0; i < NodesList.Count; i++)
            {
                if (NodesList[i].IsLeaf)
                {
                    string[] name1 = NodesList[i].name.Split('|');
                    if(MetadataDic.ContainsKey(NodesList[i].name))
                    {
                        NodesList[i].accessionID = NodesList[i].name;
                        NodesList[i].CollectionDate = MetadataDic[NodesList[i].name].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = MetadataDic[NodesList[i].name].Lineage;
                        NodesList[i].NextstrainClade = MetadataDic[NodesList[i].name].NextstrainClade;
                        NodesList[i].Location = MetadataDic[NodesList[i].name].Location;
                    }
                    else
                        if (MetadataDic.ContainsKey(name1[1]))//搜索accession id
                    {
                        NodesList[i].accessionID = name1[1];
                        NodesList[i].CollectionDate = MetadataDic[name1[1]].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = MetadataDic[name1[1]].Lineage;
                        NodesList[i].NextstrainClade = MetadataDic[name1[1]].NextstrainClade;
                        NodesList[i].Location = MetadataDic[name1[1]].Location;
                    }
                    else
                        if (MetadataDic.ContainsKey(name1[0]))//找不到看看文件名有没有
                    {
                        NodesList[i].accessionID = "Seq_Name";
                        NodesList[i].CollectionDate = MetadataDic[name1[0]].CollectionDate;
                        if (NodesList[i].CollectionDate.Length == 4)
                            NodesList[i].CollectionDate += "-99-99";
                        if (NodesList[i].CollectionDate.Length == 7)
                            NodesList[i].CollectionDate += "-99";
                        NodesList[i].Lineage = MetadataDic[name1[0]].Lineage;
                        NodesList[i].NextstrainClade = MetadataDic[name1[0]].NextstrainClade;
                        NodesList[i].Location = MetadataDic[name1[0]].Location;
                    }
                    else//确实找不到，看看文件名里写没写
                    {
                        string[] name11 = name1[0].Split('/');
                        if (name11.Length == 3)
                        {
                            NodesList[i].CollectionDate = name1[name1.Length - 1];
                            if (NodesList[i].CollectionDate.Length == 4)
                                NodesList[i].CollectionDate += "-99-99";
                            if (NodesList[i].CollectionDate.Length == 7)
                                NodesList[i].CollectionDate += "-99";
                            NodesList[i].Location = name11[0];
                        }
                    }
                }
            }
            return;
        }
        static void TreeBuild(string jsonLine, int fatherInd)
        {
            Node newNode = new Node();
            int i, j, k;
            int kuohao = 0;
            //Console.WriteLine(jsonLine);

            //---------分割序列---------
            List<int> SplitPos = new List<int>();//直系逗号的位置
            List<int> ChildSplitPos = new List<int>();//孩子分割位置
            SplitPos.Add(0);
            for (i = 0; i < jsonLine.Length; i++)
            {
                if (kuohao == 0 && jsonLine[i] == ',')
                    SplitPos.Add(i + 1);
                if (jsonLine[i] == '{' || jsonLine[i] == '[') kuohao++;
                if (jsonLine[i] == '}' || jsonLine[i] == ']') kuohao--;
            }
            for (i = 0; i < SplitPos.Count; i++)//切割序列，找突变和子代
            {
                j = SplitPos[i] + 1;
                while (jsonLine[j] != '\"') j++;
                string tag = jsonLine.Substring(SplitPos[i] + 1, j - SplitPos[i] - 1);
                if (tag == "branch_attrs")
                {
                    string lineba1 = jsonLine.Substring(SplitPos[i] + 1, SplitPos[i + 1] - SplitPos[i] - 2);
                    string lineba2 = lineba1.Replace('}', '{');
                    string[] lineba3 = lineba2.Split('{');
                    if (lineba3[4] != "\"nuc\":[]")
                    {
                        string lineba4 = lineba3[4].Substring(7, lineba3[4].Length - 8);
                        string[] lineba5 = lineba4.Split(',');
                        for (k = 0; k < lineba5.Length; k++)
                        {
                            newNode.Nuc_mutation.Add(lineba5[k].Substring(1, lineba5[k].Length - 2));
                        }
                    }
                }
                if (tag == "children")
                {
                    ChildSplitPos.Add(SplitPos[i] + 11);
                    kuohao = 0;
                    for (k = SplitPos[i] + 12; k < SplitPos[i + 1]; k++)
                    {
                        if (kuohao == 0 && jsonLine[k] == ',')
                            ChildSplitPos.Add(k);
                        if (jsonLine[k] == '{' || jsonLine[k] == '[') kuohao++;
                        if (jsonLine[k] == '}' || jsonLine[k] == ']') kuohao--;
                    }
                    ChildSplitPos.Add(SplitPos[i + 1] - 2);
                }
                if (tag == "node_attrs")
                {

                }
                if (tag == "name")
                {
                    string lineba1 = jsonLine.Substring(SplitPos[i] + 1, SplitPos[i + 1] - SplitPos[i] - 2);
                    newNode.name = lineba1.Substring(7, lineba1.Length - 8);
                }
            }

            //---------从父节点继承突变----------
            newNode.Total_Nuc_mutation = new List<string>(newNode.Nuc_mutation);
            if (fatherInd != -1)
            {
                for (i = 0; i < NodesList[fatherInd].Total_Nuc_mutation.Count; i++)
                    newNode.Total_Nuc_mutation.Add(NodesList[fatherInd].Total_Nuc_mutation[i]);
            }

            //--------递归-----------
            if (ChildSplitPos.Count <= 2)//没有子代
            {
                newNode.IsLeaf = true;
                NodesList[fatherInd].ChildernIndList.Add(NodesList.Count);
                NodesList.Add(newNode);
            }
            else//递归子代
            {
                newNode.IsLeaf = false;
                if (fatherInd != -1)
                    NodesList[fatherInd].ChildernIndList.Add(NodesList.Count);
                NodesList.Add(newNode);
                for (k = 0; k < ChildSplitPos.Count - 1; k++)
                {
                    //Console.WriteLine(jsonLine.Substring(ChildSplitPos[k] + 1, ChildSplitPos[k + 1] - ChildSplitPos[k] - 1));
                    TreeBuild(jsonLine.Substring(ChildSplitPos[k] + 2, ChildSplitPos[k + 1] - ChildSplitPos[k] - 3), NodesList.Count - 1);
                }
            }
            return;
        }
        static int TreeBuild2(int lineId,int fatherInd)
        {
            int i, j, k;
            int thisNode = -1;
            int kuohao = 0;
            Node newnode = new Node();
            //继承父节点突变
            if(fatherInd!=-1)
                for (k = 0; k < NodesList[fatherInd].Total_Nuc_mutation.Count; k++)
                    newnode.Total_Nuc_mutation.Add(NodesList[fatherInd].Total_Nuc_mutation[k]);
            while (lineId < JsonLine.Count())
            {
                if (JsonLine[lineId].IndexOf('{') != -1) kuohao++;
                if (JsonLine[lineId].IndexOf('[') != -1) kuohao++;
                if (JsonLine[lineId].IndexOf('}') != -1) kuohao--;
                if (JsonLine[lineId].IndexOf(']') != -1) kuohao--;
                if (kuohao == 0) 
                    break;
                //------------节点名称-----------
                if (JsonLine[lineId].Contains("nuc mutations"))
                {
                    string[] label = JsonLine[lineId].Split(':');
                    newnode.name = label[2].Substring(1, label[2].Length - 3);
                }
                //------------节点突变-----------
                if (JsonLine[lineId] == "\"nuc\":[")
                {
                    lineId++; kuohao--;
                    string[] mut1 = JsonLine[lineId].Split('\"');
                    for (i = 0; i < mut1.Count(); i++)
                        if (mut1[i].Length > 2)
                        {
                            newnode.Nuc_mutation.Add(mut1[i]);
                            newnode.Total_Nuc_mutation.Add(mut1[i]);
                        }
                }
                //------------递归子代-----------
                if (JsonLine[lineId] == "\"children\":[")
                {
                    newnode.IsLeaf = false;
                    NodesList.Add(newnode);
                    thisNode = NodesList.Count() - 1;
                    lineId++; 
                    while (JsonLine[lineId]!="]")
                        lineId = TreeBuild2(lineId, thisNode);
                    if (JsonLine[lineId] == "]") kuohao--;
                }
                lineId++;
            }
            //--------------和父节点建立联系--------------
            if(newnode.IsLeaf)
            {
                NodesList.Add(newnode);
                thisNode = NodesList.Count() - 1;
            }
            if (fatherInd != -1)
            {
                //for (k = 0; k < NodesList[fatherInd].Total_Nuc_mutation.Count; k++)
                //    NodesList[thisNode].Total_Nuc_mutation.Add(NodesList[fatherInd].Total_Nuc_mutation[k]);
                NodesList[fatherInd].ChildernIndList.Add(thisNode);
            }

            return lineId+1;
        }
        static void MutationEventCal(string jsfile)
        {
            //-----------------------读取json转换输出----------------------
            int i, j, k;
            FileStream fs = new FileStream(jsfile, FileMode.Open);
            TextReader read = new StreamReader(fs);
            StreamWriter writel = new StreamWriter(jsfile + ".line");
            var clen = 1024 * 1024;
            var buffer = new Char[clen];
            var count = read.Read(buffer, 0, clen);
            while (count > 0)
            {
                var str = new string(buffer, 0, count);
                StringBuilder strBuilder = new StringBuilder();
                strBuilder.Append(str);
                strBuilder.Replace("{", "{\n");
                strBuilder.Replace("[", "[\n");
                strBuilder.Replace("}", "}\n");
                strBuilder.Replace("]", "]\n");
                strBuilder.Replace(",", "");
                writel.Write(strBuilder.ToString());
                count = read.Read(buffer, 0, clen);
            }
            read.Close();
            writel.Close();
            read.Dispose();
            Console.WriteLine("Json Reform Done.");

            //------------------------读取新格式json----------------
            StreamReader readline = new StreamReader(jsfile + ".line");
            string line = readline.ReadLine();
            while(line!=null)
            {
                JsonLine.Add(line);
                line = readline.ReadLine();
            }
            Console.WriteLine("Read Json Line Done.");
            readline.Close();
            //------------------------构建树结构---------------------
            k = 0;
            while (JsonLine[k] != "\"tree\":{") k++;
            k++;
            TreeBuild2(k, -1);
            Console.WriteLine("Tree Build Done.");

            //-------------------------------------------------
            /*int kuohao = 0;//处于括号的层数
            int start = 0;
            int end = line.Length - 1;
            while (start != line.Length)//locate the start pos of the tree
            {
                if (line.Substring(start, 7) == "\"tree\":" && kuohao == 1)
                    break;
                if (line[start] == '{') kuohao++;
                if (line[start] == '}') kuohao--;
                start++;
            }
            kuohao = 0;
            while (end != 0)//locate the ending position
            {
                if (end < line.Length - 7 && line.Substring(end, 7) == "\"name\":" && kuohao == 2)
                    break;
                if (line[end] == '{') kuohao--;
                if (line[end] == '}') kuohao++;
                end--;
            }
            start += 21;
            end -= 4;
            TreeBuild(line.Substring(start, end - start + 1), -1);*/


            //---------------------------------------------------------
            LeafFindMetadata();//给叶子节点找metadata
            Console.WriteLine("Finding Metadata.");

            AssignDateToNodes();//给树中间的节点分配采样时间
            Console.WriteLine("Assign Sampling Date.");


            //-------------------------输出结果------------------------
            Console.WriteLine("Writing Result ......");
            StreamWriter write = new StreamWriter(jsfile + ".mutevent");
            write.Write("Sample\tAccessionID\tOriginalMut\tNewMut\tLineage\tCollectionDate\tLocation\tSameSeqNumber\tNextstrainClade\n");
            for (i = 0; i < NodesList.Count; i++)
            {
                if ((NodesList[i].IsLeaf || NodesList[i].SameSequenceNumber > 0) && NodesList[i].Nuc_mutation.Count == 1)
                {
                    string output = NodesList[i].name;
                    output += "\t" + NodesList[i].accessionID;
                    string mutation = "";
                    for (j = 0; j < NodesList[i].Total_Nuc_mutation.Count; j++)
                        mutation += NodesList[i].Total_Nuc_mutation[j] + ",";
                    if (mutation != "")
                        mutation = mutation.Substring(0, mutation.Length - 1);
                    output += "\t" + mutation;
                    output += "\t" + NodesList[i].Nuc_mutation[0];
                    output += "\t" + NodesList[i].Lineage;
                    output += "\t" + NodesList[i].CollectionDate;
                    output += "\t" + NodesList[i].Location;
                    output += "\t" + Convert.ToString(NodesList[i].SameSequenceNumber);
                    output += "\t" + NodesList[i].NextstrainClade;
                    write.WriteLine(output);
                }
                else
                    if ((NodesList[i].IsLeaf || NodesList[i].SameSequenceNumber > 0) && NodesList[i].Nuc_mutation.Count == 2)
                {
                    string output = NodesList[i].name;
                    output += "\t" + NodesList[i].accessionID;
                    string mutation = "";
                    for (j = 0; j < NodesList[i].Total_Nuc_mutation.Count; j++)
                        mutation += NodesList[i].Total_Nuc_mutation[j] + ",";
                    if (mutation != "")
                        mutation = mutation.Substring(0, mutation.Length - 1);
                    output += "\t" + mutation;
                    string output1 = output, output2 = output;
                    output1 += "\t" + NodesList[i].Nuc_mutation[0];
                    output1 += "\t" + NodesList[i].Lineage;
                    output1 += "\t" + NodesList[i].CollectionDate;
                    output1 += "\t" + NodesList[i].Location;
                    output1 += "\t" + Convert.ToString(NodesList[i].SameSequenceNumber);
                    output1 += "\t" + NodesList[i].NextstrainClade;
                    write.WriteLine(output1);
                    output2 += "\t" + NodesList[i].Nuc_mutation[1];
                    output2 += "\t" + NodesList[i].Lineage;
                    output2 += "\t" + NodesList[i].CollectionDate;
                    output2 += "\t" + NodesList[i].Location;
                    output2 += "\t" + Convert.ToString(NodesList[i].SameSequenceNumber);
                    output2 += "\t" + NodesList[i].NextstrainClade;
                    write.WriteLine(output2);
                }
            }
            write.Close();
            return;
        }
        static void ReadMetadata()
        {
            int i, j, k;
            string file = "/metadata.tsv";
            StreamReader read;
            string line;

            /*read = new StreamReader(file);
            line = read.ReadLine();
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Metadata newm = new Metadata();
                newm.name = line1[0];
                newm.Lineage = line1[4];
                newm.CollectionDate = line1[10];
                string[] line2 = line1[11].Split('/');
                k = line2[0].Length - 1;
                while (line2[0][k] == ' ') k--;
                newm.Location = line2[0].Substring(0, k + 1);
                if (line1[1] != "null" && line1[1] != "" && !MetadataDic.ContainsKey(line1[1]))
                    MetadataDic.Add(line1[1], newm);
                if (line1[3] != "null" && line1[3] != "" && line1[3] != " " && !MetadataDic.ContainsKey(line1[3]))
                    MetadataDic.Add(line1[3], newm);
                if(!MetadataDic.ContainsKey(line1[0]))
                    MetadataDic.Add(line1[0], newm);
                line = read.ReadLine();
            }
            read.Close();*/
            read = new StreamReader(Workfold + "/Usher1128/public-latest.metadata.tsv");
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                Metadata newm = new Metadata();
                newm.name = line1[0];
                newm.Lineage = line1[8];
                newm.NextstrainClade = line1[7];
                newm.CollectionDate = line1[2];
                newm.Location = line1[3];
                if (!MetadataDic.ContainsKey(line1[0])) 
                    MetadataDic.Add(line1[0], newm);
                if( !MetadataDic.ContainsKey(line1[0]))
                    MetadataDic.Add(line1[1], newm);
                line = read.ReadLine();
            }
            Console.WriteLine("Read Metadata Done");
        }
        static void WriteNodeInfo(string jsfile)
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(jsfile + ".nodeinfo");
            write.WriteLine("NodeName\tEarliestDate\tEarliestLocation\tLineage");
            for(i=0;i<NodesList.Count;i++)
            {
                string output = NodesList[i].name;
                output += "\t" + NodesList[i].CollectionDate;
                output += "\t" + NodesList[i].Location;
                output += "\t" + NodesList[i].Lineage;
                write.WriteLine(output);
            }
            write.Close();
            return;
        }
        static void Main(string[] args)
        {
            ReadMetadata();
            MutationEventCal(Workfold + "/Usher1128/MAT1128.json");
            WriteNodeInfo(Workfold + "/Usher1128/MAT1128.json");
            return;
        }
    }
}
