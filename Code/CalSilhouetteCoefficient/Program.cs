using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 计算轮廓系数_评价聚类好坏
{
    public class Sample
    {
        public string Name;
        public string Group;
        public List<double> FeaturesList = new List<double>();
    }
    class Program
    {
        static string file;
        static List<Sample> sampleList = new List<Sample>();
        static double EuclideanDistanceCalculator(List<double> listA, List<double> listB)//计算欧式距离
        {
            int i, j;
            double distance = 0;
            for(i=0;i<listA.Count;i++)
            {
                distance += (listA[i] - listB[i]) * (listA[i] - listB[i]);
            }
            return Math.Sqrt(distance);
        }
        static double CalculatorSilhouetteCoefficient()
        {
            int i, j,k;
            List<string> groupList = new List<string>();
            for (i = 0; i < sampleList.Count; i++)
                if (!groupList.Contains(sampleList[i].Group))
                    groupList.Add(sampleList[i].Group);
            List<double> samplesc = new List<double>();
            for(i=0;i<sampleList.Count;i++)
            {
                List<double> inclusterD = new List<double>();
                List<List<double>> betweenclusterD = new List<List<double>> ();
                for(j=0;j<groupList.Count;j++)
                {
                    List<double> newd = new List<double>();
                    betweenclusterD.Add(newd);
                }
                for(j=0;j<sampleList.Count;j++)
                {
                    if(i!=j)
                    {
                        if(sampleList[i].Group==sampleList[j].Group)//组内
                        {
                            inclusterD.Add(EuclideanDistanceCalculator(sampleList[i].FeaturesList, sampleList[j].FeaturesList));
                        }
                        else//组间
                        {
                            for(k=0;k<groupList.Count;k++)
                                if(groupList[k]==sampleList[j].Group)
                                {
                                    betweenclusterD[k].Add(EuclideanDistanceCalculator(sampleList[i].FeaturesList, sampleList[j].FeaturesList));
                                }
                        }
                    }
                }
                if (inclusterD.Count == 0)
                    samplesc.Add(0);
                else
                {
                    double a=0, b=9999999;
                    for (j = 0; j < inclusterD.Count; j++)
                        a += inclusterD[j];
                    a /= inclusterD.Count();

                    for(k=0;k<betweenclusterD.Count;k++)
                    {
                        double sumsc = 0;
                        for (j = 0; j < betweenclusterD[k].Count; j++)
                            sumsc += betweenclusterD[k][j];
                        sumsc /= betweenclusterD[k].Count;
                        if (sumsc < b) b = sumsc;
                    }

                    if (a > b)
                        samplesc.Add((b - a) / a);
                    else
                        samplesc.Add((b - a) / b);
                }
            }
            double mean = 0;
            for (i = 0; i < samplesc.Count; i++)
                mean += samplesc[i];
            mean /= samplesc.Count;
            
            StreamWriter write = new StreamWriter(file + ".SilhouetteCoefficient");
            write.WriteLine(Convert.ToString(mean));
            for(k=0;k<groupList.Count;k++)
            {
                mean = 0; int n = 0;
                for (i = 0; i < samplesc.Count; i++)
                    if (sampleList[i].Group == groupList[k])
                    {
                        n++;
                        mean += samplesc[i];
                    }
                write.WriteLine(groupList[k] + "\t" + Convert.ToString(mean/n));
            }
            write.Close();
            return mean;
        }
        static void Main(string[] args) //输入一个矩阵 行是样本 列是特征 第一列是分组信息 第二列是样本名称
        {
            /*input is a matrix 
            row as the sample 
            column as the feature. 
            The first column is the grouping information. 
            The second column is the sample name*/
            file = "/Silhouette_Coefficient/1128CladeHeatmap.txt";
            int i, j, k;
            StreamReader read = new StreamReader(file);
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                Sample news = new Sample();
                news.Group = line1[0];
                news.Name = line1[1];
                for (i = 2; i < line1.Length; i++)
                    news.FeaturesList.Add(Convert.ToDouble(line1[i]));
                sampleList.Add(news);
                line = read.ReadLine();
            }
            read.Close();
            CalculatorSilhouetteCoefficient();
            return;
        }
    }
}
