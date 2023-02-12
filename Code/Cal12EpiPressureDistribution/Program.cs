using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace 根据抗体中和数据计算12类抗体的压力谱
{
    class Antibody
    {
        public string AntibodyName;
        public string EpitopeGroup;
        public string Source;
        public Dictionary<string, double> IC50dic = new Dictionary<string, double>();
    }
    class Source
    {
        public string SourceName;
        public List<Antibody> AntibodayList = new List<Antibody>();
    }
    class AntibodyPressureSpectrum
    {
        public string Source; // WT convalescents, WT vaccinees
        public string VOC;// WT
        public List<int> AntibodyCount = new List<int>();//抗体个数分布
        public List<double> AntibodyPressure = new List<double>();//抗体个数*中和力分布
    }
    class Program
    {
        static Dictionary<string, AntibodyPressureSpectrum> SpectrumList = new Dictionary<string, AntibodyPressureSpectrum>();
        static Dictionary<string, Source> SourceList = new Dictionary<string, Source>();
        static string[] Group12 = "A,B,C,D1,D2,E1,E2.1,E2.2,E3,F1,F2,F3".Split(',');
        static List<string> GroupList = new List<string> (Group12.ToArray());
        static string workfold = "***/Usher1128";
        static void ReadinNeutral()
        {
            int i, j, k;
            StreamReader read = new StreamReader(workfold + "/Data/NeutralWTBA125_Cross.txt");
            string line = read.ReadLine();
            List<string> title = new List<string>(line.Split('\t'));
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                string[] line2 = line1[3].Split(' ');
                if(!SourceList.ContainsKey(line1[3]))
                {
                    Source newa = new Source();
                    newa.SourceName = line1[3];
                    SourceList.Add(line1[3], newa);
                }
                Antibody newb = new Antibody();
                newb.AntibodyName = line1[1];
                newb.EpitopeGroup = line1[2];
                newb.Source = line1[3];
                for(i=4;i<=8;i++)
                {
                    string IC50S = line1[i];
                    if (IC50S != "--")
                    {
                        double IC50 = Convert.ToDouble(IC50S);
                        if (IC50 > 1) IC50 = 1;
                        if (IC50 < 0.0005) IC50 = 0.0005;
                        double Neu = -1 * Math.Log10(IC50);
                        newb.IC50dic.Add(title[i], Neu);
                    }
                }
                SourceList[line1[3]].AntibodayList.Add(newb);
                line = read.ReadLine();
            }
            return;
        }
        static void CalculateSpectrum(List<string> sources, string clade)//抗体的来源表，和中和针对的毒株
        {
            int i, j, k;
            AntibodyPressureSpectrum newa = new AntibodyPressureSpectrum();
            for (i = 0; i < sources.Count; i++)
                newa.Source += sources[i] + ",";
            newa.Source = newa.Source.Substring(0, newa.Source.Length - 1);
            newa.VOC = clade;
            for(i=0;i<12;i++)
            {
                newa.AntibodyCount.Add(0);
                newa.AntibodyPressure.Add(0);
            }
            foreach(string val in SourceList.Keys)
            {
                if(sources.Contains(val))
                {
                    for(i=0;i<SourceList[val].AntibodayList.Count;i++)
                    {
                        newa.AntibodyCount[GroupList.IndexOf(SourceList[val].AntibodayList[i].EpitopeGroup)]++;
                        if(SourceList[val].AntibodayList[i].IC50dic.ContainsKey(clade))
                            newa.AntibodyPressure[GroupList.IndexOf(SourceList[val].AntibodayList[i].EpitopeGroup)] += SourceList[val].AntibodayList[i].IC50dic[clade];
                    }
                }
            }
            SpectrumList.Add(newa.Source + "2" + newa.VOC, newa);
        }
        static void WriteResult()
        {
            StreamWriter writeCount = new StreamWriter(workfold + "/Data/AntibodySpectrum12_Count.txt");
            StreamWriter writePressure = new StreamWriter(workfold + "/Data/AntibodySpectrum12_logIC50.txt");
            int i, j, k;
            writeCount.WriteLine("Clade\tA\tB\tC\tD1\tD2\tE1\tE2.1\tE2.2\tE3\tF1\tF2\tF3");
            writePressure.WriteLine("Clade\tA\tB\tC\tD1\tD2\tE1\tE2.1\tE2.2\tE3\tF1\tF2\tF3");
            foreach (string val in SpectrumList.Keys)
            {
                string output = val;
                for (i = 0; i < 12; i++)
                    output += "\t" + SpectrumList[val].AntibodyCount[i];
                writeCount.WriteLine(output);

                output = val;
                for (i = 0; i < 12; i++)
                    output += "\t" + SpectrumList[val].AntibodyPressure[i];
                writePressure.WriteLine(output);
            }
            writeCount.Close();
            writePressure.Close();
        }
        static void Main(string[] args)
        {
            ReadinNeutral();

            //-------------------Calculate----------------------

            CalculateSpectrum(new List<string>("WT convalescents+WT vaccinees".Split('+').ToArray()), "WT");
            CalculateSpectrum(new List<string>("BA.1 convalescents".Split('+').ToArray()), "BA.1");
            CalculateSpectrum(new List<string>("BA.2 convalescents".Split('+').ToArray()), "BA.2");
            CalculateSpectrum(new List<string>("BA.5 convalescents".Split('+').ToArray()), "BA.5");


            WriteResult();
            return;
        }
    }
}
