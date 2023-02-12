#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <pthread.h>
#include <mutex>
#include <map>
#include <unordered_map>
#include <time.h>
#include <math.h>
#include <algorithm>

using namespace std;

class Haplotype
{
public:
    string AAmut="";
    string AAseq = "";
    string Location = "";
    string DateStart = "";
    string DateEnd = "";
    string DateLater = "";
    string SeqProp = "";
    string StringencyInd = "";
    string Re = "";
    string SC = "";
    string PropLater = "";
    string VaccinePerMillion;
    double InfectionPerMillion = 0;
    double Pressure=0;
    double PressureCoeffcient=0;
};

class VOC
{
public:
    string Name;
    double count = 0;
    unordered_map<string, int> DefineMuts;
};

vector<vector<Haplotype>> haplotypeList;

static string refRBD_331_531 = "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST";
string AA20 = "ACDEFGHIKLMNPQRSTVWY";

unordered_map<string, double> posPressure_map;

//string split函数
vector<string> splitStr(string str, char delimiter)
{
    vector<string> r;
    string tmpstr;
    int i, j, k;
    vector<int> pointList;
    pointList.emplace_back(-1);
    for (i = 0; i < str.length(); i++)
        if (str[i] == delimiter)pointList.emplace_back(i);
    pointList.emplace_back(i);
    for (i = 1; i < pointList.size(); i++)
        r.emplace_back(str.substr(pointList[i - 1] + 1, (pointList[i] - pointList[i - 1] - 1)));
    return r;
}

int ReadinPressure(string file)
{
    ifstream read(file);
    if (!read.is_open())
    {
        cout << "/// ERROR CBPF: Can Not Open " << file << endl;
        return 1;
    }
    string line;
    int i, j, k;
    getline(read, line);
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        for (i = 3; i < line1.size(); i++)
        {
            string tmps = line1[0];
            tmps += AA20[i - 3];
            posPressure_map.insert(pair<string, double>(tmps, stof(line1[i])));
        }
    }
    read.close();
    return 0;
}

int ReadinFitness(string Fitfile, string Infectfile)
{
    ifstream readInf(Infectfile);
    ifstream read(Fitfile);
    if (!read.is_open() || !readInf.is_open())
    {
        cout << "/// ERROR CBPF: Can Not Open " << Fitfile << endl;
        cout << "/// ERROR CBPF: Can Not Open " << Infectfile << endl;
        return 1;
    }
    string line;
    int i, j, k;

    //--------------readin vaccination and infection------------
    string haploVac = "0";
    double infectPerMillion = 0;
    
    getline(readInf, line);
    getline(readInf, line);
    vector<string> lineinf1 = splitStr(line, '\t');
    if (lineinf1[1] != "")
        haploVac = lineinf1[1];
    while (getline(readInf, line))
    {
        lineinf1 = splitStr(line, '\t');
        infectPerMillion += stof(lineinf1[1]);
    }
    //----------------------------------------------------------

    k = haplotypeList.size() - 1;
    double averagePressure = 0;
    getline(read, line);
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        Haplotype newh;
        newh.AAmut = line1[0];
        newh.AAseq = refRBD_331_531;
        newh.Pressure = 0;
        newh.PressureCoeffcient = 0;
        newh.Location = line1[1];
        newh.DateStart = line1[2];
        newh.DateEnd = line1[3];
        newh.SeqProp = line1[4];
        newh.StringencyInd = line1[5];
        newh.Re = line1[6];
        newh.SC = line1[7];
        newh.DateLater = line1[8];
        newh.PropLater = line1[9];
        newh.VaccinePerMillion = haploVac;
        newh.InfectionPerMillion = infectPerMillion;

        //Replace mutation amino acids
        vector<string> mutations = splitStr(line1[0], ' ');
        if (mutations[0] != "")
        {
            vector<char> tmpAAseq;
            for (i = 0; i < refRBD_331_531.size(); i++)
                tmpAAseq.push_back(refRBD_331_531[i]);
            for (i = 0; i < mutations.size(); i++)
            {
                string pos = mutations[i].substr(1, 3);
                if(stoi(pos) >=331 && stoi(pos) <=531)
                    tmpAAseq[stoi(pos) - 331] = mutations[i][4];
            }
            newh.AAseq = "";
            for (i = 0; i < refRBD_331_531.size(); i++)
                newh.AAseq += tmpAAseq[i];
        }
        //Add up pressure from each pos
        for (i = 0; i < newh.AAseq.size(); i++)
        {
            string tmps = to_string(i + 331);
            tmps += newh.AAseq[i];
            newh.Pressure += posPressure_map[tmps];
            averagePressure += posPressure_map[tmps];
        }
        haplotypeList[k].push_back(newh);
    }

    //Calculate pressure coeffcient
    for (i = 0; i < haplotypeList[k].size(); i++)
    {
        if (haplotypeList[k].size() == 1)
            haplotypeList[k][i].PressureCoeffcient = 0;
        else
        {
            haplotypeList[k][i].PressureCoeffcient = haplotypeList[k][i].Pressure / ((averagePressure - haplotypeList[k][i].Pressure) / (haplotypeList[k].size() - 1));
            haplotypeList[k][i].PressureCoeffcient -= 1;
        }
    }
    read.close();
    return 0;
}

//based on mutation decide which voc it belongs to
//covariants.org/shared-mutations
string VOCdecide(string AAmut)
{
    vector<string> mutation = splitStr(AAmut, ' ');
    if (mutation[0] == "")
        return "WT";
    
    vector<VOC> vocList;
    
    VOC alpha;
    alpha.Name = "Alpha";
    alpha.DefineMuts.insert(pair<string, int>("N501Y", 1));
    vocList.push_back(alpha);

    VOC Beta;
    Beta.Name = "Beta";
    Beta.DefineMuts.insert(pair<string, int>("N501Y", 1));
    Beta.DefineMuts.insert(pair<string, int>("E484K", 1));
    Beta.DefineMuts.insert(pair<string, int>("K417N", 1));
    vocList.push_back(Beta);

    VOC Gamma;
    Gamma.Name = "Gamma";
    Gamma.DefineMuts.insert(pair<string, int>("N501Y", 1));
    Gamma.DefineMuts.insert(pair<string, int>("E484K", 1));
    Gamma.DefineMuts.insert(pair<string, int>("K417T", 1));
    vocList.push_back(Gamma);

    VOC Delta;
    Delta.Name = "Delta";
    Delta.DefineMuts.insert(pair<string, int>("T478K", 1));
    Delta.DefineMuts.insert(pair<string, int>("L452R", 1));
    vocList.push_back(Delta);

    VOC Kappa;
    Kappa.Name = "Kappa";
    Kappa.DefineMuts.insert(pair<string, int>("E484Q", 1));
    Kappa.DefineMuts.insert(pair<string, int>("L452R", 1));
    vocList.push_back(Kappa);

    VOC Eta;
    Eta.Name = "Eta";
    Eta.DefineMuts.insert(pair<string, int>("E484K", 1));
    vocList.push_back(Eta);

    VOC Lambda;
    Lambda.Name = "Lambda";
    Lambda.DefineMuts.insert(pair<string, int>("L452Q", 1));
    vocList.push_back(Lambda);

    VOC Omicron;
    Omicron.Name = "Omicron";
    Omicron.DefineMuts.insert(pair<string, int>("Y505H", 1));
    Omicron.DefineMuts.insert(pair<string, int>("N501Y", 1));
    Omicron.DefineMuts.insert(pair<string, int>("Q498R", 1));
    Omicron.DefineMuts.insert(pair<string, int>("E484A", 1));
    Omicron.DefineMuts.insert(pair<string, int>("T478K", 1));
    Omicron.DefineMuts.insert(pair<string, int>("S477N", 1));
    Omicron.DefineMuts.insert(pair<string, int>("N440K", 1));
    Omicron.DefineMuts.insert(pair<string, int>("K417N", 1));
    Omicron.DefineMuts.insert(pair<string, int>("S375F", 1));
    Omicron.DefineMuts.insert(pair<string, int>("S373P", 1));
    vocList.push_back(Omicron);

    int i, j, k;
    VOC tmpv;
    for (i = 0; i < mutation.size(); i++)
    {
        for (j = 0; j < vocList.size(); j++)
        {
            if (vocList[j].DefineMuts.find(mutation[i]) != vocList[j].DefineMuts.end())
                vocList[j].count++;
        }
    }
    for (i = 0; i < vocList.size(); i++)
        vocList[i].count += (double)vocList[i].count / vocList[i].DefineMuts.size() -0.01;

    for (i = 0; i < vocList.size(); i++)
    {
        for (j = i+1; j < vocList.size(); j++)
        {
            if (vocList[i].count < vocList[j].count)
            {
                tmpv = vocList[i]; vocList[i] = vocList[j]; vocList[j] = tmpv;
            }
        }
    }
    if (vocList[0].count <= 0 || vocList[0].count < vocList[0].DefineMuts.size() * 0.8)
        return "Undetermined";
    else
        return vocList[0].Name;
}

void CorrelationBuild(string fileTask, string pathPressure, string pathFitness)
{
    int i, j, k;
    string line;
    ifstream read(fileTask);
    if (!read.is_open())
    {
        cout << "/// ERROR CBPF: Can Not Open " << fileTask << endl;
        return;
    }

    string outfile = fileTask.substr(0, fileTask.size() - 4);
    outfile += ".ResultTable";
    ofstream write(outfile);
    write << "Haplo\tLocation\tDateStart\tDateEnd\tSeqProp\tStringencyInd\tRe\tSQRTSR\tSC\tDateLater\tPropLater\tInfectionPerMillion\tVaccinePerMillion\tAntibodyPressure\tPressureCoeffcient\tVOCtype\tMonth" << endl;

    getline(read, line);
    while (getline(read, line))
    {
        vector<Haplotype> tmpvh;
        haplotypeList.push_back(tmpvh);
        unordered_map<string, double> tmpmap;
        posPressure_map.swap(tmpmap);
        for (i = 0; i < line.size(); i++)
        {
            if (line[i] == ' ')line[i] = '+';
        }
        vector<string> line1 = splitStr(line, '\t');
        string DateS = line1[1].substr(0, 4) + line1[1].substr(5, 2) + line1[1].substr(8, 2);
        string DateE = line1[2].substr(0, 4) + line1[2].substr(5, 2) + line1[2].substr(8, 2);
        string DateL = line1[3].substr(0, 4) + line1[3].substr(5, 2) + line1[3].substr(8, 2);
        string Prefile = line1[0] + "_20191201_" + DateE + ".AntibodyPressure";
        string Fitfile = line1[0] + "_" + DateS + "_" + DateE + "_" + DateL + ".HaploFitness";
        string Infectfile = line1[0] + "_20191201_" + DateE  + ".InfectionList";
        //Make sure two files existed
        if (ReadinPressure(pathPressure + "/" + Prefile))
            continue;
        if (ReadinFitness(pathFitness + "/" + Fitfile, pathPressure + "/" + Infectfile))
            continue;
        //Output the result
        k = haplotypeList.size() - 1;
        for (i = 0; i < haplotypeList[k].size(); i++)
        {
            string output = haplotypeList[k][i].AAmut;
            output += "\t" + haplotypeList[k][i].Location;
            output += "\t" + haplotypeList[k][i].DateStart;
            output += "\t" + haplotypeList[k][i].DateEnd;
            output += "\t" + haplotypeList[k][i].SeqProp;
            output += "\t" + haplotypeList[k][i].StringencyInd;
            output += "\t" + haplotypeList[k][i].Re;
            output += "\t" + to_string(sqrt(stof(haplotypeList[k][i].StringencyInd) * stof(haplotypeList[k][i].Re)));
            output += "\t" + haplotypeList[k][i].SC;
            output += "\t" + haplotypeList[k][i].DateLater;
            output += "\t" + haplotypeList[k][i].PropLater;
            output += "\t" + haplotypeList[k][i].VaccinePerMillion;
            output += "\t" + to_string(haplotypeList[k][i].InfectionPerMillion);
            output += "\t" + to_string(haplotypeList[k][i].Pressure);
            output += "\t" + to_string(haplotypeList[k][i].PressureCoeffcient);
            output += "\t" + VOCdecide(haplotypeList[k][i].AAmut);
            output += "\t" + haplotypeList[k][i].DateEnd.substr(0, 6);
            write << output << endl;
        }
        cout << Infectfile << endl;
    }
    write.close();
    return;
}

/*
1. tasklist.tsv
2. path to AntibodySpectrum/Result
3. path to SelectiveCoEfficient/Result
The output will be at the same fold as the tasklist.tsv
*/
int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        cout << "/// ERROR: CorrelationBuild_PressureFitness Need 3 Input, 1 The Tasklist.tsv, 2 The path to AntibodySpectrum/Result, 3 The path to SelectiveCoEfficient/Result" << endl;
        return 1;
    }
    
    CorrelationBuild(argv[1], argv[2], argv[3]);
    return 0;
}