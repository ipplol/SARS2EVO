#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mutex>
#include <map>
#include <unordered_map>
#include <time.h>
#include <math.h>
#include <algorithm>

using namespace std;

class PressureMatrix
{
public:
    string VOCname;
    double casePerMillion;
    unordered_map<string, double> matrix;
};
vector<PressureMatrix> pressureMatrixList;
unordered_map<string, double> finalPressureMatrix;
unordered_map<string, int> VOCmap;
unordered_map<string, string> lineage2VocMap;//deside which spectrum is used for which lineage

class Antibody
{
public:
    string name;
    string source;//WT,SARS,BA1
    string type;//12 epitope
    double neutral = 0;
    double escapeScoreSum = 0;
    vector<double> posAverage;
    unordered_map<string, double> escapeScore;
};
vector<Antibody> antibodyList;
unordered_map<string, int> antibodyMap;

string AA20 = "ACDEFGHIKLMNPQRSTVWY";
string RefRBD331_531 = "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST";

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

//input the historical infection situation
void ReadinInfectionList(string file)
{
    int i, j, k;
    ifstream read(file);
    if (!read.is_open())
    {
        printf("/// ERROR CalAntibodyPresuure: Can NOT open %s\n", file);
        return;
    }
    string line;
    getline(read, line);//header
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        PressureMatrix tmppm;
        //Replace lineage name with spectrum name, 用指定的抗体谱对应变异株
        //Unable spectrum and vaccined replaced with Wild type 
        if (lineage2VocMap.find(line1[0]) != lineage2VocMap.end())
            tmppm.VOCname = lineage2VocMap[line1[0]];
        else
            tmppm.VOCname = "WT";

        tmppm.casePerMillion = stof(line1[1]);
        for (i = 331; i <= 531; i++)
        {
            for (j = 0; j < 20; j++)
            {
                string tmps = to_string(i) + AA20[j];
                tmppm.matrix.insert(pair<string, double>(tmps, 0));
            }
        }
        VOCmap.insert(pair<string, int>(line1[0], pressureMatrixList.size()));
        pressureMatrixList.push_back(tmppm);
    }
    for (i = 331; i <= 531; i++)
    {
        for (j = 0; j < 20; j++)
        {
            string tmps = to_string(i) + AA20[j];
            finalPressureMatrix.insert(pair<string, double>(tmps, 0));
        }
    }
    read.close();
    return;
}

//Readin antibody
void ReadinAntibody(string FileFold)
{
    string NeutralFile = FileFold + "/Neutral_info.csv";
    string EscapeFile = FileFold + "/use_abs_res.csv";
    string VariantTable = FileFold + "/VariantTable.tsv";
    int i, j, k;
    ifstream read(NeutralFile);
    if (!read.is_open())
    {
        printf("/// ERROR CalAntibodyPressure: Can NOT open %s\n", NeutralFile);
        return;
    }
    string line;
    getline(read, line);
    vector<string> strainList = splitStr(line, ',');
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, ',');
        Antibody tmpa;
        tmpa.name = line1[0];
        tmpa.source = line1[line1.size() - 2];
        tmpa.type = line1[line1.size() - 1];
        for (i = 0; i < strainList.size(); i++)
            if (strainList[i] == tmpa.source)
                break;
        if (line1[i] != "NA")
        {
            tmpa.neutral = (-1 * log10(stof(line1[i])) + 1) / 5;
            //tmpa.neutral = 10 / stof(line1[i]);
        }
        antibodyMap.insert(pair<string, int>(line1[0], antibodyList.size()));
        antibodyList.push_back(tmpa);
    }
    read.close();

    //Escape score
    ifstream readE(EscapeFile);
    if (!readE.is_open())
    {
        printf("/// ERROR CalAntibodyPressure: Can NOT open %s\n", EscapeFile);
        return;
    }
    getline(readE, line);
    while (getline(readE, line))
    {
        vector<string> line1 = splitStr(line, ',');
        string mut = line1[1] + line1[2];
        if (antibodyMap.find(line1[0]) != antibodyMap.end())
        {
            antibodyList[antibodyMap[line1[0]]].escapeScore[mut] = stof(line1[3]);
            //antibodyList[antibodyMap[line1[0]]].escapeScore[mut] = 1;
        }
    }
    readE.close();

    //感染lineage和抗体谱的对应关系
    ifstream readVT(VariantTable);
    if (!readVT.is_open())
    {
        printf("/// ERROR CalAntibodyPressure: Can NOT open %s\n", VariantTable);
        return;
    }
    while (getline(readVT, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        lineage2VocMap.insert(pair<string, string>(line1[0], line1[1]));
    }
    readVT.close();
    return;
}

void CalPressureMatrix(string outputPath)
{
    int i, j, k;
    k = outputPath.size() - 1;
    while (outputPath[k] != '.')k--;
    outputPath = outputPath.substr(0, k);
    outputPath += ".AntibodyPressure";

    //filter the antibody escape, 
    //to a pos, sum all amino acids
    //to an antibody, average all pos
    //only keep pos that 2 fold high than average
    for (k = 0; k < antibodyList.size(); k++)
    {
        double allPosAverage = 0;
        for (i = 331; i <= 531; i++)
        {
            double tmpd = 0;
            for (j = 0; j < 20; j++)
            {
                string tmps = to_string(i) + AA20[j];
                tmpd += antibodyList[k].escapeScore[tmps];
            }
            allPosAverage += tmpd;
            antibodyList[k].posAverage.push_back(tmpd);
        }
        antibodyList[k].escapeScoreSum = allPosAverage;
        allPosAverage /= 200;
        for (i = 331; i <= 531; i++)
        {
            if (antibodyList[k].posAverage[i - 331] < allPosAverage * 2)
            {
                for (j = 0; j < 20; j++)
                {
                    string tmps = to_string(i) + AA20[j];
                    antibodyList[k].escapeScore[tmps] = 1;
                }
            }
        }
    }

    //calculate pressure matrix for each strain
    //by adding up all antibodys of that strain
    for (i = 0; i < pressureMatrixList.size(); i++)
    {
        k = 0;
        for (j = 0; j < antibodyList.size(); j++)
        {
            if (antibodyList[j].source == pressureMatrixList[i].VOCname)
            {
                k++;
                for (auto val = antibodyList[j].escapeScore.begin(); val != antibodyList[j].escapeScore.end(); val++)
                {
                    //assign neutral force to each pos by escape score
                    string pos = (val->first).substr(0, 3);
                    double percent = antibodyList[j].posAverage[stoi(pos) - 331] / antibodyList[j].escapeScoreSum;
                    if (antibodyList[j].escapeScoreSum == 0) percent = 0;

                    //(1-Escape) * Neutral * percenage
                    pressureMatrixList[i].matrix[val->first] += (1 - val->second) * antibodyList[j].neutral * percent;
                }
            }
        }
        cout << pressureMatrixList[i].VOCname << "  Antibody number: " << k << endl;
    }

    //add up all strain, with infection percentage, to obtain the final pressure matrix
    for (i = 0; i < pressureMatrixList.size(); i++)
    {
        for (auto val = pressureMatrixList[i].matrix.begin(); val != pressureMatrixList[i].matrix.end(); val++)
        {
            finalPressureMatrix[val->first] += pressureMatrixList[i].casePerMillion * val->second / 1000000;
            //finalPressureMatrix[val->first] += val->second;
        }
    }

    ofstream write(outputPath);
    if (!write.is_open())
    {
        printf("/// ERROR CalAntibodyPresuure: Can NOT open %s\n", outputPath);
        return;
    }
    string output = "AntibodyPressureMatrix\tRefNuc\tRef";
    for (j = 0; j < 20; j++)
    {
        string tmps = "";
        tmps += AA20[j];
        output += "\t" + tmps;
    }
    write << output << endl;
    for (i = 331; i <= 531; i++)
    {
        output = to_string(i);
        string tmpsr = output + RefRBD331_531[i - 331];
        output += "\t";
        output += RefRBD331_531[i - 331];
        output += "\t" + to_string(finalPressureMatrix[tmpsr]);
        for (j = 0; j < 20; j++)
        {
            string tmps = to_string(i) + AA20[j];
            output += "\t" + to_string(finalPressureMatrix[tmps]);
        }
        write << output << endl;
    }
    write.close();
}

int main(int args, char* argv[])//3个输入 InfectionList, AntibodyNeutralFile, EscapeScoreFile
{
    if (args != 3)
    {
        printf("/// ERROR CalAntibodyPressure: 2 input needed. The InfectionList, The path to the Data fold which contains AntibodyNeutralFile, EscapeScoreFile\n");
        printf("/// The output file will be at the same fold as the InfectionList\n");
        return 1;
    }

    ReadinAntibody(argv[2]);
    ReadinInfectionList(argv[1]);
    CalPressureMatrix(argv[1]);

    return 0;
}
