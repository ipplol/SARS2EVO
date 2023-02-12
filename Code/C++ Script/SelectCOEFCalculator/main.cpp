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
    string Haplo;
    string Location;
    string DateStart;
    string DateEnd;
    string PropNow;
    double StringencyIndex;
    double ReproductionNumber = 0;
    double R_wt = 0;
    double SelectionCoefficient;
    string DateLater;
    string PropLater;
};
vector<Haplotype> HaplotypeList;

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

void Readin(string file)
{
    //obtain parameter from file name
    vector<string> line1 = splitStr(file, '.');
    vector<string> line2 = splitStr(line1[0], '_');
    int k = line2[0].size()-1;
    while (line2[0][k] != '/')k--;
    string K_Location = line2[0].substr(k+1,line2[0].size() - k);
    string K_DateStart = line2[1];
    string K_DateEnd = line2[2];
    string K_DateLater = line2[3];

    //read in haplotype
    string input2 = line1[0] + ".DailyCases.RE";
    ifstream read(file);
    if (!read.is_open())
    {
        printf("/// ERROR SelectCOEFCalculator: Can NOT open .LaterProp file");
        return;
    }
    string line;
    getline(read, line);
    while (getline(read,line))
    {
        vector<string> linex = splitStr(line, '\t');
        Haplotype newa;
        newa.Haplo = linex[1];
        newa.Location = K_Location;
        newa.DateStart = K_DateStart;
        newa.DateEnd = K_DateEnd;
        newa.StringencyIndex = stof(linex[2]);
        newa.PropNow = linex[3];
        newa.DateLater = K_DateLater;
        newa.PropLater = linex[4];
        HaplotypeList.push_back(newa);
    }
    read.close();

    //read in R
    ifstream read2(input2);
    if (!read2.is_open())
    {
        printf("/// ERROR SelectCOEFCalculator: Can NOT open .DailyCase.RE file");
        return;
    }
    int day = 0,i,j;
    getline(read2, line);
    getline(read2, line);
    getline(read2, line);
    getline(read2, line);//Discard the first 3 result of R
    while (getline(read2, line))
    {
        vector<string> linex = splitStr(line, '\t');
        day++;
        for (k = 1; k < linex.size(); k += 2)
        {
            HaplotypeList[(k - 1) / 2].ReproductionNumber += stof(linex[k]);
            HaplotypeList[(k - 1) / 2].R_wt += stof(linex[k+1]);
        }
    }
    read2.close();

    //Calculate Co-efficient and sort
    for (i = 0; i < HaplotypeList.size(); i++)
    {
        HaplotypeList[i].ReproductionNumber /= day;
        HaplotypeList[i].R_wt /= day;
        HaplotypeList[i].SelectionCoefficient = (HaplotypeList[i].ReproductionNumber / HaplotypeList[i].R_wt) - 1;
    }
    Haplotype tmpH;
    for (i = 0; i < HaplotypeList.size(); i++)
    {
        for (j = i + 1; j < HaplotypeList.size(); j++)
        {
            if (HaplotypeList[i].SelectionCoefficient < HaplotypeList[j].SelectionCoefficient)
            {
                tmpH = HaplotypeList[i];
                HaplotypeList[i] = HaplotypeList[j];
                HaplotypeList[j] = tmpH;
            }
        }
    }

    //output
    ofstream write(line1[0] + ".HaploFitness");
    write << "Haplo\tLocation\tDateStart\tDateEnd\tSeqProp\tStringencyInd\tRe\tSC\tDateLater\tPropLater" << endl;
    for (i = 0; i < HaplotypeList.size(); i++)
    {
        string output = HaplotypeList[i].Haplo;
        output += "\t" + HaplotypeList[i].Location;
        output += "\t" + HaplotypeList[i].DateStart;
        output += "\t" + HaplotypeList[i].DateEnd;
        output += "\t" + HaplotypeList[i].PropNow;
        output += "\t" + to_string(HaplotypeList[i].StringencyIndex);
        output += "\t" + to_string(HaplotypeList[i].ReproductionNumber);
        output += "\t" + to_string(HaplotypeList[i].SelectionCoefficient);
        output += "\t" + HaplotypeList[i].DateLater;
        output += "\t" + HaplotypeList[i].PropLater;
        write << output << endl;
    }
    write.close();
    return;
}

int main(int args, char*argv[])//1个输入 LaterProp
{
    if (args != 2)
    {
        printf("/// ERROR SelectCOEFCalculator: 1 input needed. The .LaterProp file");
        return 1;
    }
    Readin(argv[1]);
    return 0;
}