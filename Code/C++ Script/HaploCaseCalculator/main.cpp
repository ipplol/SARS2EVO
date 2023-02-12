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
    int SeqNumber = 0;
    int AfterSeqNumber = 0;
    unordered_map<string, int> DateSeqNumber;
    vector<double> DailyCases;
};

class SamplingDate
{
public:
    string Date;
    int Date_int;
    int DailyCases;
    int SeqNumber=0;
    double StringencyIndex;
};

vector<Haplotype> haploList;//all haplotypes
unordered_map<string, int> haploListMap;
vector<SamplingDate> dateList;//All days between collection start and end, also plus the +-3days of after collection
unordered_map<string, int> dateListMap;//from date to pos in the list

int CollectionStartDate;
int CollectionEndDate;
int AfterCollectionDate;
string CollectionLocation;
int CollectSeqNumber = 0;
int AfterSeqNumber = 0;
double AverageStrigencyIndex = 0;
int CollectionDayNumber = 0;

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

//read in dependence files
void Readin(string fold)
{
    string line;
    //Readin stringency index from OWID and figure out how many days
    ifstream read(fold + "/owid-covid-data.csv");
    while (getline(read,line))
    {
        vector<string> line1 = splitStr(line, ',');
        if (line1[2] == CollectionLocation && line1[3].size() == 10)
        {
            int seqDate = stoi(line1[3].substr(0, 4)) * 10000 + stoi(line1[3].substr(5, 2)) * 100 + stoi(line1[3].substr(8, 2));
            if (seqDate >= CollectionStartDate && seqDate <= CollectionEndDate)//采样时间
            {
                SamplingDate newDate;
                newDate.Date = line1[3];
                newDate.Date_int = seqDate;
                newDate.DailyCases = 0;
                newDate.StringencyIndex = stof(line1[47]);
                AverageStrigencyIndex += stof(line1[47]);
                CollectionDayNumber++;
                dateList.push_back(newDate);
                dateListMap.insert(pair<string, int>(line1[3], dateList.size() - 1));
            }
            else
                if (seqDate >= AfterCollectionDate - 3 && seqDate <= AfterCollectionDate + 3)//采样后时间，研究频率变化
                {
                    SamplingDate newDate;
                    newDate.Date = line1[3];
                    newDate.Date_int = seqDate;
                    newDate.DailyCases = 0;
                    newDate.StringencyIndex = stof(line1[47]);
                    dateList.push_back(newDate);
                    dateListMap.insert(pair<string, int>(line1[3], dateList.size() - 1));
                }
        }
    }
    read.close();
    AverageStrigencyIndex /= CollectionDayNumber;

    //Readin dailyCases from lancet estim
    int i;
    for (i = 2020; i <= 2022; i++)
    {
        ifstream read2(fold + "/data_download_file_reference_" + to_string(i) + ".csv");
        while (getline(read2, line))
        {
            vector<string> line1 = splitStr(line, ',');
            if (line1[3].find(CollectionLocation)!=-1 && line1[1].size() == 10)
            {
                int seqDate = stoi(line1[1].substr(0, 4)) * 10000 + stoi(line1[1].substr(5, 2)) * 100 + stoi(line1[1].substr(8, 2));
                if ((seqDate >= CollectionStartDate && seqDate <= CollectionEndDate)|| (seqDate >= AfterCollectionDate - 3 && seqDate <= AfterCollectionDate + 3))//采样时间
                {
                    if(line1[4]!="")
                        dateList[dateListMap[line1[1]]].DailyCases += stof(line1[4]);
                }
            }
        }
        read2.close();
    }
    
    //Readin haplotypes
    ifstream read3(fold + "/sample_mut_loc_time.tsv.AAmut");
    while (getline(read3, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        if (line1[2] == CollectionLocation && line1[3].size() == 10)
        {
            int seqDate = stoi(line1[3].substr(0, 4)) * 10000 + stoi(line1[3].substr(5, 2)) * 100 + stoi(line1[3].substr(8, 2));
            if ((seqDate >= CollectionStartDate && seqDate <= CollectionEndDate) || (seqDate >= AfterCollectionDate - 3 && seqDate <= AfterCollectionDate + 3))//采样时间
            {
                if (haploListMap.find(line1[1]) != haploListMap.end())//+1
                {
                    if (seqDate <= CollectionEndDate)
                    {
                        haploList[haploListMap[line1[1]]].SeqNumber++;
                        CollectSeqNumber++;
                    }
                    else
                    {
                        haploList[haploListMap[line1[1]]].AfterSeqNumber++;
                        AfterSeqNumber++;
                    }
                    haploList[haploListMap[line1[1]]].DateSeqNumber[line1[3]]++;
                }
                else//add new haplo
                {
                    Haplotype newa;
                    newa.Haplo = line1[1];
                    newa.DateSeqNumber.insert(pair<string, int>(line1[3], 1));
                    if (seqDate <= CollectionEndDate)
                    {
                        newa.SeqNumber++;
                        CollectSeqNumber++;
                    }
                    else
                    {
                        newa.AfterSeqNumber++;
                        AfterSeqNumber++;
                    }
                    haploList.push_back(newa);
                    haploListMap.insert(pair<string, int>(line1[1], haploList.size() - 1));
                }
                dateList[dateListMap[line1[3]]].SeqNumber++;
            }
        }
    }
    read3.close();
    return;
}

//daily cases * sequence prop = haplo cases
void HaploCaseEstim(string fold)
{
    string filename;
    if (CollectionLocation.find(" ") != -1)
        filename = CollectionLocation.replace(CollectionLocation.find(" "), 1, "+");
    else
        filename = CollectionLocation;
    filename += "_" + to_string(CollectionStartDate) + "_" + to_string(CollectionEndDate) + "_" + to_string(AfterCollectionDate);
    //only choose haplo that 0.05% prop and 50 seq and cover 80% of the collection time
    vector<int> filteredHaplo;
    int i, j, k;
    for (i = 0; i < haploList.size(); i++)
    {
            if(((double)haploList[i].SeqNumber / CollectSeqNumber >= 0.0005))
                if(haploList[i].DateSeqNumber.size()>=(dateList.size()-5)*0.80)
                {
                    filteredHaplo.push_back(i);
                }
    }

    //Calculate after prop
    ofstream write(fold + "/" + filename + ".LaterProp");
    write << "HaploNum\tHaplo\tStrigencyIndex\tPropNow\tProp" << to_string(AfterCollectionDate) << endl;
    for (i = 0; i < filteredHaplo.size(); i++)
    {
        string output = to_string(i) + "\t";
        output += haploList[filteredHaplo[i]].Haplo + "\t";
        output += to_string(AverageStrigencyIndex) + "\t";
        output += to_string((double)haploList[i].SeqNumber / CollectSeqNumber) + "\t";
        output += to_string((double)haploList[filteredHaplo[i]].AfterSeqNumber / AfterSeqNumber);
        write << output << endl;
    }
    write.close();

    //Calculate daily cases
    ofstream write2(fold + "/" + filename + ".DailyCases");
    string output = "Date";
    for (i = 0; i < filteredHaplo.size(); i++)
    {
        output += "\tHaplo" + to_string(i);
        output += "\tWT" + to_string(i);
    }
    write2 << output << endl;
    for (j = 0; j < dateList.size() && dateList[j].Date_int <= CollectionEndDate; j++)
    {
        output = dateList[j].Date;
        for (i = 0; i < filteredHaplo.size(); i++)
        {
            output += "\t";
            if (haploList[filteredHaplo[i]].DateSeqNumber.find(dateList[j].Date) != haploList[filteredHaplo[i]].DateSeqNumber.end())
            {
                output += to_string((double)haploList[filteredHaplo[i]].DateSeqNumber[dateList[j].Date] / dateList[j].SeqNumber * dateList[j].DailyCases);
                output += "\t" + to_string(dateList[j].DailyCases - ((double)haploList[filteredHaplo[i]].DateSeqNumber[dateList[j].Date] / dateList[j].SeqNumber * dateList[j].DailyCases));
            }
            else
            {
                output += "0";
                output += "\t" + to_string(dateList[j].DailyCases);
            }
        }
        write2 << output << endl;
    }
    write2.close();
    return;
}

/*7个输入 
1.input path (RBDmetadata, lancet daily cases)
2.output file path
3.collection date start [yyyy-mm-dd]
4.collection date end [yyyy-mm-dd]
5.collection date later [yyyy-mm-dd]
6.country [United States]
*/
int main(int argc, char*argv[])
{
    if (argc != 7)
    {
        printf("/// ERROR HaploCaseCalculator: 6 input needed. Please follow the order, 1:input fold, 2:output path, 3:collection date start, 4:end, 5:later for prop, 6:country");
        return 1;
    }
    CollectionLocation = argv[3];
    if(CollectionLocation.find("_")!=-1)
        CollectionLocation = CollectionLocation.replace(CollectionLocation.find("_"), 1, " ");//United_States -> United States
    string para = argv[4];
    CollectionStartDate = stoi(para.substr(0, 4)) * 10000 + stoi(para.substr(5, 2)) * 100 + stoi(para.substr(8, 2));
    para = argv[5];
    CollectionEndDate = stoi(para.substr(0, 4)) * 10000 + stoi(para.substr(5, 2)) * 100 + stoi(para.substr(8, 2));
    para = argv[6];
    AfterCollectionDate = stoi(para.substr(0, 4)) * 10000 + stoi(para.substr(5, 2)) * 100 + stoi(para.substr(8, 2));

    Readin(argv[1]);//input files
    HaploCaseEstim(argv[2]);//calculate and output

    return 0;
}