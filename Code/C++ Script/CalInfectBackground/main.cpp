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

class Variant
{
public:
    string Name;//BA.1.1
    string VOC;//Omicron
    int SeqNumber = 0;
    unordered_map<string, int> DateSeqNumber;
};

class SamplingDate
{
public:
    string Date;
    int Date_int;
    double DailyCases;
    int SeqNumber = 0;
    unordered_map<string, int> VariantSeqNum;
    unordered_map<string, double> VariantCaseNum;
};

vector<Variant> variantList;//historical infected variant
unordered_map<string, int> variantListMap;
vector<SamplingDate> dateList;//All days between collection start and end, also plus the +-3days of after collection
unordered_map<string, int> dateListMap;//from date to pos in the list

int CollectionStartDate;
int CollectionEndDate;
string CollectionLocation;
int CollectionDayNumber = 0;
int people_vaccinated = 0;

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

    //Readin dailyCases and vac from lancet estim
    int i;
    ifstream read2(fold + "/LancetCasePerMillion_Vac");
    while (getline(read2, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        if (line1[0] == CollectionLocation && line1[1].size() == 10)
        {
            int seqDate = stoi(line1[1].substr(0, 4)) * 10000 + stoi(line1[1].substr(5, 2)) * 100 + stoi(line1[1].substr(8, 2));
            if ((seqDate >= CollectionStartDate && seqDate <= CollectionEndDate))//采样时间
            {
                SamplingDate newDate;
                newDate.Date = line1[1];
                newDate.Date_int = seqDate;
                newDate.DailyCases = 0;
                dateList.push_back(newDate);
                dateListMap.insert(pair<string, int>(line1[1], dateList.size() - 1));

                if (line1[2] != "")
                    dateList[dateListMap[line1[1]]].DailyCases += stof(line1[2]);
                if (line1[3] != "")
                    if (stoi(line1[3]) > people_vaccinated)
                        people_vaccinated = stoi(line1[3]);
            }
        }
    }
    read2.close();

    //Readin strain
    ifstream read3(fold + "/metadata.tsv");
    getline(read3, line);
    while (getline(read3, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        if (line1[11].find(CollectionLocation)!=-1 && line1[10].size() == 10)
        {
            int seqDate = stoi(line1[10].substr(0, 4)) * 10000 + stoi(line1[10].substr(5, 2)) * 100 + stoi(line1[10].substr(8, 2));
            if ((seqDate >= CollectionStartDate && seqDate <= CollectionEndDate))//采样时间
            {
                if (variantListMap.find(line1[4]) != variantListMap.end())//+1
                {
                    
                    variantList[variantListMap[line1[4]]].SeqNumber++;
                    variantList[variantListMap[line1[4]]].DateSeqNumber[line1[10]]++;
                }
                else//add new lineage
                {
                    Variant newa;
                    newa.Name = line1[4];
                    newa.VOC = "TBD";
                    newa.DateSeqNumber.insert(pair<string, int>(line1[10], 1));
                    newa.SeqNumber++;
                    
                    variantList.push_back(newa);
                    variantListMap.insert(pair<string, int>(line1[4], variantList.size() - 1));
                }
                dateList[dateListMap[line1[10]]].SeqNumber++;
                if (dateList[dateListMap[line1[10]]].VariantSeqNum.find(line1[4]) != dateList[dateListMap[line1[10]]].VariantSeqNum.end())
                    dateList[dateListMap[line1[10]]].VariantSeqNum[line1[4]]++;
                else
                    dateList[dateListMap[line1[10]]].VariantSeqNum.insert(pair<string, int>(line1[4], 1));
            }
        }
    }
    read3.close();
    return;
}

//daily cases * sequence prop = variant cases
void VariantCaseEstim(string fold)
{
    string filename;
    if (CollectionLocation.find(" ") != -1)
        filename = CollectionLocation.replace(CollectionLocation.find(" "), 1, "+");
    else
        filename = CollectionLocation;
    filename += "_" + to_string(CollectionStartDate) + "_" + to_string(CollectionEndDate);
    ofstream write(fold + "/" + filename + ".InfectionList");

    vector<int> filteredHaplo;
    int i, j, k;
    k = 0; 
    for (i = 0; i < dateList.size(); i++)//if data is missing for someday, fill it with formal days
    {
        if (dateList[i].SeqNumber != 0)
            k = i;
        else
        {
            dateList[i].SeqNumber = dateList[k].SeqNumber;
            dateList[i].VariantSeqNum = dateList[k].VariantSeqNum;
        }

        //multiply seq percent with daily cases
        for (auto val = dateList[i].VariantSeqNum.begin(); val != dateList[i].VariantSeqNum.end(); val++)
        {
            dateList[i].VariantCaseNum.insert(pair<string, double>(val->first, ((double)val->second / dateList[i].SeqNumber * dateList[i].DailyCases)));
        }
    }

    //Calculate historical cases
    string output;
    write << "Strain\tCasePerMillion" << endl;
    write << "people_vaccinated\t" << people_vaccinated << endl;
    for (i = 0; i < variantList.size(); i++)
    {
        output = variantList[i].Name + "\t";
        double totalCasesPerMillion = 0;
        for (j = 0; j < dateList.size(); j++)
        {
            if (dateList[j].VariantCaseNum.find(variantList[i].Name) != dateList[j].VariantCaseNum.end())
                totalCasesPerMillion += dateList[j].VariantCaseNum[variantList[i].Name];
        }
        output += to_string(totalCasesPerMillion);
        write << output << endl;
    }
    write.close();
    return;
}

/*5个输入
1.input path (RBDmetadata, lancet daily cases)
2.output file path
3.country [United_States]
4.collection date start [yyyy-mm-dd]
5.collection date end [yyyy-mm-dd]
*/
int main(int argc, char* argv[])
{
    if (argc < 6)
    {
        printf("/// ERROR CalInfectBackground: 5 input needed. Please follow the order, 1:input fold, 2:output path, 3: country 4:collection date start, 5:end");
        return 1;
    }
    CollectionLocation = argv[3];
    if (CollectionLocation.find("_") != -1)
        CollectionLocation = CollectionLocation.replace(CollectionLocation.find("_"), 1, " ");//United_States -> United States
    string para = argv[4];
    //CollectionStartDate = stoi(para.substr(0, 4)) * 10000 + stoi(para.substr(5, 2)) * 100 + stoi(para.substr(8, 2));
    CollectionStartDate = 20191201;
    para = argv[5];
    CollectionEndDate = stoi(para.substr(0, 4)) * 10000 + stoi(para.substr(5, 2)) * 100 + stoi(para.substr(8, 2));

    Readin(argv[1]);//input files
    VariantCaseEstim(argv[2]);//calculate and output

    return 0;
}