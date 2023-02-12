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

class Day
{
public:
    string Date;
    double NewCases;
    double TotalVaccined = 0;
};

class Country
{
public:
    string Name;
    double Population;
    unordered_map<string, int> dateMap;//date and pos in the list
    vector<Day> dateList;
};

unordered_map<string, Country> countryMap;

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

void Readin(string fold)
{
    string line;
    int i, j, k;
    ifstream read(fold + "/owid-covid-data.csv");
    getline(read, line);
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, ',');
        if (countryMap.find(line1[2]) != countryMap.end())// country existed
        {
            Day newd;
            newd.Date = line1[3];
            if(line1[35]!="")
                newd.TotalVaccined = stof(line1[35]);
            else
            {
                if (countryMap[line1[2]].dateList.size() != 0)
                    newd.TotalVaccined = countryMap[line1[2]].dateList[countryMap[line1[2]].dateList.size() - 1].TotalVaccined;
            }
            countryMap[line1[2]].dateMap.insert(pair<string, int>(line1[3], countryMap[line1[2]].dateList.size()));
            countryMap[line1[2]].dateList.push_back(newd);
        }
        else// new 
        {
            Country newa;
            newa.Name = line1[2];
            if (line1[48] != "")
                newa.Population = stof(line1[48]);
            else
                continue;
            Day newd;
            newd.Date = line1[3];
            if (line1[35] != "")
                newd.TotalVaccined = stof(line1[35]);
            newa.dateMap.insert(pair<string, int>(line1[3], 0));
            newa.dateList.push_back(newd);
            countryMap.insert(pair<string,Country>(line1[2], newa));
        }
    }
    read.close();

    for (i = 2020; i <= 2022; i++)
    {
        ifstream read2(fold + "/data_download_file_reference_"+to_string(i)+".csv");
        getline(read2, line);
        while (getline(read2, line))
        {
            vector<string> line1 = splitStr(line, ',');
            string loc = line1[3];
            if (loc == "United States of America")
                loc = "United States";
            if (countryMap[loc].dateMap.find(line1[1]) != countryMap[loc].dateMap.end())
            {
                countryMap[loc].dateList[countryMap[loc].dateMap[line1[1]]].NewCases = stof(line1[4]);
            }
        }
        read2.close();
    }
    return;
}

void Output(string fold)
{
    ofstream write(fold + "/LancetCasePerMillion_Vac");
    write << "Country\tDate\tNewInfectedPM\tTotalVacPeoplePM" << endl;
    int i, j, k;
    for (auto val = countryMap.begin(); val != countryMap.end(); val++)
    {
        for (i = 0; i < val->second.dateList.size(); i++)
        {
            string output = val->first;
            output += "\t" + val->second.dateList[i].Date;
            output += "\t" + to_string(val->second.dateList[i].NewCases * 1000000 / val->second.Population);
            output += "\t" + to_string(val->second.dateList[i].TotalVaccined * 1000000 / val->second.Population);
            write << output << endl;
        }
    }
}

int main(int argc, char* argv[])//1个输入 datafold的目录
{
    if (argc != 2)
    {
        printf("/// ERROR CasePerMillion: 1 input needed. The path of data fold.\n");
        return 1;
    }
    Readin(argv[1]);
    Output(argv[1]);
    return 0;
}