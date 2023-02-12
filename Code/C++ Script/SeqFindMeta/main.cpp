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

class Seq
{
public:
    string SeqID;//Accession Id
    string Mutation;//Nucleotide mutation
    string Location;//Collection location
    string CollectionDate;//Collection date
    bool QCPass = false;//If "human, complete, high quality, GISAID source"
};

unordered_map<string, Seq> SeqMap;//Hash map for accession id and Seq information

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

bool Readin(string fold)
{
    int i, j, k;
    ifstream reads(fold + "/sample_mutlist.tsv");
    ifstream readm(fold + "/metadata.tsv");
    if (!reads.is_open())
    {
        cout << "/// ERROR SeqFindMeta: Can NOT open sample_mutlist.tsv ///";
        return false;
    }
    if (!readm.is_open())
    {
        cout << "/// ERROR SeqFindMeta: Can NOT open metadata.tsv ///";
        return false;
    }

    string line;
    while (getline(reads,line))//read seq id and mutation data from sample_mutlist
    {
        vector<string> line1 = splitStr(line, '\t');
        if (SeqMap.find(line1[0]) == SeqMap.end())
        {
            Seq newseq;
            newseq.SeqID = line1[0];
            newseq.Mutation = line1[1];
            SeqMap.insert(pair<string, Seq>(line1[0], newseq));
        }
    }

    getline(readm, line);
    while (getline(readm,line))//read seq metadata
    {
        vector<string> line1 = splitStr(line, '\t');

        //only keep seq that have GISAID id
        string accession = "$*TBD*$";
        if (line1[1].length() > 8 && line1[1].substr(0, 7) == "EPI_ISL")
            accession = line1[1];
        else
            if (line1[2].length() > 8 && line1[2].substr(0, 7) == "EPI_ISL")
                accession = line1[2];
        
        //the seq exist in our hash table
        if (accession != "$*TBD*$")
        {
            if (SeqMap.find(line1[1]) != SeqMap.end())
            {
                if (line1[5] == "Complete" && line1[7] == "High" && line1[9] == "Homo sapiens")//sequence fliter
                {
                    SeqMap[line1[1]].SeqID = accession;//replace with GISAID id if any
                    SeqMap[line1[1]].QCPass = true;
                    SeqMap[line1[1]].CollectionDate = line1[10];
                    //Location only keep the country level
                    vector<string> line2 = splitStr(line1[11], '/');
                    k = line2[0].size() - 1;
                    while (line2[0][k] == ' ')k--;
                    SeqMap[line1[1]].Location = line2[0].substr(0, k + 1);
                }
            }
            else
            if (SeqMap.find(line1[2]) != SeqMap.end())
            {
                if (line1[5] == "Complete" && line1[7] == "High" && line1[9] == "Homo sapiens")//sequence fliter
                {
                    SeqMap[line1[2]].SeqID = accession;
                    SeqMap[line1[2]].QCPass = true;
                    SeqMap[line1[2]].CollectionDate = line1[10];
                    //Location only keep the country level
                    vector<string> line2 = splitStr(line1[11], '/');
                    k = line2[0].size() - 1;
                    while (line2[0][k] == ' ')k--;
                    SeqMap[line1[1]].Location = line2[0].substr(0, k + 1);
                }
            }
        }
    }

    reads.close();
    readm.close();
    return true;
}

bool Output(string fold)
{
    ofstream write(fold + "/sample_mut_loc_time.tsv");
    if (!write.is_open())
    {
        cout << "/// ERROR SeqFindMeta: Can NOT create output file at fold: "<< fold <<"///";
        return false;
    }
    string output = "";
    for (auto val = SeqMap.begin(); val != SeqMap.end(); val++)
    {
        if (val->second.QCPass == true)
        {
            output = val->second.SeqID + "\t";
            if(val->second.Mutation[val->second.Mutation.size()-1]!='\r')
                output += val->second.Mutation + "\t";
            else
                output += val->second.Mutation.substr(0, val->second.Mutation.size() - 1) + "\t";
            output += val->second.Location + "\t";
            output += val->second.CollectionDate;
            write << output << endl;
        }
    }
    write.close();
}

int main(int argc, char* argv[])//1个输入，目标文件夹; 
{
    if (argc != 2)
    {
        printf("/// ERROR SeqFindMeta: 1 Input Needed. Please specify the fold that contains sample_mutlist and metadata");
        return 1;
    }

    string targetfold = argv[1];
    //string targetfold = "/mnt/g/VariationMutation/病毒逃逸预测/Fitness/SelectiveCoEfficient/data";

    if (!Readin(targetfold))
    {
        return 1;
    }
    if (!Output(targetfold))
    {
        return 1;
    }
    return 0;
}