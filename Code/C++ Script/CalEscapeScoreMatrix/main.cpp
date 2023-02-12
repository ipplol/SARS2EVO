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

string AA20 = "ACDEFGHIKLMNPQRSTVWY";

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

int Readin(string file, string path)
{
    ifstream read(file);
    if (!read.is_open())
    {
        printf("Can not open file %s\n", file);
        return 1;
    }
    int i, j, k;
    string line;
    getline(read, line);
    vector<vector<string>> ESSMatrix;
    vector<string> tmpv;
    tmpv.push_back("ESSM");
    for (j = 0; j < AA20.size(); j++)
    {
        string tmps = "";
        tmps += AA20[j];
        tmpv.push_back(tmps);
    }
    ESSMatrix.push_back(tmpv);
    for (i = 331; i <= 531; i++)
    {
        tmpv[0] = to_string(i);
        ESSMatrix.push_back(tmpv);
    }

    string antibodyName = "$TBD$";
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, ',');
        if (line1[0] != antibodyName)
        {
            if (antibodyName != "$TBD$")
            {
                //output
                ofstream write(path + "/" + antibodyName + ".EscapeMatrix");
                for (i = 0; i < ESSMatrix.size(); i++)
                {
                    string output = ESSMatrix[i][0];
                    for (j = 1; j < ESSMatrix[i].size(); j++)
                    {
                        if (i!=0 && AA20.find(ESSMatrix[i][j]) != -1)
                            output += "\tNA";
                        else
                            output += "\t" + ESSMatrix[i][j];
                    }
                    write << output << endl;
                }
                write.close();
            }
            //reset the matrix
            antibodyName = line1[0];
            for (i = 331; i <= 531; i++)
            {
                tmpv[0] = to_string(i);
                ESSMatrix[i-330] = tmpv;
            }
        }
        for (i = 331; i <= 531; i++)
        {
            if (line1[1] == ESSMatrix[i - 330][0])
            {
                ESSMatrix[i - 330][AA20.find(line1[2])+1] = line1[3];
                break;
            }
        }
    }
    ofstream write(path + "/" + antibodyName + ".EscapeMatrix");
    for (i = 0; i < ESSMatrix.size(); i++)
    {
        string output = ESSMatrix[i][0];
        for (j = 1; j < ESSMatrix[i].size(); j++)
        {
            if (i!=0&&AA20.find(ESSMatrix[i][j]) != -1)
                output += "\tNA";
            else
                output += "\t" + ESSMatrix[i][j];
        }
        write << output << endl;
    }
    write.close();

    return 0;
}

int main(int args, char* argv[])//两个输入，抗体的逃逸得分文件路径，输出文件夹路径
{
    printf("%s\n", "CalEscapeScoreMatrix");
    if (args != 3)
    {
        printf("/// ERROR CalEscapeScoreMatrix: 2 input needed. The use_abs_res.csv file & The output path");
        return 1;
    }
    return Readin(argv[1], argv[2]);
}