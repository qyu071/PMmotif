
#include <bits/stdc++.h>
#include "ReadFile.h"

ReadFile::ReadFile()
{
}

ReadFile::~ReadFile()
{
}

void ReadFile::readSequence(const char *fileName, vector<string>& strings)
{
    if (read(fileName, strings))
    {
        //不作处理
    }
    else
    {
        cerr << "文件读取出错" << endl;
    }
}
//按行读取文件
bool ReadFile::read(const char *fileName, vector<string>& strings)
{
    ifstream myfile(fileName);
    string data;

    if (myfile.is_open())
    {
        string line;

        while (myfile.good())
        {
            //按行读取字符串
            getline(myfile, line);
            //去除空格
            trim(line);
            if (line[0] == '>')
            {
                if (data.size() > 0)
                {
                    strings.push_back(data);
                    data.clear();
                }
            }
            else
            {
                data.append(line);
            }
        }

        myfile.close();

        if (data.size() > 0)
        {
            //存储最后一个数据源
            strings.push_back(data);
        }
        return true;
    }
    else
    {
        return false;
    }
}

//清除字符串的两头
void ReadFile::trim(string& s)
{
    string delim = " \n\r\t";

    s.erase(s.find_last_not_of(delim) + 1);
    s.erase(0, s.find_first_not_of(delim));
}
