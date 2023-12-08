#include <iostream>
#include <bits/stdc++.h>
#include<time.h>

#include "stdafx.h"
#include "GetFilePathName.h"
#include "ReadFile.h"
#include "Tools.h"
#include "PMSMotif.h"
using namespace std;

int main()
{
   	string ALPHABETSET = "ACGT";//字母表
   	int beta = 3; //  beta表示β

    //获取directory路径下的文件名称
	string directory = ".\\data\\";
    GetFilePathName GFP;
    vector<string> filesPath;
    GFP.getFiles(directory, filesPath);
    for (int i = 0 ; i<filesPath.size() ; ++i)//输出目录下的文件名
    {
        cout << filesPath[i] <<endl;
    }

    //获得当前时间
    string currenttime ="";
    time_t t = time(0);
    char tmp[64];
    strftime(tmp,sizeof(tmp),"%Y%m%d_%H%M%S",localtime(&t));
    puts(tmp);
    for(int i=0;i<strlen(tmp); ++i)
    {
        currenttime += tmp[i];
    }


    //输出文件（路径+名称）
    fstream f(".\\data\\output_"+currenttime+".txt",ios::out);

    //找出每个文件中的(l,d)模体
    for (int i = 0; i < filesPath.size(); ++i)  // 遍历directory路径下的文件名称
    {
        if(filesPath[i].find("data")!=filesPath[i].npos)//找出其中存储DNA序列数据集的文件
        {
            double start = clock();//程序运行开始时间
            cout << "Test File " << filesPath[i] << "..." << endl; //输出当前搜索的DNA序列数据集的文件名称

            //从文件名称中获取输入参数t,n,l,d,h（h初始的公共前缀所在的层）
            ParaConfig mc = getParameter(filesPath[i]);
            mc.H = mc.H - beta;

            cout<<mc.l<<' '<<mc.H<<endl; //输出当前l和初始公共前缀层

            //读取DNA序列(存入seqTI)
            ReadFile RF;
            vector<string> seqTI;
            string filePath = directory + filesPath[i];
            RF.readSequence(filePath.c_str(), seqTI);
            vector<string> res;

            //验证候选模体前缀
            PMSMotif PM;
            PM.getMotif(seqTI, mc, res);


            double end1 = clock(); // 程序运行结束时间
            cout << "take time :" << (end1-start)/CLOCKS_PER_SEC << endl<<endl;//输出程序运行时间

            /*
            向输出文件中写入数据
            包括：
            1、文件名称
            2、初始公共前缀层mc.h
            3、找到的（l,d）模体个数
            4、（l,d）模体
            5、程序运行时间
            */
            f << "Test File " << filesPath[i] << "..." << endl;
            f << "mc.H: " << mc.H << endl;
            f << res.size() <<endl;


            /*
            为查看结果提供方便最多输出50个（l,d）模体
            */
            int temp;
            if (res.size()<50)
            {
                temp = res.size();
            }
            else
            {
                temp = res.size();
            }
            for(int i = 0; i<temp; ++i)
            {
                f<< res[i] <<endl;
            }
            f << "take time :" <<(end1-start)/CLOCKS_PER_SEC <<endl<<endl;
        }
    }
    f.close();
}
