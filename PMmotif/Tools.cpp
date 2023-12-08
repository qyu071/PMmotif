
#include <bits/stdc++.h>
#include "Tools.h"


ParaConfig getParameter(string fileName)
{
    ParaConfig mc;

    //样本序列条数
    mc.t = atof(fileName.substr(0,fileName.find("t")).c_str());

    //样本序列长
    mc.n = atof(fileName.substr(fileName.find("t")+2,fileName.find("n")-fileName.find("t")-2).c_str());

    //模体长度
    mc.l = atof(fileName.substr(fileName.find("n")+2,fileName.find("l")-fileName.find("n")-2).c_str());

    //开始验证层
    mc.H = atof(fileName.substr(fileName.find("l")+2,fileName.find("h")-fileName.find("l")-2).c_str());

    //最大容许突变值
    mc.d = atof(fileName.substr(fileName.find("h")+2,fileName.find("d")-fileName.find("h")-2).c_str());

    //q指比例值
    mc.q = atof(fileName.substr(fileName.find("d")+2,fileName.find("q")-fileName.find("d")-2).c_str())/mc.t;

    //β
    mc.hAdd = 2;
    //随机插入模体时突变的范围(生成数据写 0-d的范围)
    mc.MINSTISFYSEQ = mc.q * mc.t;
    return mc;
}

