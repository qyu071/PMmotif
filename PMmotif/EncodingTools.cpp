#include "EncodingTools.h"
#include <typeinfo>
#include <bits/stdc++.h>

EncodingTools::EncodingTools(){}
EncodingTools::~EncodingTools(){}

#if 1
//反向编码(最大支持l＝32的长度)
string EncodingTools::deEncodingLmer(unsigned long long lmerInt, int len, string alphabet)
{
    string res ;
    int end2Bit = 0;
    int flag = len;
    while(flag!=0)
    {
        //取末尾的两位
        end2Bit = lmerInt & 3;
        char c = alphabet[end2Bit];
        res.insert(res.begin(), c);
        lmerInt = lmerInt >> 2;
        flag--;
    }
    return res;
}



//序列集编码(每一个字符串长度为l)
void EncodingTools::encodingSeqs(vector<string> &seqTI, vector<vector<unsigned int> > &bitTI, ParaConfig mc)
{
    for (int i = 0; i < seqTI.size(); ++i)
    {
        vector<unsigned int> temp;
        encodingSeq(seqTI[i], mc, temp);
        bitTI.push_back(temp);
    }
}
//重载
void EncodingTools::encodingSeqs(vector<string> &seqTI, vector<vector<unsigned long long> > &bitTI, ParaConfig mc)
{
    for (int i = 0; i < seqTI.size(); ++i)
    {
        vector<unsigned long long> temp;
        encodingSeq(seqTI[i], mc, temp);
        bitTI.push_back(temp);
    }
}

//单条序列编码
void EncodingTools::encodingSeq(string& seqI, ParaConfig mc,vector<unsigned int>& res)
{
    int len = mc.l;
    for (unsigned int i = 0; i < seqI.size()- len + 1; ++i)
    {
        unsigned int encoding = (unsigned int)encodingLmer(seqI.substr(i,len),ALPHABETSET);
        res.push_back(encoding);
    }
}
//重载
void EncodingTools::encodingSeq(string& seqI, ParaConfig mc,vector<unsigned long long>& res)
{
    int len =  mc.l;
    for (unsigned int i = 0; i < seqI.size()- len + 1; ++i)
    {

        unsigned long long encoding = encodingLmer(seqI.substr(i,len),ALPHABETSET);
        res.push_back(encoding);
    }
}


//长为l字符串的编码
unsigned long long EncodingTools::encodingLmer(string str, string alphabet)
{
    unsigned long long res = 0;
    int flag = 0;
    res = find(alphabet,str[flag]);
    for (int i = 1; i < str.size(); ++i)
    {
		flag = find(alphabet, str[i]);
        res = res << 2;
        res = res + flag;
    }
    return res;
}


//求字符c在str的位置
int EncodingTools::find(string str, char c)
{
	for (unsigned int i = 0; i < str.size(); ++i)
	{
		if (str[i] == c)
		{
			return i;
		}
	}
	return -1;
}
#endif
