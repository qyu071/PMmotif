
#include "stdafx.h"

class EncodingTools
{
public:
    EncodingTools();
    ~EncodingTools();

    public:
        /**
         编码序列集A－00；C－01；G－10；T－11；
         @Param
            输入序列seqTI
            输出序列bitTI
            参数mc（当前验证层为H）
         @return
            长度为H+2字符串长度的字符串编码
         */

        void encodingSeqs(vector<string>& seqTI, vector<vector<unsigned int> >& bitTI, ParaConfig mc);
        void encodingSeqs(vector<string>& seqTI, vector<vector<unsigned long long> >& bitTI, ParaConfig mc);

        /**
         反向编码
         @Param
            编码整数
            字符串长度
            字符集合
         @return
            编码整数对应的字符串
         */
        string deEncodingLmer(unsigned long long lmerInt, int len, string alphabet);

    private:
        //定义常量
        string ALPHABETSET = "ACGT";

        //编码单条序列
        void encodingSeq(string& seqI, ParaConfig mc, vector<unsigned int>& res);
        void encodingSeq(string& seqI, ParaConfig mc,vector<unsigned long long>& res);

        //编码lmer字符串
        unsigned long long encodingLmer(string str, string alphabet);
//        unsigned int encodingLmer(string str, string alphabet);

        //差值支付c在str的位置
        int find(string str, char c);

};
