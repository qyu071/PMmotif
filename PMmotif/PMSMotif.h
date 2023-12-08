#include "stdafx.h"
#include <bits/stdc++.h>

#include "EncodingTools.h"
//#include "HamDistance.h"
//AATAGCCGT

class PMSMotif
{
public:
    PMSMotif();
    ~PMSMotif();

    public:
         /**
          求解植入模体搜索算法
          @Param
            输入序列集合
            参数mc
            最终结果集res
          @return
            结果集res
          */
        void getMotif(vector<string>& seqTI, ParaConfig mc, vector<string>& res);


    private:



        void checkDownRight(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len);
        void checkDownRight(unsigned long long blockFunElement, vector<vector<unsigned long long> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len);

        void checkDownThree(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len);
        void checkDownThree(unsigned long long blockFunElement, vector<vector<unsigned long long> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len);

        void checkDownThreeLocation(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len);
        void checkDownThreeLocation(unsigned long long int blockFunElement, vector<vector<unsigned long long int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len);


        void checkDown(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc,  vector<string>& result, int cur_len);
        void checkDown(unsigned long long blockFunElement, vector<vector<unsigned long long> >& bitTI, ParaConfig mc,  vector<string>& result, int cur_len);




};
