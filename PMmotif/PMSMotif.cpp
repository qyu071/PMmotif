#include "PMSMotif.h"
#include "HamDistance.h"

#include <bits/stdc++.h>
#include <string>
#include <sstream>
#include <iterator>

PMSMotif::PMSMotif(){}
PMSMotif::~PMSMotif(){}

string ALPHABETSET = "ACGT";

int mer[16] = {4383, 8751, 17487, 34959, 4593, 8946, 17652, 35064, 7953, 12066, 20292, 36744, 61713, 61986, 62532, 63624};

unsigned long long int mer30[64]={1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
            65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216,
            33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648, 4294967296,
            8589934592, 17179869184, 34359738368, 68719476736, 137438953472, 274877906944, 549755813888,
            1099511627776, 2199023255552, 4398046511104, 8796093022208, 17592186044416, 35184372088832,
            70368744177664, 140737488355328, 281474976710656, 562949953421312, 1125899906842624, 2251799813685248,
            4503599627370496, 9007199254740992, 18014398509481984, 36028797018963968, 72057594037927936,
            144115188075855872,288230376151711744, 576460752303423488, 1152921504606846976, 2305843009213693952, 4611686018427387904,9223372036854775808};

unsigned long long int mer31[64] = {281479271747871, 562958543495727, 1125917086991439, 2251834173982863, 4503668347900401, 9007336695800562,
18014673391600884, 36029346783201528, 72058693566340881, 144117387132677922, 288234774265352004, 576469548530700168,
 1152939097061388561, 2305878194122715682, 4611756388245369924, 9223512776490678408, 281479558922241, 562959116861442,
 1125918232739844, 2251836464496648, 4503668647854096, 9007337279979552, 18014674544230464, 36029349072732288,
 72058694070763776, 144117387889869312, 288234775528080384, 576469550804502528, 1152939100837318656, 2305878197648105472,
 4611756391269679104, 9223512778512826368, 300299818434561, 600535212359682, 1201006000209924, 2401947575910408,
  4523326413209616, 9045622034268192, 18090213276385344, 36179395760619648, 72091751929610496, 144167011184804352,
  288317529695192064, 576618566715967488, 1153186560192024576, 2306109237593382912, 4611954592396099584, 9223645302001532928,
  1233704827217838081, 2463187529785016322, 4922152934919372804, 9840083745188085768, 1292814636752568336,
  2518075279094579232, 4968596563778601024, 9869639133146644608, 2238571589308252416, 3396279268047585792,
  5711694625526252544, 10342525340483586048, 17370682830199197696, 17447543091295690752, 17601263613488676864,
   17908704657874649088};

unsigned long long int mer32[64] = {1233723648051773439, 2463225107027329023, 4922228024978440191, 9840233860880662527, 1292834295117905919,
 2518113565017309183, 4968672104816115711, 9869789184413728767, 2238604648176025599, 3396328892856991743,
  5711777382218924031, 10342674360942788607, 17370930297105940479, 17447774138291912703, 17601461820663857151,
  17908837185407746047, 1233723652059369759, 2463225110748668463, 4922228028127265871, 9840233862884460687,
  1292834299111739889, 2518113568725869298, 4968672107954128116, 9869789186410645752, 2238604651949661969,
  3396328896361082658, 5711777385183924036, 10342674362829606792, 17370930297356415249, 17447774138524496418,
  17601461820860658756, 17908837185532983432, 1233986293891797279, 2463468992450208303, 4922434389567030351,
   9840365183800674447, 1293096039014797809, 2518356612921565938, 4968877760735102196, 9869920056362174712,
   2238851960982806289, 3396558540463288098, 5711971699424251716, 10342798017346178952, 17370946712470941969,
    17447789381130842658, 17601474718450644036, 17908845393090246792, 18446481423861747999, 18446500184565162543,
     18446537705971991631, 18446612748785649807, 18446482325818642929, 18446501022096564978, 18446538414652409076,
      18446613199764097272, 18446496757128961809, 18446514422599003938, 18446549753539088196,
18446620415419256712, 18446727658094063889, 18446728830638027298, 18446731175725954116, 18446735865901807752};


/**
 （1）确定验证层H
 （2）在该层使用向下向右的加速验证手段（实际上找到的是一组((H＋2),d）模体）
 （3）如果H＋2等于l那么退出，如果H＋2小于l那么进行（4）
 （4）继续向下探位验证（不能使用向右跳跃式验证，因为第H＋2层节点不具备向右跳跃式验证的条件）

 */

#if 1
//大框架
void PMSMotif::getMotif(vector<string>& seqTI, ParaConfig mc, vector<string>& result)
{
    EncodingTools ET;

    //根据mc.l，分两种处理方法
    if(mc.l > 0 && mc.l <= 16)
    {

        double start1 = clock();//编码计时开始时间
        vector<vector<unsigned int> > bitTI;//存储编码后的DNA序列数据集

        //对所有长度为l的字符串进行编码（A-00，C-01，G=10，T=11）
        ET.encodingSeqs(seqTI, bitTI, mc);

        double end1 = clock();//编码计时结束时间
        cout << "encodingSeqs take time :" << (end1-start1)/CLOCKS_PER_SEC << endl;//输出编码时间


        //MAX表示模式树mc.h层的节点数量,即初始公共前缀的数量
        unsigned int MAX = pow(4, mc.H) - 1;

        double start_totaltime = clock();//验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的开始时间
        double end_totaltime;//验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的结束时间

        //根据mc.q,分为两种情况
        if(mc.q>0 & mc.q<=0.4)
        {
            //公共前缀层mc.H必须小于mc.l-3
            if(mc.l - mc.H > 3)
            {
                for(unsigned int i = 0; i <= MAX; ++i)//遍历公共前缀
                {
                    //找出公共前缀下方的（l，d）模体
                    checkDownThreeLocation(i, bitTI, mc, result, mc.H);

                    //输出验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的运行时间
                    if(i%1000==0 && i != 0)
                    {
                        end_totaltime = clock();
                        cout<<"deep Search "<<i-1000<<"-"<<i<<" node take time:"<<(end_totaltime - start_totaltime)/CLOCKS_PER_SEC <<", resINT length:"<<result.size()<< endl;
                        start_totaltime = clock();
                    }

                }
            }
            else
            {
                cout<<"Please enter the correct start layer, which should be less than l"<<endl;
            }
        }
        else if(mc.q>0.4 & mc.q<=1)
        {
            if(mc.l - mc.H >= 3)
            {
                for(unsigned int i = 0; i <= MAX; ++i) //遍历公共前缀
                {
                    //找出公共前缀下方的（l，d）模体
                    checkDownThree(i, bitTI, mc, result, mc.H);

                    //输出验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的运行时间
                    if(i%1000==0 && i!= 0)
                    {
                        end_totaltime = clock();
                        cout<<"deep Search "<<i-1000<<"-"<<i<<" node take time:"<<(end_totaltime - start_totaltime)/CLOCKS_PER_SEC <<", resINT length:"<<result.size()<< endl;
                        start_totaltime = clock();
                    }
                }
            }
            else
            {
                cout<<"Please enter the correct start layer, which should be less than l"<<endl;
            }


        }

    }
    else if(mc.l >= 17 && mc.l <= 32)
    {
        double start1 = clock();//编码计时开始时间
        vector<vector<unsigned long long> > bitTI; //存储编码后的DNA序列数据集

        // 对所有长度为l的字符串进行编码
        ET.encodingSeqs(seqTI, bitTI, mc);
        double end1 = clock();//编码计时结束时间
        cout << "encodingSeqs take time :" << (end1-start1)/CLOCKS_PER_SEC << endl;//输出编码时间

        //MAX表示模式树mc.h层的节点数量,即初始公共前缀的数量
        unsigned long long MAX = pow(4, mc.H) - 1;


        double start_totaltime = clock();//验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的开始时间
        double end_totaltime;//验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的结束时间


        if(mc.q>0 & mc.q<=0.4)
        {
            if(mc.l-mc.H > 3)
            {
                for(unsigned int i = 0; i <= MAX; ++i)
                {
                    //找出公共前缀下方的（l，d）模体
                    checkDownThreeLocation(i, bitTI, mc, result, mc.H);

                    //输出验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的运行时间
                    if(i%1000==0 && i!= 0)
                    {
                        end_totaltime = clock();
                        cout<<"deep Search "<<i-1000<<"-"<<i<<" node take time:"<<(end_totaltime - start_totaltime)/CLOCKS_PER_SEC <<", resINT length:"<<result.size()<< endl;
                        start_totaltime = clock();
                    }

                }
            }
            else
            {
                cout<<"Please enter the correct start layer, which should be less than l"<<endl;
            }
        }
        else if(mc.q>0.4 & mc.q<=1)
        {
            if(mc.l-mc.H >= 3)
            {
                for(unsigned int i = 0; i <= MAX; ++i)
                {
                    //找出公共前缀下方的（l，d）模体
                    checkDownThree(i, bitTI, mc, result, mc.H);

                    //输出验证1000个公共前缀下方的候选模体前缀是否是（l，d）模体的运行时间
                    if(i%1000==0 && i!= 0)
                    {
                        end_totaltime = clock();
                        cout<<"deep Search "<<i-1000<<"-"<<i<<" node take time:"<<(end_totaltime - start_totaltime)/CLOCKS_PER_SEC <<", resINT length:"<<result.size()<< endl;
                        start_totaltime = clock();
                    }
                }
            }
            else
            {
                cout<<"Please enter the correct start layer, which should be less than or equal to l"<<endl;
            }

        }

    }
    else
    {
        cout << "Please input appropriate l, 0 < l ≤ 32 !" << endl;
    }
}

/*
函数：checkDownThreeLocation()
作用：验证公共前缀下方64个候选模体前缀的有效性，并记录64个候选模体前缀在序列中的出现位置（模体实例）
针对情况：1、mc.l-mc.h大于3
          2、mc.q小于等于0.4
          3、mc.l小于等于16
*/
void PMSMotif::checkDownThreeLocation(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;

    vector<unsigned int> res_temp;//res_temp临时存放合适的候选模体

    array<unsigned int,64> EPMTs_test = {};//EPMTs_test存储64个待验证候选模体前缀在数据集中出现的次数（以序列为单位计数，即在序列中存在模体实例，计数加1）

    vector<vector<vector<unsigned int> > > location(64,vector< vector<unsigned int> >(mc.t));//location存储64个待验证候选模体前缀在数据集中出现位置
    unsigned long long int EPMTs_temp_test;//EPMTs_temp_test记录64个待验证候选模体前缀在某条序列中的出现情况
    unsigned int temp_bitTI;
    unsigned int strLittle2;//strLittle2表示每个长为l的字符串中和公共前缀等长的部分（从头开始截取）
    int datumHam;//datumHam表示strLittle2和blockFunElement的汉明距离
    unsigned int resOne;
    unsigned int blockFunElement_new;//blockFunElement_new表示新的公共前缀
    int pre_len = cur_len+3;//cur_len表示公共前缀长度，pre_len表示当前待验证候选模体前缀长度
    int max_EPMT;
    int first;//first表示最后一个字符的二进制编码（针对待验证候选模体前缀）
    int second;//second表示倒数第二个字符的二进制编码（针对待验证候选模体前缀）
    int third;//third表示倒数第三个字符的二进制编码（针对待验证候选模体前缀）
    int moveloc = 2*(mc.l - pre_len);//由长为l的字符串变化为和待验证候选模体前缀等长字符串需要移动的位数（操作对象为二进制编码）


    for (unsigned int j = 0; j < bitTI.size(); ++j)//遍历序列
    {

        EPMTs_temp_test = 0;//初始化

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)//遍历序列i中的长为l的字符串
        {

            temp_bitTI = bitTI[j][k] >> moveloc;//从字符串中截取和待验证候选模体前缀长度相同的字符串
            strLittle2 = temp_bitTI >> 6;//strLittle2表示每个长为l的字符串中和公共前缀等长的部分（从头开始截取）


            //计算blockFunElement和strLittle2的汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
                //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer30[(temp_bitTI & 63)];//将最后3个字符和temp_bitTI最后3个字符相同的候选模体前缀的EPMTs_temp_test对应位置置为1
                location[(temp_bitTI & 63)][j].push_back(bitTI[j][k]);//记录候选模体前缀的实例信息

            }
            else if(datumHam == mc.d - 1)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer31[temp_bitTI & 63];
                first = temp_bitTI & 3;
                second = (temp_bitTI & 15) >> 2;
                third = (temp_bitTI & 63) >> 4;

                //记录候选模体前缀的实例信息
                for(int f1 = 0; f1 < 4 ; ++f1)
                {
                    //改变first
                    location[f1 + second*4 + third*16][j].push_back(bitTI[j][k]);
                    //改变second
                    location[first + f1*4 + third*16][j].push_back(bitTI[j][k]);
                    //改变third
                    location[first + second*4 + f1*16][j].push_back(bitTI[j][k]);

                }

            }
            else if(datumHam == mc.d - 2)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer32[temp_bitTI & 63];

                first = temp_bitTI & 3;
                second = (temp_bitTI & 15) >> 2;
                third = (temp_bitTI & 63) >> 4;

                for(int s1 = 0; s1 < 4 ; ++s1)
                {
                    for(int t1 = 0; t1 < 4; ++t1)
                    {
                        //first固定
                        location[first + s1*4 + t1*16][j].push_back(bitTI[j][k]);

                        //second固定
                        location[s1 + second*4 + t1*16][j].push_back(bitTI[j][k]);

                        //third固定
                        location[s1 +t1*4 + third*16][j].push_back(bitTI[j][k]);
                     }
                }
            }
            else if(datumHam <= mc.d - 3)
            {
                EPMTs_temp_test = 18446744073709551615;//64位二进制数均为1，对应的十进制数

                for(int s1 = 0; s1 < 64 ; ++s1)
                {
                    location[s1][j].push_back(bitTI[j][k]);
                }

            }
            else
            {
                cout<<"程序出错"<<endl;
            }

        }



        max_EPMT = 0;

        for(int i = 0; i < 64; ++i)
        {
            EPMTs_test[i] += ((EPMTs_temp_test >> i) & 1 );//计算完公共前缀和序列si中位点的汉明距离，更新EPMTs_test[i]（nk）

            //获得EPMTs_test中数值最大的元素
            if(max_EPMT < EPMTs_test[i])
            {
                max_EPMT = EPMTs_test[i];
            }

            //验证的候选模体前缀有效
            if(EPMTs_test[i] >= mc.MINSTISFYSEQ)
            {
                //将候选模体前缀添加到res_temp数组中
                resOne = (blockFunElement << 6) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }
        }
        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }

    }


    //获得新的待验证的候选模体前缀
    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        for(int i=0; i<res_temp.size(); ++i)
        {
            blockFunElement_new = res_temp[i];

            if(mc.l - pre_len >= 3)
            {
                checkDownThree(blockFunElement_new, location[blockFunElement_new&63], mc, result, pre_len);
            }
            else if(mc.l - pre_len == 2)
            {
                checkDownRight(blockFunElement_new, location[blockFunElement_new&63], mc, result, pre_len);
            }
            else
            {
                checkDown(blockFunElement_new,location[blockFunElement_new&63], mc, result, pre_len);
            }
        }
    }
}

/*
重载函数：checkDownThreeLocation()
作用：验证公共前缀下方64个候选模体前缀的有效性，并记录64个候选模体前缀在序列中的出现位置（模体实例）
针对情况：1、mc.l-mc.h大于3
          2、mc.q小于等于0.4
          3、mc.l大于16
*/
void PMSMotif::checkDownThreeLocation(unsigned long long int blockFunElement, vector<vector<unsigned long long int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned long long int> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,64> EPMTs_test = {};
    vector<vector<vector<unsigned long long int> > > location(64,vector< vector<unsigned long long int> >(mc.t));
    unsigned long long int EPMTs_temp_test;
    unsigned long long int temp_bitTI;
    unsigned long long int strLittle2;
    int datumHam;
    unsigned long long int resOne;
    unsigned long long int blockFunElement_new;
    int pre_len = cur_len+3;
    int max_EPMT;
    int moveloc = 2*(mc.l - pre_len);

    int first;
    int second;
    int third;
        for (unsigned int j = 0; j < bitTI.size(); ++j)
    {

       // 临时存放扩展表
        EPMTs_temp_test = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {


            temp_bitTI = bitTI[j][k] >> moveloc;
             //bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 6;


            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


           // 如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
               // 取后两位字符
                EPMTs_temp_test = EPMTs_temp_test | mer30[(temp_bitTI & 63)];
                location[(temp_bitTI & 63)][j].push_back(bitTI[j][k]);

            }
            else if(datumHam == mc.d - 1)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer31[temp_bitTI & 63];
                first = temp_bitTI & 3;
                second = (temp_bitTI & 15) >> 2;
                third = (temp_bitTI & 63) >> 4;

                for(int f1 = 0; f1 < 4 ; ++f1)
                {
                    //改变first
                    location[f1 + second*4 + third*16][j].push_back(bitTI[j][k]);
                    //改变second
                    location[first + f1*4 + third*16][j].push_back(bitTI[j][k]);
                    //改变third
                    location[first + second*4 + f1*16][j].push_back(bitTI[j][k]);
//                    cout<<ET.deEncodingLmer(f1 + second*4 + third*16, 3, ALPHABETSET)<<" "<<ET.deEncodingLmer(first + f1*4 + third*16 , 3, ALPHABETSET) <<" "<<ET.deEncodingLmer(first + second*4 + f1*16 , 3, ALPHABETSET)<<endl;

                }

            }
            else if(datumHam == mc.d - 2)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer32[temp_bitTI & 63];

                first = temp_bitTI & 3;
                second = (temp_bitTI & 15) >> 2;
                third = (temp_bitTI & 63) >> 4;

                for(int s1 = 0; s1 < 4 ; ++s1)
                {
                    for(int t1 = 0; t1 < 4; ++t1)
                    {
                        //first固定
                        location[first + s1*4 + t1*16][j].push_back(bitTI[j][k]);

                        //second固定
                        location[s1 + second*4 + t1*16][j].push_back(bitTI[j][k]);

                        //third固定
                        location[s1 +t1*4 + third*16][j].push_back(bitTI[j][k]);
                        //cout<<ET.deEncodingLmer(first + s1*4 + t1*16, 3, ALPHABETSET)<<" "<<ET.deEncodingLmer(s1 + second*4 + t1*16 , 3, ALPHABETSET) <<" "<<ET.deEncodingLmer(s1 +t1*4 + third*16 , 3, ALPHABETSET)<<endl;

                    }
                }
            }
            else if(datumHam <= mc.d - 3)
            {
                EPMTs_temp_test = 18446744073709551615;

                for(int s1 = 0; s1 < 64 ; ++s1)
                {
                    location[s1][j].push_back(bitTI[j][k]);
                }

            }
            else
            {
                cout<<"程序出错"<<endl;
            }

        }



        max_EPMT = 0;

        for(int i = 0; i < 64; ++i)
        {
            EPMTs_test[i] += ((EPMTs_temp_test >> i) & 1 );
            if(max_EPMT < EPMTs_test[i])
            {
                max_EPMT = EPMTs_test[i];
            }
            if(EPMTs_test[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+3，d）模体
                resOne = (blockFunElement << 6) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }
        }
        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }
    }

    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        for(int i=0; i<res_temp.size(); ++i)
        {
            blockFunElement_new = res_temp[i];

            if(mc.l - pre_len >= 3)
            {
                checkDownThree(blockFunElement_new, location[blockFunElement_new&63], mc, result, pre_len);
            }
            else if(mc.l - pre_len == 2)
            {
                checkDownRight(blockFunElement_new, location[blockFunElement_new&63], mc, result, pre_len);
            }
            else
            {
                checkDown(blockFunElement_new,location[blockFunElement_new&63], mc, result, pre_len);
            }
        }
    }
}

/*
函数：checkDownThree()
作用：验证公共前缀下方64个候选模体前缀的有效性
针对情况：1、mc.l-mc.h大于等于3
          2、mc.l小于等于16
*/
void PMSMotif::checkDownThree(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned int> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,64> EPMTs_test = {};
    unsigned long long int EPMTs_temp_test;
    unsigned int temp_bitTI;
    unsigned int strLittle2;
    int datumHam;
    unsigned int resOne;
    unsigned int blockFunElement_new;
    int pre_len = cur_len+3;
    int max_EPMT;
    int count38 = 0;

    int moveloc = 2*(mc.l - pre_len);

    for (unsigned int j = 0; j < bitTI.size(); ++j)
    {

       // 临时存放扩展表
        EPMTs_temp_test = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {

            temp_bitTI = bitTI[j][k] >> moveloc;
             //bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 6;

            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


           // 如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
               // 取后两位字符
                EPMTs_temp_test = EPMTs_temp_test | mer30[(temp_bitTI & 63)];
            }
            else if(datumHam == mc.d - 1)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer31[temp_bitTI & 63];

            }
            else if(datumHam == mc.d - 2)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer32[temp_bitTI & 63];

            }
            else if(datumHam <= mc.d - 3)
            {
                EPMTs_temp_test = 18446744073709551615;
                break;
            }
            else
            {
                cout<< "Program Error!"<<endl;
            }

        }



        max_EPMT = 0;
        for(int i = 0; i < 64; ++i)
        {
            EPMTs_test[i] += (EPMTs_temp_test >> i) & 1 ;

            if(max_EPMT < EPMTs_test[i])
            {
                max_EPMT = EPMTs_test[i];
            }

            if(EPMTs_test[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+3，d）模体
                resOne = (blockFunElement << 6) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }
        }

        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }

        if(res_temp.size() >= 64)
        {
            break;
        }
    }



//    cout<<blockFunElement<<" ";
//    for(int i = 0; i < 64; ++i)
//    {
//        cout<< EPMTs_test[i]<<" ";
//    }
//    cout<<endl;

    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        for(int i=0; i<res_temp.size(); ++i)
        {
            blockFunElement_new = res_temp[i];

            if(mc.l - pre_len >= 3)
            {
                checkDownThree(blockFunElement_new, bitTI, mc, result, pre_len);
            }
            else if(mc.l - pre_len == 2)
            {
                checkDownRight(blockFunElement_new, bitTI, mc, result, pre_len);
            }
            else
            {
                checkDown(blockFunElement_new, bitTI, mc, result, pre_len);
            }
        }
    }

}

/*
重载函数：checkDownThree()
作用：验证公共前缀下方64个候选模体前缀的有效性
针对情况：1、mc.l-mc.h大于等于3
          2、mc.l大于16
*/
void PMSMotif::checkDownThree(unsigned long long int blockFunElement, vector<vector<unsigned long long int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned long long int> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,64> EPMTs_test = {};
    unsigned long long int EPMTs_temp_test;
    unsigned long long int temp_bitTI;
    unsigned long long int strLittle2;
    int datumHam;
    unsigned long long int resOne;
    unsigned long long int blockFunElement_new;
    int pre_len = cur_len+3;
    int max_EPMT;

    int moveloc = 2*(mc.l - pre_len);


    for (unsigned int j = 0; j < bitTI.size(); ++j)
    {

       // 临时存放扩展表
        EPMTs_temp_test = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {
            temp_bitTI = bitTI[j][k] >> moveloc;
             //bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 6;


            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);

           // 如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
               // 取后两位字符
                EPMTs_temp_test = EPMTs_temp_test | mer30[(temp_bitTI & 63)];

            }
            else if(datumHam == mc.d - 1)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer31[temp_bitTI & 63];

            }
            else if(datumHam == mc.d - 2)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer32[temp_bitTI & 63];

            }
            else if(datumHam <= mc.d - 3)
            {
                EPMTs_temp_test = 18446744073709551615;
                break;

            }
            else
            {
                cout<<"程序出错"<<endl;
            }

        }



        max_EPMT = 0;

        for(int i = 0; i < 64; ++i)
        {
            EPMTs_test[i] += (EPMTs_temp_test >> i) & 1 ;
            if(max_EPMT < EPMTs_test[i])
            {
                max_EPMT = EPMTs_test[i];
            }
            if(EPMTs_test[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+3，d）模体
                resOne = (blockFunElement << 6) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }
        }
        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }

        if(res_temp.size() >= 64)
        {
            break;
        }

    }

    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        for(int i=0; i<res_temp.size(); ++i)
        {
            blockFunElement_new = res_temp[i];

            if(mc.l - pre_len >= 3)
            {
                checkDownThree(blockFunElement_new, bitTI, mc, result, pre_len);
            }
            else if(mc.l - pre_len == 2)
            {
                checkDownRight(blockFunElement_new, bitTI, mc, result, pre_len);
            }
            else
            {
                checkDown(blockFunElement_new, bitTI, mc, result, pre_len);
            }
        }
    }
}

/*
函数：checkDownRight()
作用：验证公共前缀下方16个候选模体前缀的有效性
针对情况：1、mc.l-mc.h=2
          2、mc.l小于等于16
*/
void PMSMotif::checkDownRight(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned int> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,16> EPMTs_test = {};
    unsigned int EPMTs_temp_test;
    unsigned int temp_bitTI;
    unsigned int strLittle2;
    int datumHam;
    unsigned int resOne;
    unsigned int blockFunElement_new;
    int pre_len = cur_len+2;
    int max_EPMT;

    int moveloc = 2*(mc.l - pre_len);
    for (unsigned int j = 0; j < bitTI.size(); ++j)
    {

       // 临时存放扩展表
        EPMTs_temp_test = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {
            temp_bitTI = bitTI[j][k] >> moveloc;
             //bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 4;


            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


           // 如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
               // 取后两位字符
                EPMTs_temp_test = EPMTs_temp_test | (1 << (temp_bitTI & 15));
            }
            else if(datumHam == mc.d - 1)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer[temp_bitTI & 15];
            }
            else if(datumHam <= mc.d - 2)
            {
                EPMTs_temp_test = 65535;
                break;
            }
            else
            {
                cout<<"Program Error!"<<endl;
            }
        }

        max_EPMT = 0;

        for(int i = 0; i < 16; ++i)
        {
            EPMTs_test[i] += (EPMTs_temp_test >> i) & 1 ;
            if(max_EPMT < EPMTs_test[i])
            {
                max_EPMT = EPMTs_test[i];
            }
            if(EPMTs_test[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+2，d）模体
                resOne = (blockFunElement << 4) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }
        }

        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }
        if(res_temp.size() >= 16)
        {
            break;
        }

    }


    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        cout<<"Program Error!"<<endl;
    }

}

/*
重载函数：checkDownRight()
作用：验证公共前缀下方16个候选模体前缀的有效性
针对情况：1、mc.l-mc.h=2
          2、mc.l小于等于16
*/

void PMSMotif::checkDownRight(unsigned long long blockFunElement, vector<vector<unsigned long long> >& bitTI, ParaConfig mc, vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned long long> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,16> EPMTs_test = {};
    unsigned int EPMTs_temp_test;
    unsigned long long temp_bitTI;
    unsigned long long strLittle2;
    int datumHam;
    unsigned long long resOne;
    unsigned long long blockFunElement_new;
    int pre_len = cur_len+2;
    int max_EPMT;

    int moveloc = 2*(mc.l - pre_len);
    for (unsigned int j = 0; j < bitTI.size(); ++j)
    {

       // 临时存放扩展表
        EPMTs_temp_test = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {
            temp_bitTI = bitTI[j][k] >> moveloc;
             //bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 4;


            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


           // 如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
               // 取后两位字符
                EPMTs_temp_test = EPMTs_temp_test | (1 << (temp_bitTI & 15));
            }
            else if(datumHam == mc.d - 1)
            {
              //或运算 1 | 0 = 1 、1 | 1 = 1、 0 | 0 = 0
                EPMTs_temp_test = EPMTs_temp_test | mer[temp_bitTI & 15];
            }
            else if(datumHam <= mc.d - 2)
            {
                EPMTs_temp_test = 65535;
                break;
            }
            else
            {
                cout<<"程序出错"<<endl;
            }

        }

        max_EPMT = 0;

        for(int i = 0; i < 16; ++i)
        {
            EPMTs_test[i] += (EPMTs_temp_test >> i) & 1 ;

            if(max_EPMT < EPMTs_test[i])
            {
                max_EPMT = EPMTs_test[i];
            }

            if(EPMTs_test[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+2，d）模体
                resOne = (blockFunElement << 4) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }
        }

        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }
        if(res_temp.size() >= 16)
        {
            break;
        }

    }
    if (pre_len >= mc.l)
    {
        for(int i = 0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        cout<<"Program Error!"<<endl;
    }

}

/*
函数：checkDown()
作用：验证公共前缀下方4个候选模体前缀的有效性
针对情况：1、mc.l-mc.h=1
          2、mc.l小于等于4
*/
void PMSMotif::checkDown(unsigned int blockFunElement, vector<vector<unsigned int> >& bitTI, ParaConfig mc,  vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned int> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,4> EPMTs = {};
    unsigned int EPMTs_temp;
    unsigned int temp_bitTI;
    unsigned int strLittle2;
    int datumHam;
    unsigned int resOne;
    int pre_len = cur_len + 1;
    int max_EPMT;

    int moveloc = 2*(mc.l - pre_len);
    for (unsigned int j = 0; j < bitTI.size(); ++j)
    {
        //临时存放扩展表
        EPMTs_temp = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {

            temp_bitTI = bitTI[j][k] >> moveloc;
            // bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 2;


            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


            //如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
                //取后两位字符
                EPMTs_temp = EPMTs_temp | (1 << (temp_bitTI & 3));
            }
            else if(datumHam <= mc.d - 1)
            {
                EPMTs_temp = 15;
                break;
            }
            else
            {
                cout<<"程序出错"<<endl;
            }
        }

        max_EPMT = 0;
        for(int i = 0; i < 4; ++i)
        {
            EPMTs[i] = EPMTs[i] + ((EPMTs_temp >> i) & 1 );
            if(max_EPMT < EPMTs[i])
            {
                max_EPMT = EPMTs[i];
            }

            if(EPMTs[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+2，d）模体
                resOne = (blockFunElement << 2) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }

        }

        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }
        if(res_temp.size() >= 4)
        {
            break;
        }

    }


    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        cout<<"error: probe one layer!"<<endl;
    }

}

/*
重载函数：checkDown()
作用：验证公共前缀下方4个候选模体前缀的有效性，
针对情况：1、mc.l-mc.h=1
          2、mc.l大于4
*/
void PMSMotif::checkDown(unsigned long long blockFunElement, vector<vector<unsigned long long> >& bitTI, ParaConfig mc,  vector<string>& result, int cur_len)
{
    HamDistance HD;
    EncodingTools ET;
    //res_temp临时存放合适的候选模体
    vector<unsigned long long> res_temp;

    //为每一个节点的每一行维护临时的一个扩张模式匹配表
    array<unsigned int,4> EPMTs = {};
    unsigned int EPMTs_temp;
    unsigned long long temp_bitTI;
    unsigned long long strLittle2;
    int datumHam;
    unsigned long long resOne;
    int pre_len = cur_len + 1;
    int max_EPMT;

    int moveloc = 2*(mc.l - pre_len);
    for (unsigned int j = 0; j < bitTI.size(); ++j)
    {
        //临时存放扩展表
        EPMTs_temp = 0;

        for (unsigned int k = 0; k < bitTI[j].size(); ++k)
        {

            temp_bitTI = bitTI[j][k] >> moveloc;
            // bitTI[j][k]长度为H+2，与长度为H的字符串进行比较，作为计算的基准
            strLittle2 = temp_bitTI >> 2;


            //基b准汉明距离
            datumHam = HD.getHamDistance(blockFunElement, strLittle2, cur_len);


            //如果候选模体与位点之间的距离大于d，意味着改组与该位点都没有必要进行计算，直接跳跃16个位置，进行下一个位点的比对
            if(datumHam > mc.d)
            {
                continue;
            }
            else if(datumHam == mc.d)
            {
                //取后两位字符
                EPMTs_temp = EPMTs_temp | (1 << (temp_bitTI & 3));
            }
            else if(datumHam <= mc.d - 1)
            {
                EPMTs_temp = 15;
                break;
            }
            else
            {
                cout<<"程序出错"<<endl;
            }
        }

        max_EPMT = 0;
        for(int i = 0; i < 4; ++i)
        {
            EPMTs[i] = EPMTs[i] + ((EPMTs_temp >> i) & 1 );
            if(max_EPMT < EPMTs[i])
            {
                max_EPMT = EPMTs[i];
            }

            if(EPMTs[i] >= mc.MINSTISFYSEQ)
            {
                //添加新的（H+2，d）模体
                resOne = (blockFunElement << 2) + i;
                if(find(res_temp.begin(), res_temp.end(), resOne) == res_temp.end())
                {
                    res_temp.push_back(resOne);
                }
            }

        }

        if((max_EPMT + mc.t - j - 1) < mc.MINSTISFYSEQ)
        {
            break;
        }
        if(res_temp.size() >= 4)
        {
            break;
        }

    }


    if (pre_len >= mc.l)
    {
        for(int i =0; i < res_temp.size(); ++i)
        {
            result.push_back(ET.deEncodingLmer(res_temp[i], mc.l, ALPHABETSET));
        }
    }
    else
    {
        cout<<"error: probe one layer!"<<endl;
    }
}


#endif
