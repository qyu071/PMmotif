#include "stdafx.h"
#include <math.h>

class HamDistance
{
public:
    HamDistance();
    ~HamDistance();

public:
    /**
     计算汉明距离（推介使用）
     @Param
        编码后的int1
        编码后的int2
     @return
        汉明距离
     */
    int getHamDistance(unsigned int int1, unsigned int int2, int l);

    //计算汉明距离
    int getHamDistance_2(unsigned long long int1, unsigned long long int2, int l);
private:
    /**
     汉明距离查询表
     */

};
