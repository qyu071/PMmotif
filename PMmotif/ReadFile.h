
#include "stdafx.h"

#include <fstream>

class ReadFile
{

public:
    ReadFile();
    ~ReadFile();

public:

    //按顺序返回字符集
    void readSequence(const char* fileName, vector<string>& strings);

private:

    //按行读取文件
    bool read(const char* fileName, vector<string>& strings);

    //清除字符串的两头
    void trim(string& s);
};
