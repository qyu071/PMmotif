//
//  GetFilePathName.h
//  PTMotif1.1
//
//  Created by apple on 18/6/4.
//  Copyright © 2018年 star. All rights reserved.
//

#include "stdafx.h"

#include<io.h>

class GetFilePathName
{

public:

    GetFilePathName(void);

    ~GetFilePathName(void);

    /*
     @Param
     文件夹路径：directory
     文件名路径：fileVector
     @return
     文件夹里所有的文件名(包含路径)：fileVectorle
     */
    void getFiles(string directory,vector<string> &fileVector);
};
