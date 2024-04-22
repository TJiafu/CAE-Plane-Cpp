#ifndef MASHDATA_H
#define MASHDATA_H



#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//#pragma pack(4)
using namespace std;
//定义节点
struct Node
{
    int tag;
    double x, y, z;
};

//定义单元，用结构图便于扩展
struct EleNode
{
    int tag, node;
};

struct EleLine
{
    int tag;
    int vecn[2];
};

struct ElePlane4N
{
    int tag;
    int vecn[4];
};


class MashData
{
public:
    MashData(const string& fname);
    bool readFile();

    void printMash();

    vector<Node> vecN;
    vector<EleNode> vecEN;
    vector<EleLine> vecEL;
    vector<ElePlane4N> vecP4N;

private:
    int countNode;
    int countEle;
    const string filename;

    vector<string> split(string& str, const string& pattern);
};

#endif // !MASHDATA_H
