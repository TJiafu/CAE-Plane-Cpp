#include "mashdata.h"


MashData::MashData(const string& fname)
    :countNode(0),countEle(0),filename(fname)
{}

bool MashData::readFile()
{
    ifstream fin(filename);
    if (fin.is_open() == false)
    {
        cout << "open file: " << filename << " failed" << endl;
        return false;
    }

    //过滤前面无用数据
    string buf;
    getline(fin, buf);
    while (buf != "$Nodes")
    {
        getline(fin, buf);
    }
    getline(fin, buf);
    countNode = stoi(buf);

    //读取结点数据
    vector<string> vecs;
    Node newN{0,0.0,0.0,0.0};
    vecN.push_back(newN);		//结点号与索引号一致
    getline(fin, buf);
    while (buf != "$EndNodes")
    {
        vecs = split(buf, " ");
        newN.tag = stoi(vecs[0]);
        newN.x = stod(vecs[1]);
        newN.y = stod(vecs[2]);
        newN.z = stod(vecs[3]);

        vecN.push_back(newN);
        getline(fin, buf);
    }

    //读取单元数据
    getline(fin, buf);			//跳过“$Elements”
    getline(fin, buf);			//读取"单元数"
    countEle = stoi(buf);

    getline(fin, buf);

    EleNode newEN{0,0};
    EleLine newEL{0,{0,0}};
    ElePlane4N newEP{0,{0,0,0,0}};
    vecEN.push_back(newEN);		//保证单元号与索引号一致
    vecEL.push_back(newEL);
    vecP4N.push_back(newEP);
    while (buf != "$EndElements")
    {
        vecs = split(buf, " ");

        if (vecs[1] == "3")
        {
            newEP.tag = stoi(vecs[0]);
            newEP.vecn[0] = stoi(vecs[5]);
            newEP.vecn[1] = stoi(vecs[6]);
            newEP.vecn[2] = stoi(vecs[7]);
            newEP.vecn[3] = stoi(vecs[8]);
            vecP4N.push_back(newEP);
            getline(fin, buf);
        }
        else if (vecs[1] == "1")
        {
            newEL.tag = stoi(vecs[0]);
            newEL.vecn[0] = stoi(vecs[5]);
            newEL.vecn[1] = stoi(vecs[6]);
            vecEL.push_back(newEL);
            getline(fin, buf);
        }
        else if (vecs[1] == "15")
        {
            newEN.tag = stoi(vecs[0]);
            newEN.node = stoi(vecs[5]);
            vecEN.push_back(newEN);
            getline(fin, buf);
        }
        else
        {
            return false;
        }
    }
    fin.close();
    return true;
}

void MashData::printMash()
{
    if (countEle!=0&&countNode!=0)
    {
        vector<Node> myNode = vecN;
        for (int i = 0; i < myNode.size(); i++)
        {
            cout << myNode[i].tag << '-' << myNode[i].x << '-' << myNode[i].y << '-' << myNode[i].z << endl;
        }

        vector<EleNode> myEN = vecEN;
        for (int i = 0; i < myEN.size(); i++)
        {
            cout << myEN[i].tag << "--" << myEN[i].node << endl;
        }

        vector<EleLine> myEL = vecEL;
        for (int i = 0; i < myEL.size(); i++)
        {
            cout << myEL[i].tag << "--" << myEL[i].vecn[0] << "--" << myEL[i].vecn[1] << endl;
        }

        vector<ElePlane4N> myP4N = vecP4N;
        for (int i = 0; i < myP4N.size(); i++)
        {
            cout << myP4N[i].tag << "--" << myP4N[i].vecn[0] << "--" << myP4N[i].vecn[1]
                 << "--" << myP4N[i].vecn[2] << "--" << myP4N[i].vecn[3] << endl;
        }
    }
    else
    {
        cout << "读取文件" << filename << "失败" << endl;
    }
}


vector<string> MashData::split(string& str, const string& pattern)
{
    string::size_type pos;
    vector<string> result;
    str += pattern;				//扩展字符串以方便操作
    int size = str.size();
    for (int i = 0; i < size; i++)
    {
        pos = str.find(pattern, i);
        if (pos < size)
        {
            string s = str.substr(i, pos - i);
            result.push_back(s);
            i = pos + pattern.size() - 1;
        }
    }
    return result;
}
