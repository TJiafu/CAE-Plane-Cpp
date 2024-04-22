#include "mashdata.h"
#include "solution.h"


/******* 前处理（参数输入） *******/
const double E = 200*pow(10,9);	//材料.弹模(pa)
const double prxy = 0.3;			//材料.泊松比
const double dens = 7850;			//材料.体密度(Kg/m^3)

const double t = 0.01;				//实常数.板厚(m)

const string filename = "beam1.msh";//网格划分

const double g = 9.8;				//(体力)重力加速度(N/Kg)

const BFI L1_Qx{ 'x',116,215,000000,000000 };
const BFI L1_Qy{ 'y',116,215,-1000000,-1000000 }; //边界面力,{direction,fyi, fyj, ELi(Pa), ELj(Pa)}

const BDI D1{ 1, true,true,0,0 };
const BDI D2{ 3, false,true,0,0 };	//简支


int main()
{
//    const string filename = "test5.msh";//网格划分
    MashData md(filename);
    md.readFile();
//    md.printMash();

    /******* D矩阵，材料本构关系（平面应力问题） *******/
    const Matrix3d D = (Matrix3d() << 1, prxy, 0,
                        prxy, 1, 0,
                        0, 0, (1 - prxy) / 2).finished() * E / (1 - pow(prxy, 2));

    /******* 力边界 *******/
    vector<BFI> vecBFI;			//将边界面力整合传参
    vecBFI.push_back(L1_Qx);
    vecBFI.push_back(L1_Qy);
    /******* 位移边界 *******/
    vector<BDI> vecBDI;
    vecBDI.push_back(D1);
    vecBDI.push_back(D2);

    Solution s1(md,E,prxy,dens,t,g,vecBFI,vecBDI);
    s1.element_K(1);

    //cout << s1.getBandwidth() << endl;
    //s1.global_K();
    //s1.KtoExcel();

    s1.solve();

    return 0;
}




