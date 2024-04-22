#pragma once
#ifndef SOLUTION_H
#define SOLUTION_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mashdata.h"
#include <math.h>
using namespace Eigen;

/******* 存放力边界条件的输入参数 *******/
typedef struct BoundaryForceInput
{

    char dir;
    int ELi, ELj;				//最小精度：单元边上线性变化
    double fi, fj;
} BFI;

/******* 存放位移边界条件的输入参数(结点位移约束) *******/
typedef struct BoundaryDisplacementInput
{
    int node;
    bool style[2];		//[1,0] -> 1:有约束；0：没有约束
    double disp[2];		//[Dx,Dy] (m)
} BDI;

/******* Guess Point(2P) *******/
const double GPcoor[2] = {-1/pow(3,0.5),1 / pow(3,0.5) };
const double GPweight[2] = { 1,1 };


/******* 2维四边形四节点单元求解器 *******/
class Solution		//静力求解
{
public:
    Solution(const MashData& md,
             const double E,				//材料.弹模(pa)
             const double prxy,			//材料.泊松比
             const double dens,			//材料.体密度(Kg/m^3)
             const double t,				//实常数.板厚(m)
             const double g,				//(体力)重力加速度(N/Kg)
             const vector<BFI>& LQxy,
             const vector<BDI>& vBDI);	//力边界（面力）（pa）

    //求解器计算接口
    void solve();					//返回位移向量

private:
    //输入参数存储
    const MashData _md;
    const double _E;
    const double _prxy;
    const double _dens;
    const double _t;
    const double _g;
    const vector<BFI> _LQxy;
    const vector<BDI> _vBDI;
    Matrix3d _D;

    //临时变量
    MatrixXd _Ke;
    SparseMatrix<double> _K;
    MatrixXd _Pbe;
    SparseMatrix<double> _Fs;
    MatrixXd _F;
    MatrixXd NodeDis;
public:
    //内部操作函数


    void element_K(int Nom);
    void global_K();
    void body_Pe(int ENom);
    void boundary_Fs();
    void global_F();

    int getBandwidth();
    void KtoExcel();

};

#endif
