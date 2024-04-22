#include "solution.h"


Solution::Solution(const MashData& md,
                   const double E,				//材料.弹模(pa)
                   const double prxy,			//材料.泊松比
                   const double dens,			//材料.体密度(Kg/m^3)
                   const double t,				//实常数.板厚(m)
                   const double g,				//(体力)重力加速度(N/Kg)
                   const vector<BFI>& LQxy,	//力边界（面力）（pa）
                   const vector<BDI>& vBDI)	//位移边界
    :_md(md), _E(E),_prxy(prxy),_dens(dens),_t(t),_g(g),_LQxy(LQxy),_vBDI(vBDI)
{
    /******* D矩阵，材料本构关系（平面应力问题） *******/
    _D = (Matrix3d() << 1, prxy, 0,
          prxy, 1, 0,
          0, 0, (1 - prxy) / 2).finished() * E / (1 - pow(prxy, 2));
}

void Solution::element_K(int ENom)  //ENom = 0、1、2、3、4....
{
    _Ke = MatrixXd::Zero(8, 8);
    Matrix<double, 4, 2> E_coor;	//每个单元的结点坐标
    Matrix<double, 2, 4> dNab;		//单元中每个高斯点处，插值函数矩阵的导数的值
    Matrix2d J;						//雅可比矩阵（某高斯点处）
    Matrix2d invJ;					//雅可比矩阵的逆矩阵（某高斯点处）
    Vector2d dNxyi;					//插值函数对物理坐标的导数的值
    Matrix<double, 3, 8> B;			//应变矩阵（某高斯点处）
    Matrix<double, 8, 8> Fab;		//高斯积分被积函数矩阵，（某高斯点处）

    for (int i = 0; i < 4; i++)
    {
        E_coor(i, 0) = _md.vecN[_md.vecP4N[ENom].vecn[i]].x;
        E_coor(i, 1) = _md.vecN[_md.vecP4N[ENom].vecn[i]].y;
    }
    //cout << E_coor << endl;
    //cout << _md.vecP4N[1].tag << endl;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            dNab = (Matrix<double, 2, 4>() << GPcoor[j] - 1, 1 - GPcoor[j], GPcoor[j] + 1, -GPcoor[j] - 1,
                    GPcoor[i] - 1, -GPcoor[i] - 1, GPcoor[i] + 1, 1 - GPcoor[i]).finished() / 4;
            J = (dNab * E_coor);
            invJ = J.inverse();
            //cout << "dNab = " << endl;
            //cout << dNab << endl;
            //cout << "E_coor = " << endl;
            //cout << E_coor << endl;
            //cout << "J = " << endl;
            //cout << J << endl;
            for (int k = 0; k < 4; k++)
            {
                dNxyi = invJ * dNab.col(k);
                B(all, seq(k * 2, k * 2 + 1)) = (Matrix<double, 3, 2>() <<
                                                     dNxyi[0], 0, 0, dNxyi[1], dNxyi[1], dNxyi[0]).finished();
            }
            //cout << "B = " << endl;
            //cout << B << endl;
            Fab = B.transpose() * _D * B * _t * J.determinant();
            _Ke += GPweight[i] * GPweight[j] * Fab;
        }
    }
//    cout<<_md.vecN[_md.vecEL[0].vecn[0]].x;
    //cout << "_Ke = " << endl;
    //cout << _Ke << endl;
}


void Solution::global_K()
{
    int rows = (_md.vecN.size() - 1) * 2;
    int Bandwidth = getBandwidth();
    int size = _md.vecP4N.size();
    _K.resize(rows,rows);
    _K.reserve(VectorXi::Constant(rows, Bandwidth));
    ElePlane4N thisElem;

    for (int i = 1; i < size; i++)
    {
        element_K(i);
        thisElem = _md.vecP4N[i];
        int row, column;
        Matrix2d value;
        for (int p = 0; p < 4; p++)
        {
            for (int q = 0; q < 4; q++)
            {
                row = thisElem.vecn[p];
                column = thisElem.vecn[q];
                value = _Ke.block<2, 2>(p*2, q*2);
                for(int m=0;m<2;m++)
                {
                    for (int n = 0; n < 2; n++) {
                        _K.coeffRef((row - 1) * 2 + m, (column - 1) * 2 + n) += value(m, n);
                    }
                }
            }
        }
    }
    _K.makeCompressed();
}


void Solution::body_Pe(int ENom)
{
    _Pbe = MatrixXd::Zero(8, 1);
    Matrix<double, 4, 2> E_coor;	//每个单元的结点坐标
    Matrix<double, 8, 2> G;			//变换矩阵
    Matrix3d E_coorforA1;			//将单元分成两个三角形，求其中一个面积的方阵
    Matrix3d E_coorforA2;			//求另一个面积的方阵
    double Ae;						//单元面积
    Vector2d gxy{ {0},{-_dens * _g} };	//体力（(gx,gy) (N/m^3)）

    for (int i = 0; i < 4; i++)
    {
        E_coor(i, 0) = _md.vecN[_md.vecP4N[ENom].vecn[i]].x;
        E_coor(i, 1) = _md.vecN[_md.vecP4N[ENom].vecn[i]].y;
    }
    G << 1,0, 0,1, 1,0, 0,1, 1,0, 0,1, 1,0, 0,1;
    E_coorforA1 << E_coor.block<3,2>(0, 0), Vector3d::Constant(1);
    E_coorforA2 << E_coor.block<3,2>(1, 0), Vector3d::Constant(1);
    //cout << "E_coorforA1 = " << endl;
    //cout <<	E_coorforA1 << endl;
    //cout << "E_coorforA2 = " << endl;
    //cout <<	E_coorforA2 << endl;
    Ae = (E_coorforA1.determinant() + E_coorforA2.determinant()) / 2;
    _Pbe =  _t * Ae * G * gxy/4;
    //cout << "Ae=" << Ae << endl;
    //cout << "G=" << endl;
    //cout << G << endl;
    //cout << "gxy=" << endl;
    //cout << gxy << endl;
    //cout << "_Pbe=" << endl;
    //cout << _Pbe << endl;
}


void Solution::boundary_Fs()
{
    int NoN = _md.vecN.size() - 1;
    _Fs.resize(NoN * 2, 1);
    _Fs.reserve(VectorXi::Constant(NoN*2, 1));
    //_Fs = MatrixXd::Zero(NoN*2, 1);

    for (int i = 0; i < _LQxy.size();i++)
    {
        double Qadd = (_LQxy[i].fj - _LQxy[i].fi) / (_LQxy[i].ELj - _LQxy[i].ELi + 1); // 面力沿单元增长的步长
        Vector2d Qex = Vector2d::Zero();
        Vector2d Fes = Vector2d::Zero();
        Matrix2d le_coor = Matrix2d::Zero();
        double ls;
        for (int j = _LQxy[i].ELi; j <= _LQxy[i].ELj; j++)							//第i个x方向面力的第j个受力单元
        {
            Qex(0) = _LQxy[i].fi + Qadd * (j - _LQxy[i].ELi);
            Qex(1) = _LQxy[i].fi + Qadd * (j - _LQxy[i].ELi + 1);
            le_coor(0, 0) = _md.vecN[_md.vecEL[j].vecn[0]].x;
            le_coor(0, 1) = _md.vecN[_md.vecEL[j].vecn[0]].y;
            le_coor(1, 0) = _md.vecN[_md.vecEL[j].vecn[1]].x;
            le_coor(1, 1) = _md.vecN[_md.vecEL[j].vecn[1]].y;
            ls = pow(pow((le_coor(1, 1) - le_coor(0, 1)), 2) + pow((le_coor(1, 0) - le_coor(0, 0)), 2), 0.5);

            Fes(0) = ls * _t * (Qex(0) / 2 + (Qex(1) - Qex(0)) / 6);
            Fes(1) = ls * _t * (Qex(0) / 2 + (Qex(1) - Qex(0)) / 3);

            if (_LQxy[i].dir == 'x')
            {
                _Fs.coeffRef((_md.vecEL[j].vecn[0] - 1) * 2,0) += Fes(0);
                _Fs.coeffRef((_md.vecEL[j].vecn[1] - 1) * 2,0) += Fes(1);
            }
            else
            {
                _Fs.coeffRef((_md.vecEL[j].vecn[0] - 1) * 2+1, 0) += Fes(0);
                _Fs.coeffRef((_md.vecEL[j].vecn[1] - 1) * 2+1, 0) += Fes(1);
            }
        }
    }
    _Fs.makeCompressed();
    //cout << _Fs << endl;
}


void Solution::global_F()
{

    int NoN = _md.vecN.size() - 1;		//结点数
    int NoE = _md.vecP4N.size();		//单元数+1
    _F = VectorXd::Zero(NoN*2);
    ElePlane4N thisElem;
    int rowe = 0;
    for (int i = 1; i < NoE; i++)		//组装Fbe
    {
        body_Pe(i);
        thisElem = _md.vecP4N[i];
        for (int j = 0; j < 4; j++)
        {
            rowe = thisElem.vecn[j];
            _F((rowe - 1) * 2) = _Pbe(j * 2);
            _F((rowe - 1) * 2 + 1) = _Pbe(j * 2 + 1);
        }
    }
    //cout << _F << endl;
    boundary_Fs();
    _F += _Fs;
    //cout << "after add _Fs:" << endl;
    //cout << _F << endl;
}


void Solution::solve()
{
    //乘大数法引入位移边界
    int NoN = _md.vecN.size() - 1;
    //NodeDis = VectorXd::Zero(NoN * 2);

    global_K();
    KtoExcel();
    global_F();
    //cout << _F << endl;
    SparseMatrix<double> _KK = _K;
    MatrixXd _FF = _F;

    double maxval = 0, minval = 0,alpha=0;
    int Nindex;

    for (int i = 0; i < _K.outerSize(); i++)
    {
        for (SparseMatrix<double>::InnerIterator it(_K,i);it;++it)
        {
            if (it.value() > maxval)
                maxval = it.value();
            if (it.value() < minval)
                minval = it.value();
        }
    }
    alpha = (maxval > -minval ? maxval : -minval) * 10000;
    for (int i = 0; i < _vBDI.size(); i++)
    {
        Nindex = _vBDI[i].node;
        for (int j = 0; j < 2; j++)
        {
            if (_vBDI[i].style[j])
            {
                _FF((Nindex - 1) * 2 + j) = alpha * _vBDI[i].disp[j] * _KK.coeff((Nindex - 1) * 2 + j, (Nindex - 1) * 2 + j);
                _KK.coeffRef((Nindex - 1) * 2 + j, (Nindex - 1) * 2 + j) =
                    alpha * _K.coeff((Nindex - 1) * 2 + j, (Nindex - 1) * 2 + j);
            }
        }
    }
    ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
    NodeDis = solver.compute(_KK).solve(_FF);
    cout << NodeDis << endl;        //输出变形（结点位移）
}


int Solution::getBandwidth()
{
    int maxtag, mintag,thistag;
    int size = _md.vecP4N.size();
    int Bandwidth = 0;
    for (int i = 1; i < size; i++)
    {
        maxtag = _md.vecP4N[i].vecn[0];
        mintag = _md.vecP4N[i].vecn[0];
        for (int j = 1; j < 4; j++)
        {
            thistag = _md.vecP4N[i].vecn[j];
            if (maxtag < thistag)
                maxtag = thistag;

            if (mintag > thistag)
                mintag = thistag;
        }
        //cout << maxtag << "--" << mintag << endl;

        if (Bandwidth < (maxtag - mintag + 1) * 2)
            Bandwidth = (maxtag - mintag + 1) * 2;
    }
    return Bandwidth;
}


void Solution::KtoExcel()
{
    ofstream outfile("Matrix_K.xls");
    for (int i = 0; i < _K.rows(); i++)
    {
        for (int j = 0; j < _K.cols(); j++)
        {
            outfile << _K.coeff(i,j) <<'\t';
        }
        outfile << endl;
    }
    outfile.close();
}










