//
// Created by wasku on 19-11-12.
//

#ifndef EIGENSAMPLE_JFUNC_H
#define EIGENSAMPLE_JFUNC_H
#include "MdisFunc.h"

#include <eigen3/Eigen/Eigen>
#include <vector>

using namespace Eigen;

struct LinePairData{
    int nIndex1;
    int nIndex2;
    double distance;
};


/**
 *
 *（L12为T1、T2的点间距 ...）
 * (其中m_x,m_y 为点P1,P2等在近截面(Z=1)上的偏移量,  P1_z=P2_z=1)
 * f12 = (S1*(P1_x+m_x) -S2*(P2_x+m_x))² +(S1*(P1_y+m_y) -S2*(P2_y+m_y))² +(S1*P1_z -S2*P2_z)² - L12² = 0;
 *
 *
 *
 * J(S1,S2,...,Sn,m_x,m_y) = f12²+f13+f14+...+f1n+
 *                             f23+f24+...+f2n+
 *                                 f34+...+f3n
 *                                       .
 *                                       .
 *                                       .
 *                                       = 0
 * */

class Jfunc {
    typedef Vector2d point2d_;
public:
    void show();
    Jfunc(){}
    Jfunc(const std::vector<Vector2d>& points,const std::vector<LinePairData>& vecLPdata);
    double operator()(const VectorXd& x,VectorXd& grad);

private:
    double calc_derivative_Sn(int n,const VectorXd &x);

private:
    std::vector<MdisFunc> m_vec_funcs;
};


#endif //EIGENSAMPLE_JFUNC_H
