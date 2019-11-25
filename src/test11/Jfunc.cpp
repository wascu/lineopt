//
// Created by wasku on 19-11-12.
//

#include "../../include/test11/Jfunc.h"


#include <iostream>

using namespace std;
/**
 *
 *（L12为T1、T2的点间距 ...）
 * (其中m_x,m_y 为点P1,P2等在近截面(Z=1)上的偏移量,  P1_z=P2_z=1)
 * f12 = (S1*(P1_x+m_x) -S2*(P2_x+m_x))² +(S1*(P1_y+m_y) -S2*(P2_y+m_y))² +(S1*P1_z -S2*P2_z)² - L12² = 0;
 *
 *
 *
 * J(S1,S2,...,Sn,m_x,m_y) = f12²+f13+f14+...+f1n+
 *                                f23+f24+...+f2n+
 *                                    f34+...+f3n
 *                                       .
 *                                       .
 *                                       .
 *                                       = 0
 * */
void Jfunc::show() {
    cout << "show from Jfunc...\n";
}

Jfunc::Jfunc(const std::vector<Vector2d> &points, const std::vector<LinePairData> &vecLPdata) {
    LinePairData tmpData;
    for (int i = 0; i < vecLPdata.size(); ++i) {
        tmpData = vecLPdata[i];
        m_vec_funcs.push_back(
                MdisFunc(points[tmpData.nIndex1], points[tmpData.nIndex2], tmpData.nIndex1, tmpData.nIndex2,
                         tmpData.distance));
    }

    std::cout<<"the size of vector MdisFunc is: "<<m_vec_funcs.size()<<std::endl;

}

double Jfunc::operator()(const VectorXd &x, VectorXd &grad) {
    double error = 0;
    double tmp = 0;
    int nVars = m_vec_funcs.size();
    for (int i = 0; i < nVars-2; ++i) {
        tmp = m_vec_funcs[i].calc_func_value(x);
        error+= tmp*tmp;

        grad[i] = calc_derivative_Sn(i,x);
    }
    grad[nVars-2] = 0;
    grad[nVars-1] = 0;
    std::cout<<"the current error is: "<<error<<std::endl;
    return error;
}

double Jfunc::calc_derivative_Sn(int n,const VectorXd &x) {
    double r = 0;
    for(int i = 0;i<m_vec_funcs.size();++i){
        if(n==m_vec_funcs[i].getS1()){
            r+= 2*m_vec_funcs[i].calc_func_value(x)
                    *m_vec_funcs[i].calc_derivative_s1(x);
        } else if(n==m_vec_funcs[i].getS2()){
            r+= 2*m_vec_funcs[i].calc_func_value(x)
                *m_vec_funcs[i].calc_derivative_s2(x);
        }
    }
    return r;
}
