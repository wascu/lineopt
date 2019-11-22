//
// Created by wasku on 19-11-11.
//

#include "test11/MdisFunc.h"

#include <iostream>

using namespace std;

/**
 *
 *（L12为T1、T2的点间距 ...）
 * (其中m_x,m_y 为点P1,P2等在近截面(Z=1)上的偏移量,  P1_z=P2_z=1)
 * f12 = (S1*(P1_x+m_x) -S2*(P2_x+m_x))² +(S1*(P1_y+m_y) -S2*(P2_y+m_y))² +(S1*(P1_z) -S2*(P2_z)² - L12² = 0;
 *
 *
 * f12 = (S1*(P1_x+m_x) -S2*(P2_x+m_x))² +(S1*(P1_y+m_y) -S2*(P2_y+m_y))² +(S1*(P1_z+m_z) -S2*(P2_z+m_z))² - L12² = 0;
 *
 * */

void MdisFunc::show() {
    cout << "show from MdisFunc...\n";
}

MdisFunc::MdisFunc(const Vector2d &pnt1, const Vector2d &pnt2, int nIdx1, int nIdx2, double distance) {
    m_P1 = pnt1;
    m_P2 = pnt2;
    m_index_S1 = nIdx1;
    m_index_S2 = nIdx2;
    m_distance = distance;
}

//double MdisFunc::calc_derivative_s1(const VectorXd &x) {
//    updateParams(x);
//    double part_1 = s1 * (m_P1[0] + m_x) - s2 * (m_P2[0] + m_x);
//    double part_2 = s1 * (m_P1[1] + m_y) - s2 * (m_P2[1] + m_y);
//    double part_3 = s1 * (1 + m_z) - s2 * (1 + m_z);
//    double r = 2 * part_1 * (m_P1[0] + m_x);
//    r += 2 * part_2 * (m_P1[1] + m_y);
//    r += 2 * part_3 * (1 + m_z);
//    return r;
//}

double MdisFunc::calc_derivative_s1(const VectorXd &x) {
    updateParams(x);
    double part_1 = s1 * (m_P1[0] + m_x) - s2 * (m_P2[0] + m_x);
    double part_2 = s1 * (m_P1[1] + m_y) - s2 * (m_P2[1] + m_y);
    double part_3 = s1 - s2 ;
    double r = 2 * part_1 * (m_P1[0] + m_x);
    r += 2 * part_2 * (m_P1[1] + m_y);
    r += 2 * part_3;
    return r;
}

double MdisFunc::calc_derivative_s2(const VectorXd &x) {
    updateParams(x);
    double part_1 = s1 * (m_P1[0] + m_x) - s2 * (m_P2[0] + m_x);
    double part_2 = s1 * (m_P1[1] + m_y) - s2 * (m_P2[1] + m_y);
    double part_3 = s1 - s2 ;
    double r = -2 * part_1 * (m_P2[0] + m_x);
    r -= 2 * part_2 * (m_P2[1] + m_y);
    r -= 2 * part_3;
    return r;
}

double MdisFunc::calc_func_value(const VectorXd &x) {
    updateParams(x);
    //S1*(P1_x+m_x) -S2*(P2_x+m_x)
    double part_1 = s1 * (m_P1[0] + m_x) - s2 * (m_P2[0] + m_x);
    double part_2 = s1 * (m_P1[1] + m_y) - s2 * (m_P2[1] + m_y);
    double part_3 = s1  - s2 ;
    double square_part_1 = part_1 * part_1;
    double square_part_2 = part_2 * part_2;
    double square_part_3 = part_3 * part_3;
    double sum = square_part_1 + square_part_2 + square_part_3 - m_distance * m_distance;

    return sum;
}

double MdisFunc::calc_derivative_m_x(const VectorXd &x) {
    int length = x.size();
    double m_x = x[length - 3];
    double s1 = x[m_index_S1];
    double s2 = x[m_index_S2];
    //S1*(P1_x+m_x) -S2*(P2_x+m_x)
    double part_1 = s1 * (m_P1[0] + m_x) - s2 * (m_P2[0] + m_x);
    double r = 2 * part_1 * (s1 - s2);

    return r;
}

double MdisFunc::calc_derivative_m_y(const VectorXd &x) {
    int length = x.size();
    double m_y = x[length - 2];
    double s1 = x[m_index_S1];
    double s2 = x[m_index_S2];
    double part_2 = s1 * (m_P1[1] + m_y) - s2 * (m_P2[1] + m_y);
    double r = 2 * part_2 * (s1 - s2);
    return r;
}

//double MdisFunc::calc_derivative_m_z(const VectorXd &x) {
//    int length = x.size();
//    double m_z = x[length - 1];
//    double s1 = x[m_index_S1];
//    double s2 = x[m_index_S2];
//    double part_2 = s1 * (m_P1[1] + m_y) - s2 * (m_P2[1] + m_y);
//    double r = 2 * part_2 * (s1 - s2);
//    return r;
//}