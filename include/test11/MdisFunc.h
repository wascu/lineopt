//
// Created by wasku on 19-11-11.
//

#ifndef EIGENSAMPLE_MDISFUNC_H
#define EIGENSAMPLE_MDISFUNC_H


#include <eigen3/Eigen/Eigen>



using namespace Eigen;
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
class MdisFunc {
public:
    MdisFunc(const Vector2d &pnt1, const Vector2d &pnt2, int nIdx1, int nIdx2, double distance);

    void show();


public:
    int getS1() {
        return m_index_S1;
    }

    int getS2() {
        return m_index_S2;
    }

    bool containderivative(int n) {
        return n == m_index_S1 || n == m_index_S2 || false;
    }

    double calc_func_value(const VectorXd &x);

    double calc_derivative_s1(const VectorXd &x);

    double calc_derivative_s2(const VectorXd &x);

    double calc_derivative_m_x(const VectorXd &x);

    double calc_derivative_m_y(const VectorXd &x);

    double calc_derivative_m_z(const VectorXd &x);

private:
    void updateParams(const VectorXd &x) {
        int length = x.size();
        m_x = x[length - 2];
        m_y = x[length - 1];
//        m_z = x[length - 1];
        s1 = x[m_index_S1];
        s2 = x[m_index_S2];
    }

private:
    Vector2d m_P1,m_P2;
    int m_index_S1 = 0;
    int m_index_S2 = 0;
    double m_distance = 0;

    double m_x = 0;
    double m_y = 0;
//    double m_z = 0;
    double s1 = 0;
    double s2 = 0;


};


#endif //EIGENSAMPLE_MDISFUNC_H
