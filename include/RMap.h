//
// Created by wasku on 19-8-26.
//

#ifndef EIGENSAMPLE_RMAP_H
#define EIGENSAMPLE_RMAP_H
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>


/**
 *
 * f12 = [m1(t1_x - p0_x)-m2(t2_x - p0_x)]^² +
 *       [m1(t1_y - p0_y)-m2(t2_y - p0_y)]^² +
 *       [m1(t1_z - p0_z)-m2(t2_z - p0_z)]^² +
 *       - L12^² = 0
 *
 * d(p0_x)f12 = 2[m1(t1_x - p0_x)-m2(t2_x - p0_x)]*(m2-m1) ;
 * d(p0_y)f12 = 2[m1(t1_y - p0_y)-m2(t2_y - p0_y)]*(m2-m1) ;
 * d(p0_z)f12 = 2[m1(t1_z - p0_z)-m2(t2_z - p0_z)]*(m2-m1) ;
 *
 * d(m1)f12   = 2[m1(t1_x - p0_x)-m2(t2_x - p0_x)]*(t1_x - p0_x) +
 *              2[m1(t1_y - p0_y)-m2(t2_y - p0_y)]*(t1_y - p0_y) +
 *              2[m1(t1_z - p0_z)-m2(t2_z - p0_z)]*(t1_z - p0_z) ;
 *
 * 注意这里有个负号
 * -d(m2)f12   = 2[m1(t1_x - p0_x)-m2(t2_x - p0_x)]*(t2_x - p0_x) +
 *               2[m1(t1_y - p0_y)-m2(t2_y - p0_y)]*(t2_y - p0_y) +
 *               2[m1(t1_z - p0_z)-m2(t2_z - p0_z)]*(t2_z - p0_z) ;
 *
 *
 * J(m,p) = f12^² + f23^² + f13^² +...=0;
 *
 *dJ(m1) = 2*f12*d(m1)f12 + 2*f13*d(m1)f13 ;
 *
 * dJ(p0_x) = 2f12*d(p0_x)f12 + 2*f23*d(p0_x)f23 + ... ;
 *
 *
 *
 *
 *
 * */


class Dfun{
public:
    typedef double (*tmFunc)(int argc,double* argv);

    //m1(t1_x - p0_x)-m2(t2_x - p0_x)类似的形式
private:
    static double tfunc(int argc,double* argv){
        assert(argc==5);
        double m1_ = argv[0];
        double m2_ = argv[1];
        double t1_dim_i = argv[2];
        double t2_dim_i = argv[3];
        double p0_dim_i = argv[4];

        return m1_*(t1_dim_i-p0_dim_i) - m2_*(t2_dim_i - p0_dim_i);
    }

private:
    double m1,m2;
    Eigen::Vector3d t1,t2,p0;
    tmFunc m_func;

public:
    Dfun(double m1_,double m2_,Eigen::Vector3d t1_,Eigen::Vector3d t2_):m1(m1_),
    m2(m2_),
    t1(t1_),
    t2(t2_){
//        tmFunc = &Dfun::tfunc;
        m_func = &tfunc;
    }

    void setOriPos(Eigen::Vector3d p0_){
        p0 = p0_;
    }

    double dm1();
    double dm2();

    double dp0_x(){
        double a[5] ={m1,m2,t1(0),t2(0),p0(0)};
        double r = m_func(5,&a[0]);
        return 2*r*(m2-m1);
    }
    double dp0_y(){
        double a[5] ={m1,m2,t1(1),t2(1),p0(1)};
        double r = m_func(5,a);
        return 2*r*(m2-m1);
    }
    double dp0_z(){
        double a[5] ={m1,m2,t1(1),t2(1),p0(1)};
        double r = m_func(5,a);
        return 2*r*(m2-m1);
    }



};


class RMap {
private:
    Eigen::Vector3d t1,t2;

public:
    RMap(Eigen::Vector3d t1_,Eigen::Vector3d t2_):t1(t1_),t2(t2_){}
};


#endif //EIGENSAMPLE_RMAP_H
