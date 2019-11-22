//
// Created by wasku on 19-11-20.
//
#include <LBFGS.h>
#include <iostream>
#include "test11/Jfunc.h"

using namespace LBFGSpp;
using namespace std;


typedef Eigen::Vector3d Geo3d;

Eigen::Vector3d remapPoint(Eigen::Vector3d& pos,Eigen::Matrix4d& compMat) {
    Eigen::Vector4d vec4d(pos[0],pos[1],pos[2],1);
    vec4d = compMat*vec4d;
    return vec4d.head<3>();
}

void test_9(){
    /**
     * 1. 构建测试数据
     * 1.1 构建测试用的旋转矩阵，平移矩阵，缩放矩阵
     * 1.2 根据1.1将一个三角形顶点信息经过矩阵映射(缩放、旋转、平移)得到新的三角形边长数据
     * 2. 根据前后的三角形信息恢复矩阵变换
     *
     * */



    double rx = M_PI/2;
    double ry = M_PI/5;
    double rz = M_PI/3;


    double scale = 2;
    Vector3d p_a (0,0,0);
    Vector3d p_b (4,0,0);
    Vector3d p_c (4,3,0);
    Vector3d p_d (0,3,0);
    Vector4d p_test (4,5,0,1);

    /**************************************************************************/

    Eigen::Matrix4d rotateMat,scaleMat,offsetMat,compMat;

    Eigen::Vector3d offsetVec(30,12,15);
    offsetMat.setIdentity();
    offsetMat.block<3,1>(0,3) = offsetVec;

    scaleMat.setIdentity();
    scaleMat*=scale;
    scaleMat(3,3) = 1;

    Eigen::AngleAxisd::QuaternionType quat = Eigen::AngleAxisd(rz,Geo3d::UnitZ())
                                             *Eigen::AngleAxisd(ry,Geo3d::UnitY())
                                             *Eigen::AngleAxisd(rx,Geo3d::UnitX());
    auto rotMat = quat.matrix();
    rotateMat.setIdentity();
    rotateMat.block<3,3>(0,0) = rotMat;

    compMat = offsetMat*rotateMat*scaleMat;


    Vector3d p_A = remapPoint(p_a,compMat);
    Vector3d p_B = remapPoint(p_b,compMat);
    Vector3d p_C = remapPoint(p_c,compMat);
    Vector3d p_D = remapPoint(p_d,compMat);

    std::cout<<"the p_A-p_D is:\n"<<p_A<<std::endl<<std::endl
             <<p_B<<std::endl<<std::endl
             <<p_C<<std::endl<<std::endl
             <<p_D<<std::endl<<std::endl;




    Vector3d v_AB = p_B-p_A;
    Vector3d v_AC = p_C-p_A;
    Vector3d v_AD = p_D-p_A;
    Vector3d v_BC = p_C-p_B;
    Vector3d v_BD = p_D-p_B;
    Vector3d v_CD = p_D-p_C;
    double l_AB = v_AB.norm();
    double l_AC = v_AC.norm();
    double l_AD = v_AD.norm();
    double l_BC = v_BC.norm();
    double l_BD = v_BD.norm();
    double l_CD = v_CD.norm();

    std::cout<<"the length is:\n"
             <<l_AB<<std::endl
             <<l_AC<<std::endl
             <<l_AD<<std::endl
             <<l_BC<<std::endl
             <<l_BD<<std::endl
             <<l_CD<<std::endl;

    Vector3d pro_A = p_A/p_A[2];
    Vector3d pro_B = p_B/p_B[2];
    Vector3d pro_C = p_C/p_C[2];
    Vector3d pro_D = p_D/p_D[2];

    std::cout<<"the pro_A-pro_D is:\n"<<pro_A<<std::endl<<std::endl
             <<pro_B<<std::endl<<std::endl
             <<pro_C<<std::endl<<std::endl
             <<pro_D<<std::endl<<std::endl;



    std::cout<<"the compMat is:\n"<<compMat<<std::endl;
    std::cout<<"the compMat*p_test is:\n"<<compMat*p_test<<std::endl;

    std::cout<<"the rotateMat4 is:\n"<<rotateMat<<std::endl;


    Vector3d pc_ = rotMat*p_c;

    std::cout<<pc_<<std::endl;

    std::cout<<rotMat<<std::endl;
}


void test10(){
    //calibrate front camera

    //left
    /*
       237.5328  515.4148       1
       155.24864 466.6159       0
       263.27258 390.5183       3
       359.22452 405.14447      2
       */

    //right
    /*
       1097.5 503.5
       1011.5 541.5
       903.0 425.0
       999.5 419.5
       */
    std::vector<LinePairData> vecLPdata;
    std::vector<Vector2d> points;
    points.emplace_back(2,0.8);
    points.emplace_back(3.22752,1.70961);
    points.emplace_back(2.30992,1.36348);
    points.emplace_back(1.59984,0.758242);

    std::vector<Vector2d> distancePoints;
    distancePoints.emplace_back(155.24864,720 - 466.6159);
    distancePoints.emplace_back(237.5328,720 - 515.4148);
    distancePoints.emplace_back(359.22452,720 - 405.14447);
    distancePoints.emplace_back(263.27258,720 - 390.5183);

    std::vector<double> vec_distances;
    vec_distances.push_back(8);
    vec_distances.push_back(10);
    vec_distances.push_back(6);
    vec_distances.push_back(6);
    vec_distances.push_back(10);
    vec_distances.push_back(8);

    int tmp = 0;
    for(int i=0;i<distancePoints.size();++i){
        for(int j=i+1;j<distancePoints.size();++j){
            LinePairData linePairData{};
            linePairData.nIndex1 = i;
            linePairData.nIndex2 =j;
            linePairData.distance = vec_distances[tmp++];
            vecLPdata.push_back(linePairData);
        }
    }


/********************************************************************/

    Jfunc jfunc(points,vecLPdata);
    const int n = 6;
    LBFGSParam<double> param;
    LBFGSSolver<double> solver(param);

    VectorXd x = VectorXd::Zero(n);
    x[0]=15;
    x[1]=10;
    x[2]=15;
    x[3]=20;
    x[4]=0;
    x[5]=0;
    double fx;
    int niter = solver.minimize(jfunc, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cout<<"end of test10...\n";
}


void test10_1(){
    //calibrate front camera

    //left
    /*
       237.5328  515.4148       1
       155.24864 466.6159       0
       263.27258 390.5183       3
       359.22452 405.14447      2
       */

    //right
    /*
       1097.5 503.5
       1011.5 541.5
       903.0 425.0
       999.5 419.5
       */
    std::vector<LinePairData> vecLPdata;
    std::vector<Vector2d> points;
    points.emplace_back(2,0.8);
    points.emplace_back(3.22752,1.70961);
    points.emplace_back(2.30992,1.36348);
    points.emplace_back(1.59984,0.758242);

    for(int i=0;i<points.size();++i){
        points[i] *=1.000019;
    }

    std::vector<Vector2d> distancePoints;
    distancePoints.emplace_back(155.24864,720 - 466.6159);
    distancePoints.emplace_back(237.5328,720 - 515.4148);
    distancePoints.emplace_back(359.22452,720 - 405.14447);
    distancePoints.emplace_back(263.27258,720 - 390.5183);

    std::vector<double> vec_distances;
    vec_distances.push_back(8);
    vec_distances.push_back(10);
    vec_distances.push_back(6);
    vec_distances.push_back(6);
    vec_distances.push_back(10);
    vec_distances.push_back(8);

    int tmp = 0;
    for(int i=0;i<distancePoints.size();++i){
        for(int j=i+1;j<distancePoints.size();++j){
            LinePairData linePairData{};
            linePairData.nIndex1 = i;
            linePairData.nIndex2 =j;
            linePairData.distance = vec_distances[tmp++];
            vecLPdata.push_back(linePairData);
        }
    }


/********************************************************************/

    Jfunc jfunc(points,vecLPdata);
    const int n = 6;
    LBFGSParam<double> param;
    LBFGSSolver<double> solver(param);

    VectorXd x = VectorXd::Zero(n);
    x[0]=15;
    x[1]=10;
    x[2]=15;
    x[3]=20;
    x[4]=0;
    x[5]=0;
    double fx;
    int niter = solver.minimize(jfunc, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cout<<"end of test10_1...\n";
}

int main(){
    test_9();
    cout<<"begin test10...\n";
    test10();
    cout<<"begin test10_1...\n";
    test10_1();
    cout<<"end...\n";
    return 0;
}