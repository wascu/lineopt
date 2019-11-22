#include <iostream>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <LBFGS.h>

#include "test/TCallback.h"

#include "test/TDelegate.h"
#include "../../include/test11/Jfunc.h"
//#include "include/aws.h"

#include <string>
#include <cmath>


using namespace Eigen;
using namespace std;


void test0() {
    int n = 2;
    int m = 4;
    Eigen::VectorXd x(n), b(n), b_d(m);
    SparseMatrix<double> A(n, n);
    SparseMatrix<double> D(m, n);

    std::vector<Eigen::Triplet<double> > triplets_D;
    triplets_D.push_back({0, 0, 15});
    triplets_D.push_back({0, 1, 0});

    triplets_D.push_back({1, 0, 10});
    triplets_D.push_back({1, 1, 0});

    triplets_D.push_back({2, 0, 0});
    triplets_D.push_back({2, 1, 9});

    triplets_D.push_back({3, 0, 0});
    triplets_D.push_back({3, 1, 4.5});

    b_d[0] = 7;
    b_d[1] = 4.5;
    b_d[2] = 27;
    b_d[3] = 17;


    D.setFromTriplets(triplets_D.cbegin(), triplets_D.cend());
    std::cout << "D is:\n " << D << std::endl;

    SparseMatrix<double> _D = D.transpose();
    std::cout << "_D is:\n " << _D << std::endl;
    A = _D * D;
    std::cout << "A is:\n " << A << std::endl;


    Eigen::VectorXd vb = _D * b_d;



//    std::vector<Eigen::Triplet<double> > triplets;
//    triplets.push_back({0,0,4});
//    triplets.push_back({0,1,1});
//
//    triplets.push_back({1,0,1});
//    triplets.push_back({1,1,3});

//    A.setFromTriplets(triplets.cbegin(),triplets.cend());
    b[0] = 1;
    b[1] = 2;


    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;

    cg.compute(A);
    x = cg.solve(vb);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;

    std::cout << "x is: " << x << std::endl;


}

void test0_0() {
    int n = 2;
    int m = 6;
    Eigen::VectorXd x(n), b(n), b_d(m);
    SparseMatrix<double> A(n, n);
    SparseMatrix<double> D(m, n);

    std::vector<Eigen::Triplet<double> > triplets_D;
    triplets_D.push_back({0, 0, 15});
//    triplets_D.push_back({0,1,0});

    triplets_D.push_back({1, 0, 10});
//    triplets_D.push_back({1,1,0});

//    triplets_D.push_back({2,0,0});
    triplets_D.push_back({2, 1, 9});

//    triplets_D.push_back({3,0,0});
    triplets_D.push_back({3, 1, 4.5});

    triplets_D.push_back({4, 0, 24.7});
    triplets_D.push_back({4, 1, -4.1});

    triplets_D.push_back({5, 0, 17.8});
    triplets_D.push_back({5, 1, -3});


    b_d[0] = 7;
    b_d[1] = 4.5;
    b_d[2] = 27;
    b_d[3] = 17;
    b_d[4] = 0;
    b_d[5] = 0;


    D.setFromTriplets(triplets_D.cbegin(), triplets_D.cend());
    Eigen::SparseMatrix<double> blk = D.block(4, 0, 2, 2);
    std::cout << "the blk is:\n" << blk << std::endl;


    std::cout << "D is:\n " << D << std::endl;

    SparseMatrix<double> _D = D.transpose();
    std::cout << "_D is:\n " << _D << std::endl;
    A = _D * D;
    std::cout << "A is:\n " << A << std::endl;


    Eigen::VectorXd vb = _D * b_d;



//    std::vector<Eigen::Triplet<double> > triplets;
//    triplets.push_back({0,0,4});
//    triplets.push_back({0,1,1});
//
//    triplets.push_back({1,0,1});
//    triplets.push_back({1,1,3});

//    A.setFromTriplets(triplets.cbegin(),triplets.cend());
    b[0] = 1;
    b[1] = 2;


    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;

    cg.compute(A);
    x = cg.solve(vb);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;

    std::cout << "x is: " << x << std::endl;


}


void test1() {
    int n = 4;
    Eigen::VectorXd x(n), b(n);
    SparseMatrix<double> A(n, n);

    std::vector<Eigen::Triplet<double> > triplets;
    triplets.push_back({0, 0, 15});
    triplets.push_back({0, 1, 1e-12});
    triplets.push_back({0, 2, 1e-12});
    triplets.push_back({0, 3, 1e-12});

    triplets.push_back({1, 0, 10});
    triplets.push_back({1, 1, 1e-12});
    triplets.push_back({1, 2, 1e-12});
    triplets.push_back({1, 3, 1e-12});

    triplets.push_back({2, 0, 1e-12});
    triplets.push_back({2, 1, 1e-12});
    triplets.push_back({2, 2, 9});
    triplets.push_back({2, 3, 1e-12});

    triplets.push_back({3, 0, 1e-12});
    triplets.push_back({3, 1, 1e-12});
    triplets.push_back({3, 2, 6.5});
    triplets.push_back({3, 3, 1e-12});

    A.setFromTriplets(triplets.cbegin(), triplets.cend());
    b[0] = 7;
    b[1] = 4.5;
    b[2] = 27;
    b[3] = 17;


    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;

    cg.compute(A);
    x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;

    std::cout << "x is: " << x << std::endl;
// update b, and solve again
    x = cg.solve(b);
    std::cout << "22  #iterations:     " << cg.iterations() << std::endl;
    std::cout << "22  estimated error: " << cg.error() << std::endl;

    std::cout << "22  x is: " << x << std::endl;

//    A.setFromTriplets()
}

void test2() {
    int n = 2;
    Eigen::VectorXd x(n), b(n);
    SparseMatrix<double> A(n, n);

    std::vector<Eigen::Triplet<double> > triplets;
    triplets.push_back({0, 0, 4});
    triplets.push_back({0, 1, 1});

    triplets.push_back({1, 0, 1});
    triplets.push_back({1, 1, 3});

    A.setFromTriplets(triplets.cbegin(), triplets.cend());
    b[0] = 1;
    b[1] = 2;


    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;

    cg.compute(A);
    x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;

    std::cout << "x is: " << x << std::endl;

}

void test3() {
    double x1 = 818;
    double y1 = 159;

    double x2 = 735;
    double y2 = 160;

    double r = pow((x2 - x1), 2) + pow((y2 - y1), 2);
    r = sqrt(r);
    std::cout << "r is:" << r << std::endl;
}


void createOriData(std::vector<MatrixXd> &vec_r) {
    int row = 2;
    int col = 4;
    MatrixXd lm(row, col);
    MatrixXd fm(row, col);
    MatrixXd nm(row, col);
    MatrixXd rm(row, col);
    lm << 71.2609, 57.4591, 70.0207, 57.5016,
            76.542, 66.3955, 75.4821, 65.821;
    fm << 79.006, 62.000, 82.006, 62.07,
            83.006, 65.122, 85.006, 65.03;
    nm << 75.9519, 65.7169, 75.1111, 66.222,
            73.239, 58.3631, 72.7181, 57.384;
    rm << 56.3097, 71.7912, 56.374, 71.702,
            54.8468, 70.0997, 56.2216, 69.9497;

    vec_r.clear();
    vec_r.push_back(lm);
    vec_r.push_back(fm);
    vec_r.push_back(nm);
    vec_r.push_back(rm);


//    MatrixXd m1(row,col);
////    m1.setConstant(0);
//    m1.row(0)<<lm.row(0);
//    std::cout<<m1<<std::endl;

}

vector<string> split(const string &str, const string &pattern)
{
    vector<string> res;
    if(str == "")
        return res;
    //在字符串末尾也加入分隔符，方便截取最后一段
    string strs = str + pattern;
    size_t pos = strs.find(pattern);

    while(pos != strs.npos)
    {
        string temp = strs.substr(0, pos);
        res.push_back(temp);
        //去掉已分割的字符串,在剩下的字符串中进行分割
        strs = strs.substr(pos+1, strs.size());
        pos = strs.find(pattern);
    }

    return res;
}

SparseMatrix<double>  createOriDataBystr() {
    std::string str_vals="71.2609, 57.4591, 70.0207, 57.5016,"\
                         "76.542, 66.3955, 75.4821, 65.821,"\
                         "79.006, 62.000, 82.006, 62.07,"\
                         "83.006, 65.122, 85.006, 65.03,"\
                         "75.9519, 65.7169, 75.1111, 66.222,"\
                         "73.239, 58.3631, 72.7181, 57.384,"\
                         "56.3097, 71.7912, 56.374, 71.702,"\
                         "54.8468, 70.0997, 56.2216, 69.9497";

    vector<string> vec_str = split(str_vals,",");


    int row = 8;
    int col = 4;
    std::vector<Eigen::Triplet<double> > triplets_D;
    MatrixXd md(row ,col);
    string tm;
    for(int i=0;i<row;++i){
        for(int j=0;j<col;++j){
            tm =vec_str.at(col*i+j);
            triplets_D.push_back({i,j,atof(tm.c_str())});
        }
    }

    SparseMatrix<double> D(row, col);





    D.setFromTriplets(triplets_D.begin(),triplets_D.end());
    return D;
}



void remodelingData(std::vector<MatrixXd> &vecMats, std::vector<MatrixXd> &vec_rmol) {
    int row = 2;
    int col = 4;
    int v_size = vecMats.size();
    vec_rmol.clear();


    for (int i = 1; i < v_size; ++i) {
        MatrixXd t(row, col);
        t.row(0) = vecMats[i - 1].row(1);
        t.row(1) = -1 * vecMats[i].row(0);
        vec_rmol.push_back(t);
    }

    MatrixXd t(row, col);
    t.row(0) = vecMats[v_size - 1].row(1);
    t.row(1) = -1 * vecMats[0].row(0);
    vec_rmol.push_back(t);
}

void createSparseMat(std::vector<MatrixXd> &vecMatrs, SparseMatrix<double> &sp_mat) {
    SparseMatrix<double> tmp_sp;
    MatrixXd mat_xd(16, 4);

    int v_size = vecMatrs.size();
    for (int i = 0; i < v_size - 1; ++i) {
        tmp_sp = vecMatrs[i].sparseView();
        tmp_sp = tmp_sp.transpose();
        mat_xd.block(4 * i, i, 4, 2) = tmp_sp;
    }
    tmp_sp = vecMatrs[v_size - 1].sparseView();
    tmp_sp = tmp_sp.transpose();
    mat_xd.block(12, 3, 4, 1) = tmp_sp.col(0);
    mat_xd.block(12, 0, 4, 1) = tmp_sp.col(1);
    sp_mat = mat_xd.sparseView();
}

void getData() {
    int row = 4;
    int col = 2;

    SparseMatrix<double> l_mat;
    SparseMatrix<double> f_mat;
    SparseMatrix<double> n_mat;
    SparseMatrix<double> r_mat;

    MatrixXd t(col, row);


    t << 76.542, 66.3955, 75.4821, 65.821,
            79.006, 62.000, 82.006, 62.07;
    t.row(0) *= -1;

    SparseMatrix<double, ColMajor> csp = t.sparseView();
    std::cout << "the csp is:\n" << csp << std::endl;
    SparseMatrix<double> ctsp = csp.transpose();
    std::cout << "the ctsp is:\n" << ctsp << std::endl;
}

void test4_0(){
    SparseMatrix<double> nsp = createOriDataBystr();
    SparseMatrix<double> nsp_T = nsp.transpose();


    //构造16X4矩阵
    SparseMatrix<double> spMatr(16,4);

    MatrixXd mmxd(16,4);
    mmxd.setZero();

    for(int i=0;i<3;++i){
        mmxd.block(4*i,i,4,1)=nsp_T.block(0,2*i+1,4,1);
        mmxd.block(4*i,i+1,4,1)=-1* nsp_T.block(0,2*i+2,4,1);
    }
    mmxd.block(4*3,3,4,1)=nsp_T.block(0,2*3+1,4,1);
    mmxd.block(4*3,0,4,1)=-1*nsp_T.block(0,0,4,1);
    std::cout<<"the mmxd is:\n"<<mmxd<<std::endl;

    MatrixXd A, b;
    MatrixXd left_mat = mmxd.block(0, 1, 16, 3);
    std::cout<<"the sp_mat is:"<<mmxd<<std::endl;

    MatrixXd left_mat_t = left_mat.transpose();

    MatrixXd right_mat = mmxd.col(0);
    right_mat *= -1;

    ConjugateGradient<MatrixXd, Lower | Upper> cg;
    A = left_mat_t * left_mat;
    b = left_mat_t * right_mat;
    cg.compute(A);
    Eigen::VectorXd x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "x is: " << x << std::endl;
}

void test4() {
    test4_0();

    std::vector<MatrixXd> vec_ori, vec_rmd;
    SparseMatrix<double> sp_mat;
    createOriData(vec_ori);
    remodelingData(vec_ori, vec_rmd);
    createSparseMat(vec_rmd, sp_mat);

    SparseMatrix<double> A, b;
    SparseMatrix<double> left_mat = sp_mat.block(0, 1, 16, 3);
    std::cout<<"the sp_mat is:"<<sp_mat<<std::endl;

    SparseMatrix<double> left_mat_t = left_mat.transpose();

    SparseMatrix<double> right_mat = sp_mat.col(0);
    right_mat *= -1;

    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
    A = left_mat_t * left_mat;
    b = left_mat_t * right_mat;
    cg.compute(A);
    Eigen::VectorXd x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "x is: " << x << std::endl;


    std::cout << "the sp mat is:\n" << sp_mat << std::endl;
}
//广义逆矩阵  （利用SVD矩阵分解求）
Eigen::MatrixXd pinv(Eigen::MatrixXd  A)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);//M=USV*
    double  pinvtoler = 1.e-8; //tolerance
    int row = A.rows();
    int col = A.cols();
    int k = std::min(row,col);
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(col,row);
    Eigen::MatrixXd singularValues_inv = svd.singularValues();//奇异值
    Eigen::MatrixXd singularValues_inv_mat = Eigen::MatrixXd::Zero(col, row);
    for (long i = 0; i<k; ++i) {
        if (singularValues_inv(i) > pinvtoler)
            singularValues_inv(i) = 1.0 / singularValues_inv(i);
        else singularValues_inv(i) = 0;
    }
    for (long i = 0; i < k; ++i)
    {
        singularValues_inv_mat(i, i) = singularValues_inv(i);
    }
    X=(svd.matrixV())*(singularValues_inv_mat)*(svd.matrixU().transpose());//X=VS+U*

    return X;

}

void test5(){


    //构建一个3X2的矩阵
    Eigen::MatrixXd mat(3,2);
    mat<<23,54,
    324,56,
    456,1;
    Eigen::MatrixXd matId=Eigen::MatrixXd::Identity(3,3);
    Eigen::MatrixXd matRes = pinv(mat);
    std::cout<<"the mat is:\n"<<mat<<std::endl;
    std::cout<<"the matInv is:\n"<<matRes<<std::endl;
}


MatrixXd getPerspectiveMatrix(){
    double w =800;
    double h = 600;
    double n = 0.5;
    double f = 1000;

    MatrixXd m(4,4);
    m.setZero();
    m(0,0) = 2*n/w;
    m(1,1) = 2*n/h;
    m(2,2) = f/(f-n);
    m(2,3) = n*f/(n-f);
    m(3,2) = 1;
    return m;
}

void test6_0(){
    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(150,250),(771,0),(0,1023),(650,1023)
    //(0.1,0.2),(0.2,0.2),(0.1,0.1),(0.2,0.1)


    double w =800;
    double h = 600;
    double n = 0.5;
    double f = 1000;

    double a = 2*n/w;
    double b =  2*n/h;

    MatrixXd mat(14,13);
    mat.setZero();
    mat<<a,0,-0.1, 0,0,0, 0,0,0, 0,0,0, 0,
    0,b,-0.2, 0,0,0, 0,0,0, 0,0,0, 0,
    0,0,0, a,0,-0.1, 0,0,0, 0,0,0, 0,
    0,0,0, 0,b,-0.2, 0,0,0, 0,0,0, 0,
    0,0,0, 0,0,0, a,0,-0.1, 0,0,0, 0,
    0,0,0, 0,0,0, 0,b,-0.2, 0,0,0, 0,
    0,0,0, 0,0,0, 0,0,0, a,0,-0.1, 0,
    0,0,0, 0,0,0, 0,0,0, 0,b,-0.2, 0,
    -1,0,0, 1,0,0, 0,0,0, 0,0,0,   771-150,
    0,0,0, -1,0,0, 1,0,0, 0,0,0,   -771,
    0,0,0, 0,0,0, -1,0,0, 1,0,0,   650,
    0,-1,0, 0,1,0, 0,0,0, 0,0,0,   -250,
    0,0,0, 0,-1,0, 0,1,0, 0,0,0,   1023,
    0,0,0, 0,0,0, 0,-1,0, 0,1,0,   0;

    MatrixXd A = mat.block(0,0,14,12);
    MatrixXd mb = mat.block(0,12,14,1);

    MatrixXd a_T = A.transpose();

    ConjugateGradient<MatrixXd, Lower | Upper> cg;
    A = a_T * A;
    mb = a_T * mb;
    cg.compute(A);
    Eigen::VectorXd x = cg.solve(mb);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "x is: \n" << x << std::endl;
}

void test6_1(){
    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1),(0.34,0.1),(0.1,0.3),(0.34,0.3)
    double x1 =1254;
    double x2 =1362.5;
    double x3 =1199;
    double x4 =1421;

    double y1 =330.6;
    double y2 =328.5;
    double y3 =403;
    double y4 =401;


    double x1_ =0.1;
    double x2_ =0.34;
    double x3_ =0.1;
    double x4_ =0.34;

    double y1_ =0.1;
    double y2_ =0.1;
    double y3_ =0.3;
    double y4_ =0.3;

    double w = 60;
    double h = 60;
    double n = 10;
    double f = 1000;

    double a = 2*n/w;
    double b =  2*n/h;

    MatrixXd mat(14,13);
    mat.setZero();
    mat<<a,0,-x1_, 0,0,0, 0,0,0, 0,0,0, 0,
            0,b,-y1_, 0,0,0, 0,0,0, 0,0,0, 0,
            0,0,0, a,0,-x2_, 0,0,0, 0,0,0, 0,
            0,0,0, 0,b,-y2_, 0,0,0, 0,0,0, 0,
            0,0,0, 0,0,0, a,0,-x3_, 0,0,0, 0,
            0,0,0, 0,0,0, 0,b,-y3_, 0,0,0, 0,
            0,0,0, 0,0,0, 0,0,0, a,0,-x4_, 0,
            0,0,0, 0,0,0, 0,0,0, 0,b,-y4_, 0,
            -1,0,0, 1,0,0, 0,0,0, 0,0,0,   x2-x1,
            0,0,0, -1,0,0, 1,0,0, 0,0,0,   x3-x2,
            0,0,0, 0,0,0, -1,0,0, 1,0,0,   x4-x3,
            0,-1,0, 0,1,0, 0,0,0, 0,0,0,   y2-y1,
            0,0,0, 0,-1,0, 0,1,0, 0,0,0,   y3-y2,
            0,0,0, 0,0,0, 0,-1,0, 0,1,0,   y4-y3;

    MatrixXd A = mat.block(0,0,14,12);
    MatrixXd mb = mat.block(0,12,14,1);

    MatrixXd a_T = A.transpose();

    ConjugateGradient<MatrixXd, Lower | Upper> cg;
    A = a_T * A;
    mb = a_T * mb;
    cg.compute(A);
    Eigen::VectorXd x = cg.solve(mb);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "x is: \n" << x << std::endl;
}

void test6_2(){
    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1),(0.34,0.1),(0.1,0.3),(0.34,0.3)
    double x1 =1254;
    double x2 =1362.5;
    double x3 =1199;
    double x4 =1421;

    double y1 =330.6;
    double y2 =328.5;
    double y3 =403;
    double y4 =401;


    double x1_ =0.1;
    double x2_ =0.34;
    double x3_ =0.1;
    double x4_ =0.34;

    double y1_ =0.1;
    double y2_ =0.1;
    double y3_ =0.3;
    double y4_ =0.3;

    double w =2560;
    double h = 680;
    double n = 10;
    double f = 1000;

    double a = 2*n/w;
    double b =  2*n/h;

    a = 1/a;
    b = 1/b;

    MatrixXd mat(8,10);
    mat.setZero();
    mat<<x1,y1,0, 0,0,0, -x1_*a*x1,-x1_*a*y1,-x1_*a, 0,
            x2,y2,0, 0,0,0, -x2_*a*x2,-x2_*a*y2,-x2_*a, 0,
            x3,y3,0, 0,0,0, -x3_*a*x3,-x3_*a*y3,-x3_*a, 0,
            x4,y4,0, 0,0,0, -x4_*a*x4,-x4_*a*y4,-x4_*a, 0,

            0,0,0, x1,y1,0, -y1_*b*x1,-y1_*b*y1,-y1_*b, 0,
            0,0,0, x2,y2,0, -y2_*b*x2,-y2_*b*y2,-y2_*b, 0,
            0,0,0, x3,y3,0, -y3_*b*x3,-y3_*b*y3,-y3_*b, 0,
            0,0,0, x4,y4,0, -y4_*b*x4,-y4_*b*y4,-y4_*b, 0;


    MatrixXd A = mat.block(0,0,8,8);
    MatrixXd mb = mat.block(0,8,8,1);
    mb = -1*mb;
    std::cout<<"the mb is:"<<mb<<std::endl;

    MatrixXd a_T = A.transpose();

    ConjugateGradient<MatrixXd, Lower | Upper> cg;
    A = a_T * A;
    mb = a_T * mb;
    cg.compute(A);
    Eigen::VectorXd x = cg.solve(mb);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "x is: \n" << x << std::endl;
}

void test6_3(){
    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1),(0.34,0.1),(0.1,0.3),(0.34,0.3)
    double x1 =1254;
    double x2 =1362.5;
    double x3 =1199;
    double x4 =1421;

    double y1 =330.6;
    double y2 =328.5;
    double y3 =403;
    double y4 =401;


    double x1_ =0.1;
    double x2_ =0.34;
    double x3_ =0.1;
    double x4_ =0.34;

    double y1_ =0.1;
    double y2_ =0.1;
    double y3_ =0.3;
    double y4_ =0.3;

    double w = 60;
    double h = 60;
    double n = 10;
    double f = 1000;

    double a = 2*n/w;
    double b =  2*n/h;

    MatrixXd mat(14,13);
    mat.setZero();
    mat<<a,0,-x1_, 0,0,0, 0,0,0, 0,0,0, 0,
            0,b,-y1_, 0,0,0, 0,0,0, 0,0,0, 0,
            0,0,0, a,0,-x2_, 0,0,0, 0,0,0, 0,
            0,0,0, 0,b,-y2_, 0,0,0, 0,0,0, 0,
            0,0,0, 0,0,0, a,0,-x3_, 0,0,0, 0,
            0,0,0, 0,0,0, 0,b,-y3_, 0,0,0, 0,
            0,0,0, 0,0,0, 0,0,0, a,0,-x4_, 0,
            0,0,0, 0,0,0, 0,0,0, 0,b,-y4_, 0,
            -1,0,0, 1,0,0, 0,0,0, 0,0,0,   x2-x1,
            0,0,0, -1,0,0, 1,0,0, 0,0,0,   x3-x2,
            0,0,0, 0,0,0, -1,0,0, 1,0,0,   x4-x3,
            0,-1,0, 0,1,0, 0,0,0, 0,0,0,   y2-y1,
            0,0,0, 0,-1,0, 0,1,0, 0,0,0,   y3-y2,
            0,0,0, 0,0,0, 0,-1,0, 0,1,0,   y4-y3;

    MatrixXd A = mat.block(0,0,14,12);
    MatrixXd mb = mat.block(0,12,14,1);

    MatrixXd a_T = A.transpose();

    ConjugateGradient<MatrixXd, Lower | Upper> cg;
    A = a_T * A;
    mb = a_T * mb;
    cg.compute(A);
    Eigen::VectorXd x = cg.solve(mb);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "x is: \n" << x << std::endl;
}

void test6_4(){
    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1,1),(0.34,0.1,1),(0.1,0.3,1),(0.34,0.3,1)
    Vector3d v1(0.1,0.1,1);
    Vector3d v2(0.34,0.1,1);
    Vector3d v3(0.1,0.3,1);


    double x1 =1254;
    double x2 =1362.5;
    double x3 =1199;
    double x4 =1421;

    double y1 =330.6;
    double y2 =328.5;
    double y3 =403;
    double y4 =401;


    double x1_ =0.1;
    double x2_ =0.34;
    double x3_ =0.1;
    double x4_ =0.34;

    double y1_ =0.1;
    double y2_ =0.1;
    double y3_ =0.3;
    double y4_ =0.3;

    double z1_ =1;
    double z2_ =1;
    double z3_ =1;
    double z4_ =1;

    double w = 60;
    double h = 60;
    double n = 10;
    double f = 1000;

    double a = 2*n/w;
    double b =  2*n/h;

    double cos_12 =v1.dot(v2)/v1.norm()/v2.norm();
    double cos_23 =v3.dot(v2)/v3.norm()/v2.norm();
    double cos_31 =v1.dot(v3)/v1.norm()/v3.norm();


    double k11 = x1_*x1_+y1_*y1_+z1_*z1_;
    double k12 = x2_*x2_+y2_*y2_+z2_*z2_;
    double k13 = -2*(x1_*x2_+y1_*y2_+z1_*z2_);
    double k14 = (1362.5-1254)*(1362.5-1254) + (330.6-328.5)*(330.6-328.5);

    double k21 =  x2_*x2_+y2_*y2_+z2_*z2_;
    double k22 = x3_*x3_+y3_*y3_+z3_*z3_;
    double k23 = -2*(x3_*x2_+y3_*y2_+z3_*z2_);
    double k24 = (1362.5-1199)*(1362.5-1199) + (403-328.5)*(403-328.5);

    double k31 = x3_*x3_+y3_*y3_+z3_*z3_;
    double k32 = x1_*x1_+y1_*y1_+z1_*z1_;
    double k33 = -2*(x1_*x3_+y1_*y3_+z1_*z3_);
    double k34 = (1199-1254)*(1199-1254) + (330.6-403)*(330.6-403);
    ////////////////////////////////////////////////////////////////////////////////////
    double a1 = k24*k11;
    double b1 = k24*k13/(2*a1);
    double c1 = k24*k12-a1*b1*b1;

    double a2 = k14*k22;
    double b2 = k14*k23/(2*a2);
    double c2 = k14*k21-a2*b2*b2;

    double a3 = k34*k21;
    double b3 = k34*k23/(2*a3);
    double c3 = k34*k22-a3*b3*b3;

    double a4 = k24*k32;
    double b4 = k24*k33/(2*a3);
    double c4 = k24*k31-a4*b4*b4;

    double a5 = k34*k12;
    double b5 = k34*k13/(2*a5);
    double c5 = k34*k11-a5*b5*b5;

    double a6 = k14*k31;
    double b6 = k14*k33/(2*a6);
    double c6 = k14*k32-a6*b6*b6;
    ////////////////////////////////////////////////////////////////////////////////////




}


void convertParams(const double params[],double r[]){
    r[0] = params[0]+params[9];
    r[1] = params[1]+params[4];
    r[2] = params[5]+params[8];
    r[3] = params[2];
    r[4] = params[6];
    r[5] = params[10];
    r[6] = -(params[3]+params[7]+params[11]);
}

void test6_6(double vals[]);
//牛顿迭代法
void test6_5(){
    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1,1),(0.34,0.1,1),(0.1,0.3,1),(0.34,0.3,1)
    Vector3d v1(0.1,0.1,1);
    Vector3d v2(0.34,0.1,1);
    Vector3d v3(0.1,0.3,1);


    double x1 =1254;
    double x2 =1362.5;
    double x3 =1199;
    double x4 =1421;

    double y1 =330.6;
    double y2 =328.5;
    double y3 =403;
    double y4 =401;


    double x1_ =0.1;
    double x2_ =0.34;
    double x3_ =0.1;
    double x4_ =0.34;

    double y1_ =0.1;
    double y2_ =0.1;
    double y3_ =0.3;
    double y4_ =0.3;

    double z1_ =1;
    double z2_ =1;
    double z3_ =1;
    double z4_ =1;

    double w = 60;
    double h = 60;
    double n = 10;
    double f = 1000;

    double a = 2*n/w;
    double b =  2*n/h;

    double cos_12 =v1.dot(v2)/v1.norm()/v2.norm();
    double cos_23 =v3.dot(v2)/v3.norm()/v2.norm();
    double cos_31 =v1.dot(v3)/v1.norm()/v3.norm();


    double k11 = x1_*x1_+y1_*y1_+z1_*z1_;
    double k12 = x2_*x2_+y2_*y2_+z2_*z2_;
    double k13 = -2*(x1_*x2_+y1_*y2_+z1_*z2_);
    double k14 = (1362.5-1254)*(1362.5-1254) + (330.6-328.5)*(330.6-328.5);

    double k21 =  x2_*x2_+y2_*y2_+z2_*z2_;
    double k22 = x3_*x3_+y3_*y3_+z3_*z3_;
    double k23 = -2*(x3_*x2_+y3_*y2_+z3_*z2_);
    double k24 = (1362.5-1199)*(1362.5-1199) + (403-328.5)*(403-328.5);

    double k31 = x3_*x3_+y3_*y3_+z3_*z3_;
    double k32 = x1_*x1_+y1_*y1_+z1_*z1_;
    double k33 = -2*(x1_*x3_+y1_*y3_+z1_*z3_);
    double k34 = (1199-1254)*(1199-1254) + (330.6-403)*(330.6-403);

    ////////////////////////////////////////////////////////////////////////////////////
    double pas[12]={k11,k12,k13,k14,
                  k21,k22,k23,k24,
                  k31,k32,k33,k34};
    double arr[7];
    convertParams(pas,arr);
    test6_6(arr);
    ////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////
    //Y1=k11*m1*m1+k12*m2*m2-k13*m1*m2-k14=0
    //Y2=k21*m2*m2+k22*m3*m3-k23*m2*m3-k24=0
    //Y3=k31*m3*m3+k32*m1*m1-k33*m3*m1-k34=0

#define N 3
#define eps 0.0000001
#define MaxLoop 5000
    const int N2 = N*N;

    using namespace std;
    double m[N],y[N],h1[N],es,esmax,jacobi[N][N],aij;
    int i,j,k,it=0,iter=0;
    //给h0【N】赋值
    m[0]=1;
    m[1]=1;
    m[1]=1;

    do{
        it++;
        for(i=0;i<N;i++){
            //计算雅克比矩阵
            jacobi[0][0]=2*k11*m[0]+k13*m[1];
            jacobi[0][1]=2*k12*m[1]+k13*m[0];
            jacobi[0][2]=0;

            jacobi[1][0]=0;
            jacobi[1][1]=2*k21*m[1]+k23*m[2];
            jacobi[1][2]=2*k22*m[2]+k23*m[1];

            jacobi[2][0]=2*k32*m[0]+k33*m[2];
            jacobi[2][1]=0;
            jacobi[2][2]=2*k31*m[2]+k33*m[0];

            aij=0;
            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    aij=aij+jacobi[i][j]*jacobi[i][j];
                    if(aij>1){
                        cout<<"Sorry, aij>1...\n";
                        k=rand()%3;
                        m[k]=((double)rand())/RAND_MAX;
                        cout<<"由计算机第"<<it<<" 次随机初值：\n";
                        for(i=0;i<N;i++){
                            cout<<m[i]<<"\t \n";
                        }
                    }
                }
            }
        }
    }
    while (aij>1);
    //保证aij<1

    for(i=0;i<N;i++){
        h1[i]=m[i];
    }

    do{
        iter++;

        y[0] = k11*h1[0]*h1[0]+k12*h1[1]*h1[1]+k13*h1[0]*h1[1]-k14;
        y[1] = k21*h1[1]*h1[1]+k22*h1[2]*h1[2]+k23*h1[1]*h1[2]-k24;
        y[2] = k31*h1[2]*h1[2]+k32*h1[0]*h1[0]+k33*h1[2]*h1[0]-k34;

        esmax=0.0;
        for(i=0;i<N;i++){
            es=y[i]-h1[i];
            if(fabs(es)>fabs(esmax))
                esmax=es;
        }
        if(fabs(esmax)<eps){
            cout<<"方程组的解为：\n";
            for(i=0;i<N;i++){
                cout<<h1[i]<<" \t";
            }
            cout<<endl;
            break;
        }

        for(i=0;i<N;i++){
            h1[i]=y[i];
        }
    }
    while (iter<MaxLoop);

    std::cout<<"end...\n";
    ////////////////////////////////////////////////////////////////////////////////////

}


typedef double (*cal_fx)(double params[],double xs[]);
typedef double (*cal_dfx)(double params[],double xs[],int index);
typedef double (*cal_ddfx)(double params[],double xs[],int index);



//Y = a1*x1*x1+a2*x2*x2+a3*x3*x3+a4*x1*x2+a5*x2*x3+a6*x3*x1+a7 = 0
double calculate_fx_back(double params[],double xs[],int index){
    double a1 = params[0]+params[9];
    double a2 = params[1]+params[4];
    double a3 = params[5]+params[8];
    double a4 = params[2];
    double a5 = params[6];
    double a6 = params[10];
    double a7 = -(params[3]+params[7]+params[11]);

    double sum = a1*xs[0]*xs[0]+
            a2*xs[1]*xs[1]+
            a3*xs[2]*xs[2]+
            a4*xs[0]*xs[1]+
            a5*xs[1]*xs[2]+
            a6*xs[2]*xs[0]+
            a7;
    return sum;
}
double calculate_fx(double params[],double xs[]){
    double sum = params[0]*xs[0]*xs[0]+
            params[1]*xs[1]*xs[1]+
            params[2]*xs[2]*xs[2]+
            params[3]*xs[0]*xs[1]+
            params[4]*xs[1]*xs[2]+
            params[5]*xs[2]*xs[0]+
            params[6];
    return sum;
}
enum PFLAG{
    D=0,
    P1=1,
    P2=2,
    P3=4
};

int bflag(int p1,int p2){
    if(p1&p2){
        return 1;
    } else{
        return 0;
    }
}
double calculate_dfx(double params[], double xs_c[],int index){
    double xs[3];
    for(int i=0;i<3;++i){
        xs[i] = xs_c[i];
    }

    int flag = 1;
    flag=flag<<index;

    double d11 = 2*params[0]*xs[0]*bflag(PFLAG::P1,flag);
    double d12 = 2*params[1]*xs[1]*bflag(PFLAG::P2,flag);
    double d13 = 2*params[2]*xs[2]*bflag(PFLAG::P3,flag);
    double tem_d1 = d11+d12+d13;

    double d1 = 2*params[0]*xs[0]*bflag(PFLAG::P1,flag)+
                 2*params[1]*xs[1]*bflag(PFLAG::P2,flag)+
                 2*params[2]*xs[2]*bflag(PFLAG::P3,flag);

    xs[index] = 1;

    double d21 = params[3]*xs[0]*xs[1]*bflag(PFLAG::P1|PFLAG::P2,flag);
    double d22 = params[4]*xs[1]*xs[2]*bflag(PFLAG::P2|PFLAG::P3,flag);
    double d23 = params[5]*xs[2]*xs[0]*bflag(PFLAG::P3|PFLAG::P1,flag);
    double tem_d2 = d21+d22+d23;

    double d2 = params[3]*xs[0]*xs[1]*bflag(PFLAG::P1|PFLAG::P2,flag)+
            params[4]*xs[1]*xs[2]*bflag(PFLAG::P2|PFLAG::P3,flag)+
            params[5]*xs[2]*xs[0]*bflag(PFLAG::P3|PFLAG::P1,flag);
    return d1+d2;
}

double calculate_ddfx(double params[], double xs_c[],int index){
    double xs[3];
    for(int i=0;i<3;i++){
        xs[i]=0;
    }

    xs[index]=1;

    int flag = 1;
    flag=flag<<index;
    double d1 = 2*params[0]*xs[0]*bflag(PFLAG::P1,flag)+
                2*params[1]*xs[1]*bflag(PFLAG::P2,flag)+
                2*params[2]*xs[2]*bflag(PFLAG::P3,flag);
    return d1;
}

//vals[]数组为a1-a7
void test6_6(double vals[]){

    cal_fx func_fx = &calculate_fx;
    cal_dfx func_dfx = &calculate_dfx;
    cal_ddfx func_ddfx = &calculate_ddfx;
    double x[3];
    x[0]=50;
    x[1]=30;
    x[2]=80;
    double err=1.0e-10;
    double delta_fx=1.0e-6;
    double fx;
    double dfx;
    double ddfx;
    int idx =0;
    int loopCnt=0;
    //更新x【0】
    while(loopCnt<4000){
        loopCnt++;
        fx= func_fx(vals,x);
        if(abs(fx)<err){
            break;
        }
        dfx = func_dfx (vals,x,idx);
        if(abs(dfx)<delta_fx){
            if(idx<2){
                idx++;
            } else{
                idx = 0;
            }
            continue;
        }
        ddfx = func_ddfx(vals,x,idx);
        x[idx]=x[idx]-fx/dfx;
    }

    std::cout<<"the err is:"<<fx<<std::endl;
    std::cout<<"loop cnt is:"<<loopCnt<<" .\n"<<
    "the x[] is:\n";
    for(int mmm=0;mmm<3;++mmm){
        std::cout<<" "<<x[mmm];
    }
    std::cout<<std::endl;
}
void test8_lbfgs_t(Eigen::MatrixXd& k_mat);
void test6_7(){
//    double x1 =1254;
//    double x2 =1362.5;
//    double x3 =1199;
//    double x4 =1421;
//
//    double y1 =330.6;
//    double y2 =328.5;
//    double y3 =403;
//    double y4 =401;


//    double x1_ =100;
//    double x2_ =200;
//    double x3_ =150;
//    double x4_ =0.34;
//
//    double y1_ =200;
//    double y2_ =400;
//    double y3_ =700;
//    double y4_ =0.3;

    double x1 =249;
    double x2 =238;
    double x3 =461;
    double x4 =1421;

    double y1 =700-352;
    double y2 =700-425;
    double y3 =700-424;
    double y4 =401;

    double x1_ =-50;
    double x2_ =-50;
    double x3_ =70;
    double x4_ =0.34;

    double y1_ =50;
    double y2_ =-50;
    double y3_ =-50;
    double y4_ =0.3;

    double z1_ =1;
    double z2_ =1;
    double z3_ =1;
    double z4_ =1;

    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1,1),(0.34,0.1,1),(0.1,0.3,1),(0.34,0.3,1)
    Vector3d v1(100,200,1);
    Vector3d v2(200,400,1);
    Vector3d v3(150,700,1);

    Vector3d v1_(x1_,y1_,z1_);
    Vector3d v2_(x2_,y2_,z2_);
    Vector3d v3_(x3_,y3_,z3_);

    Vector2d p1(x1,y1);
    Vector2d p2(x2,y2);
    Vector2d p3(x3,y3);


    double w = 60;
    double h = 60;
    double n = 10;
    double f = 1000;

    double cos_12 =v1.dot(v2)/v1.norm()/v2.norm();
    double cos_23 =v3.dot(v2)/v3.norm()/v2.norm();
    double cos_31 =v1.dot(v3)/v1.norm()/v3.norm();


    double k11 = x1_*x1_+y1_*y1_+z1_*z1_;
    double k12 = x2_*x2_+y2_*y2_+z2_*z2_;
    double k13 = -2*(x1_*x2_+y1_*y2_+z1_*z2_);
    double k14 = (1362.5-1254)*(1362.5-1254) + (330.6-328.5)*(330.6-328.5);
    k11 = v1_.squaredNorm();
    k12 = v2_.squaredNorm();
    k13 = -2*v1_.dot(v2_);
    k14 = (p1-p2).squaredNorm();//模长的平方

    double k21 =  x2_*x2_+y2_*y2_+z2_*z2_;
    double k22 = x3_*x3_+y3_*y3_+z3_*z3_;
    double k23 = -2*(x3_*x2_+y3_*y2_+z3_*z2_);
    double k24 = (1362.5-1199)*(1362.5-1199) + (403-328.5)*(403-328.5);
    k21 = v2_.squaredNorm();
    k22 = v3_.squaredNorm();
    k23 = -2*v2_.dot(v3_);
    k24 = (p3-p2).squaredNorm();//模长的平方

    double k31 = x3_*x3_+y3_*y3_+z3_*z3_;
    double k32 = x1_*x1_+y1_*y1_+z1_*z1_;
    double k33 = -2*(x1_*x3_+y1_*y3_+z1_*z3_);
    double k34 = (1199-1254)*(1199-1254) + (330.6-403)*(330.6-403);
    k31 = v1_.squaredNorm();
    k32 = v3_.squaredNorm();
    k33 = -2*v1_.dot(v3_);
    k34 = (p1-p3).squaredNorm();//模长的平方

    ////////////////////////////////////////////////////////////////////////////////////
    MatrixXd k_mat(3,4);
    k_mat<<k11,k12,k13,k14,
            k21,k22,k23,k24,
            k31,k32,k33,k34;
    test8_lbfgs_t(k_mat);
    ///////////////////////////////////////////////////////////////////////////////////
    //验证
    Eigen::Vector3d v_x={0.983655 ,0.249174 ,0.398666};
    Eigen::Vector3d vx1 = v1*v_x[0];
    Eigen::Vector3d vx2 = v2*v_x[1];
    Eigen::Vector3d vx3 = v3*v_x[2];

    Eigen::Vector3d vx12 = vx1 - vx2;
    Eigen::Vector3d vx23 = vx3 - vx2;
    Eigen::Vector3d vx13 = vx1 - vx3;
    std::cout<<"vx12 length is:"<<vx12.norm()<<std::endl;
    std::cout<<"vx23 length is:"<<vx23.norm()<<std::endl;
    std::cout<<"vx13 length is:"<<vx13.norm()<<std::endl;

    std::cout<<"k14 length is:"<<sqrt(k14)<<std::endl;
    std::cout<<"k24 length is:"<<sqrt(k24)<<std::endl;
    std::cout<<"k34 length is:"<<sqrt(k34)<<std::endl;


}

void test6_8(){
//    double x1 =1254;
//    double x2 =1362.5;
//    double x3 =1199;
//    double x4 =1421;
//
//    double y1 =330.6;
//    double y2 =328.5;
//    double y3 =403;
//    double y4 =401;


//    double x1_ =100;
//    double x2_ =200;
//    double x3_ =150;
//    double x4_ =0.34;
//
//    double y1_ =200;
//    double y2_ =400;
//    double y3_ =700;
//    double y4_ =0.3;

    double x1 =249;
    double x2 =238;
    double x3 =461;
    double x4 =1421;

    double y1 =700-352;
    double y2 =700-425;
    double y3 =700-424;
    double y4 =401;

    double x1_ =-50;
    double x2_ =-50;
    double x3_ =70;
    double x4_ =0.34;

    double y1_ =50;
    double y2_ =-50;
    double y3_ =-50;
    double y4_ =0.3;

    double z1_ =1;
    double z2_ =1;
    double z3_ =1;
    double z4_ =1;

    //(x1,y1)    (x2,y2) (x3,y3) (x4,y4)
    //(1254,330.6),(1362.5,328.5),(1199,403),(1421,401)
    //(0.1,0.1,1),(0.34,0.1,1),(0.1,0.3,1),(0.34,0.3,1)
    Vector3d v1(100,200,1);
    Vector3d v2(200,400,1);
    Vector3d v3(150,700,1);

    Vector3d v1_(x1_,y1_,z1_);
    Vector3d v2_(x2_,y2_,z2_);
    Vector3d v3_(x3_,y3_,z3_);

    Vector2d p1(x1,y1);
    Vector2d p2(x2,y2);
    Vector2d p3(x3,y3);



    double k11 =0;
    double k12 =0;
    double k13 = 0;
    double k14 = 0;
    k11 = v1_.squaredNorm();
    k12 = v2_.squaredNorm();
    k13 = -2*v1_.dot(v2_);
    k14 = (p1-p2).squaredNorm();//模长的平方

    double k21 = 0;
    double k22 =0;
    double k23 = 0;
    double k24 = 0;
    k21 = v2_.squaredNorm();
    k22 = v3_.squaredNorm();
    k23 = -2*v2_.dot(v3_);
    k24 = (p3-p2).squaredNorm();//模长的平方

    double k31 =0;
    double k32 = 0;
    double k33 = 0;
    double k34 = 0;
    k31 = v1_.squaredNorm();
    k32 = v3_.squaredNorm();
    k33 = -2*v1_.dot(v3_);
    k34 = (p1-p3).squaredNorm();//模长的平方

    ////////////////////////////////////////////////////////////////////////////////////
    MatrixXd k_mat(3,4);
    k_mat<<k11,k12,k13,k14,
            k21,k22,k23,k24,
            k31,k32,k33,k34;
    test8_lbfgs_t(k_mat);
    ///////////////////////////////////////////////////////////////////////////////////
    //验证
    Eigen::Vector3d v_x={0.641743, 0.823504 , 2.39353};
    Eigen::Vector3d vx1 = v1*v_x[0];
    Eigen::Vector3d vx2 = v2*v_x[1];
    Eigen::Vector3d vx3 = v3*v_x[2];

    Eigen::Vector3d vx12 = vx1 - vx2;
    Eigen::Vector3d vx23 = vx3 - vx2;
    Eigen::Vector3d vx13 = vx1 - vx3;
    std::cout<<"vx12 length is:"<<vx12.norm()<<std::endl;
    std::cout<<"vx23 length is:"<<vx23.norm()<<std::endl;
    std::cout<<"vx13 length is:"<<vx13.norm()<<std::endl;

    std::cout<<"k14 length is:"<<sqrt(k14)<<std::endl;
    std::cout<<"k24 length is:"<<sqrt(k24)<<std::endl;
    std::cout<<"k34 length is:"<<sqrt(k34)<<std::endl;
}

void test6(){}

typedef Eigen::Vector3d Geo3d;
typedef Eigen::Vector4d Geo4d;
typedef Eigen::Matrix4d GeoMat4d;


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

    std::cout<<"the length is:\n"<<l_AB<<std::endl<<std::endl
             <<l_AC<<std::endl<<std::endl
             <<l_AD<<std::endl<<std::endl
             <<l_BC<<std::endl<<std::endl
             <<l_BD<<std::endl<<std::endl
             <<l_CD<<std::endl<<std::endl;

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


void test7(){
    //Y = x*x -3 = 0
    // dfx = 2*x
    double x =-1;
    double y = x*x-3;
    double esp=1.0e-10;
    while (abs(y)>esp){
        x = x-y/(2*x);
        std::cout<<"x is:"<<x<<std::endl;
        y = x*x-3;
    }

    std::cout<<"finished,the x is:"<<x<<std::endl;
}

void test8_lbfgs(){


}

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace LBFGSpp;


//f1 = k11*x1*x1 + k12*x2*x2 + k13*x1*x2 -k14 =0
//二次型
//             |k11       k13/2     0|
// [x1,x2,x3]* |k13/2     k12       0| *[x1,x2,x3]T
//             |0         0         0|

//f2 = k21*x2*x2 + k22*x3*x3 + k23*x2*x3 -k24 =0
//二次型
//             |0     0           0|
// [x1,x2,x3]* |0     k21     k23/2| *[x1,x2,x3]T
//             |0     k23/2     k22|


//f3 = k31*x3*x3 + k32*x1*x1 + k33*x3*x1 -k34 =0
//二次型
//             |k32       0     k33/2|
// [x1,x2,x3]* |0         0         0| *[x1,x2,x3]T
//             |k33/2     0       k31|


//J(x1,x2,x3)=f1*f1 + f2*f2 + f3*f3 =0
//判断J是否有解，可以通过判断J的凸凹性来得到
//这里为了快速得到结果，先假设J是有解的，直接调用L-BFGS来求
//df1(x1) = 2*k11*x1 + k13 *1*x2;
//df1(x2) = 2*k12*x2 + k13 *x1*1;
//df1(x3) = 0

//df2(x1) = 0
//df2(x2) = 2*k21*x2 + k23*1*x3;
//df2(x3) = 2*k22*x3 + k23*1*x2;

//df3(x1) = 2*k32*x1 + k33*x3*1;
//df3(x2) = 0;
//df3(x3) = 2*k31*x3 + k33*1*x1;


//dJ(x1) = 2*f1*df1 + 2*f2*df2 + 2*f3*df3;

class Rosenbrock
{
private:
    int n;
    Eigen::MatrixXd params;//3X4的参数矩阵

    Eigen::MatrixXd f1_quad_mat;//f1 二次型 矩阵
    Eigen::MatrixXd f2_quad_mat;
    Eigen::MatrixXd f3_quad_mat;

public:
    Rosenbrock(int n_,Eigen::MatrixXd& params_) : n(n_),params(params_) {
        f1_quad_mat = Eigen::MatrixXd(3,3);
        f2_quad_mat = Eigen::MatrixXd(3,3);
        f3_quad_mat = Eigen::MatrixXd(3,3);
        f1_quad_mat.setZero();
        f2_quad_mat.setZero();
        f3_quad_mat.setZero();

        f1_quad_mat(0,0) = params(0,0);
        f1_quad_mat(0,1) = params(0,2)/2;
        f1_quad_mat(1,0) = params(0,2)/2;
        f1_quad_mat(1,1) = params(0,1);

        f2_quad_mat(1,1) = params(1,0);
        f2_quad_mat(1,2) = params(1,2)/2;
        f2_quad_mat(2,1) = params(1,2)/2;
        f2_quad_mat(2,2) = params(1,1);

        f3_quad_mat(0,0) = params(2,1);
        f3_quad_mat(0,2) = params(2,2)/2;
        f3_quad_mat(2,0) = params(2,2)/2;
        f3_quad_mat(2,2) = params(2,0);
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = cal_J(x);

        grad[0] = cal_J_x1(x);
        grad[1] = cal_J_x2(x);
        grad[2] = cal_J_x3(x);
        std::cout<<"cur err is:"<<fx<<std::endl;
        return fx;
    }

private:
    double cal_f1(const VectorXd& x){
        return x.transpose()*f1_quad_mat*x-params(0,3);
    }
    double cal_f2(const VectorXd& x){
        return x.transpose()*f2_quad_mat*x-params(1,3);
    }
    double cal_f3(const VectorXd& x){
        return x.transpose()*f3_quad_mat*x-params(2,3);
    }

    double cal_df1_x1(const VectorXd& x){
        return 2*params(0,0)*x[0]+params(0,2)*x[1];
    }
    double cal_df1_x2(const VectorXd& x){
        return 2*params(0,1)*x[1]+params(0,2)*x[0];
    }
    double cal_df1_x3(const VectorXd& x){
        return 0;
    }
    double cal_df2_x1(const VectorXd& x){
        return 0;
    }
    double cal_df2_x2(const VectorXd& x){
        return 2*params(1,0)*x[1]+params(1,2)*x[2];
    }
    double cal_df2_x3(const VectorXd& x){
        return 2*params(1,1)*x[2]+params(1,2)*x[1];
    }

    double cal_df3_x1(const VectorXd& x){
        return 2*params(2,1)*x[0]+params(2,2)*x[2];
    }
    double cal_df3_x2(const VectorXd& x){
        return 0;
    }
    double cal_df3_x3(const VectorXd& x){
        return 2*params(2,0)*x[2]+params(2,2)*x[0];
    }
    double cal_J_x1(const VectorXd& x){
        return 2*cal_f1(x)*cal_df1_x1(x) + 2*cal_f2(x)*cal_df2_x1(x)+ 2*cal_f3(x)*cal_df3_x1(x);
    }
    double cal_J_x2(const VectorXd& x){
        return 2*cal_f1(x)*cal_df1_x2(x) + 2*cal_f2(x)*cal_df2_x2(x)+ 2*cal_f3(x)*cal_df3_x2(x);
    }
    double cal_J_x3(const VectorXd& x){
        return 2*cal_f1(x)*cal_df1_x3(x) + 2*cal_f2(x)*cal_df2_x3(x)+ 2*cal_f3(x)*cal_df3_x3(x);
    }

    double cal_J(const VectorXd& x){
        double f1 = cal_f1(x);
        double f2 = cal_f2(x);
        double f3 = cal_f3(x);
        return f1*f1+f2*f2+f3*f3;
    }
};










/**
 * 近截面z=-1;
 * 近截面显示投影映射的点（P1,P2,P3...）
 * 约束：已知近截面的P1等坐标，和目标点T1等 点间距
 * 目标：求目标点T1的坐标
 *
 * S1*P1 = T1
 *
 * （L12为T1、T2的点间距 ...）
 * f12 = (S1*P1_x -S2*P2_x)² +(S1*P1_y -S2*P2_y)² +(S1*P1_z -S2*P2_z)² - L12² = 0;
 * f13 = (S1*P1_x -S3*P3_x)² +(S1*P1_y -S3*P3_y)² +(S1*P3_z -S2*P3_z)² - L13² = 0;
 *
 * J(S1,S2,S3...) = f12² + f13²+ f14²+...+
 *                         f23²+ f24²+...+
 *                               f34²+...+
 *                               ... = 0
 *
 *
 * */
class mrBrock{

};

void test8_lbfgs_t(Eigen::MatrixXd& k_mat){
    const int n = 3;
    LBFGSParam<double> param;
    LBFGSSolver<double> solver(param);
    Rosenbrock fun(n,k_mat);

    VectorXd x = VectorXd::Zero(n);
    x[0]=10;
    x[1]=10;
    x[2]=10;
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

}

void pHello1(const CB_Flag f) {
    if (f == CB_Flag::RF_SCALE)
        std::cout << "Hello1 scale...\n";
}

void pHello2(const CB_Flag f) {
    if (f == CB_Flag::RF_FLAT)
        std::cout << "Hello2 flat...\n";
}

void testCallback() {
    CallBackManager cbm;
    CB_Flag f(CB_Flag::RF_FLAT);
    cbm.push_back(pHello1);
    cbm.push_back(pHello2);
    cbm.run(f);
    std::cout << "end of cbm run...\n";
    cbm.remove(pHello1);
    cbm.run(f);
    std::cout << "end of cbm run...\n";
}



struct A
{
    void Fun(int i){std::cout<<i<<std::endl;}
    void Fun1(int i, double j){std::cout<<i+j<<std::endl;}
};

void Fun2(int i, double j){std::cout<<"this is:"<<i+i+j<<std::endl;}

void testDelegate(){
//    A a;
//    auto d = CreateDelegate(&a, &A::Fun); //创建委托
//    d(1); //调用委托，将输出1
//    auto d1 = CreateDelegate(&a, &A::Fun1); //创建委托
//    d1(1, 2.5); //调用委托，将输出3.5

    TDelegate<void,int,double> d3 = CreateDelegate(Fun2); //创建委托
    d3(1,2); //调用委托，将输出3.5
}






int main() {
//    test0();
//    std::cout<<"\nbegin 0_0:"<<std::endl;
//    test0_0();
//    getData();


//    test4();
//    test5();
    test_9();
//    test10();
//    test6_8();
//    testDelegate();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}