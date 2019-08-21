//
// Created by wasku on 19-7-27.
//

#ifndef EIGENSAMPLE_TCALLBACK_H
#define EIGENSAMPLE_TCALLBACK_H

#include <vector>

enum CB_Flag{
    DEFAULT,
    RF_FLAT,
    RF_SCALE,
    RF_ROTATION
};

typedef void (*glCallBack)(const CB_Flag);

class CallBackManager{
private:
    typedef std::vector<glCallBack >::iterator t_iter;

public:
    void push_back(glCallBack pFunc){
        m_vec.push_back(pFunc);
    }
    void remove(glCallBack pFunc){
        t_iter i;
        for( i = m_vec.begin();i!=m_vec.end();++i){
            if((*i) == pFunc){
                m_vec.erase(i);
                break;
            }
        }
    }

    void run(const CB_Flag flag){

        for(glCallBack item :m_vec){
            (*item)(flag);
        }
    }

    void clear(){
        m_vec.clear();
    }

private:
    std::vector<glCallBack> m_vec;
};



#endif //EIGENSAMPLE_TCALLBACK_H
