#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <ctime>
#include <cstdlib>
using namespace std;
#define RING_BUFFER_LEN 10//ring buffer的长度，我们假设突触延迟都相同，都是这个（ms）.多个文件需要同时修改！ (多了可以放入parameter.h)
class Poisson_generator
{

public:
    int get_rand=0;
    int les_i_time=4;
    double rate;
    double weight;
        // 左闭右闭区间
    int getRand(int min, int max) {
        return ( rand() % (max - min + 1) ) + min ;
    }
    
    Poisson_generator(double rate_,int weight_)
    {
        // srand(time(0));
        // cout<<"time(0) is "<<time(0)/100<<endl;
        rate=rate_;
        weight=weight_;
        if(rate>10000)
        {
            cout<<"rate_toOOOOOOO_large!more than 10000"<<endl;
        }
    }

    void Get_poisson_num(double re_spikes[RING_BUFFER_LEN])
    {
    double res[RING_BUFFER_LEN];
    double add_table[RING_BUFFER_LEN]={0};
    for (size_t i = 0; i < RING_BUFFER_LEN; i++)
    {
        srand(1042710+get_rand*100);
        // srand(time(0)+get_rand*100);
        //cout<<"time(0) is "<<time(0)+get_rand*100<<endl;
        get_rand=getRand(0,99);
        if(get_rand<rate/100){
            res[i]=weight;
        }
        else
        {
            res[i]=0;
        }
    }
    for (size_t i = 0; i < RING_BUFFER_LEN; i++)
    {
        for (size_t j = 0; j < les_i_time; j++)
        {
            if (i+j+1<RING_BUFFER_LEN)
            {
                add_table[i+j+1]=res[i]*(1/(2+j));
            }
        }    
    }
    for (size_t i = 0; i < RING_BUFFER_LEN; i++)
    {
        res[i]+=add_table[i];
        re_spikes[i]=res[i];
    }
    }


    // void Get_poisson_num(double re_spikes[RING_BUFFER_LEN])
    // {
    // double step_I_in=weight*(rate/10000.);
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::poisson_distribution<> d(step_I_in);
    // double res[RING_BUFFER_LEN];
    // double add_table[RING_BUFFER_LEN]={0};
    // for (size_t i = 0; i < RING_BUFFER_LEN; i++)
    // {
    //     res[i]=d(gen);
    // }
    // for (size_t i = 0; i < RING_BUFFER_LEN; i++)
    // {
    //     for (size_t j = 0; j < les_i_time; j++)
    //     {
    //         if (i+j+1<RING_BUFFER_LEN)
    //         {
    //             add_table[i+j+1]=res[i]*(1/(2+j));
    //         }
    //     }    
    // }
    // for (size_t i = 0; i < RING_BUFFER_LEN; i++)
    // {
    //     res[i]+=add_table[i];
    //     re_spikes[i]=res[i];
    // }
    

    // }
};





// int main()


// {
//     // Poisson_generator ps(5);
//     // double re_spikes[RING_BUFFER_LEN];
//     // ps.Get_poisson_num(re_spikes);
//     // for (size_t i = 0; i < RING_BUFFER_LEN; i++)
//     // {
//     //     cout<<re_spikes[i]<<" ";
//     // }
//     // cout<<endl;



//     // std::random_device rd;
//     // std::mt19937 gen(rd());
//     // // for (int i = 0; i < 100; ++i)
//     // // {
//     // //     std::cout<<"gen is "<<gen()<<std::endl;
//     // // }

//     // // 若平均一分钟出现 4 次事件
//     // // 则在一分钟内出现 n 次的频率如何？
//     // std::poisson_distribution<> d(4);
 
//     // std::map<int, int> hist;
//     // for(int n=0; n<10000; ++n) {
//     //     ++hist[d(gen)];
//     //     cout<<d(gen)<<endl;
//     // }
//     // for(auto p : hist) {
//     //     std::cout << p.first <<
//     //             ' ' << std::string(p.second/100, '*') << '\n';
//     // }
// }