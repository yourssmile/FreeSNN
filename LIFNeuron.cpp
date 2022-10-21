#include <vector>
// #include <Event.cpp>
#include<deque>
#include <string>
#include"Histentry.h"
#include <cmath>
#include <iostream>
#include "Poisson_generator.cpp"
#include <sys/time.h>
using namespace std;
#define RING_BUFFER_LEN 10//ring buffer的长度，我们假设突触延迟都相同，都是这个(ms).多个文件需要同时修改！
#define less_i 4 
#define rate_record_step 10//记录步数
#define print_spike_and_recv 0//是否开启脉冲接受跟发射时打印信息


// struct timeval
// {
// #ifdef __USE_TIME_BITS64
//   __time64_t tv_sec;		/* Seconds.  */
//   __suseconds64_t tv_usec;	/* Microseconds.  */
// #else
//   __time_t tv_sec;		/* Seconds.  */
//   __suseconds_t tv_usec;	/* Microseconds.  */
// #endif
// };


//#define Record 1//是否记录Vm跟spike，是1，否0
class LIFNeuron
{
private:

    
public:
    //STDP突触需要的参数!!
    int open_STDP=1;//开启STDP
    int N_syn = 0;//定义神经元的入度
    double tau=10.0;//时间常数，超参数
    double t_lastspike=0.0;//神经元上次发放的时间
    //t_sp是脉冲发生的事件，K_minus是t_sp时刻的K_,counter_sp计算脉冲信息被突触访问的次数
    deque<Histentry> history;//每一个神经元都维护这一个脉冲历史记录发放表
    double Kminus=0.0;//抑制因子
    //STDP突触需要的参数!!

    
    int spike_time=0;
    int spike_count = 0;//神经元的计数器
    vector<int> rate_record;//脉冲发射；频率记录表
    bool on_record = false;//是否开启记录
    
    int debug=1;//测试，可以删除
    //Simulation config (may not all be needed!!)
    double dt=0.1;//没啥用，因为dt不是在这儿设置的
    //double t_rest=0;
    
    //记录信息 
    int Record=1;//是否记录Vm，true=1，false=0
    int Poisson_switch=0;//是否开启神经元的Poisson发生器（注意，Poisson发生器的逻辑是镶嵌于神经元内部的，这样在分布式中占有优势）
    double Recording_Vm[RING_BUFFER_LEN];
    int Recording_spikes[RING_BUFFER_LEN];

    double Ring_buffer[RING_BUFFER_LEN];
    double Poisson_rate=0;
    double Poisson_weight=0;
    int global_num=1;//神经元编号

    double Vm=-65;//初始化膜电位
    double Vth=-50;//发射阈值
    double V_reset=-65;//恢复变量
    // double time;
    // double spikes;
    //LIF Properties
    double t=0;
    double Rm=1;
    double Cm=250;
    // double Cm=10;
    double tau_m=Rm*Cm;
    //int less_i=3;//收到脉冲后，对后几步进行影响
    double less_i_table[less_i];
    int tau_ref=15; //不应期时间//10的话可以固定
    int ref_less=0;//剩下的不应期时间
    //double V_spike=1;
    string type="Leaky Integrate and Fire";

    // double Vm=0;//初始化膜电位
    // double Vth=50;//发射阈值
    // double V_reset=0;//恢复变量
    // // double time;
    // // double spikes;
    // //LIF Properties
    // double t=0;
    // double Rm=1;
    // double Cm=10;
    // double tau_m=Rm*Cm;
    // //int less_i=3;//收到脉冲后，对后几步进行影响
    // double less_i_table[less_i];
    // int tau_ref=2; //不应期时间
    // int ref_less=0;//剩下的不应期时间
    // //double V_spike=1;
    // string type="Leaky Integrate and Fire";

    //STDP函数
        //记录神经元的脉冲发放历史
    void set_spiketime(double t_spike);
        //返回（t1，t2]时刻的历史发放记录
    void get_history(double t1,double t2,deque<Histentry>::iterator*start, deque<Histentry>::iterator* end);
        //返回t处的Kminus（突触轨迹）值
    double get_K_value(double t);
    //STDP函数

    LIFNeuron(int global_num_=-1)
    {
        global_num=global_num_;
        for (size_t i = 0; i < RING_BUFFER_LEN; i++)//初始化Ring_buffer跟Recording_spikes
        {
            Ring_buffer[i]=0;
            Recording_spikes[i]=0;
        }
        for (size_t i = 0; i < less_i; i++)//初始化残余电流
        {
            less_i_table[i]=0;
        }
        
        
        
    }
    int Get_globel_num()
    {
        return global_num;
    }

    void Clear_Ring_buffer()
    {
        for (size_t i = 0; i < RING_BUFFER_LEN; i++)//初始化Ring_buffer
        {
            Ring_buffer[i]=0;
        }
    }

    void Open_Poisson(double rate,double weight)//此处的rate是指每个step的脉冲强度
    {
        Poisson_rate=rate;
        Poisson_weight=weight;
        Poisson_switch=1;
    }
    void Closs_Poisson()
    {
        Poisson_switch=0;
    }

    void Poisson_inject()//给每步Ring_buffer添加Poisson概率的rate
    {
        Poisson_generator ps(Poisson_rate,Poisson_weight);
        double re_spikes[RING_BUFFER_LEN];
        ps.Get_poisson_num(re_spikes);
        for (size_t i = 0; i < RING_BUFFER_LEN; i++)
        {
            Ring_buffer[i]+=re_spikes[i];
        }
        // if(debug)
        // {
        //     cout<<"Ring_buffer is :";
        //     for (size_t i = 0; i < RING_BUFFER_LEN; i++)
        //     {
        //         cout<<" "<<Ring_buffer[i];
        //     }
        //     cout<<endl;
            
        // }
    }
    // void spike(int time)
    // {
    //     //发射脉冲之后需要获取脉冲发射时间，然后记录进ring_buffer
    // }
    // void updata(double I_syn[RING_BUFFER_LEN])  //测试版，可自由输入spikes
    void updata(int loop_step)
    {

        if(on_record)
        {
            if(loop_step % rate_record_step == 0)
            {
                double spike_rate = spike_count / (rate_record_step * 0.001);
                rate_record.push_back(spike_rate);
                spike_count  = 0;
            }
        }
        set_less_i_to_ring_buffer();
        if(Poisson_switch==1)
        {
            Poisson_inject();
        }
       for (size_t i = 0; i < RING_BUFFER_LEN; i++)
       
        {
            if(ref_less==0)
            {    
                if(Record)
                {
                    Recording_Vm[i]=Vm;
                }
                //Recording_spikes[i]=0;
                Vm=Vm+((Rm*Ring_buffer[i]-(Vm-V_reset))/tau_m)*dt;
                Recording_spikes[i]=0;//记录脉冲，作为该神经元的输出
            }
            else
            {   
                if(Record)
                {
                    Recording_Vm[i]=Vm;                   
                }
                Recording_spikes[i]=0;
                ref_less--;
                // spike(i);
            }
            if(Vm>Vth)
            {
                spike_time++;
                if(Record)
                {
                    Recording_Vm[i]=Vm;
                }
                ref_less=tau_ref;
                Vm=V_reset;
                Recording_spikes[i]=1;//记录下，发射了脉冲
                double spike_time=double(loop_step)+0.1*i;
                
                if(print_spike_and_recv){
                cout<<"n : "<<global_num<<" spiked!!"<<" time: "<<spike_time<<endl;
                }

                if(on_record){
                spike_count += 1;
                }
                //STDP
                if(open_STDP)
                {
                   set_spiketime(spike_time); 
                }
                //STDP

            }
           
        }
        Clear_Ring_buffer();
        clear_less_i_table();
    }

    void set_less_i_to_ring_buffer()//输入上一轮留下的残余电流
    {
        for (size_t i = 0; i < less_i; i++)
        {
            Ring_buffer[i]+=less_i_table[i];
        }
        
    }

    void clear_less_i_table()
    {
        for (size_t i = 0; i < less_i; i++)
        {
            less_i_table[i]=0;
        }
        
    }

    // double Get_Vm()
    // {
    //     return Vm;
    // }
    double Set_Vm(double Vm_)
    {
        Vm=Vm_;
        return Vm;
    }

    double Get_Vth()
    {
        return Vth;
    }
    double Set_Vth(double Vth_)
    {
        Vth=Vth_;
        return Vth;
    }

    double Get_V_reset()
    {
        return V_reset;
    }
    double Set_V_reset(double V_reset_)
    {
        V_reset=V_reset_;
        return V_reset;
    }

    void Get_Recording_Vm(double Recording_Vm_[RING_BUFFER_LEN])
    {
        for (size_t i = 0; i < RING_BUFFER_LEN; i++)
        {
            Recording_Vm_[i]=Recording_Vm[i];
        }  
    }  

    void Get_Recording_spikes(int Recording_spikes_[RING_BUFFER_LEN])
    {
        for (size_t i = 0; i < RING_BUFFER_LEN; i++)
        {
            Recording_spikes_[i]=Recording_spikes[i];
        }
    } 

    void Set_Ring_buffer(double inject_spikes[RING_BUFFER_LEN])//注入脉冲到Ring_buffer
    {
        for (size_t i = 0; i < RING_BUFFER_LEN; i++)
        {
            Ring_buffer[i]+=inject_spikes[i];
            if(inject_spikes[i]>0)//设置一个残余电流机制
            {
                if (i<(RING_BUFFER_LEN-less_i))
                {
                    for (size_t j = 0; j < less_i; j++)
                    {
                        Ring_buffer[i+j+1]+=inject_spikes[i]*(1/((j+1)*2));
                    }
                    
                }
                else
                {
                    int less_in=0;
                    for (size_t j = 0; j < less_i; j++)
                    {
                        if ((i+j+1)<RING_BUFFER_LEN)
                        {
                            Ring_buffer[i+j+1]+=inject_spikes[i]*(1/((j+1)*2));
                        }
                        else
                        {
                            less_i_table[less_in++]+=inject_spikes[i]*(1/((j+1)*2));
                        }
                    }
                    
                }
                
            }
        }
        
    }

    void Set_one_Ring_buffer(int mini_step,double weight,int loop_step)//第一个是注入的时间步，第二个是权重
    {
        if(print_spike_and_recv)
        {
        cout<<"n : "<<global_num<<" recved spike! "<<"time : "<<loop_step<<"."<<mini_step<<endl;
        }

        Ring_buffer[mini_step]+=weight;
        
        if (mini_step<(RING_BUFFER_LEN-less_i))//残余电流机制
        {
            for (size_t j = 0; j < less_i; j++)
            {
                Ring_buffer[mini_step+j+1]+=weight*(1/((j+1)*2));
            }
            
        }
        else
        {
            int less_in=0;
            for (size_t j = 0; j < less_i; j++)
            {
                if ((mini_step+j+1)<RING_BUFFER_LEN)
                {
                    Ring_buffer[mini_step+j+1]+=weight*(1/((j+1)*2));
                }
                else
                {
                    less_i_table[less_in++]+=weight*(1/((j+1)*2));
                }
            }
            
        }
    }

    void Set_Record(int record_)//是否记录Vm的开关
    {
        Record=record_;
    }

};
//STDP使用的函数
void LIFNeuron::set_spiketime(double t_spike)
{
    //当神经元的入度不为零时，查看历史记录的某一条记录是否被至少使用N_syn次
    //如果超出这个，则证明此条记录不会使用，放出队列
    if (N_syn > 0)
    {

        while (history.size() > 1)
        {
            if (history.front().count_sp >= N_syn)
            {
                history.pop_front();    
            }
            else
            {
                break;
            }
        }
    }
    //将此时刻的Kminus和时间记录到记录表中
    Kminus = Kminus * std::exp(-(t_spike - t_lastspike) / tau) + 1; 
    history.push_back(Histentry(t_spike, Kminus, 0));
    t_lastspike = t_spike;
}

double LIFNeuron::get_K_value(double t)
{
    double K_value;
    //当记录表为零时，返回0
    if (history.empty())
    {
        K_value = 0.0;
        return K_value;
    }
    //从队列的队尾开始，加快访问速度
    deque<Histentry>::reverse_iterator it;
    for (it = history.rbegin(); it != history.rend(); it++)
    {
        if (t - it->t_sp > 0)
        {
            K_value = it->Kminus * std::exp((it->t_sp - t) / tau);
            return K_value;
        }
    }
    //如果遍历表没有的话，默认返回0
    K_value = 0.0;
    return K_value;
}
// double LIFNeuron::get_K_value(double t)
// {
//     double K_value;
//     //当记录表为零时，返回0
//     if (history.empty())
//     {
//         K_value = 0.0;
//         return K_value;
//     }
//     //从队列的队尾开始，加快访问速度
//     for (int i = history.size() - 1; i >=0 ; i--)
//     {
//         if (t - history[i].t_sp > 0)
//         {
//             K_value = Kminus * std::exp((history[i].t_sp - t) / tau);
//             return K_value;
//         }
//     }
//     //如果遍历表没有的话，默认返回0
//     K_value = 0.0;
//     return K_value;
// }

void LIFNeuron::get_history(double t1, double t2, deque<Histentry>::iterator* start, deque<Histentry>::iterator* end)
{
    *end = history.end();
    //当队列为空，返回空指针
    while (history.empty())
    {
        *start = *end;
        return;
    }

    deque<Histentry>::reverse_iterator it = history.rbegin();
    //寻找（t1，t2]时刻的记录，返回指向这两个位置的指针
    while (it != history.rend() && it->t_sp > t2)
    //while (it != history.rend() && it->t_sp < t2)
    {
        it++;
    }
    *end = it.base();
    while (it != history.rend() && it->t_sp > t1)
    {
        it->count_sp += 1;
        it++;
    }
    *start = it.base();
}

//STDP使用的函数




// #include <vector>
// // #include <Event.cpp>
// #include<deque>
// #include <string>
// #include"Histentry.h"
// #include <cmath>
// #include <iostream>
// #include "Poisson_generator.cpp"
// using namespace std;
// #define RING_BUFFER_LEN 10//ring buffer的长度，我们假设突触延迟都相同，都是这个(ms).多个文件需要同时修改！
// #define less_i 4 
// #define rate_record_step 10//记录步数
// //#define Record 1//是否记录Vm跟spike，是1，否0
// class LIFNeuron
// {
// private:


// public:
//     int spike_count = 0;//神经元的计数器
//     vector<int> rate_record;//记录表
//     bool on_record = false;//是否开启记录

//     //STDP突触需要的参数!!
//     int N_syn = 0;//定义神经元的入度
//     double tau=10.0;//时间常数，超参数
//     double t_lastspike=0.0;//神经元上次发放的时间
//     //t_sp是脉冲发生的事件，K_minus是t_sp时刻的K_,counter_sp计算脉冲信息被突触访问的次数
//     deque<Histentry> history;//每一个神经元都维护这一个脉冲历史记录发放表
//     double Kminus=0.0;//抑制因子
//     //STDP突触需要的参数!!

//     int debug=1;//测试，可以删除
//     //Simulation config (may not all be needed!!)
//     double dt=0.1;//没啥用，因为dt不是在这儿设置的
//     //double t_rest=0;
    
//     //记录信息 
//     int Record=1;//是否记录Vm，true=1，false=0
//     int Poisson_switch=0;//是否开启神经元的Poisson发生器（注意，Poisson发生器的逻辑是镶嵌于神经元内部的，这样在分布式中占有优势）
//     double Recording_Vm[RING_BUFFER_LEN];
//     int Recording_spikes[RING_BUFFER_LEN];

//     double Ring_buffer[RING_BUFFER_LEN];
//     double Poisson_rate=0;
//     int global_num=1;//神经元编号
//     double Vm=-65;//初始化膜电位
//     double Vth=-50;//发射阈值
//     double V_reset=-65;//恢复变量
//     // double time;
//     // double spikes;
//     //LIF Properties
//     double t=0;
//     double Rm=1;
//     double Cm=250;
//     // double Cm=10;
//     double tau_m=Rm*Cm;
//     //int less_i=3;//收到脉冲后，对后几步进行影响
//     double less_i_table[less_i];
//     int tau_ref=20; //不应期时间
//     int ref_less=0;//剩下的不应期时间
//     //double V_spike=1;
//     string type="Leaky Integrate and Fire";

//     //STDP函数
//         //记录神经元的脉冲发放历史
//     void set_spiketime(double t_spike);
//         //返回（t1，t2]时刻的历史发放记录
//     void get_history(double t1,double t2,deque<Histentry>::iterator*start, deque<Histentry>::iterator* end);
//         //返回t处的Kminus（突触轨迹）值
//     double get_K_value(double t);
//     //STDP函数

//     LIFNeuron(int global_num_=-1)
//     {
//         global_num=global_num_;
//         for (size_t i = 0; i < RING_BUFFER_LEN; i++)//初始化Ring_buffer跟Recording_spikes
//         {
//             Ring_buffer[i]=0;
//             Recording_spikes[i]=0;
//         }
//         for (size_t i = 0; i < less_i; i++)//初始化残余电流
//         {
//             less_i_table[i]=0;
//         }
        
        
        
//     }
//     int Get_globel_num()
//     {
//         return global_num;
//     }

//     void Clear_Ring_buffer()
//     {
//         for (size_t i = 0; i < RING_BUFFER_LEN; i++)//初始化Ring_buffer
//         {
//             Ring_buffer[i]=0;
//         }
//     }

//     void Open_Poisson(double rate)//此处的rate是指每个step的脉冲强度
//     {
//         Poisson_rate=rate;
//         Poisson_switch=1;
//     }
//     void Closs_Poisson()
//     {
//         Poisson_switch=0;
//     }

//     void Poisson_inject()//给每步Ring_buffer添加Poisson概率的rate
//     {
//         Poisson_generator ps(Poisson_rate);
//         double re_spikes[RING_BUFFER_LEN];
//         ps.Get_poisson_num(re_spikes);
//         for (size_t i = 0; i < RING_BUFFER_LEN; i++)
//         {
//             Ring_buffer[i]+=re_spikes[i];
//         }
//         // if(debug)
//         // {
//         //     cout<<"Ring_buffer is :";
//         //     for (size_t i = 0; i < RING_BUFFER_LEN; i++)
//         //     {
//         //         cout<<" "<<Ring_buffer[i];
//         //     }
//         //     cout<<endl;
            
//         // }
//     }
//     // void spike(int time)
//     // {
//     //     //发射脉冲之后需要获取脉冲发射时间，然后记录进ring_buffer
//     // }
//     // void updata(double I_syn[RING_BUFFER_LEN])  //测试版，可自由输入spikes
//     void updata(int loop_step)
//     {
//         if(on_record)
//         {
//             if(loop_step % rate_record_step == 0)
//             {
//                 double spike_rate = spike_count / (rate_record_step * 0.001);
//                 rate_record.push_back(spike_rate);
//                 spike_count  = 0;
//             }
//         }

//         set_less_i_to_ring_buffer();
//         if(Poisson_switch==1)
//         {
//             Poisson_inject();
//         }
//        for (size_t i = 0; i < RING_BUFFER_LEN; i++)
       
//         {
//             if(ref_less==0)
//             {    
//                 if(Record)
//                 {
//                     Recording_Vm[i]=Vm;
//                 }
//                 //Recording_spikes[i]=0;
//                 Vm=Vm+((Rm*Ring_buffer[i]-(Vm-V_reset))/tau_m)*dt;
//                 Recording_spikes[i]=0;//记录脉冲，作为该神经元的输出
//             }
//             else
//             {   
//                 if(Record)
//                 {
//                     Recording_Vm[i]=Vm;                   
//                 }
//                 Recording_spikes[i]=0;
//                 ref_less--;
//                 // spike(i);
//             }
//             if(Vm>Vth)
//             {
//                 if(Record)
//                 {
//                     Recording_Vm[i]=Vm;
//                 }
//                 ref_less=tau_ref;
//                 Vm=V_reset;
//                 Recording_spikes[i]=1;//记录下，发射了脉冲
//                 cout<<"n : "<<global_num<<" spiked!!"<<" time: "<<loop_step<<"."<<i<<endl;
//             }
           
//         }
//         Clear_Ring_buffer();
//         clear_less_i_table();
//     }

//     void set_less_i_to_ring_buffer()//输入上一轮留下的残余电流
//     {
//         for (size_t i = 0; i < less_i; i++)
//         {
//             Ring_buffer[i]+=less_i_table[i];
//         }
        
//     }

//     void clear_less_i_table()
//     {
//         for (size_t i = 0; i < less_i; i++)
//         {
//             less_i_table[i]=0;
//         }
        
//     }

//     // double Get_Vm()
//     // {
//     //     return Vm;
//     // }
//     double Set_Vm(double Vm_)
//     {
//         Vm=Vm_;
//         return Vm;
//     }

//     double Get_Vth()
//     {
//         return Vth;
//     }
//     double Set_Vth(double Vth_)
//     {
//         Vth=Vth_;
//         return Vth;
//     }

//     double Get_V_reset()
//     {
//         return V_reset;
//     }
//     double Set_V_reset(double V_reset_)
//     {
//         V_reset=V_reset_;
//         return V_reset;
//     }

//     void Get_Recording_Vm(double Recording_Vm_[RING_BUFFER_LEN])
//     {
//         for (size_t i = 0; i < RING_BUFFER_LEN; i++)
//         {
//             Recording_Vm_[i]=Recording_Vm[i];
//         }  
//     }  

//     void Get_Recording_spikes(int Recording_spikes_[RING_BUFFER_LEN])
//     {
//         for (size_t i = 0; i < RING_BUFFER_LEN; i++)
//         {
//             Recording_spikes_[i]=Recording_spikes[i];
//         }
//     } 

//     void Set_Ring_buffer(double inject_spikes[RING_BUFFER_LEN])//注入脉冲到Ring_buffer
//     {
//         for (size_t i = 0; i < RING_BUFFER_LEN; i++)
//         {
//             Ring_buffer[i]+=inject_spikes[i];
//             if(inject_spikes[i]>0)//设置一个残余电流机制
//             {
//                 if (i<(RING_BUFFER_LEN-less_i))
//                 {
//                     for (size_t j = 0; j < less_i; j++)
//                     {
//                         Ring_buffer[i+j+1]+=inject_spikes[i]*(1/((j+1)*2));
//                     }
                    
//                 }
//                 else
//                 {
//                     int less_in=0;
//                     for (size_t j = 0; j < less_i; j++)
//                     {
//                         if ((i+j+1)<RING_BUFFER_LEN)
//                         {
//                             Ring_buffer[i+j+1]+=inject_spikes[i]*(1/((j+1)*2));
//                         }
//                         else
//                         {
//                             less_i_table[less_in++]+=inject_spikes[i]*(1/((j+1)*2));
//                         }
//                     }
                    
//                 }
                
//             }
//         }
        
//     }

//     void Set_one_Ring_buffer(int mini_step,double weight,int loop_step)//第一个是注入的时间步，第二个是权重
//     {
//         cout<<"n : "<<global_num<<" recved spike! "<<"time : "<<loop_step<<"."<<mini_step<<endl;

//         Ring_buffer[mini_step]+=weight;
        
//         if (mini_step<(RING_BUFFER_LEN-less_i))//残余电流机制
//         {
//             for (size_t j = 0; j < less_i; j++)
//             {
//                 Ring_buffer[mini_step+j+1]+=weight*(1/((j+1)*2));
//             }
            
//         }
//         else
//         {
//             int less_in=0;
//             for (size_t j = 0; j < less_i; j++)
//             {
//                 if ((mini_step+j+1)<RING_BUFFER_LEN)
//                 {
//                     Ring_buffer[mini_step+j+1]+=weight*(1/((j+1)*2));
//                 }
//                 else
//                 {
//                     less_i_table[less_in++]+=weight*(1/((j+1)*2));
//                 }
//             }
            
//         }
//     }

//     void Set_Record(int record_)//是否记录Vm的开关
//     {
//         Record=record_;
//     }

// };
// //STDP使用的函数
// void LIFNeuron::set_spiketime(double t_spike)
// {
//     //当神经元的入度不为零时，查看历史记录的某一条记录是否被至少使用N_syn次
//     //如果超出这个，则证明此条记录不会使用，放出队列
//     if (N_syn > 0)
//     {

//         while (history.size() > 1)
//         {
//             if (history.front().count_sp >= N_syn)
//             {
//                 history.pop_front();    
//             }
//             else
//             {
//                 break;
//             }
//         }
//     }
//     //将此时刻的Kminus和时间记录到记录表中
//     Kminus = Kminus * std::exp(-(t_spike - t_lastspike) / tau) + 1; 
//     history.push_back(Histentry(t_spike, Kminus, 0));
//     t_lastspike = t_spike;
// }

// double LIFNeuron::get_K_value(double t)
// {
//     double K_value;
//     //当记录表为零时，返回0
//     if (history.empty())
//     {
//         K_value = 0.0;
//         return K_value;
//     }
//     //从队列的队尾开始，加快访问速度
//     for (int i = history.size() - 1; i >=0 ; i--)
//     {
//         if (t - history[i].t_sp > 0)
//         {
//             K_value = Kminus * std::exp((history[i].t_sp - t) / tau);
//             return K_value;
//         }
//     }
//     //如果遍历表没有的话，默认返回0
//     K_value = 0.0;
//     return K_value;
// }

// void LIFNeuron::get_history(double t1, double t2, deque<Histentry>::iterator* start, deque<Histentry>::iterator* end)
// {
//     *end = history.end();
//     //当队列为空，返回空指针
//     while (history.empty())
//     {
//         *start = *end;
//         return;
//     }

//     deque<Histentry>::reverse_iterator it = history.rbegin();
//     //寻找（t1，t2]时刻的记录，返回指向这两个位置的指针
//     while (it != history.rend() && it->t_sp > t2)
//     //while (it != history.rend() && it->t_sp < t2)
//     {
//         it++;
//     }
//     *end = it.base();
//     while (it != history.rend() && it->t_sp > t1)
//     {
//         it->count_sp += 1;
//         it++;
//     }
//     *start = it.base();
// }
// //STDP使用的函数

// // int main(int argc, char const *argv[])
// // {
// //  /* code */
// //     double spike_in[RING_BUFFER_LEN]={0,2,1000,70,3,1,7,1,1,1};
// //     double Recording_Vm[RING_BUFFER_LEN];
// //     int Recording_spikes[RING_BUFFER_LEN];
// //     LIFNeuron neuron;
    
// //         neuron.Set_Vm(0.2);
// //         neuron.updata(spike_in);
// //         neuron.Get_Recording_Vm(Recording_Vm);
// //         neuron.Get_Recording_spikes(Recording_spikes);
// //         for (size_t i = 0; i < RING_BUFFER_LEN; i++)
// //         {
// //             cout<<i<<" Vm is "<<Recording_Vm[i]<<endl;
// //             cout<<i<< " spike is"<<Recording_spikes[i]<<endl;
// //         }
        

    
    
// //     cout<<"finished"<<endl;
// //  return 0;
// // }