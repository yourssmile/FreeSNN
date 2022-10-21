#pragma once
//#include"LIFNeuron.cpp"
#include <cmath>
#include <iostream>
#include<deque>
#include "Histentry.h"
using namespace std;
class STDP
{
public:
	STDP();

    double facilitate(double w, double Kplus);
    double depresss(double w, double Kminus);
    void update_weight(double t_spike,LIFNeuron target);
    void set_weight(double weight_);
    void set_source_Gid(int Gid);
    void set_target_Gid(int Gid);
    double weight;//STDP权重
    int source_Gid=-1;
    int target_Gid=-1;

private:
    //以下为STDP公式的超参数
    double lamda;
    double alpha;
    double mu;
    double Wmax;
    double tau;


    double Kplus;//增强因子，每一个STDP突触都有

    double t_lastspike;//该突触上次发放脉冲的时间(该突触上次接收到source神经元脉冲的时间)


};
STDP::STDP()
    : weight(10.0),tau(20.0),lamda(0.1),alpha(0.5),Wmax(2000.0),t_lastspike(0.0),Kplus(0.0),mu(1.0)
{

};
void STDP::set_weight(double weight_)
{
    weight=weight_;
}
void STDP::set_source_Gid(int Gid)
{
    source_Gid=Gid;
}
void STDP::set_target_Gid(int Gid)
{
    target_Gid=Gid;
}
//w代表当前权重，Kplus时增强程度，返回增强后的权重
double STDP::facilitate(double w, double Kplus)
{
    double norm_w = (w / Wmax) + (lamda * pow(1.0 - (w / Wmax), mu) * Kplus);
    return norm_w < 1.0 ? norm_w * Wmax : Wmax;
    //return norm_w < Wmax ? norm_w * Wmax : Wmax;
};
//w代表当前权重，Kminus时抑制程度，返回抑制后的权重
double STDP::depresss(double w, double Kminus)
{
    double norm_w = (w / Wmax) - (alpha * lamda * pow(w / Wmax, mu) * Kminus);
    return norm_w > 0.0 ? norm_w * Wmax : 0.0;
};
//当脉冲传递到后神经元时，调用此函数，t_spike代表源神经元的发放时间，
//target代表后神经元
void STDP::update_weight(double t_spike, LIFNeuron target)
{
    double delay = 1.0;
    std::deque<Histentry>::iterator start;
    std::deque<Histentry>::iterator end;
    //从突触后神经元获取相关范围（t1、t2]内的尖峰历史
    //历史记录（t_lastspike - delay，…，t_spike-delay]
    target.get_history(t_lastspike - delay, t_spike - delay, &start, &end);

    //在前神经元的脉冲发放时间后后神经元发放对突触权重起到增强作用
    double dt;
    while (start != end)
    {
        dt = t_lastspike - (start->t_sp + delay);
        ++start;
        weight = facilitate(weight, Kplus * exp(dt / tau));
    }

    //突触后神经元脉冲发作时的抑制应仅取决于自最近一次突触前尖峰以来的时间。
    double Kmius = target.get_K_value(t_spike - delay);
    weight = depresss(weight, Kmius);
    //更新Kplus值
    Kplus = Kplus * exp((t_lastspike - t_spike) / tau) + 1;
    //更新突触脉冲发放值
    t_lastspike = t_spike;

    //std::cout << weight << endl;
};