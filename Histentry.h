#pragma once
class Histentry
{
public:
	double t_sp; //脉冲发生的时间
	double Kminus;//此时刻对应的K_
	int count_sp;//计数器
	Histentry(double t_sp,double Kminus,int count_sp);

};

Histentry::Histentry(double t_sp, double Kminus,int count_sp)
	:t_sp(t_sp),
	Kminus(Kminus),
	count_sp(count_sp)
{
}