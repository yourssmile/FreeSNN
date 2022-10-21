#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <cstdlib>
#include <ctime>
class Event
{
private:
    int Source_G_id=-1;
    double time=-1;
    int mini_step=-1;//为了更新设置，STDP不用理会
    double weight=-1;
public:
    Event(/* args */)
    {

    }

    void set_Source_G_id(int Source_G_id_)
    {
        Source_G_id=Source_G_id_;
    }
    int get_Source_G_id()
    {
        return Source_G_id;
    }

    void set_time(double time_)
    {
        time=time_;
    }
    double get_time()
    {
        return time;
    }

    void set_mini_step(int step)
    {
        mini_step=step;
    }
    int get_mini_step()
    {
        return mini_step;
    }

    void set_weight(double weight_)
    {
        weight=weight_;
    }
    double get_weight()
    {
        return weight;
    }

    ~Event()
    {
        
    }

};

