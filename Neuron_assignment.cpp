#include<vector>
#include<iostream>
using namespace std;
class Neuron_assignment
{
private:
    /* data */
public:
    Neuron_assignment(/* args */);
    
    //传入参数分别是，mpi进程数，静态连接表，STDP连接表。返回一个二维vector，第一维是进程，第二维是该进程分配到的神经元
    vector<vector<int>> rotation_assignment(int MPI_size, vector<vector<int>> adject_table,vector<vector<int>> STDP_adj_table);
    
    ~Neuron_assignment();

};

Neuron_assignment::Neuron_assignment(/* args */)
{
}

Neuron_assignment::~Neuron_assignment()
{
}

vector<vector<int>> Neuron_assignment::rotation_assignment(int MPI_size, vector<vector<int>> adject_table,vector<vector<int>> STDP_adj_table)
{
    if(MPI_size<1)
    {
        cout<<"MPI_size<1 ERROR!"<<endl;
    }
    int neuron_size=adject_table.size();
    int Rank_index=0;
    vector<vector<int>> res;
    res.resize(MPI_size);
    for (size_t i = 0; i < neuron_size; i++)
    {
        if(Rank_index<MPI_size)
        {
            res[Rank_index].push_back(i);
            Rank_index++;
        }
        else
        {
            Rank_index=0;
            res[Rank_index].push_back(i);
            Rank_index++;   
        }
    }

    return res;
    

}
