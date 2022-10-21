#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

vector<vector<int>> read_STDP_adj_table(string STDP_adj_table_txt)
{
    vector<vector<int>> temp_STDP_adj_table;
    int line_count=0;
    int line_index=0;
    int i=0;
    int datalen=0;
    int num;
    ifstream file(STDP_adj_table_txt,ios::in);
    while( ! file.eof() ){
    file>>num;
    if (num==-1)
    {
        vector<int> adj;
        temp_STDP_adj_table.push_back(adj);
        line_count++;
    }
    if (i!=0&&num==-1)
    {
        line_index++;
    }
    temp_STDP_adj_table[line_index].push_back(num); 
    i++;
    //cout<<"line_index is : "<<line_index<<endl;
    //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
    }
    file.close();


    return temp_STDP_adj_table;          
}

vector<vector<int>> read_STDP_in_table(string STDP_adj_table_txt)
{
    cout<<"read_STDP_in_1"<<endl;
    vector<vector<int>> STDP_out_adj_table=read_STDP_adj_table(STDP_adj_table_txt);
    cout<<"read_STDP_in_2"<<endl;
    vector<vector<int>> STDP_in_adj_table;
    for (size_t i = 0; i < STDP_out_adj_table.size(); i++)
    {
        vector<int> adj={-1};
        STDP_in_adj_table.push_back(adj);
    }
    cout<<"STDP_in_adj_table size : "<<STDP_in_adj_table.size()<<endl;
    cout<<"read_STDP_in_3"<<endl;
    for (size_t i = 0; i < STDP_out_adj_table.size(); i++)
    {
        for (size_t j = 1; j < STDP_out_adj_table[i].size(); j++)
        {
            cout<<"STDP_out_adj_table.size() is :"<<STDP_out_adj_table.size()<<endl;
            cout<<"STDP_out_adj_table[i].size() is :"<<STDP_out_adj_table[i+1].size()<<endl;
            cout<<"STDP_out_adj_table["<<i<<"]["<<j<<"] is : "<<STDP_out_adj_table[i][j]<<endl;
            STDP_in_adj_table[STDP_out_adj_table[i][j]].push_back(i);
        }
        
    }
    cout<<"read_STDP_in_4"<<endl;
    return STDP_in_adj_table;

}



vector<vector<int>> read_adj_table(string adj_table_txt)
{
    vector<vector<int>> temp_adj_table;
    int line_count=0;
    int line_index=0;
    int i=0;
    int datalen=0;
    int num;
    ifstream file(adj_table_txt,ios::in);
    while( ! file.eof() ){
    file>>num;
    if (num==-1)
    {
        vector<int> adj;
        temp_adj_table.push_back(adj);
        line_count++;
    }
    if (i!=0&&num==-1)
    {
        line_index++;
    }
    temp_adj_table[line_index].push_back(num);   
    i++;
    //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
    }
    file.close();


    return temp_adj_table;         
}

void record_res_hypergraph(vector<vector<double>> res_hypergraph,string adders)//记录adjtable到txt
{
    ofstream ofs;
    ofs.open(adders,ios::out);
    for (size_t i = 0; i < res_hypergraph.size(); i++)
    {
        for (size_t j = 0; j < res_hypergraph[i].size(); j++)
        {
            ofs<<res_hypergraph[i][j]<<" ";
        }
        ofs<<endl;
        
    }
    ofs.close();
}

int main(int argc, char const *argv[])
{
    // string adj_table_txt="adj_table.txt";
    // string STDP_adj_table_txt="STDP_adj_table.txt";
    string adj_table_txt="/home/xc/Documents/SNN_c++_lif_1.11/part_8_record/adj.txt";
    string STDP_adj_table_txt="/home/xc/Documents/SNN_c++_lif_1.11/part_8_record/STDP_adj.txt";
    vector<vector<int>> adj_table=read_adj_table(adj_table_txt);//读取静态连接
    vector<vector<int>> STDP_adj_table=read_STDP_adj_table(STDP_adj_table_txt);//读取动态连接
    vector<vector<int>> STDP_in_table=read_STDP_in_table(STDP_adj_table_txt);//读取动态连接入度表
    vector<double> Neuron_weight;//每个神经元的权重
    //vector<double> Neuron_spikes;
    double pre_spikes=75;
    double STDP_1_time_1_cal_weight=0.0000081;
    double Neuron_1_weight=0.0584;
    // double pre_spikes=375;
    // double STDP_1_time_1_cal_weight=0.0000280673;
    // double Neuron_1_weight=0.2920;
    vector<vector<double>> total_adj_table;//总连接表
    for (size_t i = 0; i < adj_table.size(); i++)
    {
        vector<double> adj;
        total_adj_table.push_back(adj);
    }

    for (int i = 0; i < adj_table.size(); i++)
    {
        total_adj_table[i].push_back(i+1);//把本神经元编号放入(!!!!注意，神经元编号在运行超图算法的时候需要从1开始!!!!)
        for (int j = 1; j < adj_table[i].size(); j++)
        {
            total_adj_table[i].push_back(adj_table[i][j]+1);
        }

        for (int j = 1; j < STDP_adj_table[i].size(); j++)
        {
            total_adj_table[i].push_back(STDP_adj_table[i][j]+1);
        }
    }
    
    // for (int i = 0; i < STDP_in_table.size(); i++)
    // {
    //     int STDP_in_degree=(STDP_in_table[i].size()-1);
    //     double Node_weight=STDP_in_degree*pre_spikes*STDP_1_time_1_cal_weight+Neuron_1_weight;
    //     Neuron_weight.push_back(Node_weight*100);

    // }
    // cout<<"Neuron_weight size() : "<<Neuron_weight.size()<<endl;
    cout<<"total_adj_table size() : "<<total_adj_table.size()<<endl;
 

    vector<vector<double>> res_hypergraph;//需要返回的超图
    vector<double> begin_vec;//超图的划分信息；
    begin_vec.push_back(total_adj_table.size());
    begin_vec.push_back(total_adj_table.size());
    //begin_vec.push_back(int(10));
    res_hypergraph.push_back(begin_vec);
    for (size_t i = 0; i < total_adj_table.size(); i++)
    {
        res_hypergraph.push_back(total_adj_table[i]);
    }

    // for (size_t i = 0; i < Neuron_weight.size(); i++)
    // {
    //     vector<double> adj={Neuron_weight[i]};
    //     res_hypergraph.push_back(adj);
    // }

    cout<<"res_hypegraph size : "<<res_hypergraph.size()<<endl;
    string hypergraph_2="hyper2_noweight.txt";
    record_res_hypergraph(res_hypergraph,hypergraph_2);

    
    

    

    return 0;
}
