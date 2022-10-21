#include <vector>
#include <string>
#include <iostream>
#include "Mpi_Manage_9_27.cpp"
using namespace std;
int main(int argc, char const *argv[])
{

//MPI 测试
    Mpi_Manage ma;
    // cout<<ma.MPI_rank<<" "<<ma.MPI_size<<endl;
    double pop_de_rate=1;
    int pop_num=0;
    int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};

    
    if(ma.MPI_rank==0)//这些是在一个进程上做
    {
    
    double syn_weight=3.14;
    
    // int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};
    double connect_pr[8][8]={
        {0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0., 0.0076, 0.},
             {0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0., 0.0042, 0.},
             {0.0077, 0.0059, 0.0497, 0.135, 0.0067, 0.0003, 0.0453, 0.},
             {0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0., 0.1057, 0.},
             {0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.},
             {0.0548, 0.0269, 0.0257, 0.0022, 0.06, 0.3158, 0.0086, 0.},
             {0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252},
             {0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443}
    };
    // double connect_pr[8][8]={
    //     {0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0., 0.0076, 0.},
    //          {0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0., 0.0042, 0.},
    //          {0.0077, 0.0059, 0.0497, 0.135, 0.0067, 0.0003, 0.0453, 0.},
    //          {0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0., 0.1057, 0.},
    //          {0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.},
    //          {0.0548, 0.0269, 0.0257, 0.0022, 0.06, 0.3158, 0.0086, 0.},
    //          {0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252},
    //          {0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443}
    // };

    int L23E=ma.create(1,pop_nums[0]);
    int L23I=ma.create(1,pop_nums[1]);
    int L4E=ma.create(1,pop_nums[2]);
    int L4I=ma.create(1,pop_nums[3]);
    int L5E=ma.create(1,pop_nums[4]);
    int L5I=ma.create(1,pop_nums[5]);
    int L6E=ma.create(1,pop_nums[6]);
    int L6I=ma.create(1,pop_nums[7]);
    int populatins[8]={L23E,L23I,L4E,L4I,L5E,L5I,L6E,L6I};
    double mu[8][8]={
        {277,-1110,555,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110},
        {277,-1110,277,-1110,277,-1110,277,-1110}
    };
    double sigma[8][8]={
        {27,111,55,111,27,111,27,111},
        {27,111,27,111,27,111,27,111},
        {27,111,27,111,27,111,27,111},
        {27,111,27,111,27,111,27,111},
        {27,111,27,111,27,111,27,111},
        {27,111,27,111,27,111,27,111},
        {27,111,27,111,27,111,27,111},
        {27,111,27,111,27,111,27,111}
    };

    int conn_weight=0;
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            conn_weight=gaussrand()*sigma[i][j]*syn_weight+mu[i][j]*syn_weight;
            //ma.FIX_random_connect_bad_4(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
            ma.FIX_random_connect_STDP(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);            
            // if(conn_weight>0)
            // ma.STDP_random_connect(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
            // else
            // ma.random_connect(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
        }
        
    }
    ma.neuron_assignment(1);//根据进程以及神经元分配方案进行神经元分配，单进程
    }
    cout<<endl;
    // ma.prepare_msg();
    if(ma.MPI_rank==0){
    string Gid_r_txt="part_4_record/Gid_to_rank.txt";
    string adj_txt="part_4_record/adj.txt";
    string STDP_adj_txt="part_4_record/STDP_adj.txt";
    string w_txt="part_4_record/w.txt";
    string STDP_w_txt="part_4_record/STDP_w.txt";
    if(ma.MPI_size==4)
    {
    string Gid_r_txt="part_4_record/Gid_to_rank.txt";
    string adj_txt="part_4_record/adj.txt";
    string STDP_adj_txt="part_4_record/STDP_adj.txt";
    string w_txt="part_4_record/w.txt";
    string STDP_w_txt="part_4_record/STDP_w.txt";
    }
    else if(ma.MPI_size==8)
    {
    Gid_r_txt="part_8_record/Gid_to_rank.txt";
    adj_txt="part_8_record/adj.txt";
    STDP_adj_txt="part_8_record/STDP_adj.txt";
    w_txt="part_8_record/w.txt";
    STDP_w_txt="part_8_record/STDP_w.txt";
    }
    else if(ma.MPI_size==16)
    {
    Gid_r_txt="part_16_record/Gid_to_rank.txt";
    adj_txt="part_16_record/adj.txt";
    STDP_adj_txt="part_16_record/STDP_adj.txt";
    w_txt="part_16_record/w.txt";
    STDP_w_txt="part_16_record/STDP_w.txt";
    }
    else{
    string Gid_r_txt="else_record/Gid_to_rank.txt";
    string adj_txt="else_record/adj.txt";
    string STDP_adj_txt="else_record/STDP_adj.txt";
    string w_txt="else_record/w.txt";
    string STDP_w_txt="else_record/STDP_w.txt";        
    }
    ma.record_Gid_to_rank(Gid_r_txt);
    ma.record_adjtable(adj_txt);
    ma.record_STDP_adj_table(STDP_adj_txt);
    ma.record_weight_table(w_txt);
    ma.record_STDP_w_table(STDP_w_txt);
    cout<<"finished_record"<<endl;
    }
   


}
