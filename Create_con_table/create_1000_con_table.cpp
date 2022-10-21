#include <vector>
#include <string>
#include <iostream>
#include "../Mpi_Manage_STDP_9_15.cpp"
int main(int argc, char const *argv[])
{

//MPI 测试
    Mpi_Manage ma;
    // cout<<ma.MPI_rank<<" "<<ma.MPI_size<<endl;
    double pop_de_rate=0.1;
    int pop_num=0;
    int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};
    // for (size_t i = 0; i < 8; i++)
    // {
    //     pop_num+=pop_nums[i];
    //     cout<<pop_nums[i]<<"pop_num: "<<pop_num<<endl;
    // }
    
    if(ma.MPI_rank==0)//这些是在一个进程上做
    {
    
    // double syn_weight=3.14;
    // int pop1,pop2;
    // // int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};
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
    // // pop1=ma.create(1,2);//创建神经元（此时没有创建实体，只是构建表格）单进程
    // // pop2=ma.create(1,2);//创建神经元（此时没有创建实体，只是构建表格）单进程
    // // cout<<"pop1 is "<<pop1<<" pop2 is "<<pop2<<endl;
    // // ma.full_connect(pop1,pop2,100);//创建连接（此时没有创建实体，只是构建表格）单进程
    // int L23E=ma.create(1,pop_nums[0]);
    // int L23I=ma.create(1,pop_nums[1]);
    // int L4E=ma.create(1,pop_nums[2]);
    // int L4I=ma.create(1,pop_nums[3]);
    // int L5E=ma.create(1,pop_nums[4]);
    // int L5I=ma.create(1,pop_nums[5]);
    // int L6E=ma.create(1,pop_nums[6]);
    // int L6I=ma.create(1,pop_nums[7]);
    // int populatins[8]={L23E,L23I,L4E,L4I,L5E,L5I,L6E,L6I};
    // double mu[8][8]={
    //     {277,-1110,555,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110},
    //     {277,-1110,277,-1110,277,-1110,277,-1110}
    // };
    // double sigma[8][8]={
    //     {27,111,55,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111},
    //     {27,111,27,111,27,111,27,111}
    // };

    // int conn_weight=0;
    // for (int i = 0; i < 8; i++)
    // {
    //     for (int j = 0; j < 8; j++)
    //     {
    //         conn_weight=gaussrand()*sigma[i][j]*syn_weight+mu[i][j]*syn_weight;

    //         ma.FIX_random_connect_bad_4(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);            
    //         // if(conn_weight>0)
    //         // ma.STDP_random_connect(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
    //         // else
    //         // ma.random_connect(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
    //     }
        
    // }
    // ma.neuron_assignment(1);//根据进程以及神经元分配方案进行神经元分配，单进程
    }
    // ma.prepare_msg();
    if(ma.MPI_rank==0){
    // ma.record_Gid_to_rank();
    // ma.record_adjtable();
    // ma.record_STDP_adj_table();
    // ma.record_weight_table();
    // ma.record_STDP_w_table();
    }
    ma.prepare_read_msg();
    //possion
    double k_ex[8]={1600, 1500, 2100, 1900, 2000, 1900, 2900, 2100};
    double g=8.;
    double Poisson_weight=277;
    double possion_fixed=0.0001;
    double possion_weight=878;
    int Gid_index=0;
    for (int i = 0; i < 8; i++)
    {
        ma.Poisson_full_connect(k_ex[i]*g*0.05,278,Gid_index,Gid_index+pop_nums[i]);
        Gid_index=Gid_index+pop_nums[i];
    }
    
    
    // ma.Open_rate_record();
    // ma.print_n_rate();
    //possion
    struct timeval t1,t2;
    gettimeofday(&t1, nullptr);
    ma.simulate(60000);
    gettimeofday(&t2, nullptr);

    //ma.print_n_rate();
    double usetime_s = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000.0;
    cout<<"rank: "<<ma.MPI_rank<<" use time: "<<usetime_s<<endl;
    cout<<"rank: "<<ma.MPI_rank<<" STDP_count: "<<ma.STDP_count<<endl;        
    cout<<"N : "<<ma.neuron_set[1].global_num<<" spike_times : "<<ma.neuron_set[1].spike_time<<endl;
    cout<<"N : "<<ma.neuron_set[100].global_num<<" spike_times : "<<ma.neuron_set[100].spike_time<<endl;
    // for (size_t i = 0; i < 10; i++)
    // {
    //     cout<<ma.neuron_set[0].Recording_Vm[i]<<" ";
    // }
    // cout<<endl;

    // for (size_t i = 0; i < 10; i++)
    // {
    //     cout<<ma.neuron_set[0].Recording_Vm[i]<<" ";
    // }
    cout<<endl;



}

