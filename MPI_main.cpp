#include <vector>
#include <string>
#include <iostream>
#include "Mpi_Manage_9_27.cpp"
int main(int argc, char const *argv[])
{

//MPI 测试
    Mpi_Manage ma;
    // cout<<ma.MPI_rank<<" "<<ma.MPI_size<<endl;
    double pop_de_rate=1;
    int pop_num=0;
    int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};

    //string G_t_r_txt="part_4_record/hyper2_noweight.txt.part.4.txt";
    // string G_t_r_txt="part_4_record/hyper2.txt.part.4.txt";
    string G_t_r_txt="part_4_record/Gid_to_rank.txt";
    string adj_table_txt="part_4_record/adj.txt";
    string adj_w_table_txt="part_4_record/w.txt";
    string STDP_adj_table_txt="part_4_record/STDP_adj.txt";
    string STDP_adj_w_table_txt="part_4_record/STDP_w.txt";

    if(ma.MPI_size==8)
    {
    cout<<"MPI size is 8!8!"<<endl;
    
    //G_t_r_txt="part_8_record/hyper2_noweight.txt.part.8.txt";
    //G_t_r_txt="part_8_record/hyper2.txt.part.8.txt";
    G_t_r_txt="part_8_record/Gid_to_rank.txt";
    //string G_t_r_txt="part_8_record/Gid_to_rank.txt";
    adj_table_txt="part_8_record/adj.txt";
    adj_w_table_txt="part_8_record/w.txt";
    STDP_adj_table_txt="part_8_record/STDP_adj.txt";
    STDP_adj_w_table_txt="part_8_record/STDP_w.txt";        
    }
    // //string G_t_r_txt="922_record/Gid_to_rank922.txt";
    // string G_t_r_txt="922_record/hyper2.txt.part.8.txt";
    // //string G_t_r_txt="922_record/hyper2_noweight.txt.part.8.txt";
    // //string G_t_r_txt="922_record/hyper2.txt.part5.txt";
    // string adj_table_txt="922_record/adj_922.txt";
    // string adj_w_table_txt="922_record/w.txt";
    // string STDP_adj_table_txt="922_record/STDP_adj.txt";
    // string STDP_adj_w_table_txt="922_record/STDP_w.txt";

    // string G_t_r_txt="hyper2_noweight.txt";
    // string adj_table_txt="adj_table.txt";
    // string adj_w_table_txt="weight_table.txt";
    // string STDP_adj_table_txt="STDP_adj_table.txt";
    // string STDP_adj_w_table_txt="STDP_w_table.txt";
    cout<<"start_prepare_"<<endl;
    ma.prepare_read_msg(G_t_r_txt,adj_table_txt,adj_w_table_txt,STDP_adj_table_txt,STDP_adj_w_table_txt);
    cout<<"end_prepare_"<<endl;
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
    int sime_time=500;
    struct timeval t1,t2;
    gettimeofday(&t1, nullptr);
    ma.simulate(sime_time);
    gettimeofday(&t2, nullptr);

    //ma.print_n_rate();
    double usetime_s = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000.0;
    cout<<"rank: "<<ma.MPI_rank<<" use time: "<<usetime_s<<endl;
    cout<<"rank: "<<ma.MPI_rank<<" STDP_count: "<<ma.STDP_count<<" STDP_t/count is : "<<ma.STDP_time_/ma.STDP_count<<endl;        
    cout<<"rank: "<<ma.MPI_rank<<" neuron_N is "<<ma.neuron_set.size()<<" per_n_ms is : "<<ma.update_time_/ma.neuron_set.size()<<endl;
    cout<<"N : "<<ma.neuron_set[1].global_num<<" spike_times : "<<ma.neuron_set[1].spike_time<<endl;

    cout<<endl;



}
