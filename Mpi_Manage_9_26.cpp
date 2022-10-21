#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <ctime>
#include "Event.cpp"
#include <mpi.h>
#include "LIFNeuron.cpp"
#include "Neuron_assignment.cpp"
#include "MPI_ass_fun.cpp"
#include <cmath>
#include<deque>
//#include "Histentry.h"
#include "STDP.h"
using namespace std;
#define random(x) rand()%(x)
#define RING_BUFFER_LEN 10//ring buffer的长度，我们假设突触延迟都相同，都是这个（ms）.多个文件需要同时修改！ (多了可以放入parameter.h)
#define PUT_BUF_LEN 60000//MPI通信的缓冲区
#define RECV_EVENT_LEN 60000//MPI接收脉冲缓冲区
#define SEND_E_LEN 2 //每个event给别的rank发送的信息量
#define num_len_G_t_r 9000 //神经元数量
#define adjtable_len 2000000 //突触连接数量
//#define 

#define Event_debug 0//
#define get_Event_debug 0
#define static_debug 0
#define send_self_debug 0
#define send_out_debug 0
#define debug_STDP 0

class Mpi_Manage
{
    public:
        //record time
        double get_handle_time_;
        double handle_STDP_time_;
        double STDP_time_;
        double update_time_;
        double send_event_time_;
        double mini_send_time_;
        int STDP_count=0;
        //record time

        //MPI
        // MPI_Status status;
        int MPI_size,MPI_rank;
        int *target_rank_num;
        vector<vector<int>> neuron_assi_table;//第一维是进程,第二维是该进程的神经元G_id
        vector<LIFNeuron> neuron_set;//存放当前进程的神经元
        vector<STDP> STDP_set;//存储STDP突触
        map<int,int> G_id_to_N_set_index;//为多线程准备，G_id标映射到neuron_set的下标
        //double Put_buf[PUT_BUF_LEN];
        map<int,int> Gid_to_rank;//神经元到进程的映射表。
        map<int,map<int,int>> source_target_to_STDP_index;//source神经元到target神经元的STDP_set的下标

        MPI_Win adj_win;//共享内存窗口adj
        MPI_Win weight_win;
        MPI_Win stdp_adj_win;//共享内存窗口adj
        MPI_Win stdp_w_win;
        MPI_Win Gid_rank_map_win;
        //recv_from_rank_win
        vector<MPI_Win> r_f_win;

        vector<vector<int>> send_e_buf;//发送脉冲的buf 第一维是rank，第二维是发送的信息
        int recv_neuron_as[PUT_BUF_LEN];//接收当前进程的神经元连接表
        int recv_Gid_to_rank[PUT_BUF_LEN]={0};//接收当前进程的Gid_to_rank表
        int send_table[PUT_BUF_LEN]={0};
        double double_send_table[PUT_BUF_LEN]={0};
        int recv_adj_table[PUT_BUF_LEN]={0};//接收神经元连接表，一个是数字是神经元数量N（vector行数），接下来N个数字是每行的数量，然后是adjtable的平铺
        double recv_weight_table[PUT_BUF_LEN]={0};//接收神经元连接表，一个是数字是神经元数量N（vector行数），接下来N个数字是每行的数量
        int recv_stdp_adj_table[PUT_BUF_LEN]={0};//接收神经元连接表，一个是数字是神经元数量N（vector行数），接下来N个数字是每行的数量
        double recv_stdp_w_table[PUT_BUF_LEN]={0};//接收神经元连接表，一个是数字是神经元数量N（vector行数），接下来N个数字是每行的数量

        //send_to_rank
        int **send_to_rank;//二维数组，第一维度是rank，第二维是送到该rank的数据

        //send_to_rank

        //recv_from_rank
        int **recv_from_rank;
        //recv_from_rank
        //MPI

        int debug=1;//测试使用，可以删除
        int Record=1;//是否记录Vm，关闭可加快计算.打开(1),关闭(0)
        int Used_Neuron_assignment=0;//一次只能使用一种神经元分配算法，使用过之后变成1

        int print_vm_switch=0;//是否输出Vm，打开(1),关闭(0)，需要输入要打印的神经元的G_id
        int print_Vm_G_id=0;//需要输出Vm的神经元的Gid；
        
        int global_num=0;//初始化神经元编号,神经元编号从0开始
        
        vector<Event> Event_set;//存放脉冲事件
        map<int,int> pop_map;//第一个int表示神经元种群的起始位置(开始的global_num),第二个是群落中神经元的数量

        vector<vector<int>> adjacency_table;//第一维是神经元的编号，第二维是该神经元的突触后神经元
        vector<vector<double>> weight_table;//第一维是神经元编号，第二维度是连接到该神经元的突触后神经元的权重(shape跟adjacency_table一致)

        vector<vector<int>> STDP_adj_table;//第一维是神经元的编号，第二维是该神经元的突触后神经元
        vector<vector<double>> STDP_w_table;//第一维是神经元编号，第二维度是连接到该神经元的突触后神经元的权重(shape跟adjacency_table一致)
        Mpi_Manage()
        {

            srand((int)time(0));//产生随机种子
            MPI_Init(NULL, NULL);
            MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);//the size of mpiranks;
            MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);//the rank of this thread;
            init_send_to_rank();
            // target_rank_num=new int[MPI_size];
            //cout<<"recv_from_rank[i][0] is !!!"<<recv_from_rank[0][0]<<endl;
            MPI_Win_create(&recv_adj_table[0],sizeof(int)*PUT_BUF_LEN,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&adj_win);//创建内存共享窗口
            MPI_Win_create(&recv_weight_table[0],sizeof(double)*PUT_BUF_LEN,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&weight_win);//创建内存共享窗口
            MPI_Win_create(&recv_stdp_adj_table[0],sizeof(int)*PUT_BUF_LEN,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&stdp_adj_win);//创建内存共享窗口
            MPI_Win_create(&recv_stdp_w_table[0],sizeof(double)*PUT_BUF_LEN,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&stdp_w_win);//创建内存共享窗口
            MPI_Win_create(&recv_Gid_to_rank[0],sizeof(int)*PUT_BUF_LEN,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&Gid_rank_map_win);//创建内存共享窗口
      
            
            init_rank_win();
            for (size_t i = 0; i < MPI_size; i++)
            {
                MPI_Win_create(&recv_from_rank[i][0] ,sizeof(int)*RECV_EVENT_LEN,sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&(r_f_win[i]));
            }
            
         
            
        }
        ~Mpi_Manage()
        {
            MPI_Win_free(&adj_win);
            MPI_Win_free(&weight_win);
            MPI_Win_free(&stdp_adj_win);
            MPI_Win_free(&stdp_w_win);
            MPI_Win_free(&Gid_rank_map_win);
       
            for (size_t i = 0; i < MPI_size; i++)
            {
                MPI_Win_free(&(r_f_win[i]));
            }
            
          
            MPI_Finalize();
        }
        void init_send_to_rank()
        {
            send_to_rank=new int *[MPI_size];
            recv_from_rank=new int *[MPI_size];
            target_rank_num=new int[MPI_size];
            for (size_t i = 0; i < MPI_size; i++)
            {
                send_to_rank[i]=new int[RECV_EVENT_LEN];
                send_to_rank[i][0]=-1;
                
                recv_from_rank[i]=new int[RECV_EVENT_LEN];
                recv_from_rank[i][0]=-1;
            }
            // for (size_t i = 0; i < MPI_size; i++)
            // {
            //     for (size_t j = 0; j < RECV_EVENT_LEN; j++)
            //     {
            //         cout<<send_to_rank[i][j]<<" ";
            //     }
            //     cout<<endl;
            // }
            
        }

        void init_rank_win()
        {
            for (size_t i = 0; i < MPI_size; i++)
            {
                MPI_Win win;
                r_f_win.push_back(win);
            }
            

        }
        
        int create(int neuron_type=1,int neuron_num=1)//目前之支持1号神经元类型，第二个为神经元数量，此函数会创建一个神经元种群
        {
            int start_num=global_num;
            pop_map[start_num]=neuron_num;
            if(neuron_type==1)
            {
                for (size_t i = 0; i < neuron_num; i++)
                {
                    global_num++;
                    // LIFNeuron lif(global_num++);
                    // lif.Set_Record(Record);
                    // neuron_set.push_back(lif);

                    vector<int> adj_tab={-1};//往adjacency_table中添加一个vector，占个位置，为之后的connect作准备；
                    vector<double> wei_tab={-1};//同上，只不过是权重矩阵
                    adjacency_table.push_back(adj_tab);//占位
                    weight_table.push_back(wei_tab);//占位
                    STDP_adj_table.push_back(adj_tab);//占位
                    STDP_w_table.push_back(wei_tab);//占位
                }
                
            }
            else
                cout<<"laozi only support neuron_type==1!"<<endl;
            return start_num;//返回一个creat创建的第一个神经元的编号，结合pop_map就能查询群落 
        }
        
        void Set_Vth(int pop, double Vth)//给pop种群的神经元设置Vth的阈值
        {
            for (size_t i = 0; i < pop_map[pop]; i++)
            {
                neuron_set[pop+i].Set_Vth(Vth);
            }
            
        }

        void full_connect(int left_pop,int right_pop,double weight)//全连接
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        adjacency_table[left_pop_point].push_back(right_pop_point);
                        weight_table[left_pop_point].push_back(weight); 
                    }
                }
                
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }

        }


        void random_connect(int left_pop,int right_pop,double weight,double rate)//全连接
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        double random_num=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                        if(rate>=random_num)//按照rate的大小进行概率连接
                        {
                        adjacency_table[left_pop_point].push_back(right_pop_point);
                        weight_table[left_pop_point].push_back(weight); 
                        }
                    }
                    
                }
                
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }

        }


        void Poisson_full_connect(double poisson_rate,double weight,int Left,int Right)//输入需要的poisson发射率以及需要激活的神经元的Gid范围[Left,Right)
        {
            //int start_neu_num=Right_pop;
            for (size_t i = 0; i < neuron_set.size(); i++)
            {
                if (neuron_set[i].Get_globel_num()>=Left&&neuron_set[i].Get_globel_num()<Right)
                {
                    neuron_set[i].Open_Poisson(poisson_rate,weight);//开启该神经元的自发生Poisson开关，现在他可以自己给自己poisson_rate大小的poisson激励啦！（这么设计是为了分布式考虑）
                }            
            }    
        }


        void STDP_full_connect(int left_pop,int right_pop,double weight)
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        STDP_adj_table[left_pop_point].push_back(right_pop_point);
                        STDP_w_table[left_pop_point].push_back(weight); 
                    }
                }
                
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }
        }


        void STDP_random_connect(int left_pop,int right_pop,double weight,double rate)//全连接
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        double random_num=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                        if(rate>=random_num)//按照rate的大小进行概率连接
                        {
                        STDP_adj_table[left_pop_point].push_back(right_pop_point);
                        STDP_w_table[left_pop_point].push_back(weight); 
                        }
                    }
                }  
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }
        }

        void FIX_random_connect(int left_pop,int right_pop,double weight,double rate)//全连接
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        double random_num=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                        if(rate>=random_num)//按照rate的大小进行概率连接
                        {
                            if(weight>0)
                            {
                                double random_num_stdp=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                                if(0.7>=random_num_stdp)
                                {
                                    STDP_adj_table[left_pop_point].push_back(right_pop_point);
                                    STDP_w_table[left_pop_point].push_back(weight);  
                                }
                                else
                                {
                                    adjacency_table[left_pop_point].push_back(right_pop_point);
                                    weight_table[left_pop_point].push_back(weight);                                    
                                }
                            }
                            else
                            {
                                adjacency_table[left_pop_point].push_back(right_pop_point);
                                weight_table[left_pop_point].push_back(weight);    
                            }
                        //STDP_adj_table[left_pop_point].push_back(right_pop_point);
                        //STDP_w_table[left_pop_point].push_back(weight); 
                        }
                    }
                }  
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }
        }

        void FIX_random_connect_bad_4(int left_pop,int right_pop,double weight,double rate)//全连接
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        double random_num=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                        if(rate>=random_num)//按照rate的大小进行概率连接
                        {
                            if(weight>0)
                            {
                                //double random_num_stdp=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                                if(right_pop_point%4==0)
                                {
                                    STDP_adj_table[left_pop_point].push_back(right_pop_point);
                                    STDP_w_table[left_pop_point].push_back(weight);  
                                }
                                else
                                {
                                    adjacency_table[left_pop_point].push_back(right_pop_point);
                                    weight_table[left_pop_point].push_back(weight);                                    
                                }
                            }
                            else
                            {
                                adjacency_table[left_pop_point].push_back(right_pop_point);
                                weight_table[left_pop_point].push_back(weight);    
                            }
                        //STDP_adj_table[left_pop_point].push_back(right_pop_point);
                        //STDP_w_table[left_pop_point].push_back(weight); 
                        }
                    }
                }  
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }
        }

        void FIX_random_connect_STDP(int left_pop,int right_pop,double weight,double rate)//全连接
        {
            //全连接
            if(pop_map.count(left_pop)==1 && pop_map.count(right_pop)==1)//左右种群都存在
            {
                for (int left_pop_point = left_pop; left_pop_point < left_pop+pop_map[left_pop]; left_pop_point++)
                {
                    for (int right_pop_point = right_pop; right_pop_point < right_pop+pop_map[right_pop]; right_pop_point++)
                    {
                        double random_num=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                        if(rate>=random_num)//按照rate的大小进行概率连接
                        {
                            if(weight>0)
                            {
                                double random_num_stdp=double(random(1000))/double(1000);//生成一个小数点后三位的小数
                                if(0.7>random_num_stdp)
                                {
                                    STDP_adj_table[left_pop_point].push_back(right_pop_point);
                                    STDP_w_table[left_pop_point].push_back(weight);  
                                }
                                else
                                {
                                    adjacency_table[left_pop_point].push_back(right_pop_point);
                                    weight_table[left_pop_point].push_back(weight);                                    
                                }
                            }
                            else
                            {
                                adjacency_table[left_pop_point].push_back(right_pop_point);
                                weight_table[left_pop_point].push_back(weight);    
                            }
                        //STDP_adj_table[left_pop_point].push_back(right_pop_point);
                        //STDP_w_table[left_pop_point].push_back(weight); 
                        }
                    }
                }  
            }
            else
            {
                cout<<"left_pop or right_pop not exist!"<<endl;
            }
        }



        //1是轮询，2是只根据神经元权重图划分，3是考虑突触图划分
        void neuron_assignment(int assignment_type)//1是轮询，2是只根据神经元权重图划分，3是考虑突触图划分
        {
            Neuron_assignment N_ass;
            if(assignment_type==1)
            {
                neuron_assi_table=N_ass.rotation_assignment(MPI_size,adjacency_table,STDP_adj_table);
            }
            else
            {
                cout<<"assignment_type not exsist !"<<endl;
            }    
            for (size_t i = 0; i < MPI_size ; i++)
            {
                for (size_t j = 0; j < neuron_assi_table[i].size(); j++)
                {
                    Gid_to_rank[neuron_assi_table[i][j]]=i;
                }
            }
            //cout<<"in send adj_table"<<endl;
                  
        }

        void send_Gid_to_rank()
        {
            MPI_Win_fence(0, Gid_rank_map_win);
            if(MPI_rank==0)
            {
            vector<int> send_Gid_to_rank;
            send_Gid_to_rank.push_back(-1);
            map<int,int>::iterator it;
            for (it=Gid_to_rank.begin(); it!=Gid_to_rank.end(); it++)
            {
                send_Gid_to_rank.push_back(it->first);
                send_Gid_to_rank.push_back(it->second);
            }
            send_Gid_to_rank[0]=send_Gid_to_rank.size();
            vector_to_int_set(send_Gid_to_rank,send_table);
            int send_count=send_Gid_to_rank.size();
            cout<<"send_count is "<<send_count<<endl;
            for (int i = 1; i < MPI_size; i++)
            {
                MPI_Put(&send_table[0],send_count , MPI_INT, i,0, send_count, MPI_INT, Gid_rank_map_win); 

            }
            }
            MPI_Win_fence(0, Gid_rank_map_win);
            
            if(MPI_rank!=0)
            {
                int neuron_rec_num=recv_Gid_to_rank[0];
                vector<int> rec_v=int_set_to_vector(&recv_Gid_to_rank[0],neuron_rec_num);
                map<int,int> rec_map;
                for (size_t i = 1; i < neuron_rec_num; i+=2)
                {
                    rec_map[rec_v[i]]=rec_v[i+1];
                }
                Gid_to_rank=rec_map;
            //map<int,int>::iterator it;
            // for (it=Gid_to_rank.begin(); it!=Gid_to_rank.end(); it++)
            // {
            //     cout<<"it->first is "<<it->first;
            //     cout<<" it->second is "<<it->second<<endl;
            // }
            }
        }

        void send_adj_table()//把表从0号进程传送去其他进程，进行连接表同步
        {

            MPI_Win_fence(0, adj_win);
            if(MPI_rank==0)
            {
                vector<int> send_table_v=create_send_adj_table(adjacency_table);
                vector_to_int_set(send_table_v,send_table);
                
                int send_count=send_table_v.size();
                cout<<"send_count is "<<send_count<<endl;
                for (int i = 1; i < MPI_size; i++)
                {
                   MPI_Put(&send_table[0],send_count , MPI_INT, i,0, send_count, MPI_INT, adj_win); 

                }
                
            }
            MPI_Win_fence(0, adj_win);  
            if(MPI_rank!=0)
            {
                int neuron_rec_num=recv_adj_table[0];
                int rec_count=0;
                for (int i = 1; i < neuron_rec_num+1; ++i)
                    {
                        rec_count+=recv_adj_table[i];
                    }
                    rec_count+=1;
                    rec_count+=recv_adj_table[0];    
                cout<<"this is rank "<<MPI_rank<<" rec_count is "<<rec_count<<endl;
                vector<int> rec_v=int_set_to_vector(&recv_adj_table[0],rec_count);
                adjacency_table=one_v_to_2_v(rec_v);
            }

        }

        void send_weight_table()//把表从一号进程传送去其他进程，进行连接表同步
        {

            MPI_Win_fence(0, weight_win);
            if(MPI_rank==0)
            {
                vector<double> send_table_v=create_send_adj_table(weight_table);
                vector_to_int_set(send_table_v,double_send_table);
                
                int send_count=send_table_v.size();
                cout<<"weight send_count is "<<send_count<<endl;
                for (int i = 1; i < MPI_size; i++)
                {
                   MPI_Put(&double_send_table[0],send_count , MPI_DOUBLE, i,0, send_count, MPI_DOUBLE, weight_win); 

                }
                
            }
            MPI_Win_fence(0, weight_win);  
            if(MPI_rank!=0)
            {
                int neuron_rec_num=recv_weight_table[0];
                int rec_count=0;
                for (int i = 1; i < neuron_rec_num+1; ++i)
                    {
                        rec_count+=recv_weight_table[i];
                    }
                    rec_count+=1;
                    rec_count+=recv_weight_table[0];    
                cout<<"this is rank "<<MPI_rank<<" rec_count is "<<rec_count<<endl;
                vector<double> rec_v=int_set_to_vector(&recv_weight_table[0],rec_count);
                weight_table=one_v_to_2_v(rec_v);
            }

        }

        void send_stdp_adj_table()//把表从一号进程传送去其他进程，进行连接表同步
        {

            MPI_Win_fence(0, stdp_adj_win);
            if(MPI_rank==0)
            {
                vector<int> send_table_v=create_send_adj_table(STDP_adj_table);
                vector_to_int_set(send_table_v,send_table);
                
                int send_count=send_table_v.size();
                cout<<"send_count is "<<send_count<<endl;
                for (int i = 1; i < MPI_size; i++)
                {
                   MPI_Put(&send_table[0],send_count , MPI_INT, i,0, send_count, MPI_INT, stdp_adj_win); 

                }
                
            }
            MPI_Win_fence(0, stdp_adj_win);  
            if(MPI_rank!=0)
            {
                int neuron_rec_num=recv_stdp_adj_table[0];
                int rec_count=0;
                for (int i = 1; i < neuron_rec_num+1; ++i)
                    {
                        rec_count+=recv_stdp_adj_table[i];
                    }
                    rec_count+=1;
                    rec_count+=recv_stdp_adj_table[0];    
                cout<<"this is rank "<<MPI_rank<<" rec_count is "<<rec_count<<endl;
                vector<int> rec_v=int_set_to_vector(&recv_stdp_adj_table[0],rec_count);
                STDP_adj_table=one_v_to_2_v(rec_v);
            }

        }

        void send_stdp_w_table()//把表从一号进程传送去其他进程，进行连接表同步
        {

            MPI_Win_fence(0, stdp_w_win);
            if(MPI_rank==0)
            {
                vector<double> send_table_v=create_send_adj_table(STDP_w_table);
                vector_to_int_set(send_table_v,double_send_table);
                
                int send_count=send_table_v.size();
                cout<<"weight send_count is "<<send_count<<endl;
                for (int i = 1; i < MPI_size; i++)
                {
                   MPI_Put(&double_send_table[0],send_count , MPI_DOUBLE, i,0, send_count, MPI_DOUBLE, stdp_w_win); 

                }
                
            }
            MPI_Win_fence(0, stdp_w_win);  
            if(MPI_rank!=0)
            {
                int neuron_rec_num=recv_stdp_w_table[0];
                int rec_count=0;
                for (int i = 1; i < neuron_rec_num+1; ++i)
                    {
                        rec_count+=recv_stdp_w_table[i];
                    }
                    rec_count+=1;
                    rec_count+=recv_stdp_w_table[0];    
                cout<<"this is rank "<<MPI_rank<<" rec_count is "<<rec_count<<endl;
                vector<double> rec_v=int_set_to_vector(&recv_stdp_w_table[0],rec_count);
                STDP_w_table=one_v_to_2_v(rec_v);
            }

        }

        void record_Gid_to_rank(string Gid_to_r_txt)
        {
            int neuron_num=0;
            vector<int> send_Gid_to_rank;
            send_Gid_to_rank.push_back(-1);
            map<int,int>::iterator it;
            for (it=Gid_to_rank.begin(); it!=Gid_to_rank.end(); it++)
            {
                neuron_num++;
                //send_Gid_to_rank.push_back(it->first);
                send_Gid_to_rank.push_back(it->second);
            }
            send_Gid_to_rank[0]=send_Gid_to_rank.size();
            vector_to_int_set(send_Gid_to_rank,send_table);
            int send_count=send_Gid_to_rank.size();

            ofstream ofs;
            ofs.open(Gid_to_r_txt,ios::out);
            for (size_t i = 0; i < send_Gid_to_rank.size(); i++)
            {
                    ofs<<send_Gid_to_rank[i]<<" ";
            }

            cout<<"neuron_num is "<<neuron_num<<endl;
        }

        void record_adjtable(string adj_table_txt)//记录adjtable到txt
        {
            int adj_count=0;
            ofstream ofs;
            ofs.open(adj_table_txt,ios::out);
            for (size_t i = 0; i < adjacency_table.size(); i++)
            {
                for (size_t j = 0; j < adjacency_table[i].size(); j++)
                {
                    ofs<<adjacency_table[i][j]<<" ";
                }
                ofs<<endl;
                adj_count+=adjacency_table[i].size()-1;
                
            }
            cout<<"static adj_count is "<<adj_count<<endl;
            
        }

        void record_weight_table(string weight_table_txt)//记录adjtable到txt
        {
            ofstream ofs;
            ofs.open(weight_table_txt,ios::out);
            for (size_t i = 0; i < weight_table.size(); i++)
            {
                for (size_t j = 0; j < weight_table[i].size(); j++)
                {
                    ofs<<weight_table[i][j]<<" ";
                }
                ofs<<endl;
                
            }
            
        }


        void record_STDP_adj_table(string STDP_adj_table_txt)//记录adjtable到txt
        {
            int adj_count=0;
            ofstream ofs;
            ofs.open(STDP_adj_table_txt,ios::out);
            for (size_t i = 0; i < STDP_adj_table.size(); i++)
            {
                for (size_t j = 0; j < STDP_adj_table[i].size(); j++)
                {
                    ofs<<STDP_adj_table[i][j]<<" ";
                }
                ofs<<endl;
                adj_count+=STDP_adj_table[i].size()-1;
            }
            cout<<"STDP adj_count is "<<adj_count<<endl;
            
        }

        void record_STDP_w_table(string STDP_w_table_txt)//记录adjtable到txt
        {
            ofstream ofs;
            ofs.open(STDP_w_table_txt,ios::out);
            for (size_t i = 0; i < STDP_w_table.size(); i++)
            {
                for (size_t j = 0; j < STDP_w_table[i].size(); j++)
                {
                    ofs<<STDP_w_table[i][j]<<" ";
                }
                ofs<<endl;
                
            }
            
        }


    void read_Gid_to_rank(string G_t_r_read_txt)
    {
        int i,datalen=0;
        int num[num_len_G_t_r];
        //ifstream file("hyper2_part_4_10.txt");
        ifstream file(G_t_r_read_txt);
        // ifstream file("hyper2_part4.txt");
        // ifstream file("part4_hypergraph.txt");
        while( ! file.eof() ){
        file>>num[datalen++];
        //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
        }
        file.close();

        int neuron_rec_num=num[0];
        vector<int> rec_v=int_set_to_vector(&num[0],neuron_rec_num);
        cout<<"read_Gid_to_rank_rec_v size : "<<rec_v.size()<<endl;
        map<int,int> temp_Gid_to_rank;
        for (size_t i = 1; i < rec_v.size(); i++)
        {
            temp_Gid_to_rank[i-1]=rec_v[i];
        }
        Gid_to_rank=temp_Gid_to_rank;    
    }

        void read_adj_table(string adj_table_read_txt)
        {
            vector<vector<int>> temp_adj_table;
            int line_count=0;
            int line_index=0;

            cout<<"read_adj_table"<<endl;
            int datalen=0;
            int i=0;
            int num;
            ifstream file(adj_table_read_txt,ios::in);
            while( ! file.eof() ){
                //cout<<"read_line : "<<i<<endl;
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
                datalen++;
            //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
            }
            file.close();
            cout<<"adj_table_datalen is: "<<datalen<<endl;
            adjacency_table=temp_adj_table;         
        }

        void read_weight_table(string w_table_read_txt)
        {
            vector<vector<double>> temp_weight_table;
            int line_count=0;
            int line_index=0;
            cout<<"read_w_table"<<endl;
            //cout<<"read_in"<<endl;
            int datalen=0;
            int i=0;
            double num;
            ifstream file(w_table_read_txt,ios::in);
            while( ! file.eof() ){
            file>>num;
                if (num==-1)
                {
                    vector<double> adj;
                    temp_weight_table.push_back(adj);
                    line_count++;
                }
                if (i!=0&&num==-1)
                {
                    line_index++;
                }
                temp_weight_table[line_index].push_back(num); 
                i++;
                datalen++;

            //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
            }
            file.close();
            cout<<"read_w_table len is : "<<datalen<<endl;
            weight_table=temp_weight_table;       
        }

        void read_STDP_adj_table(string STDP_adj_table_read_txt)
        {
           vector<vector<int>> temp_STDP_adj_table;
            int line_count=0;
            int line_index=0;

            cout<<"read_adj_table"<<endl;
            int datalen=0;
            int i=0;
            int num;
            ifstream file(STDP_adj_table_read_txt,ios::in);
            while( ! file.eof() ){
                //cout<<"read_line : "<<i<<endl;
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
                datalen++;
            //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
            }
            file.close();
            cout<<"STDP_adj_table_datalen is: "<<datalen<<endl;
            STDP_adj_table=temp_STDP_adj_table;       


      
        }

        void read_STDP_w_table(string STDP_w_table_read_txt)
        {
            vector<vector<double>> temp_STDP_w_table;
            int line_count=0;
            int line_index=0;
            cout<<"read_w_table"<<endl;
            //cout<<"read_in"<<endl;
            int datalen=0;
            int i=0;
            double num;
            ifstream file(STDP_w_table_read_txt,ios::in);
            while( ! file.eof() ){
            file>>num;
                if (num==-1)
                {
                    vector<double> adj;
                    temp_STDP_w_table.push_back(adj);
                    line_count++;
                }
                if (i!=0&&num==-1)
                {
                    line_index++;
                }
                temp_STDP_w_table[line_index].push_back(num); 
                i++;
                datalen++;

            //cout<<"len is :"<<datalen<<" num is : "<<num[datalen-1]<<endl;
            }
            file.close();
            cout<<"read_STDP_w_table len is : "<<datalen<<endl;
            STDP_w_table=temp_STDP_w_table;       


        
        }
        
        void create_this_rank_neuron()//创建当前进程上的神经元
        {
            map<int,int>::iterator it;
            int i=0;
            for (it = Gid_to_rank.begin(); it!=Gid_to_rank.end(); it++)
            {
                if(it->second==MPI_rank)
                {
                LIFNeuron lif(it->first);
                neuron_set.push_back(lif);
                G_id_to_N_set_index[it->first]=i;
                i++;
                }
                //N_set_index_to_G_id[i]=neuron_assi_table[MPI_rank][i];   
            
            }
        }
        

        void create_STDP()
        {
            vector<int> neuron_Gid;//存放当前进程的神经元Gid
            for (size_t i = 0; i < neuron_set.size(); i++)
            {
                neuron_Gid.push_back(neuron_set[i].global_num);
            }
            int stdp_set_index=0;
            for (size_t i = 0; i < STDP_adj_table.size(); i++)
            {
                for (size_t j = 1; j < STDP_adj_table[i].size(); j++)
                {
                    if(search(neuron_Gid,STDP_adj_table[i][j]))//如果STDP突触后神经元在该进程，则创建一个STDP实体
                    {
                        STDP stdp_;
                        stdp_.set_source_Gid(i);
                        stdp_.set_target_Gid(j);
                        stdp_.set_weight(STDP_w_table[i][j]);
                        STDP_set.push_back(stdp_);
                        neuron_set[G_id_to_N_set_index[j]].N_syn+=1;//stdp_indegree

                        source_target_to_STDP_index[i][j]=stdp_set_index;
                        stdp_set_index++;                        
                    
                    }
                }
                
            }
            
        }

        void prepare_msg()
        {
            init_send_event();//对发射事件进行初始化
            send_Gid_to_rank();//发送Gid_to_rank的数据
            send_adj_table();
            send_stdp_adj_table();
            send_weight_table();
            send_stdp_w_table();
            create_this_rank_neuron();
            create_STDP();
        }

        void prepare_read_msg(string read_gid_to_rank_txt,string read_adj_table_txt,string read_w_table_txt,string read_STDP_adj_table_txt,string read_STDP_w_adj_table_txt)
        {
            init_send_event();//对发射事件进行初始化
            read_Gid_to_rank(read_gid_to_rank_txt);
            read_adj_table(read_adj_table_txt);
            //cout<<"read_adj_table"<<endl;
            read_weight_table(read_w_table_txt);
            //cout<<"read_weight_table"<<endl;
            read_STDP_adj_table(read_STDP_adj_table_txt);
            //cout<<"read_STDP_adj_table"<<endl;
            read_STDP_w_table(read_STDP_w_adj_table_txt);
            //cout<<"read_STDP_w_table"<<endl;
            create_this_rank_neuron();
            //cout<<"create_this_rank_neuron"<<endl;
            create_STDP();
            //cout<<"create_STDP"<<endl;
        }


        void simulate(int sim_time)//开始模拟，输入模拟时间
        {
            struct timeval get_handle_event_t1,get_handle_event_t2;
            //struct timeval handle_STDP_t1,handle_STDP_t2;
            struct timeval update_t1,update_t2;
            struct timeval STDP_t1,STDP_t2;
            struct timeval send_event_t1,send_event_t2;
            double get_handle_time=0;
            double handle_STDP_time=0;
            double h_STDP_time=0;
            double send_time=0;
            double mini_send_time=0;
            double STDP_time=0;
            double update_time=0;
            double send_event_time=0;
            vector<vector<int>> merge_adj_table = merge_adj_and_stdp_table();
            // cout<<"the send buf is "<<send_e_buf[0][0]<<send_e_buf[0][1]<<endl;
            int loop_steps=sim_time;
            if(debug)
            {
                cout<<"steps is "<<loop_steps<<endl;
            }
            for (int loop_step = 0; loop_step < loop_steps; loop_step++)
            {
                gettimeofday(&get_handle_event_t1, nullptr);
                get_event();
                //gettimeofday(&get_handle_event_t2, nullptr);
                handle_event(loop_step,&h_STDP_time);//获取跟处理事件
                handle_STDP_time+=h_STDP_time;
                gettimeofday(&get_handle_event_t2, nullptr);
                get_handle_time += (get_handle_event_t2.tv_sec - get_handle_event_t1.tv_sec) + (double)(get_handle_event_t2.tv_usec - get_handle_event_t1.tv_usec) / 1000000.0;

                gettimeofday(&STDP_t1, nullptr);
                update_STDP(loop_step);//进行STDP更新
                gettimeofday(&STDP_t2, nullptr);
                STDP_time += (STDP_t2.tv_sec - STDP_t1.tv_sec) + (double)(STDP_t2.tv_usec - STDP_t1.tv_usec) / 1000000.0;
 

                Clear_event_set();//清除事件

                gettimeofday(&update_t1, nullptr);
                sim_update(loop_step);//更新神经元并且生成事件
                gettimeofday(&update_t2, nullptr);
                update_time += (update_t2.tv_sec - update_t1.tv_sec) + (double)(update_t2.tv_usec - update_t1.tv_usec) / 1000000.0;

                gettimeofday(&send_event_t1, nullptr);
                send_event(merge_adj_table,&send_time);//分发事件
                mini_send_time+=send_time;
                gettimeofday(&send_event_t2, nullptr);
                send_event_time += (send_event_t2.tv_sec - send_event_t1.tv_sec) + (double)(send_event_t2.tv_usec - send_event_t1.tv_usec) / 1000000.0;
                Clear_event_set();//清除事件
                // if(print_vm_switch==1)
                // {
                //     print_vm(print_Vm_G_id);
                // }
            }   
            cout<<"TIME:::!rank : "<<MPI_rank<<" get_handle_time : "<<get_handle_time<<" handle_STDP_time : "<<handle_STDP_time<<" STDP_time : "<<STDP_time<<" update_time : "<<update_time<<" send_event_time : "<<send_event_time<<" mini_send_event_time : "<<mini_send_time<<endl;
            get_handle_time_=get_handle_time;
            handle_STDP_time_=handle_STDP_time;
            STDP_time_=STDP_time;
            update_time_=update_time;
            send_event_time_=send_event_time;
            mini_send_time_=mini_send_time;
        }
        
        void update_STDP(int loop_step)
        {
            for (int i = 0; i < STDP_set.size(); i++)
            {
                for (int j = 0; j < Event_set.size(); j++)
                {
                    if(Event_set[j].get_Source_G_id()==STDP_set[i].source_Gid)
                    {
                        double spike_time=loop_step-1+0.1*Event_set[j].get_mini_step();
                        STDP_set[i].update_weight(spike_time,neuron_set[G_id_to_N_set_index[STDP_set[i].target_Gid]]);
                        STDP_count++;
                        if(debug_STDP)
                        {
                            cout<<"STDP: "<<i<<" weight is "<<STDP_set[i].weight<<" getEvent_time: "<<spike_time<<"+++++++++"<<endl;
                        }
                    }
                }
                
            }
            
        }

        // void prepare()//每一轮模拟开始时，进行上一轮脉冲的发放(发放到对应神经元的Ring_buffer中)
        // {
          
        // }
        
        void get_event()//每一轮模拟开始时，进行上一轮脉冲的发放(发放到对应神经元的Ring_buffer中)
        {
            for (size_t i = 0; i < MPI_size; i++)
            {
                //cout<<"recv_from_rank[i][0] is "<<recv_from_rank[i][0]<<" "<<endl;
                
                for (size_t j = 1; j < recv_from_rank[i][0]+1; j+=2)
                {
                    if(get_Event_debug){
                    cout<<"Get_event:::!!this is rank "<<MPI_rank<<" recv_from_rank["<< i<<"][0] is "<<recv_from_rank[i][0]<<" recv_from_rank["<< i<<"][1] is "<<recv_from_rank[i][1]<<" recv_from_rank["<< i<<"][2] is "<<recv_from_rank[i][2]<<endl<<endl;    
                    }
                    Event e;
                    e.set_Source_G_id(recv_from_rank[i][j]);
                    e.set_mini_step(recv_from_rank[i][j+1]);
                    Event_set.push_back(e);
                    // recv_from_rank[i][0]=-2;
                }
                
            }
            
            //void reset_recv_from_rank();
            for (size_t i = 0; i < MPI_size; i++)//重置recv_from_rank
            {
                recv_from_rank[i][0]=-1;
                // for (size_t j = 0; j < count; j++)
                // {
                //     /* code */
                // }
                
            }
            // if(get_Event_debug){
            // cout<<"Afterreset_recv_from_rank "<<MPI_rank<<" recv_from_rank[0][0] is "<<recv_from_rank[0][0]<<" recv_from_rank[0][1] is "<<recv_from_rank[0][1]<<" recv_from_rank[0][2] is "<<recv_from_rank[0][2]<<endl<<endl;    
            // }
        }

        void reset_recv_from_rank()
        {
            for (size_t i = 0; i < MPI_size; i++)
            {
                recv_from_rank[i][0]=-2;
                // for (size_t j = 0; j < count; j++)
                // {
                //     /* code */
                // }
                
            }
            
        }

        void handle_event(int loop_step,double *h_STDP_time)
        {
            struct timeval handle_STDP_t1,handle_STDP_t2;
            double temp_STDP_time=0;
            for (size_t i = 0; i < Event_set.size(); i++)
            {
                gettimeofday(&handle_STDP_t1, nullptr);
                STDP_syn_update(i,loop_step);
                gettimeofday(&handle_STDP_t2, nullptr);

                Static_syn_update(i,loop_step);
                
                temp_STDP_time+=(handle_STDP_t2.tv_sec - handle_STDP_t1.tv_sec) + (double)(handle_STDP_t2.tv_usec - handle_STDP_t1.tv_usec) / 1000000.0;
            }
            *h_STDP_time=temp_STDP_time;   
        }

        void STDP_syn_update(int i,int loop_step)
        {
            int source_n=Event_set[i].get_Source_G_id();
            int mini_step=Event_set[i].get_mini_step();
            // if(STDP_set.size()>0)
            // {

            // }
            if(STDP_adj_table[source_n].size()>1)//如果这个source神经元后面有STDP连接
            {
                for (size_t STDP_syn_point = 1; STDP_syn_point < STDP_adj_table[source_n].size(); STDP_syn_point++)
                {
                    int target_n=STDP_adj_table[source_n][STDP_syn_point];
                    if(G_id_to_N_set_index.count(target_n)==0)
                    {
                        continue;
                    }
                    //double weight=STDP_w_table[source_n][STDP_syn_point];

                    int stdp_set_index=source_target_to_STDP_index[source_n][target_n];
                    double weight=STDP_set[stdp_set_index].weight;
                    // for (int stdp_index = 0; stdp_index < STDP_set.size(); stdp_index++)//遍历stdp_index，查找对应的STDP
                    // {
                    //     if(STDP_set[stdp_index].source_Gid==source_n && STDP_set[stdp_index].target_Gid==target_n)
                    //     {
                    //         weight=STDP_set[stdp_index].weight;
                    //     }
                    //     // else
                    //     // {
                    //     //     cout<<"STDP_syn_update  NOT mached!"<<endl;
                    //     //     cout<<"STDP_set[stdp_index].source_Gid is :"<<STDP_set[stdp_index].source_Gid<<" != "<<source_n<<" STDP_set[stdp_index].target_Gid is"<<STDP_set[stdp_index].target_Gid<<" != "<<target_n<<endl;
                    //     // }
                    // }
                    
                    // double weight=STDP_w_table[source_n][STDP_syn_point];
                    neuron_set[G_id_to_N_set_index[target_n]].Set_one_Ring_buffer(mini_step,weight,loop_step);
                }
            }            
        }
        
        void Static_syn_update(int i,int loop_step)
        {
            int source_n=Event_set[i].get_Source_G_id();
            int mini_step=Event_set[i].get_mini_step();
            //cout<<"this is rank "<<MPI_rank<<" Event source id is "<<source_n<<"mini_step_is "<<mini_step<<" loop_step is "<<loop_step<<endl;
            if(adjacency_table[source_n].size()>1)//如果这个source神经元后面有连接
            {
                for (size_t Static_syn_point = 1; Static_syn_point < adjacency_table[source_n].size(); Static_syn_point++)
                {
                    int target_n=adjacency_table[source_n][Static_syn_point];
                    if(G_id_to_N_set_index.count(target_n)==0)
                    {
                        continue;
                    }
                    double weight=weight_table[source_n][Static_syn_point];
                    if(static_debug)
                    {
                    cout<<"this is rank "<<MPI_rank<<" G_id_to_N_set_index["<<target_n<<"] is"<<G_id_to_N_set_index[target_n]<<endl;
                    }
                    neuron_set[G_id_to_N_set_index[target_n]].Set_one_Ring_buffer(mini_step,weight,loop_step);
                }
            }            
        }

        void sim_update(int loop_step)//神经元进行一轮模拟。并压入事件
        {
            for (size_t i = 0; i < neuron_set.size(); i++)
            {
                neuron_set[i].updata(loop_step);
                Push_events(neuron_set[i].Recording_spikes,neuron_set[i].global_num,loop_step);
            }          
        }


        void Push_events(int Spikes_buffer[RING_BUFFER_LEN],int source_n,int loop_step)//每一轮更新后，把产生的事件放入Event_buffer
        {
            for (size_t i = 0; i < RING_BUFFER_LEN; i++)
            {
                if(Spikes_buffer[i]==1)
                {
                    if(Event_debug){
                    cout<<"this is rank "<<MPI_rank<<" generate Event source is "<<source_n<<" time: "<<loop_step<<"."<<i<<endl;
                    }
                    Event e;
                    e.set_Source_G_id(source_n);
                    e.set_time(double(loop_step)+i*0.1);
                    e.set_mini_step(i);
                    Event_set.push_back(e);
                }
            }
            
        }

        void send_event(vector<vector<int>> merge_adj_table,double * send_time)
        {  
            //cout<<"before taget"<< endl;
            for (size_t i = 0; i < MPI_size; i++)
            {
                target_rank_num[i]=0;
            } 
            //cout<<"end taget"<< endl;
            

            //cout<<"before eventsize "<<endl;
            for (size_t i = 0; i < Event_set.size(); i++)//遍历所有的事件，并分发
            {
                for (size_t j = 1; j < merge_adj_table[Event_set[i].get_Source_G_id()].size(); j++)//分发有连接的事件
                {
                    target_rank_num[Gid_to_rank[merge_adj_table[Event_set[i].get_Source_G_id()][j]]]+=1;//如果需要给突触后神经元所在的rank发送，则target_tank_num[rank]++
                }
                for (size_t target_rank = 0; target_rank < MPI_size; target_rank++)
                {
                    if (target_rank_num[target_rank]>0)//给相应的buffer添加事件信息，需要传输更多信息可以在此添加
                    {
                        send_e_buf[target_rank].push_back(Event_set[i].get_Source_G_id());
                        send_e_buf[target_rank].push_back(Event_set[i].get_mini_step());
                    }
                    
                }
                for (size_t re_index = 0; re_index < MPI_size; re_index++)//重置target_rank_num
                {
                    target_rank_num[re_index]=0;
                }     
            }
            // cout<<"end eventsize "<<endl;
            // cout<<"before 1 "<<endl;
            put_send_e_buf_to_send_f_r();//把send_e_buf的数据放入send_to_rank中
            //cout<<"end 1 "<<endl;
            reset_send_e_buf();//重置
            //cout<<"start send "<<endl;
            struct timeval s_e_t_r_t1,s_e_t_r_t2;
            double temp_setr_time=0;
            gettimeofday(&s_e_t_r_t1, nullptr);
            send_event_to_rank();//把send_to_rank数据分发到不同rank；
            gettimeofday(&s_e_t_r_t2, nullptr);
            temp_setr_time+=(s_e_t_r_t2.tv_sec - s_e_t_r_t1.tv_sec) + (double)(s_e_t_r_t2.tv_usec - s_e_t_r_t1.tv_usec) / 1000000.0;
            *send_time=temp_setr_time;
            //cout<<"end send "<<endl;
            reset_send_to_rank();//重置

        }

        vector<vector<int>> merge_adj_and_stdp_table()
        {
            vector<vector<int>> res;
            if (adjacency_table.size()!=STDP_adj_table.size())
            {
                std::cout<<"adjacency_table.size()!=STDP_adj_table.size() have some error in table create!"<<std::endl;
            }
            for (int i = 0; i < adjacency_table.size(); ++i)
            {
                vector<int> adj={-1};
                res.push_back(adj);
            }
            for (int i = 0; i < adjacency_table.size(); ++i)
            {
                for (int j = 1; j < adjacency_table[i].size(); ++j)
                {
                    res[i].push_back(adjacency_table[i][j]);
                }   
                for (int k = 1; k < STDP_adj_table[i].size(); ++k)
                {
                    res[i].push_back(STDP_adj_table[i][k]);    
                }
            }

            for (int i = 0; i < adjacency_table.size(); ++i)
            {
                sort(res[i].begin(),res[i].end());
                res[i].erase(unique(res[i].begin(),res[i].end()),res[i].end());
            }


            return res;
        }

        void put_send_e_buf_to_send_f_r()
        {
            //cout<<"send_e_buf.size() is "<<send_e_buf.size()<<endl;
            //cout<<"MPI_size is "<<MPI_size<<endl;
            for (size_t i = 0; i < send_e_buf.size(); i++)
            {  
                // if(MPI_rank==i)//如果要发送给自己，就跳过
                // {
                //     continue;
                // }
                send_to_rank[i][0]=send_e_buf[i].size()-1;//第一位是要发送到该rank的长度（减去第一个“-1”）

                if(send_e_buf[i].size()>RECV_EVENT_LEN)
                {
                    cout<<"send_e_buf[i].size()>RECV_EVENT_LEN!! please extended quantification of RECV_EVENT_LEN"<<endl;
                    exit(3);
                }

                for (size_t j = 1; j < send_e_buf[i].size(); j++)
                {   
                    send_to_rank[i][j]=send_e_buf[i][j];
                }
                
            }
               
        }

        void init_send_event()//重置
        {
            send_e_buf.clear();
            // cout<<"init_send_event send_e_buf is "<<send_e_buf.size()<<endl; 
            for (size_t i = 0; i < MPI_size; i++)
            {
                vector<int> a={-1};
                send_e_buf.push_back(a);
            }
            //cout<<"init_send_event send_e_buf is "<<send_e_buf.size()<<endl;            
        }

        void reset_send_e_buf()//重置
        {
            send_e_buf.clear();
            //send_e_buf.resize(MPI_size);
            for (size_t i = 0; i < MPI_size; i++)
            {
                vector<int> a={-1};
                send_e_buf.push_back(a);
            }  
        }

        void reset_send_to_rank()//重置
        {
            for (size_t i = 0; i < MPI_size; i++)
            {
                send_to_rank[i][0]=-1;
            }      
        }

        void send_event_to_rank()
        {
            for (int com_rank = 0; com_rank < MPI_size; ++com_rank)
            {

                MPI_Win_fence(0, r_f_win[com_rank]);

                    if (com_rank==MPI_rank)
                    {

                        for (int rank_num = 0; rank_num < MPI_size; rank_num++)
                        {

                    if(send_to_rank[rank_num][0]>0 && send_out_debug)
                    {
                     cout<<"send_to_rank["<<rank_num<<"][0]= "<<send_to_rank[rank_num][0]<<" rank is "<<MPI_rank<<" ["<<rank_num<<"][1] "<<send_to_rank[rank_num][1]<<" ["<<rank_num<<"][2] "<<send_to_rank[rank_num][2]<<endl;
                    }

                            //cout<<"in MPI_rank send "<<endl;
                            if(send_to_rank[rank_num][0]>0&&MPI_rank!=rank_num)
                            {
                            //cout<<"start put "<<endl;
                            MPI_Put(&send_to_rank[rank_num][0],send_to_rank[rank_num][0]+2,MPI_INT,rank_num ,0,send_to_rank[rank_num][0]+2,MPI_INT,r_f_win[com_rank]); 
                            //cout<<"end put "<<endl;
                            }
                            else if (send_to_rank[rank_num][0]>0&&MPI_rank==rank_num)
                            {
                                if(send_self_debug)
                                {
                                cout<<"recved!!!!!!!!!!!!!!!!!1"<<"send_to_rank["<<rank_num<<"][0] MPI_rank is "<<MPI_rank<<endl;
                                }
                                for (int i = 0; i < send_to_rank[rank_num][0]+1; i++)
                                {
                                    recv_from_rank[rank_num][i]=send_to_rank[rank_num][i];
                                }
                            } 
                        }
                    }
                    MPI_Win_fence(0, r_f_win[com_rank]);
                // if(MPI_rank==0)
                // {
                //     if (recv_from_rank[1][0]>0)
                //     {
                //         cout<<"recv_from_rank[1][0]>0"<<endl;
                //     }
                    
                // }

            }
        }




        void Clear_event_set()//清除所有的事件
        {
            Event_set.clear();
        }

        void Open_rate_record()//开启神经元记录rate
        {
            for (size_t i = 0; i < neuron_set.size(); i++)
            {
                neuron_set[i].on_record=true;
            }
            
        }

        void print_n_rate()
        {
            for (size_t i = 0; i < neuron_set.size(); i++)
            {
                cout<<"N: "<<neuron_set[i].global_num<<" rate_set : ";
                for (size_t j = 0; j < neuron_set[i].rate_record.size(); j++)
                {
                    cout<<neuron_set[i].rate_record[j]<<" ";
                }
                cout<<endl;
                
            }
            
        }


        void print_vm(int neuron_globle_id)//输出一个神经元的Vm
        {
            for (size_t i = 0; i < neuron_set.size(); i++)
            {
                for (size_t j = 0; j < RING_BUFFER_LEN; j++)
                {
                    cout<<"Vm: "<<neuron_set[neuron_globle_id].Recording_Vm[j]<<endl;
                }
            }
        }


        void Set_Record(int record_)//设置为0神经元就不记录Vm的变化，可以加快计算速度
        {
            Record=record_;
        } 

        void set_print_Vm(int print_switch)
        {
            if(print_switch==1)
            {
                print_vm_switch=1;
            }
            else
            {
                print_vm_switch=0;
            }
        }
};




// int main(int argc, char const *argv[])
// {

// //MPI 测试
//     Mpi_Manage ma;
//     // cout<<ma.MPI_rank<<" "<<ma.MPI_size<<endl;
//     double pop_de_rate=0.1;
//     int pop_num=0;
//     int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};
//     // for (size_t i = 0; i < 8; i++)
//     // {
//     //     pop_num+=pop_nums[i];
//     //     cout<<pop_nums[i]<<"pop_num: "<<pop_num<<endl;
//     // }
    
//     if(ma.MPI_rank==0)//这些是在一个进程上做
//     {
    
//     // double syn_weight=3.14;
//     // int pop1,pop2;
//     // // int pop_nums[8]={int(2068*pop_de_rate),int(583*pop_de_rate),int(2191*pop_de_rate),int(547*pop_de_rate),int(485*pop_de_rate),int(106*pop_de_rate),int(1439*pop_de_rate),int(294*pop_de_rate)};
//     // double connect_pr[8][8]={
//     //     {0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0., 0.0076, 0.},
//     //          {0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0., 0.0042, 0.},
//     //          {0.0077, 0.0059, 0.0497, 0.135, 0.0067, 0.0003, 0.0453, 0.},
//     //          {0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0., 0.1057, 0.},
//     //          {0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.},
//     //          {0.0548, 0.0269, 0.0257, 0.0022, 0.06, 0.3158, 0.0086, 0.},
//     //          {0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252},
//     //          {0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443}
//     // };
//     // // pop1=ma.create(1,2);//创建神经元（此时没有创建实体，只是构建表格）单进程
//     // // pop2=ma.create(1,2);//创建神经元（此时没有创建实体，只是构建表格）单进程
//     // // cout<<"pop1 is "<<pop1<<" pop2 is "<<pop2<<endl;
//     // // ma.full_connect(pop1,pop2,100);//创建连接（此时没有创建实体，只是构建表格）单进程
//     // int L23E=ma.create(1,pop_nums[0]);
//     // int L23I=ma.create(1,pop_nums[1]);
//     // int L4E=ma.create(1,pop_nums[2]);
//     // int L4I=ma.create(1,pop_nums[3]);
//     // int L5E=ma.create(1,pop_nums[4]);
//     // int L5I=ma.create(1,pop_nums[5]);
//     // int L6E=ma.create(1,pop_nums[6]);
//     // int L6I=ma.create(1,pop_nums[7]);
//     // int populatins[8]={L23E,L23I,L4E,L4I,L5E,L5I,L6E,L6I};
//     // double mu[8][8]={
//     //     {277,-1110,555,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110},
//     //     {277,-1110,277,-1110,277,-1110,277,-1110}
//     // };
//     // double sigma[8][8]={
//     //     {27,111,55,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111},
//     //     {27,111,27,111,27,111,27,111}
//     // };

//     // int conn_weight=0;
//     // for (int i = 0; i < 8; i++)
//     // {
//     //     for (int j = 0; j < 8; j++)
//     //     {
//     //         conn_weight=gaussrand()*sigma[i][j]*syn_weight+mu[i][j]*syn_weight;

//     //         ma.FIX_random_connect_bad_4(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);            
//     //         // if(conn_weight>0)
//     //         // ma.STDP_random_connect(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
//     //         // else
//     //         // ma.random_connect(populatins[i],populatins[j],conn_weight,connect_pr[i][j]);
//     //     }
        
//     // }
//     // ma.neuron_assignment(1);//根据进程以及神经元分配方案进行神经元分配，单进程
//     }
//     // ma.prepare_msg();
//     if(ma.MPI_rank==0){
//     // ma.record_Gid_to_rank();
//     // ma.record_adjtable();
//     // ma.record_STDP_adj_table();
//     // ma.record_weight_table();
//     // ma.record_STDP_w_table();
//     }
//     ma.prepare_read_msg();
//     //possion
//     double k_ex[8]={1600, 1500, 2100, 1900, 2000, 1900, 2900, 2100};
//     double g=8.;
//     double Poisson_weight=277;
//     double possion_fixed=0.0001;
//     double possion_weight=878;
//     int Gid_index=0;
//     for (int i = 0; i < 8; i++)
//     {
//         ma.Poisson_full_connect(k_ex[i]*g*0.05,278,Gid_index,Gid_index+pop_nums[i]);
//         Gid_index=Gid_index+pop_nums[i];
//     }
    
    
//     // ma.Open_rate_record();
//     // ma.print_n_rate();
//     //possion
//     struct timeval t1,t2;
//     gettimeofday(&t1, nullptr);
//     ma.simulate(60000);
//     gettimeofday(&t2, nullptr);

//     //ma.print_n_rate();
//     double usetime_s = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000.0;
//     cout<<"rank: "<<ma.MPI_rank<<" use time: "<<usetime_s<<endl;
//     cout<<"rank: "<<ma.MPI_rank<<" STDP_count: "<<ma.STDP_count<<endl;        
//     cout<<"N : "<<ma.neuron_set[1].global_num<<" spike_times : "<<ma.neuron_set[1].spike_time<<endl;
//     cout<<"N : "<<ma.neuron_set[100].global_num<<" spike_times : "<<ma.neuron_set[100].spike_time<<endl;
//     // for (size_t i = 0; i < 10; i++)
//     // {
//     //     cout<<ma.neuron_set[0].Recording_Vm[i]<<" ";
//     // }
//     // cout<<endl;

//     // for (size_t i = 0; i < 10; i++)
//     // {
//     //     cout<<ma.neuron_set[0].Recording_Vm[i]<<" ";
//     // }
//     cout<<endl;



// }
