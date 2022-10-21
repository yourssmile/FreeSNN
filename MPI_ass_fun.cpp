//#include<mpi.h>

#include<iostream>
#include<vector>
#include<cmath>
#include<bits/stdc++.h>
using namespace std;
#define PUT_BUF_LEN 50000//MPI通信的缓冲区

double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
 
    return X;
}

template<typename T>
inline vector<T> create_send_adj_table(vector<vector<T>> adj_table)//把二维数组变成一维，方便通信
{
    int N=adj_table.size();
    vector<T> res;
    res.push_back(N);
    for (size_t i = 0; i < N; i++)
    {
        res.push_back(adj_table[i].size());
    }
    for (size_t i = 0; i < adj_table.size(); i++)
    {
        for (size_t j = 0; j < adj_table[i].size(); j++)
        {
            res.push_back(adj_table[i][j]);
        }
        
    }
    
    
    return res;
};

template<typename T>
inline void vector_to_int_set(vector<T> table,T rec[PUT_BUF_LEN])
{
    int size=table.size();
    if (size>PUT_BUF_LEN)
    {
        std::cout<<"table size oversize than PUT_BUF_LEN! please change PUT_BUF_LEN (by XC)"<<std::endl;
    }
    
    for (size_t i = 0; i < size; i++)
    {
        rec[i]=table[i];
    }
    
}

template<typename T>
inline vector<T> int_set_to_vector(T *in_set,int in_set_size)
{
    vector<T> rec_v;
    for (int i = 0; i < in_set_size; ++i)
    {
        rec_v.push_back(in_set[i]);
    }
    return rec_v;
}

template<typename T>
inline vector<vector<T>> one_v_to_2_v(vector<T> vec_in)//一维恢复成二维
{
    vector<vector<T>> res;
    if (vec_in.size()<1)
    {
        std::cout<<"vec_in size is 0!"<<std::endl;
        return res;
    }
    int N=vec_in[0];
    int bias=N+1;
    res.resize(N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < vec_in[i+1]; ++j)
        {
            res[i].push_back(vec_in[bias]);
            bias++;
        }
    }
    return res;

}

template<typename T>
inline void print_2_v(vector<vector<T>> vec_in)
{
    for (size_t i = 0; i < vec_in.size(); i++)
    {
        for (size_t j = 0; j < vec_in[i].size(); j++)
        {
            cout<<vec_in[i][j]<<" ";
        }
        cout<<endl;
        
    }
    
};

template<typename T>
inline void print_1_v(vector<T> vec_in)
{
    for (size_t i = 0; i < vec_in.size(); i++)
    {
        cout<<vec_in[i]<<" ";
    }
    cout<<endl;
    
}

//template<typename T>
inline bool search(vector<int> nums, int target) //搜索target是否在nums中
{

  vector<int>::iterator t;
  
  t = find(nums.begin(),nums.end(),target);
  
  if(t != nums.end()){
  	return true;
  } 
  return false ;
}

inline void print_map(map<int,int> in_map)
{

    map<int,int>::iterator it;
    for(it=in_map.begin();it!=in_map.end();it++)
    {
        cout<<it->first<<" "<<it->second<<endl;
    }
 
}


// int main(int argc, char const *argv[])
// {
//     vector<vector<int>> test_in;
    
//     test_in.resize(10);
//     for (size_t i = 0; i < 10; i++)
//     {
//         test_in[i].push_back(-1);
        
//     }
//     for (size_t i = 0; i < 10; i++)
//     {
//         test_in[i].push_back(i+1);
//     }
//     test_in[3].push_back(98);
//     test_in[1].push_back(78);
//     test_in[1].push_back(78);
//     print_vec(test_in);
//     vector<int> rec=create_send_adj_table(test_in);
//     print_sing_vec(rec);
//     vector<vector<int>> res_1=send_table_to_vec(rec);
//     print_vec(res_1);
//     int rec_set[PUT_BUF_LEN]={0};
//     vector_to_int_set(rec,rec_set);
//     for (size_t i = 0; i < rec.size()+1; i++)
//     {
//         std::cout<<rec_set[i]<<endl;
//     }
    
//     return 0;
// }
