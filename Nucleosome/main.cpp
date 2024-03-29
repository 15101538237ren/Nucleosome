//
//  main.cpp
//  Methylation_Server
//
//  Created by 任红雷 on 2016/12/29.
//  Copyright (c) 2016 ___renhonglei___. All rights reserved.
//

#include<iostream>
#include<cstdio>
#include<cmath>
#include <stdlib.h>
#include<sstream>
#include<string>
#include<vector>
#include <algorithm>
#include <functional>
#include<cstdlib>
#include<ctime>
#include<time.h>
#include<regex>
#include <numeric>
#include <functional>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#define pb push_back

using namespace std;

vector<int> nucleosome_start_pos_list;
vector<int> nucleosome_end_pos_list;
string nucleosome_status;

//构建一个最大CpG数量为max_cpg_sites，位置符合超几何分布G(p)的CpG位置向量
vector<int> construct_n_cpg_sites_for_exp_distribution(int max_cpg_sites, float p)
{
    vector<int> index_pos_list;
    int first_cpg_pos = rand() % 9 + 1;
    index_pos_list.pb(first_cpg_pos);
    
    default_random_engine generator(time(NULL));
    geometric_distribution<int> distribution(p);
    int number;
    for (int i=1; i< max_cpg_sites; ++i) {
        number= distribution(generator)+2;
        first_cpg_pos=first_cpg_pos+number;
        //        cout << i << " : " << first_cpg_pos <<endl;
        index_pos_list.pb(first_cpg_pos);
    }
    
    return index_pos_list;
}

//从输入input_bed_file_path的bed文件中，获取max_cpg_sites个CpG的位置列表
vector<int> get_pos_list_from_bed_file(string input_bed_file_path,int max_cpg_sites=1000)
{
    ifstream bed_file(input_bed_file_path);
    int  buf_size=100;
    int cnt_umh=0;
    int u_count=0;
    int m_count=0;
    int h_count=0;
    
    char buffer[buf_size];
    vector<int> pos_list;
    
    if (!bed_file)
    {
        cout << "Unable to open " << input_bed_file_path << endl;
    }
    else
    {
        int pos=0;
        float methy_level=0.0;
        int cnt=0;
        
        while (! bed_file.eof() )
        {
            if (max_cpg_sites < 0 or (cnt < max_cpg_sites and max_cpg_sites > 0))
            {
                bed_file.getline(buffer,buf_size);
                sscanf(buffer,"%d %f",&pos,&methy_level);
                if (methy_level<0.2)
                {
                    u_count=u_count+1;
                }
                else if(methy_level>0.8)
                {
                    m_count=m_count+1;
                }
                else
                {
                    h_count=h_count+1;
                }
                pos_list.pb(pos);
                cnt++;
            }
            else
            {
                break;
            }
        }
        cnt_umh=m_count+h_count+u_count;
//        float u_ratio=float(u_count)/float(cnt_umh);
//        float m_ratio=float(m_count)/float(cnt_umh);
        
        bed_file.close();
        printf("get pos list from %s successfully!\n",input_bed_file_path.c_str() );
    }
    return pos_list;
}

//输入CpG位点的数量,M,U状态所占的比例,输出一条符合UMH比例分布的string链，每个状态用'U','H','M'代表
string generate_CpG_in_methylation_percent_UHM(int CpG_sites_counts,double m_ratio, double u_ratio)
{
    string CpG_str;
    
    for (int i = 0; i < CpG_sites_counts; i++)
    {
        double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
        
        if (random_num < m_ratio)
        {
            CpG_str.append(1, 'M');
        }
        else if (random_num > m_ratio and random_num < m_ratio + u_ratio)
        {
            CpG_str.append(1, 'U');
        }
        else
        {
            CpG_str.append(1, 'H');
        }
    }
    
    return CpG_str;// 返回生成的状态字符串, 长度=CpG_sites_counts
}

//生成核小体的位置文件，输入: 起始核小体在染色体上的位置, 结束位置, 输出文件路径, 每个核小体长度, 每个gap的长度, 每种状态核小体的比例(分别为AA,AU,UU三种状态的浮点列表，它们的大小比值即为比例)
void generate_nucleosome_pos_list(int start_pos, int end_pos, string output_file_path, int constant_len, int constant_gap, double ratio_list[])
{
    double sum_ratio = 0.0;
    for(int i=0; i < 3; i++){
        sum_ratio += ratio_list[i];
    }
    
    char nucleosome_status_list[] = {'A','H','U'};
    
    srand(time(0));
    
    ofstream outfile(output_file_path);
    int buffer_size = 50;
    if(!outfile){
        cout << "Unable to open file: " << output_file_path << endl;
        exit(1);
    }
    
    int now_start = start_pos;
    int now_end = now_start;
    
    int nucleosome_count = 0;
    while (now_end + constant_len < end_pos) {
        nucleosome_count += 1;
        if (nucleosome_count % 10000 == 0)
        {
            printf("Now write %d nucleosomes\n" , nucleosome_count);
        }
        
        now_end = now_start + constant_len;
        
        double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
        double tmp_sum_ratio = 0.0;
        random_num = random_num * sum_ratio;
        char n_status='U';
        
        for(int i = 0; i< 3; i++){
            tmp_sum_ratio = tmp_sum_ratio + ratio_list[i];
            if(random_num < tmp_sum_ratio){
                n_status = nucleosome_status_list[i];
                break;
            }
        }
        
        char buffer[buffer_size];
        sprintf(buffer, "%d\t%d\t%c\n",now_start , now_end, n_status);
        outfile << buffer;
        
        now_start = now_end + constant_gap;
    }
    outfile.close();
    printf("Write successful\n");
}

//产生核小体位置文件
void prepare_nucleosome_pos_file()
{
    string base_dir = "/Users/Ren/XCodeProjects/Nucleosome/Nucleosome/output";
    
    int constant_len = 160 ;
    
    int start_pos = 3001600;
    
    int end_pos = 195365800 + constant_len;
    
    string output_file_path = base_dir + "nucleosome_positions.np";
    
    int constant_gap = 5;
    
    double ratio_list[] = {0.5, 0.3, 0.2};
    
    generate_nucleosome_pos_list(start_pos , end_pos , output_file_path , constant_len , constant_gap , ratio_list);
}

vector<int> get_neucleosome_id_list_from_pos_list(vector<int> index_pos_list, string nucleosome_pos_file_path)
{
    vector<int> nucleo_id_list;
    ifstream nucleosome_pos_file(nucleosome_pos_file_path);
    
    if(!nucleosome_pos_file)
    {
        cout << "Unable to open read " << nucleosome_pos_file_path << endl;
        exit(1);
    }
    int buffer_size = 50;
    char buff[buffer_size];
    int now_nuceo_start = -1;
    int now_nuceo_end = -1;
    char status = 'N';
    
    int cpg_pos_index = 0;
    int nuceo_pos_index = 0;
    
    nucleosome_pos_file.getline(buff,buffer_size);
    sscanf(buff,"%d\t%d\t%c\n",&now_nuceo_start,&now_nuceo_end,&status);
    nucleosome_start_pos_list.push_back(now_nuceo_start);
    nucleosome_end_pos_list.push_back(now_nuceo_end);
    nucleosome_status.append(1,status);
    
    for (int i=0; i< index_pos_list.size(); i++) {
        int now_cpg_pos = index_pos_list[i];
        if ( now_cpg_pos <= now_nuceo_end)
        {
            if(now_cpg_pos >= now_nuceo_start)
            {
                //在组蛋白的范围内,则把该组蛋白的id加入到nucleo_id_list中
                nucleo_id_list.push_back(nuceo_pos_index);
            }
            else
            {
                //在gap中的CpG对应的核小体编号设置为-1
                nucleo_id_list.push_back(-1);
            }
        }
        else{
            while(now_cpg_pos > now_nuceo_end && ! nucleosome_pos_file.eof())
            {
                //当前CpG的位置>核小体的范围,继续读核小体的文件
                nucleosome_pos_file.getline(buff,buffer_size);
                sscanf(buff,"%d\t%d\t%c\n",&now_nuceo_start,&now_nuceo_end,&status);
                
                nucleosome_start_pos_list.push_back(now_nuceo_start);
                nucleosome_end_pos_list.push_back(now_nuceo_end);
                nucleosome_status.append(1,status);
                
                nuceo_pos_index += 1;
            }
            i--;
        }
    }
    nucleosome_pos_file.close();
    return nucleo_id_list;
}

void simulate(int generation, int time_step, string init_cell, string detail_file_dir, string ratio_file_dir,string nucleosome_pos_file_path, vector<double> propensity_list, vector<int> index_pos_list, int max_cells){
    
}

vector<vector<int>> get_corresponding_cpg_id_list(int cpg_id_list_size,vector<int> nucleosome_id_list)
{
    vector<vector<int>> cpg_id_list(cpg_id_list_size);
    for (int i = 0; i < nucleosome_id_list.size(); i++)
    {
        if(nucleosome_id_list[i] != -1)
        {
            cpg_id_list[nucleosome_id_list[i]].push_back(i);
        }
    }
    return cpg_id_list;
}

void start_simulation()
{
    int round_start = 1;
    int round_end = 6;
    
    int generations = 30;
    int max_cpg_sites = -1;
    string init_cell;
    double m_ratio = 0.181214;
    double u_ratio = 0.391004;
    
    vector<double> propensity_list = {0.0001,0.9899, 0.005, 0.005};
    
    vector<int> index_pos_list;
    
    //当前文件夹的全路径
    string path_dir="/Users/Ren/XCodeProjects/Nucleosome/Nucleosome/";
    
    string output_dir_name = "output/";
    string input_dir_name = "input/";
    string output_dir = path_dir + output_dir_name;
    string input_bed_file_path = path_dir + input_dir_name + "chr1.bed";
    string nucleosome_pos_file_path = output_dir + "nucleosome_positions.np";
    string ratio_file_dir;
    string detail_file_dir;
    string bed_file_dir;
    string rd_without_dir;
    
    int time_step=100;
    int d_max=1000;
    bool calc_interval=false;
    bool ignore_d=false;
    
    int max_cells = 2;
    
    bool simulation=true;
    bool calc_corr=true;
    
    bool real_chr_pos=true;
    
    if(!real_chr_pos){
        //超几何分布参数
        float geometric_p = 0.3;
        index_pos_list = construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p);
    }
    else{
        index_pos_list = get_pos_list_from_bed_file(input_bed_file_path,max_cpg_sites);
    }
    
    //根据CpG在数组index_pos_list中的下标,查找对应的核小体下标
    vector<int> neucleosome_id_list = get_neucleosome_id_list_from_pos_list(index_pos_list,nucleosome_pos_file_path);
    
    int nucleosome_list_size = (int) nucleosome_end_pos_list.size();
    vector<vector<int>> cpg_id_list = get_corresponding_cpg_id_list(nucleosome_list_size,neucleosome_id_list);
    
    init_cell = generate_CpG_in_methylation_percent_UHM(max_cpg_sites,m_ratio,u_ratio);
    ratio_file_dir = output_dir+"ratio/";//detail的生成序列
    detail_file_dir=output_dir+"detail/";//detail的生成序列
    
    if (simulation){
        for(int round_i=round_start;round_i<=round_end;round_i++)
        {
            simulate(generations,time_step,init_cell,detail_file_dir,ratio_file_dir,nucleosome_pos_file_path,propensity_list,index_pos_list,max_cells);
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}

int main(int argc, const char * argv[]) {
    
    start_simulation();
//    system("say Mission completed!");
    return 0;
}
