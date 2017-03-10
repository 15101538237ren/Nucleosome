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
using namespace std;

void generate_nucleosome_pos_list(long start_pos, long end_pos, string output_file_path, int constant_len, int constant_gap, double ratio_list[])
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
    
    long now_start = start_pos;
    long now_end = now_start;
    
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
        sprintf(buffer, "%ld\t%ld\t%c\n",now_start , now_end, n_status);
        outfile << buffer;
        
        now_start = now_end + constant_gap;
    }
    outfile.close();
    printf("Write successful\n");
}

int main(int argc, const char * argv[]) {
    string base_dir = "/Users/Ren/XCodeProjects/Nucleosome/Nucleosome/";
    long start_pos = 3001600;
    
    long end_pos = 195365800;
    
    string output_file_path = base_dir + "nucleosome_positions.np";
    
    int constant_len = 160 ;
    
    int constant_gap = 5;
    
    double ratio_list[] = {0.5, 0.3, 0.2};
    
    generate_nucleosome_pos_list(start_pos , end_pos , output_file_path , constant_len , constant_gap , ratio_list);
    return 0;
}
