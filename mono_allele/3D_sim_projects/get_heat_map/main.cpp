#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstring>
#include <cmath>
#include <chrono>
#include <random>
#include <map>
#include <set>
#include <algorithm>
#include <list>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <experimental/filesystem>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>


struct Command_content
{
    std::vector<std::string> sr_nows;
    std::vector<std::string> capture_ratios;
    std::experimental::filesystem::path out_path;
    double density;
    std::string weights;
    Command_content(int argc, char **argv)
    {
        for(int i=1; i<argc-1; ++i)
        {
            if(!strcmp(argv[i],"-sr_nows"))
                sr_nows=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-capture_ratios"))
                capture_ratios=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-out_path"))
                out_path.assign(argv[i+1]);
            if(!strcmp(argv[i],"-density"))
                density=stod(std::string(argv[i+1]));
            if(!strcmp(argv[i],"-weights"))
                weights=argv[i+1];
        }
    }
public:
    std::vector<std::string> my_split(const char* argv_cmp)
    {
        std::vector<std::string> tmp(1);
        for(const char* c=argv_cmp; *c; ++c)
        {
            if(*c==',')
                tmp.push_back(std::string());
            else
                tmp.back().push_back(*c);
        }
        return tmp;
    }
};

double cal_squ_dis(std::array<double,3>& pos1, std::array<double,3>& pos2, double length);
std::string exec(const char* cmd);
int main(int argc, char **argv) 
{
    Command_content cc(argc, argv);
    std::vector<double> square_capture_ratios;
    for(auto& str : cc.capture_ratios)
        square_capture_ratios.push_back(std::pow(stod(str),2.0));
    std::vector<int> record_1Ds;
    for(auto& sr_now : cc.sr_nows)
        record_1Ds.push_back(stoi(cc.my_split(exec(("wc -l "+sr_now).c_str()).c_str())[0]));
    int cum_record_1D=std::accumulate(record_1Ds.begin(),record_1Ds.end(),0);
    std::vector<double> weights(cum_record_1D,1.0/cum_record_1D);
    if(!cc.weights.empty())
    {
        std::ifstream fin(cc.weights);
        for(auto& weight : weights)
            fin >> weight;
        fin.close();
    }
    
    std::vector<std::array<double,3>> positions;
    std::string str_tmp;
    std::ifstream fin(cc.sr_nows.front());
    getline(fin,str_tmp);
    fin.close();
    std::istringstream iss(str_tmp);
    iss >> str_tmp;
    while(iss.good())
    {
        positions.emplace_back();
        iss >> str_tmp >> str_tmp >> str_tmp;
    }
    double length=std::pow(positions.size()/cc.density,1.0/3.0);
    std::vector<std::vector<std::vector<double>>> hics(square_capture_ratios.size(),std::vector<std::vector<double>>(positions.size(),std::vector<double>(positions.size(),0.0)));
    
    int we=0;
    for(int sr=0; sr<cc.sr_nows.size(); ++sr)
    {
        std::ifstream fin(cc.sr_nows[sr]);
        for(int re=0; re<record_1Ds[sr]; ++re)
        {
            getline(fin,str_tmp);
            std::istringstream iss(str_tmp);
            for(auto& position : positions)
            {
                iss >> str_tmp;
                position[0]=stod(str_tmp);
                iss >> str_tmp;
                position[1]=stod(str_tmp);
                iss >> str_tmp;
                position[2]=stod(str_tmp);
            }
            for(int i=0; i<positions.size(); ++i)
                for(int j=i+1; j<positions.size(); ++j)
                {
                    double seq_dis=cal_squ_dis(positions[i], positions[j], length);
                    for(int k=0; k<square_capture_ratios.size(); ++k)
                        if(seq_dis<square_capture_ratios[k])
                            hics[k][i][j]+=weights[we];
                }
            ++we;
        }
        fin.close();
    }
    for(int k=0; k<square_capture_ratios.size(); ++k)
    {
        cc.out_path/="hic_"+std::to_string(int(0.5+std::sqrt(square_capture_ratios[k])));
        std::ofstream fout(cc.out_path.string());
        cc.out_path=cc.out_path.parent_path();
        for(int i=0; i<positions.size(); ++i)
            for(int j=0; j<positions.size(); ++j)
            {
                if(j<i)
                    fout << hics[k][j][i];
                else
                    fout << hics[k][i][j];
                if(j<positions.size()-1)
                    fout << '\t';
                else
                    fout << '\n';
            }
        fout.close();
    }
    return 0;
}

double cal_squ_dis(std::array<double,3>& pos1, std::array<double,3>& pos2, double length)
{
    std::array<double,3> vec;
    for(int dim=0; dim<3; ++dim)
    {
        vec[dim]=std::fmod(pos1[dim]-pos2[dim]+length/2,length);
        if(vec[dim]<0)
            vec[dim]+=length;
        vec[dim]-=length/2;
    }
    return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
}

std::string exec(const char* cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}
