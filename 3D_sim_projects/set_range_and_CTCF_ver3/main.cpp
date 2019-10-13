#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <map>
#include <algorithm>
#include <cmath>
#include "/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/Eigen/Dense"
#include <numeric>
#include <experimental/filesystem>

const int bin=600;
const double mu=3;
const int VPE=10;
const double nipbl_ratio=0.80;

struct Command_content
{
    int start;
    int end;
    std::vector<std::string> C4s_rep1;
    std::vector<std::string> C4s_rep2;
    std::vector<std::string> CTCF_peaks;
    std::string mut_inform;
    std::string chromosome;
    std::vector<std::string> nipbls;
    
    std::vector<std::string> stiffnesses;
    std::vector<std::string> densities;
    std::vector<std::string> processivities;
    std::vector<std::string> separations;
    std::vector<std::string> D3_per_D1s;
    std::vector<std::string> capture_ratios;
    
    int record_1D;
    Command_content(int argc, char **argv)
    {
        record_1D=-1;
        for(int i=1; i<argc-1; ++i)
        {
            if(!strcmp(argv[i],"-start"))
                start=stoi(std::string(argv[i+1]));
            if(!strcmp(argv[i],"-end"))
                end=stoi(std::string(argv[i+1]));
            if(!strcmp(argv[i],"-C4s_rep1"))
                C4s_rep1=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-C4s_rep2"))
                C4s_rep2=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-CTCF_peaks"))
                CTCF_peaks=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-mut_inform"))
                mut_inform=std::string(argv[i+1]);
            if(!strcmp(argv[i],"-chromosome"))
                chromosome=std::string(argv[i+1]);
            if(!strcmp(argv[i],"-nipbls"))
                nipbls=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-stiffnesses"))
                stiffnesses=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-densities"))
                densities=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-processivities"))
                processivities=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-separations"))
                separations=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-D3_per_D1s"))
                D3_per_D1s=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-capture_ratios"))
                capture_ratios=my_split(argv[i+1]);
            if(!strcmp(argv[i],"-record_1D"))
                record_1D=stoi(std::string(argv[i+1]));
        }
    }
    std::vector<std::string> my_split(char* argv_cmp)
    {
        std::vector<std::string> tmp(1);
        for(char* c=argv_cmp; *c; ++c)
        {
            if(*c==',')
                tmp.push_back(std::string());
            else
                tmp.back().push_back(*c);
        }
        return tmp;
    }
};

double stall_pro(double CTCF_enrich);
bool isNumberC(const std::string& s);
std::tuple<std::vector<int>,std::vector<std::map<std::pair<int,int>,double>>> get_single_target(std::vector<std::vector<std::string>> C4s_reps, std::string chromosome, int start, int end, int length);
std::tuple<std::vector<double>,std::vector<double>> calculate_stall(std::vector<std::pair<std::pair<int,int>,std::pair<std::string,std::string>>>& mut_list, std::vector<std::string> CTCF_peaks, std::string chromosome, int start, int end, int length);
std::vector<double> calculate_birth(std::vector<std::pair<std::pair<int,int>,std::pair<std::string,std::string>>>& mut_list, std::vector<std::string> nipbls, std::string chromosome, int start, int end, int length);
int main(int argc, char **argv) 
{
    Command_content cc(argc,argv);
    // extend length to even multiples of bin
    if((cc.end-cc.start)%bin>0)
        cc.end=cc.start+((cc.end-cc.start)/bin+1)*bin;
    if((cc.end-cc.start)/bin%2>0)
        cc.end+=bin;
    
    int length=(cc.end-cc.start)/bin;
    
    std::vector<std::pair<std::pair<int,int>,std::pair<std::string,std::string>>> mut_list;
    if(!cc.mut_inform.empty())
    {
        std::ifstream fin(cc.mut_inform);
        std::string str_tmp;
        fin >> str_tmp;
        while(fin.good())
        {
            mut_list.push_back(std::pair<std::pair<int,int>,std::pair<std::string,std::string>>());
            mut_list.back().first.first=stoi(str_tmp);
            fin >> str_tmp;
            mut_list.back().first.second=stoi(str_tmp);
            fin >> mut_list.back().second.first;
            fin >> mut_list.back().second.second;
            for(auto& ch : mut_list.back().second.second) ch=toupper(ch);
            fin >> str_tmp;
        }
        fin.close();
    }
    
    std::ofstream fout("pos_file");
    
    if(!cc.C4s_rep1.empty())
    {
        std::vector<int> VPBs;
        std::vector<std::map<std::pair<int,int>,double>> targets;
        std::tie(VPBs,targets)=get_single_target(std::vector<std::vector<std::string>>{cc.C4s_rep1,cc.C4s_rep2}, cc.chromosome, cc.start, cc.end, length);
        
        fout << "VP\n";
        for(auto& vpb : VPBs) fout << vpb << '\n';
        
        fout << "target\n";
        std::vector<std::map<std::pair<int,int>,double>::iterator> iters;
        for(auto& target : targets) iters.push_back(target.begin());
        while(iters.front()!=targets.front().end())
        {
            double mean=0;
            for(auto& iter : iters)
                mean+=iter->second;
            mean/=targets.size();
            double variance=0;
            for(auto& iter : iters)
                variance+=std::pow(iter->second-mean,2.0);
            variance/=targets.size();
            if(mean>0)
                fout << iters.front()->first.first << '\t' << iters.front()->first.second << '\t' << mean << '\t' << variance << '\n';
            for(auto& iter : iters) ++iter;
        }
    }
    
    if(!cc.CTCF_peaks.empty())
    {
        std::vector<double> left_stall, right_stall;
        std::tie(left_stall, right_stall)=calculate_stall(mut_list, cc.CTCF_peaks, cc.chromosome, cc.start, cc.end, length);
    
        fout << "left_stall\n";
        for(auto& ls : left_stall) fout << ls << '\n';
        
        fout << "right_stall\n";
        for(auto& rs : right_stall) fout << rs << '\n';
    }
    
    if(!cc.nipbls.empty())
    {
        std::vector<double> birth=calculate_birth(mut_list, cc.nipbls, cc.chromosome, cc.start, cc.end, length);
        
        fout << "birth\n";
        for(auto& bi : birth) fout << bi << '\n';
    }
    
    if(!cc.stiffnesses.empty())
    {
        fout << "stiffness\n";
        for(auto& stiffness : cc.stiffnesses) fout << stiffness << '\n';
    }
    
    if(!cc.densities.empty())
    {
        fout << "density\n";
        for(auto& density : cc.densities) fout << density << '\n';
    }
    
    if(!cc.processivities.empty())
    {
        fout << "processivity\n";
        for(auto& processivity : cc.processivities) fout << processivity << '\n';
    }
    
    if(!cc.separations.empty())
    {
        fout << "separation\n";
        for(auto& separation : cc.separations) fout << separation << '\n';
    }
    
    if(!cc.D3_per_D1s.empty())
    {
        fout << "D3_per_D1\n";
        for(auto& D3_per_D1 : cc.D3_per_D1s) fout << D3_per_D1 << '\n';
    }
    
    if(!cc.capture_ratios.empty())
    {
        fout << "capture_ratio\n";
        for(auto& capture_ratio : cc.capture_ratios) fout << capture_ratio << '\n';
    }
    
    if(cc.record_1D>0)
    {
        fout << "record_1D\n";
        fout << cc.record_1D << '\n';
    }
    
    fout.close();
    
    return 0;
}

double stall_pro(double CTCF_enrich)
{
    return 1/(1+std::exp(-CTCF_enrich/20-mu));
}

bool isNumberC(const std::string& s)
{
    char* p;
    strtod(s.c_str(), &p);
    return *p == 0;
}

std::tuple<std::vector<int>,std::vector<std::map<std::pair<int,int>,double>>> get_single_target(std::vector<std::vector<std::string>> C4s_reps, std::string chromosome, int start, int end, int length)
{
    std::vector<int> VPBs;
    std::vector<std::vector<std::vector<double>>> bin_heightss_reps;
    for(auto& C4s : C4s_reps)
    {
        bin_heightss_reps.emplace_back();
        for(auto C4 : C4s)
        {
            bin_heightss_reps.back().emplace_back(length,0.0);
            std::ifstream fin(C4);
            std::string str_tmp;
            fin >> str_tmp;
            int VP=stoi(str_tmp);
            int VPB=(VP-start)/bin;
            if(&C4s==C4s_reps.data())
                VPBs.push_back(VPB);
            if(VPB<0 || VPB>=length)
                continue;
            std::vector<std::pair<std::pair<int,int>,double>> C4_data;
            fin >> str_tmp;
            while(fin.good())
            {
                if(str_tmp==chromosome)
                {
                    fin >> str_tmp;
                    int tmp1=stoi(str_tmp);
                    fin >> str_tmp;
                    int tmp2=stoi(str_tmp);
                    if(std::min(tmp2,end)>std::max(tmp1,start))
                    {
                        C4_data.emplace_back();
                        C4_data.back().first.first=tmp1;
                        C4_data.back().first.second=tmp2;
                        fin >> str_tmp;
                        C4_data.back().second=stod(str_tmp); 
                    }
                    else
                        getline(fin,str_tmp);
                    fin >> str_tmp;
                }
                else
                {
                    getline(fin,str_tmp);
                    fin >> str_tmp;
                }
            }
            fin.close();
            std::sort(C4_data.begin(),C4_data.end());
            
            for(int j=0; j<C4_data.size(); ++j)
            {
                int tmp1=(C4_data[j].first.first-start)/bin;
                int tmp2=(C4_data[j].first.second-start)/bin;
                for(int k=std::max(0,tmp1); k<=std::min(length-1,tmp2); ++k)
                {
                    int LB, RB;
                    if(k==tmp1)
                        LB=C4_data[j].first.first;
                    else
                        LB=start+k*bin;
                    if(k==tmp2)
                        RB=C4_data[j].first.second;
                    else
                        RB=start+(k+1)*bin;
                    bin_heightss_reps.back().back()[k]+=double(RB-LB)/bin*C4_data[j].second;
                }
            }
        }
    }
    
    int ord=bin_heightss_reps.size()*bin_heightss_reps.front().size();
    Eigen::MatrixXd A=Eigen::MatrixXd::Zero(ord,ord);
    Eigen::VectorXd b=Eigen::VectorXd::Zero(ord);
    std::vector<int> nw;
    double meab=0;
    for(auto& bin_heightss : bin_heightss_reps)
        for(auto& bin_heights : bin_heightss)
        {
            nw.push_back(0);
            int VPB=VPBs[&bin_heights-bin_heightss.data()];
            for(int j=0; j<bin_heights.size(); ++j)
            {
                int sj=std::abs(j-VPB);
                if(sj>VPE && bin_heights[j]>0)
                {
                    ++nw.back();
                    meab+=std::log(double(sj));
                }
            }
        }
    int n=accumulate(nw.begin(),nw.end(),0);
    meab/=n;
    for(auto& bin_heightss : bin_heightss_reps)
        for(auto& bin_heights : bin_heightss)
        {
            int i=(&bin_heightss-bin_heightss_reps.data())*bin_heightss.size()+(&bin_heights-bin_heightss.data());
            int VPB=VPBs[&bin_heights-bin_heightss.data()];
            for(int j=0; j<bin_heights.size(); ++j)
            {
                int sj=std::abs(j-VPB);
                if(sj>VPE && bin_heights[j]>0)
                {
                    double lsj=std::log(double(sj));
                    double tmp=lsj-meab;
                    b(0)-=2*std::log(bin_heights[j])*tmp;
                    A(0,0)+=2*lsj*tmp;
                    if(i>0)
                        A(0,i)-=2*tmp;
                }
            }
        }
    for(int w=1; w<ord; ++w)
        for(auto& bin_heightss : bin_heightss_reps)
            for(auto& bin_heights : bin_heightss)
            {
                int i=(&bin_heightss-bin_heightss_reps.data())*bin_heightss.size()+(&bin_heights-bin_heightss.data());
                int VPB=VPBs[&bin_heights-bin_heightss.data()];
                double tmp=(-(i==w)+double(nw[w])/n);
                for(int j=0; j<bin_heights.size(); ++j)
                {
                    int sj=std::abs(j-VPB);
                    if(sj>VPE && bin_heights[j]>0)
                    {
                        double lsj=std::log(double(sj));
                        b(w)-=2*std::log(bin_heights[j])*tmp;
                        A(w,0)+=2*lsj*tmp;
                        if(i>0)
                            A(w,i)-=2*tmp;
                    }
                }
            }
    Eigen::VectorXd x=A.lu().solve(b);
    for(auto& bin_heightss : bin_heightss_reps)
        for(auto& bin_heights : bin_heightss)
        {
            int i=(&bin_heightss-bin_heightss_reps.data())*bin_heightss.size()+(&bin_heights-bin_heightss.data());
            if(i>0)
                for(int j=0; j<bin_heights.size(); ++j)
                    bin_heights[j]/=std::exp(x(i));
        }
    std::vector<std::map<std::pair<int,int>,double>> targets;
    for(auto& bin_heightss : bin_heightss_reps)
    {
        targets.emplace_back();
        for(auto& bin_heights : bin_heightss) 
        {
            int VPB=VPBs[&bin_heights-bin_heightss.data()];
            for(int j=0; j<length; ++j)
                if(std::abs(j-VPB)>VPE)
                {
                    std::pair<int,int> tmp(std::min(VPB,j),std::max(VPB,j));
                    if(!targets.back().insert(std::make_pair(tmp,bin_heights[j])).second)
                        targets.back()[tmp]=(targets.back()[tmp]+bin_heights[j])/2;
                }
        }
    }
    return std::make_tuple(VPBs,targets);
}

std::tuple<std::vector<double>,std::vector<double>> calculate_stall(std::vector<std::pair<std::pair<int,int>,std::pair<std::string,std::string>>>& mut_list, std::vector<std::string> CTCF_peaks, std::string chromosome, int start, int end, int length)
{
    std::sort(mut_list.begin(), mut_list.end());
    int x=start;
    if(!mut_list.empty())
        x=std::min(x,mut_list.front().first.first-1);
    int y=x;
    std::map<int,int> mut2wt;
    auto iter=mut_list.begin();
    while(y<end)
    {
        if(iter==mut_list.end())
        {
            if(y>=start) mut2wt[x]=y;
            ++x;
            ++y;
        }
        else
        {
            if(iter->second.first=="DL")
            {
                if(y!=iter->first.first)
                {
                    if(y>=start) mut2wt[x]=y;
                    ++x;
                    ++y;
                }
                else
                {
                    y=iter->first.second;
                    ++iter;
                }
            }
            else
            {
                if(y!=iter->first.first)
                {
                    if(y>=start) mut2wt[x]=y;
                    ++x;
                    ++y;
                }
                else
                {
                    for(int k=0; k<iter->first.second; ++k)
                    {
                        if(y>=start) mut2wt[x]=y;
                        ++x;
                    }
                    ++iter;
                }
            }
        }
    }
    int xUB=x, xLB=mut2wt.begin()->first;
    std::vector<double> CTCF_enrich(length,0.0);
    std::vector<bool> CF(length,false);
    std::vector<bool> CR(length,false);
    std::string str_tmp;
    for(auto& CTCF_peak : CTCF_peaks)
    {
        std::ifstream fin(CTCF_peak);
        fin >> str_tmp;
        while(fin.good())
        {
            if(str_tmp!=chromosome)
            {
                getline(fin,str_tmp);
                fin >> str_tmp;
                continue;
            }
            fin >> str_tmp;
            int tmp=stoi(str_tmp);
            fin >> str_tmp;
            tmp=(tmp+stoi(str_tmp))/2;
            if(tmp<xLB || tmp>=xUB)
            {
                getline(fin,str_tmp);
                fin >> str_tmp;
            }
            else
            {
                int index=(mut2wt[tmp]-start)/bin;
                fin >> str_tmp >> str_tmp >> str_tmp;
                if(str_tmp=="+")
                    CF[index]=true;
                else if(str_tmp=="-")
                    CR[index]=true;
                fin >> str_tmp;
                CTCF_enrich[index]+=stod(str_tmp);
                getline(fin,str_tmp);
                fin >> str_tmp;
            }
        }
        fin.close();
    }
    for(auto& enrich : CTCF_enrich) enrich/=CTCF_peaks.size();
        
    std::vector<double> left_stall(length,0.0);
    std::vector<double> right_stall(length,0.0);
    std::map<char,int> nt2int{{'A',0},{'C',1},{'G',2},{'T',3}};
    for(int i=0; i<length; ++i)
    {
        if(CTCF_enrich[i]>0)
        {
            if(CF[i] && CR[i])
            {
                left_stall[i]=stall_pro(CTCF_enrich[i]/2);
                right_stall[i]=left_stall[i];
            }
            else if(CF[i])
                left_stall[i]=stall_pro(CTCF_enrich[i]);
            else if(CR[i])
                right_stall[i]=stall_pro(CTCF_enrich[i]);
        }
    }
    return std::make_tuple(left_stall,right_stall);
}

std::vector<double> calculate_birth(std::vector<std::pair<std::pair<int,int>,std::pair<std::string,std::string>>>& mut_list, std::vector<std::string> nipbls, std::string chromosome, int start, int end, int length)
{
    std::vector<double> nipbl_signal(length,0.0);
    std::string str_tmp;
    for(auto& nipbl : nipbls)
    {
        std::ifstream fin(nipbl);
        getline(fin,str_tmp);
        fin >> str_tmp;
        while(fin.good())
        {
            if(str_tmp!=chromosome)
            {
                getline(fin,str_tmp);
                fin >> str_tmp;
                continue;
            }
            fin >> str_tmp;
            int PS=std::max(start,stoi(str_tmp));
            fin >> str_tmp;
            int PE=std::min(end,stoi(str_tmp));
            if(PS>=PE)
            {
                getline(fin,str_tmp);
                fin >> str_tmp;
            }
            else
            {
                fin >> str_tmp;
                double cover=stod(str_tmp);
                int sindex=(PS-start)/bin;
                int eindex=(PE-1-start)/bin;
                for(int i=sindex; i<=eindex; ++i)
                {
                    int IS=std::max(start+i*bin,PS);
                    int IE=std::min(start+(i+1)*bin,PE);
                    int width=IE-IS;
                    for(auto& mut : mut_list)
                        if(mut.second.first=="DL")
                            width-=std::max(0,std::min(mut.first.second,IE)-std::max(mut.first.first,IS));
                    nipbl_signal[i]+=width*cover;
                }
                fin >> str_tmp;
            }
        }
        fin.close();
    }
    double statis_birth=std::accumulate(nipbl_signal.begin(),nipbl_signal.end(),0.0)*(1-nipbl_ratio)/nipbl_ratio/length;
    if(statis_birth==0) statis_birth=1;
    for(auto& signal : nipbl_signal) signal+=statis_birth;
    return nipbl_signal;
}
