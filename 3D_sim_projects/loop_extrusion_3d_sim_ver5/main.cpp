#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstring>
#include <queue>
#include <mutex>
#include <future>
#include <condition_variable>
#include <functional>
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

const int DEFAULT_CHAIN_LENGTH=1000;

const double kB=1.380649e-23;
const double T=300;
const double sigma=1e-9;
const double NA=6.02214076e23;
const double amu=1e-3/NA;
const double m=100*amu;
const double gam=1e10*m;
const double timestep=8e-14;
const double timeunit=std::sqrt(m*std::pow(sigma,2.0)/kB/T);
const double alpha=sigma/std::sqrt(m*kB*T);
const double sdu=timestep/timeunit;
const double ag=alpha*gam;
const double a=(2-ag*sdu)/(2+ag*sdu);
const double tb=std::sqrt(2*ag*sdu);
const double c=2*sdu/(2+ag*sdu);

const double REPe=1.5;
const double REPrmin=std::sqrt(6.0/7.0);
const double REPsigma=1.05;
const double REPemin=46656.0/823543.0;
const double coff_14=14*REPe/REPemin*std::pow(REPrmin/REPsigma,14.0);
const double coff_12=12*REPe/REPemin*std::pow(REPrmin/REPsigma,12.0);

template<typename T>
class thread_safe_queue
{
private:
    mutable std::mutex mut;
    std::queue<T> data_queue;
    std::condition_variable data_cond;
public:
    thread_safe_queue() {}
    void push(T new_value)
    {
        std::lock_guard<std::mutex> lk(mut);
        data_queue.push(std::move(new_value));
        data_cond.notify_one();
    }
    bool try_pop(T& value)
    {
        std::lock_guard<std::mutex> lk(mut);
        if(data_queue.empty())
            return false;
        value=std::move(data_queue.front());
        data_queue.pop();
        return true;
    }
    bool empty() const
    {
        std::lock_guard<std::mutex> lk(mut);
        return data_queue.empty();
    }
};

class join_threads
{
    std::vector<std::thread>& threads;
public:
    explicit join_threads(std::vector<std::thread>& threads_) : threads(threads_) {}
    ~join_threads()
    {
        for(unsigned long i=0;i<threads.size();++i)
            if(threads[i].joinable())
                threads[i].join();
    }
};

class function_wrapper
{
    struct impl_base 
    {
        virtual void call()=0;
        virtual ~impl_base() {}
    };
    std::unique_ptr<impl_base> impl;
    template<typename F>
    struct impl_type: impl_base
    {
        F f;
        impl_type(F&& f_): f(std::move(f_)) {}
        void call() { f(); }
    };
public:
    template<typename F>
    function_wrapper(F&& f):
        impl(new impl_type<F>(std::move(f))) {}
    void operator()() { impl->call(); }
    function_wrapper() = default;
    function_wrapper(function_wrapper&& other):
        impl(std::move(other.impl)) {}
    function_wrapper& operator=(function_wrapper&& other)
    {
        impl=std::move(other.impl);
        return *this;
    }
    function_wrapper(const function_wrapper&)=delete;
    function_wrapper(function_wrapper&)=delete;
    function_wrapper& operator=(const function_wrapper&)=delete;
};

class thread_pool
{
    std::atomic_bool done;
    thread_safe_queue<function_wrapper> work_queue;
    std::vector<std::thread> threads;
    join_threads joiner;
    void worker_thread()
    {
        while(!done)
        {
            function_wrapper task;
            if(work_queue.try_pop(task))
                task();
            else
                std::this_thread::yield();
        }
    }
public:
    thread_pool(int THR_MAX=std::thread::hardware_concurrency()):
        done(false),joiner(threads)
    {
        try
        {
            for(unsigned i=0;i<THR_MAX;++i)
                threads.push_back(std::thread(&thread_pool::worker_thread,this));
        }
        catch(...)
        {
            done=true;
            throw;
        }
    }
    ~thread_pool()
    {
        done=true;
    }
    template<typename FunctionType>
    std::future<typename std::result_of<FunctionType()>::type>
        submit(FunctionType f)
    {
        typedef typename std::result_of<FunctionType()>::type result_type;
        std::packaged_task<result_type()> task(std::move(f));
        std::future<result_type> res(task.get_future());
        work_queue.push(std::move(task));
        return res;
    }
};

struct Command_content
{
    std::string pos_file;
    int THR_MAX;
    Command_content(int argc, char **argv)
    {
        THR_MAX=1;
        for(int i=1; i<argc-1; ++i)
        {
            if(!strcmp(argv[i],"-pos_file"))
                pos_file=argv[i+1];
            if(!strcmp(argv[i],"-THR_MAX"))
                THR_MAX=stoi(std::string(argv[i+1]));
        }
    }
};

struct Atom
{
    std::array<double,3> position;
    std::array<double,3> velocity;
    std::array<double,3> force;
    std::array<double,3> acum_shift;
};

class Lattice
{
protected:
    std::vector<std::pair<int,int>> cohesin;
    std::vector<bool> left_occupied;
    std::vector<bool> right_occupied;
    
    double processivity;
    double separation;
    std::vector<double> left_stall;
    std::vector<double> right_stall;
    
    std::discrete_distribution<int> birth_distribution;
    std::uniform_real_distribution<double> uniform_distribution{0.0,1.0};
    
public:
    void lattice_initialize(double processivity_, double separation_, std::vector<double>& birth_, std::vector<double>& left_stall_, std::vector<double>& right_stall_, std::default_random_engine& generator)
    {
        processivity=processivity_;
        separation=separation_;
        birth_distribution=std::discrete_distribution<int>(birth_.begin(),birth_.end());
        left_stall=left_stall_;
        right_stall=right_stall_;
        left_occupied.assign(birth_.size(),false);
        right_occupied.assign(birth_.size(),false);
        cohesin.resize(int(birth_.size()/separation));
        
        for(int i=0; i<cohesin.size(); ++i)
            lattice_supply(generator, i);
    }
    
    void lattice_step(std::default_random_engine& generator)
    {
        for(int i=0; i<cohesin.size(); ++i)
        {
            if(uniform_distribution(generator)<2/processivity)
            {
                left_occupied[cohesin[i].first]=false;
                right_occupied[cohesin[i].second]=false;
                lattice_supply(generator, i);
            }
        }
        for(int i=0; i<cohesin.size(); ++i)
        {
            if(cohesin[i].first>0 && !left_occupied[cohesin[i].first-1] && !right_occupied[cohesin[i].first-1])
                if(uniform_distribution(generator)>left_stall[cohesin[i].first])
                {
                    left_occupied[cohesin[i].first]=false;
                    left_occupied[--cohesin[i].first]=true;
                }
            if(cohesin[i].second<right_occupied.size()-1 && !left_occupied[cohesin[i].second+1] && !right_occupied[cohesin[i].second+1])
                if(uniform_distribution(generator)>right_stall[cohesin[i].second])
                {
                    right_occupied[cohesin[i].second]=false;
                    right_occupied[++cohesin[i].second]=true;
                }
        }
    }
    
    void lattice_supply(std::default_random_engine& generator, int i)
    {
        while(true)
        {
            int tmp=birth_distribution(generator);
            if(!left_occupied[tmp] && !right_occupied[tmp])
            {
                if(uniform_distribution(generator)<0.5)
                {
                    cohesin[i].first=cohesin[i].second=tmp;
                    left_occupied[tmp]=right_occupied[tmp]=true;
                    return;
                }
                else if(tmp<left_occupied.size()-1 && !left_occupied[tmp+1] && !right_occupied[tmp+1])
                {
                    cohesin[i].first=tmp;
                    cohesin[i].second=tmp+1;
                    left_occupied[tmp]=right_occupied[tmp+1]=true;
                    return;
                }
            }
        }
    }
};

class Polymer : Lattice
{
    std::vector<Atom> chain;
    
    double twice_harmonic_bond;
    double twice_harmonic_cohesin;
    double distance_bond;
    double distance_cohesin;
    double stiffness;
    double density;
    
    int D3_per_D1;
    
    double cut_off;
    double nei_off;
    double nei_off2;
    double shell;
    
    std::array<double,3> lengths;
    std::array<double,3> units;
    std::vector<std::vector<std::vector<int>>> axis2index;
    
    bool rearrange;
    
    std::vector<std::pair<std::vector<std::list<Atom*>*>,std::list<Atom*>>> cell_atom_ptrs;
    std::vector<std::pair<Atom*,Atom*>> neigh_atom_pairs;
    
    std::normal_distribution<double> norm_distribution;
public:
    
    void initialize(std::vector<std::array<double,3>>& velocity_, double harmonic_bond_, double harmonic_cohesin_, double distance_bond_, double distance_cohesin_, double stiffness_, double density_, double processivity_, double separation_, std::vector<double>& birth_, std::vector<double>& left_stall_, std::vector<double>& right_stall_, int D3_per_D1_, double cut_off_, double nei_off_, std::default_random_engine& generator)
    {
        chain.resize(velocity_.size());
        for(int i=0; i<chain.size(); ++i)
            chain[i].velocity=velocity_[i];
        
        twice_harmonic_bond=2*harmonic_bond_;
        twice_harmonic_cohesin=2*harmonic_cohesin_;
        distance_bond=distance_bond_;
        distance_cohesin=distance_cohesin_;
        stiffness=stiffness_;
        density=density_;
        lattice_initialize(processivity_, separation_, birth_, left_stall_, right_stall_, generator);
        D3_per_D1=D3_per_D1_;
        cut_off=cut_off_;
        nei_off=nei_off_;
        nei_off2=nei_off*nei_off;
        shell=nei_off-cut_off;
        
        if(chain.size()%2==1)
            throw std::exception();
        lengths[0]=lengths[1]=lengths[2]=std::pow(chain.size()/density,1.0/3.0);
        std::array<int,3> THR_NUMs{int(lengths[0]/nei_off),int(lengths[1]/nei_off),int(lengths[2]/nei_off)};
        for(int i=0; i<3; ++i)
            units[i]=lengths[i]/THR_NUMs[i];
        cell_atom_ptrs.assign(THR_NUMs[0]*THR_NUMs[1]*THR_NUMs[2],std::make_pair(std::vector<std::list<Atom*>*>(),std::list<Atom*>()));
        axis2index.assign(THR_NUMs[0],std::vector<std::vector<int>>(THR_NUMs[1],std::vector<int>(THR_NUMs[2])));
        int index=0;
        for(int k=0; k<THR_NUMs[2]; ++k)
            for(int j=0; j<THR_NUMs[1]; ++j)
                for(int i=0; i<THR_NUMs[0]; ++i)
                {
                    axis2index[i][j][k]=index;
                    ++index;
                }
        for(int k=0; k<THR_NUMs[2]; ++k)
            for(int j=0; j<THR_NUMs[1]; ++j)
                for(int i=0; i<THR_NUMs[0]; ++i)
                {
                    int index_me=axis2index[i][j][k];
                    for(int ds=-1; ds<=1; ++ds)
                        for(int dc=-1; dc<=1; ++dc)
                            for(int dr=-1; dr<=1; ++dr)
                            {
                                int index_neigh=axis2index[(i+dr+THR_NUMs[0])%THR_NUMs[0]][(j+dc+THR_NUMs[1])%THR_NUMs[1]][(k+ds+THR_NUMs[2])%THR_NUMs[2]];
                                if(index_me<index_neigh)
                                    cell_atom_ptrs[index_me].first.push_back(&cell_atom_ptrs[index_neigh].second);
                            }
                }
        
        int size=int(lengths[0])-2;
        int t = size/2;
        std::vector<std::array<int,3>> a;
        for(int i=1; i<size; ++i)
            a.push_back(std::array<int,3>{t,t,i});
        for(int i=size-1; i>0; --i)
            a.push_back(std::array<int,3>{t,t-1,i});
        std::vector<std::vector<std::vector<bool>>> b(size+1,std::vector<std::vector<bool>>(size+1,std::vector<bool>(size+1,false)));
        for(auto& pos : a)
            b[pos[0]][pos[1]][pos[2]]=true;
        std::uniform_int_distribution<int> uniform_int2(1,2);
        int in=(chain.size()-a.size())/2;
        for(int i=0; i<in; ++i)
        {
            std::uniform_int_distribution<int> uniform_int1(0,a.size()-1);
            while(true)
            {
                int t=uniform_int1(generator);
                std::array<int,3> t0=a[t], t1=(t!=a.size()-1)?a[t+1]:a[0], c{std::abs(t0[0]-t1[0]),std::abs(t0[1]-t1[1]),std::abs(t0[2]-t1[2])};
                int direction=(std::max_element(c.begin(),c.end())-c.begin()+uniform_int2(generator))%3;
                int shift=(uniform_distribution(generator)>0.5)?1:-1;
                t0[direction]+=shift;
                t1[direction]+=shift;
                if(!b[t0[0]][t0[1]][t0[2]] && !b[t1[0]][t1[1]][t1[2]] && t0[direction]>0 && t1[direction]>0 && t0[direction]<size && t1[direction]<size)
                {
                    a.insert(a.begin()+t+1,t0);
                    a.insert(a.begin()+t+2,t1);
                    b[t0[0]][t0[1]][t0[2]]=true;
                    b[t1[0]][t1[1]][t1[2]]=true;
                    break;
                }
            }
        }
        for(int i=0; i<chain.size(); ++i)
            chain[i].position={double(a[i][0]),double(a[i][1]),double(a[i][2])};
        
        setup_pair();
        
        cal_force_up_velocity_position(generator);
    }
    
    void setup_pair()
    {
        for(auto& head_list : cell_atom_ptrs)
            head_list.second.clear();
        
        for(auto& atom : chain)
        {
            atom.acum_shift={0.0,0.0,0.0};
            for(int dim=0; dim<3; ++dim)
            {
                atom.position[dim]=std::fmod(atom.position[dim],lengths[dim]);
                if(atom.position[dim]<0)
                    atom.position[dim]+=lengths[dim];
            }
            int index=axis2index[int(atom.position[0]/units[0])][int(atom.position[1]/units[1])][int(atom.position[2]/units[2])];
            cell_atom_ptrs[index].second.push_back(&atom);
        }
        
        neigh_atom_pairs.clear();
        for(auto& head_list : cell_atom_ptrs)
        {
            for(auto iter1=head_list.second.begin(); iter1!=head_list.second.end(); ++iter1)
                for(auto iter2=std::next(iter1); iter2!=head_list.second.end(); ++iter2)
                {
                    std::array<double,3> vec=d_axis((*iter1)->position, (*iter2)->position);
                    if(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]<nei_off2)
                        neigh_atom_pairs.emplace_back(*iter1,*iter2);
                }
            for(auto& ptr : head_list.first)
                for(auto iter1=head_list.second.begin(); iter1!=head_list.second.end(); ++iter1)
                    for(auto iter2=ptr->begin(); iter2!=ptr->end(); ++iter2)
                    {
                        std::array<double,3> vec=d_axis((*iter1)->position, (*iter2)->position);
                        if(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]<nei_off2)
                            neigh_atom_pairs.emplace_back(*iter1,*iter2);
                    }
        }
    }
    
    void cal_force_up_velocity_position(std::default_random_engine& generator)
    {
        for(auto& atom : chain)
            atom.force={0.0,0.0,0.0};
        
        for(auto& pair : neigh_atom_pairs)
        {
            std::array<double,3> rv=d_axis(pair.second->position,pair.first->position);
            double r=std::sqrt(rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2]);
            if(r<cut_off)
            {
                double tmp=(coff_14*std::pow(r,12.0)-coff_12*std::pow(r,10.0));
                for(int j=0; j<3; ++j)
                {
                    double ft=tmp*rv[j];
                    pair.first->force[j]+=ft;
                    pair.second->force[j]-=ft;
                }
            }
        }
        
        for(auto& pair : cohesin)
            if(pair.first!=pair.second)
            {
                std::array<double,3> rv=d_axis(chain[pair.second].position,chain[pair.first].position);
                double r=std::sqrt(rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2]);
                double tmp=twice_harmonic_cohesin*(r-distance_cohesin)/r;
                for(int j=0; j<3; ++j)
                {
                    double ft=rv[j]*tmp;
                    chain[pair.first].force[j]+=ft;
                    chain[pair.second].force[j]-=ft;
                }
            }
        
        std::vector<std::array<double,3>> vec;
        std::vector<double> r;
        for(int i=0; i<chain.size()-1; ++i)
        {
            vec.push_back(d_axis(chain[i+1].position,chain[i].position));
            r.push_back(std::sqrt(vec.back()[0]*vec.back()[0]+vec.back()[1]*vec.back()[1]+vec.back()[2]*vec.back()[2]));
            double tmp=twice_harmonic_bond*(r.back()-distance_bond)/r.back();
            for(int j=0; j<3; ++j)
            {
                double ft=vec.back()[j]*tmp;
                chain[i].force[j]+=ft;
                chain[i+1].force[j]-=ft;
            }
        }
        for(int i=0; i<chain.size()-2; ++i)
        {
            double inner=vec[i][0]*vec[i+1][0]+vec[i][1]*vec[i+1][1]+vec[i][2]*vec[i+1][2];
            double tmp1=stiffness/(r[i]*r[i+1]);
            double tmp2=tmp1*inner/(r[i]*r[i]);
            for(int j=0; j<3; ++j)
            {
                double ft=tmp2*vec[i][j]-tmp1*vec[i+1][j];
                chain[i].force[j]+=ft;
                chain[i+1].force[j]-=ft;
            }
            tmp2=tmp1*inner/(r[i+1]*r[i+1]);
            for(int j=0; j<3; ++j)
            {
                double ft=tmp1*vec[i][j]-tmp2*vec[i+1][j];
                chain[i+2].force[j]+=ft;
                chain[i+1].force[j]-=ft;
            }
        }
        
        double max_dis2=0, next_dis2=0;
        for(auto& atom : chain)
        {
            for(int i=0; i<3; ++i)
            {
                atom.velocity[i]=a*atom.velocity[i]+sdu*atom.force[i]+tb*norm_distribution(generator);
                atom.position[i]+=c*atom.velocity[i];
                atom.acum_shift[i]+=c*atom.velocity[i];
            }
            double dis2=atom.acum_shift[0]*atom.acum_shift[0]+atom.acum_shift[1]*atom.acum_shift[1]+atom.acum_shift[2]*atom.acum_shift[2];
            if(dis2>next_dis2)
            {
                if(dis2>max_dis2)
                {
                    next_dis2=max_dis2;
                    max_dis2=dis2;
                }
                else
                    next_dis2=dis2;
            }
        }
        rearrange=std::sqrt(max_dis2)+std::sqrt(next_dis2)>shell;
    }
    
    void D3_step(std::default_random_engine& generator)
    {
        if(rearrange)
            setup_pair();
        
        cal_force_up_velocity_position(generator);
    }
    
    void simulation(int pre_1D, int anneal_1D, int per_1D, int record_1D, std::default_random_engine& generator)
    {
        for(int D1=0; D1<pre_1D; ++D1)
            lattice_step(generator);
        
        for(int D1=0; D1<anneal_1D; ++D1)
        {
            lattice_step(generator);
            for(int D3=0; D3<D3_per_D1; ++D3)
                D3_step(generator);
        }
        
        std::ostringstream streamObj;
        streamObj.unsetf(std::ios::floatfield);
        streamObj << "sr_now_" << twice_harmonic_bond/2 << "_" << twice_harmonic_cohesin/2 << "_" << distance_bond << "_" << distance_cohesin << "_" << stiffness << "_" << density << "_" << processivity << "_" << separation << "_" << D3_per_D1 << "_" << chain.size();
        std::ofstream sr(streamObj.str());
        
        streamObj.str("");
        streamObj << "co_now_" << twice_harmonic_bond/2 << "_" << twice_harmonic_cohesin/2 << "_" << distance_bond << "_" << distance_cohesin << "_" << stiffness << "_" << density << "_" << processivity << "_" << separation << "_" << D3_per_D1 << "_" << chain.size();
        std::ofstream co(streamObj.str());
        
        for(int rec=0; rec<record_1D; ++rec)
        {
            for(int D1=0; D1<per_1D; ++D1)
            {
                lattice_step(generator);
                for(int D3=0; D3<D3_per_D1; ++D3)
                    D3_step(generator);
            }
                
            for(int i=0; i<chain.size(); ++i)
            {
                for(int dim=0; dim<3; ++dim)
                {
                    chain[i].position[dim]=std::fmod(chain[i].position[dim],lengths[dim]);
                    if(chain[i].position[dim]<0)
                        chain[i].position[dim]+=lengths[dim];
                }
                sr << chain[i].position[0] << '\t' << chain[i].position[1] << '\t' << chain[i].position[2];
                if (i<chain.size()-1)
                    sr << '\t';
                else
                    sr << '\n';
            }
            
            for(int i=0; i<cohesin.size(); ++i)
            {
                co << cohesin[i].first << '\t' << cohesin[i].second;
                if (i<cohesin.size()-1)
                    co << '\t';
                else
                    co << '\n';
            }
        }
        
        sr.close();
        co.close();
    }
    
    std::array<double,3> d_axis(std::array<double,3>& pos1, std::array<double,3>& pos2)
    {
        std::array<double,3> vec;
        for(int dim=0; dim<3; ++dim)
        {
            vec[dim]=std::fmod(pos1[dim]-pos2[dim]+lengths[dim]/2,lengths[dim]);
            if(vec[dim]<0)
                vec[dim]+=lengths[dim];
            vec[dim]-=lengths[dim]/2;
        }
        return vec;
    }
};

bool isNumberC(const std::string& s);
void readin(std::ifstream& fin, std::string& str_tmp, std::vector<std::array<double,3>>& container);
void readin(std::ifstream& fin, std::string& str_tmp, std::vector<double>& container);
void readin(std::ifstream& fin, std::string& str_tmp, std::vector<int>& container);
void readin(std::ifstream& fin, std::string& str_tmp, double& container);
void readin(std::ifstream& fin, std::string& str_tmp, int& container);
void importdata(std::string& file, std::vector<std::array<double,3>>& velocity, std::vector<double>& harmonic_bonds, std::vector<double>& harmonic_cohesins, std::vector<double>& distance_bonds, std::vector<double>& distance_cohesins, std::vector<double>& stiffnesses, std::vector<double>& densities, std::vector<double>& processivities, std::vector<double>& separations, std::vector<double>& birth, std::vector<double>& left_stall, std::vector<double>& right_stall, std::vector<int>& D3_per_D1s, double& cut_off, double& nei_off, int& pre_1D, int& anneal_1D, int& per_1D, int& record_1D);
void one_set_paramemter(double harmonic_bond, double harmonic_cohesin, double distance_bond, double distance_cohesin, double stiffness, double density, double processivity, double separation, int D3_per_D1, std::vector<std::array<double,3>>& velocity, std::vector<double>& birth, std::vector<double>& left_stall, std::vector<double>& right_stall, double cut_off, double nei_off, int pre_1D, int anneal_1D, int per_1D, int record_1D);

int main(int argc, char **argv) 
{
    Command_content cc(argc, argv);
    
    std::vector<std::array<double,3>> velocity;
    std::vector<double> harmonic_bonds;
    std::vector<double> harmonic_cohesins;
    std::vector<double> distance_bonds;
    std::vector<double> distance_cohesins;
    std::vector<double> stiffnesses;
    std::vector<double> densities;
    std::vector<double> processivities;
    std::vector<double> separations;
    std::vector<double> birth;
    std::vector<double> left_stall;
    std::vector<double> right_stall;
    std::vector<int> D3_per_D1s;
    double cut_off, nei_off;
    int pre_1D, anneal_1D, per_1D, record_1D;
    
    importdata(cc.pos_file, velocity, harmonic_bonds, harmonic_cohesins, distance_bonds, distance_cohesins, stiffnesses, densities, processivities, separations, birth, left_stall, right_stall, D3_per_D1s, cut_off, nei_off, pre_1D, anneal_1D, per_1D, record_1D);
    
    thread_pool threads(cc.THR_MAX);
    std::vector<std::future<void>> futures;
    for(double& harmonic_bond : harmonic_bonds)
        for(double& harmonic_cohesin : harmonic_cohesins)
            for(double& distance_bond : distance_bonds)
                for(double& distance_cohesin : distance_cohesins)
                    for(double& stiffness : stiffnesses)
                        for(double& density : densities)
                            for(double& processivity : processivities)
                                for(double& separation : separations)
                                    for(int& D3_per_D1 : D3_per_D1s)
                                        futures.push_back(threads.submit(std::bind(one_set_paramemter, harmonic_bond, harmonic_cohesin, distance_bond, distance_cohesin, stiffness, density, processivity, separation, D3_per_D1, std::ref(velocity), std::ref(birth), std::ref(left_stall), std::ref(right_stall), cut_off, nei_off, pre_1D, anneal_1D, per_1D, record_1D)));
    for(auto& future : futures) future.wait();
    return 0;
}

void one_set_paramemter(double harmonic_bond, double harmonic_cohesin, double distance_bond, double distance_cohesin, double stiffness, double density, double processivity, double separation, int D3_per_D1, std::vector<std::array<double,3>>& velocity, std::vector<double>& birth, std::vector<double>& left_stall, std::vector<double>& right_stall, double cut_off, double nei_off, int pre_1D, int anneal_1D, int per_1D, int record_1D)
{
    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    Polymer polymer;
    polymer.initialize(velocity, harmonic_bond, harmonic_cohesin, distance_bond, distance_cohesin, stiffness, density, processivity, separation, birth, left_stall, right_stall, D3_per_D1, cut_off, nei_off, generator);
    polymer.simulation(pre_1D, anneal_1D, per_1D, record_1D, generator);
}

void importdata(std::string& file, std::vector<std::array<double,3>>& velocity, std::vector<double>& harmonic_bonds, std::vector<double>& harmonic_cohesins, std::vector<double>& distance_bonds, std::vector<double>& distance_cohesins, std::vector<double>& stiffnesses, std::vector<double>& densities, std::vector<double>& processivities, std::vector<double>& separations, std::vector<double>& birth, std::vector<double>& left_stall, std::vector<double>& right_stall, std::vector<int>& D3_per_D1s, double& cut_off, double& nei_off, int& pre_1D, int& anneal_1D, int& per_1D, int& record_1D)
{
    pre_1D=1e6;
    anneal_1D=2e3;
    per_1D=1;
    record_1D=8e3;
    cut_off=REPsigma;
    nei_off=3*cut_off;
    
    std::ifstream fin(file);
    std::string str_tmp;
    fin >> str_tmp;
    while(fin.good())
    {
        if(str_tmp=="velocity")
            readin(fin, str_tmp, velocity);
        else if(str_tmp=="harmonic_bond")
            readin(fin, str_tmp, harmonic_bonds);
        else if(str_tmp=="harmonic_cohesin")
            readin(fin, str_tmp, harmonic_cohesins);
        else if(str_tmp=="distance_bond")
            readin(fin, str_tmp, distance_bonds);
        else if(str_tmp=="distance_cohesin")
            readin(fin, str_tmp, distance_cohesins);
        else if(str_tmp=="stiffness")
            readin(fin, str_tmp, stiffnesses);
        else if(str_tmp=="density")
            readin(fin, str_tmp, densities);
        else if(str_tmp=="processivity")
            readin(fin, str_tmp, processivities);
        else if(str_tmp=="separation")
            readin(fin, str_tmp, separations);
        else if(str_tmp=="birth")
            readin(fin, str_tmp, birth);
        else if(str_tmp=="left_stall")
            readin(fin, str_tmp, left_stall);
        else if(str_tmp=="right_stall")
            readin(fin, str_tmp, right_stall);
        else if(str_tmp=="D3_per_D1")
            readin(fin, str_tmp, D3_per_D1s);
        else if(str_tmp=="cut_off")
            readin(fin, str_tmp, cut_off);
        else if(str_tmp=="nei_off")
            readin(fin, str_tmp, nei_off);
        else if(str_tmp=="pre_1D")
            readin(fin, str_tmp, pre_1D);
        else if(str_tmp=="anneal_1D")
            readin(fin, str_tmp, anneal_1D);
        else if(str_tmp=="per_1D")
            readin(fin, str_tmp, per_1D);
        else if(str_tmp=="record_1D")
            readin(fin, str_tmp, record_1D);
        else
            fin >> str_tmp;
    }
    fin.close();
    
    int chain_length=DEFAULT_CHAIN_LENGTH;
    if(!velocity.empty())
        chain_length=velocity.size();
    else if(!birth.empty())
        chain_length=birth.size();
    else if(!left_stall.empty())
        chain_length=left_stall.size();
    else if(!right_stall.empty())
        chain_length=right_stall.size();
    
    if(velocity.empty())
        velocity=std::vector<std::array<double,3>>(chain_length,std::array<double,3>{0,0,0});
    if(harmonic_bonds.empty())
        harmonic_bonds.push_back(100);
    if(harmonic_cohesins.empty())
        harmonic_cohesins.push_back(25);
    if(distance_bonds.empty())
        distance_bonds.push_back(1);
    if(distance_cohesins.empty())
        distance_cohesins.push_back(0.5);
    if(stiffnesses.empty())
        stiffnesses.push_back(2);
    if(densities.empty())
        densities.push_back(0.2);
    if(processivities.empty())
        processivities.push_back(200);
    if(separations.empty())
        separations.push_back(200);
    if(birth.empty())
        birth=std::vector<double>(chain_length,1.0/double(chain_length));
    if(left_stall.empty())
        left_stall=std::vector<double>(chain_length,0.0);
    if(right_stall.empty())
        right_stall=std::vector<double>(chain_length,0.0);
    if(D3_per_D1s.empty())
        D3_per_D1s.push_back(1250);
}

bool isNumberC(const std::string& s)
{
    char* p;
    strtod(s.c_str(), &p);
    return *p == 0;
}

void readin(std::ifstream& fin, std::string& str_tmp, std::vector<std::array<double,3>>& container)
{
    fin >> str_tmp;
    while(fin.good() && isNumberC(str_tmp))
    {
        container.push_back(std::array<double,3>());
        for(int i=0; i<3; ++i)
        {
            container.back()[i]=stod(str_tmp);
            fin >> str_tmp;
        }
    }
}

void readin(std::ifstream& fin, std::string& str_tmp, std::vector<double>& container)
{
    fin >> str_tmp;
    while(fin.good() && isNumberC(str_tmp))
    {
        container.push_back(stod(str_tmp));
        fin >> str_tmp;
    }
}

void readin(std::ifstream& fin, std::string& str_tmp, std::vector<int>& container)
{
    fin >> str_tmp;
    while(fin.good() && isNumberC(str_tmp))
    {
        container.push_back(stoi(str_tmp));
        fin >> str_tmp;
    }
}

void readin(std::ifstream& fin, std::string& str_tmp, double& container)
{
    fin >> str_tmp;
    container=stod(str_tmp);
    fin >> str_tmp;
}

void readin(std::ifstream& fin, std::string& str_tmp, int& container)
{
    fin >> str_tmp;
    container=stoi(str_tmp);
    fin >> str_tmp;
}
