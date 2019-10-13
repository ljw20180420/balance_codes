clc
clear

exc=10;
var_min=1e-2;
heapa='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/get_heat_map/build/get_heat_map';
parent_dic='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/';
wild='HEC-1-B_WT';
expers={'HEC-1-B_WT','HEC-1-B_F','HEC-1-B_FF','HEC-1-B_FFRR','HEC-1-B_R1','HEC-1-B_R2','HEC-1-B_RF','HEC-1-B_RRFF'};
runs={'run1','run2','run3','run4','run5','run6','run7'};

tmp=dir(fullfile(parent_dic,wild,runs{1}));
sr_now={tmp.name};
sr_now=sr_now(~cell2mat({tmp.isdir}));
sr_now(~startsWith(sr_now,'sr_now'))=[];
para=split(sr_now,'_');

chain_length=str2double(para{12});
density=str2double(para{8});
len=(chain_length/density)^(1/3);

[rcmv,capture_ratio]=load_target_CTCF(fullfile(parent_dic,wild,'pos_file'));
rcmv(rcmv(:,1)<=exc | rcmv(:,2)>=chain_length-exc,:)=[];

command_line='';
for ad=1:length(expers)
	sr_now_list=join(fullfile(parent_dic,expers{ad},runs,sr_now),',');
    sr_now_list=sr_now_list{1};
    command_line=[command_line,sprintf('%s -sr_nows %s -capture_ratios %f -out_path %s -density %f &\n',heapa,sr_now_list,capture_ratio,fullfile(parent_dic,expers{ad}),density)];
end
command_line=[command_line,sprintf('wait\n')];
system(command_line);
for ad=1:length(expers)
    movefile(fullfile(parent_dic,expers{ad},['hic_',num2str(capture_ratio)]),fullfile(parent_dic,expers{ad},'hic'));
end
    
for ad=1:length(expers)
    intersss=cell(length(runs),1);
    for rr=1:length(runs)
        positionss=load(fullfile(parent_dic,expers{ad},runs{rr},sr_now{1}));
        intersss{rr}=cell(size(positionss,1),1);
        for po=1:size(positionss,1)
            tmp=reshape(positionss(po,:),3,[]);
            intersss{rr}{po}=find(vecnorm(mod(tmp(:,rcmv(:,1)+1)-tmp(:,rcmv(:,2)+1)+len/2,len)-len/2)<capture_ratio);
        end
    end
    
    if ~strcmp(wild,expers{ad})
        lambda=load(fullfile(parent_dic,wild,'lambda'));
    else
        hic=load(fullfile(parent_dic,wild,'hic'));
        sims=hic(rcmv(:,1)+1+rcmv(:,2)*chain_length);
        medval=median(hic(:));
        
        tmp=max(0,sims-medval);
        error_fun=@(x) norm(tmp/x-rcmv(:,3),4);
        coff=fminunc(error_fun,tmp'*tmp/(tmp'*rcmv(:,3)));
        rcmv(:,3)=min(rcmv(:,3)*coff+min(sims,medval),1);
        rcmv(:,4)=max(var_min,rcmv(:,4)*coff^2);
        options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective');
        options.FunctionTolerance=0;
        lambda=fminunc(@(lambda) cal_value(lambda, intersss, rcmv),zeros(size(rcmv,1),1),options);
        save(fullfile(parent_dic,wild,'lambda'), 'lambda', '-ascii', '-double', '-tabs');
    end
    
    weights=cal_weights(lambda, intersss);
    save(fullfile(parent_dic,expers{ad},'weights'), 'weights', '-ascii', '-double', '-tabs');
    Kish_effective_sample_size=1/sum(weights.^2)
end

command_line='';
for ad=1:length(expers)
    sr_now_list=join(fullfile(parent_dic,expers{ad},runs,sr_now),',');
    sr_now_list=sr_now_list{1};
    command_line=[command_line,sprintf('%s -sr_nows %s -capture_ratios %f -out_path %s -density %f -weights %s &\n',heapa,sr_now_list,capture_ratio,fullfile(parent_dic,expers{ad}),density,fullfile(parent_dic,expers{ad},'weights'))];
end
command_line=[command_line,sprintf('wait\n')];
system(command_line);
for ad=1:length(expers)
    movefile(fullfile(parent_dic,expers{ad},['hic_',num2str(capture_ratio)]),fullfile(parent_dic,expers{ad},'hic_max_entropy'));
end

function [value,gradi,hessi]=cal_value(lambda, intersss, rcmv)
st=size(rcmv,1);
cum_record_1D=0;
for rr=1:length(intersss)
    cum_record_1D=cum_record_1D+length(intersss{rr});
end
weights=zeros(cum_record_1D,1);
we=0;
for rr=1:length(intersss)
    for po=1:length(intersss{rr})
        we=we+1;
        weights(we)=-sum(lambda(intersss{rr}{po}));
    end
end
max_weights=max(weights);
weights=exp(weights-max(weights));
value=log(sum(weights))+max_weights+lambda'*rcmv(:,3)+0.5*(lambda').^2*rcmv(:,4);
if nargout > 1
    weights=weights/sum(weights);
    gradi=zeros(st,1);
    if nargout > 2
        hessi=zeros(st);
    end
    we=0;
    for rr=1:length(intersss)
        for po=1:length(intersss{rr})
            we=we+1;
            gradi(intersss{rr}{po})=gradi(intersss{rr}{po})-weights(we);
            if nargout > 2
                hessi(intersss{rr}{po},intersss{rr}{po})=hessi(intersss{rr}{po},intersss{rr}{po})+weights(po);
            end
        end
    end
    if nargout > 2
        hessi(1:st+1:end)=hessi(1:st+1:end)'+rcmv(:,4);
        hessi=hessi-gradi*gradi';
    end
    gradi=gradi+rcmv(:,3)+lambda.*rcmv(:,4);
end
end

function weights=cal_weights(lambda, intersss)
cum_record_1D=0;
for rr=1:length(intersss)
    cum_record_1D=cum_record_1D+length(intersss{rr});
end
weights=zeros(cum_record_1D,1);
we=0;
for rr=1:length(intersss)
    for po=1:length(intersss{rr})
        we=we+1;
        weights(we)=-sum(lambda(intersss{rr}{po}));
    end
end
weights=exp(weights-max(weights));
weights=weights/sum(weights);
end