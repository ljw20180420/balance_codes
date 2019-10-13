clc
clear

exc=10;
heapa='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/get_heat_map/build/get_heat_map';
parent_dic='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/hulu/';
runs={'run1','run2','run3','run4','run5','run6','run7'};

tmp=dir(fullfile(parent_dic,runs{1}));
sr_now={tmp.name};
sr_now=sr_now(~cell2mat({tmp.isdir}));
sr_now(~startsWith(sr_now,'sr_now'))=[];
para=split(sr_now,'_');

density=str2double(para{8});

capture_ratio=[];
fid=fopen(fullfile(parent_dic,'pos_file'),'r');
line=fgetl(fid);
while ~feof(fid)
    if strcmp(line,'capture_ratio')
        while ~feof(fid)
            line=fgetl(fid);
            tmp=str2double(line);
            if(isnan(tmp))
                break;
            else
                capture_ratio=[capture_ratio,tmp];
            end
        end
    else
        line=fgetl(fid);
    end
end
fclose(fid);

sr_now_list=join(fullfile(parent_dic,runs,sr_now),',');
sr_now_list=sr_now_list{1};
command_line=sprintf('%s -sr_nows %s -capture_ratios %f -out_path %s -density %f &\nwait\n',heapa,sr_now_list,capture_ratio,parent_dic,density);
system(command_line);
movefile(fullfile(parent_dic,['hic_',num2str(capture_ratio)]),fullfile(parent_dic,'hic'));

parent_dic='~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/hulu';
shorten={'hulu'};
data_name={'sim'};
uplim=0.02;

[~,~,vps,XTLE,XTL,XT]=load_target_CTCF(fullfile(save_path,'pos_file'));

hic=load(fullfile(save_path,'hic'));
hic=max(0,hic-median(hic(:)));

figure('position',[100 100 1200 1000],'NumberTitle','off', 'Name',['sim_hulu_',num2str(vps(1:3))]);
for vv=1:3
    nor=sum(hic(vps(vv)+1,[1:vps(vv)+1-exc,vps(vv)+1+exc:end]));
    subplot(3,1,vv)
    bar(hic(vps(vv)+1,:)/nor,1)
    ax=gca;
    ax.XTickLabel=XTL;
    ax.XTick=XT;
    ylim([0,uplim]);
end
print('-painters',fullfile(parent_dic,'figures',['sim_hulu_',num2str(vps(1:3))]),'-djpeg');
print('-painters',fullfile(parent_dic,'figures',['sim_hulu_',num2str(vps(1:3))]),'-depsc');

figure('position',[100 100 1200 1000],'NumberTitle','off', 'Name',['sim_hulu_',num2str(vps(4:6))])
for vv=4:6
    nor=sum(hic(vps(vv)+1,[1:vps(vv)+1-exc,vps(vv)+1+exc:end]));
    subplot(3,1,vv)
    bar(hic(vps(vv)+1,:)/nor,1)
    ax=gca;
    ax.XTickLabel=XTL;
    ax.XTick=XT;
    ylim([0,uplim]);
end
print('-painters',fullfile(parent_dic,'figures',['sim_hulu_',num2str(vps(4:6))]),'-djpeg');
print('-painters',fullfile(parent_dic,'figures',['sim_hulu_',num2str(vps(4:6))]),'-depsc');