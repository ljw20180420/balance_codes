clc
clear
exc=10;
heapa='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/get_heat_map/build/get_heat_map';
parent_dic='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B';
fit_dic='HEC-1-B_WT_fit';
runs={'run1','run2'};
[rcmv,capture_ratios,vps,XTLE,XTL,XT]=load_target_CTCF(fullfile(parent_dic,fit_dic,'pos_file'));
chain_length=length(XTLE);
cr_list='';
for ca=1:length(capture_ratios)
    if ca<length(capture_ratios)
        cr_list=[cr_list,num2str(capture_ratios(ca)),','];
    else
        cr_list=[cr_list,num2str(capture_ratios(ca))];
    end
end
uplim=0.02:0.04:0.14;
tmp=dir(fullfile(parent_dic,fit_dic,runs{1}));
sr_nows={tmp.name};
sr_nows=sr_nows(~cell2mat({tmp.isdir}));
sr_nows(~startsWith(sr_nows,'sr_now'))=[];

% command_line='';
% for ss=1:length(sr_nows)
%     para=split(sr_nows{ss},'_');
%     tmp=cell(1,length(runs));
%     for rr=1:length(runs)
%         tmp{rr}=fullfile(parent_dic,fit_dic,runs{rr},sr_nows{ss});
%     end
%     tmp=join(tmp,',');
%     mkdir(fullfile(parent_dic,fit_dic,'tmp'),num2str(ss));
%     command_line=[command_line sprintf('%s -sr_nows %s -capture_ratios %s -out_path %s -density %s &\n',heapa,tmp{1},cr_list,fullfile(parent_dic,fit_dic,'tmp',num2str(ss)),para{8})];
% end
% command_line=[command_line,sprintf('wait\n')];
% system(command_line);

handles=[];
names=[];
min_error=inf;
all_errors=zeros(length(sr_nows),length(capture_ratios));
fid=fopen(fullfile(parent_dic,fit_dic,'fit_para'),'w');
for ss=1:length(sr_nows)
    para=split(sr_nows{ss},'_');
    for ca=1:length(capture_ratios)
        hic=load(fullfile(parent_dic,fit_dic,'tmp',num2str(ss),['hic_',num2str(capture_ratios(ca))]));
        hic=max(0,hic-median(hic(:)));
        
        name=['hic_',sr_nows{ss},'_',num2str(capture_ratios(ca))];
        names=[names,{name}];
        handles=[handles,{figure('position',[100 100 1200 1000],'NumberTitle','off', 'Name',name, 'Visible','off')}];
        he=heatmap(hic);
        he.ColorLimits=[0,uplim(ca)];
        he.Colormap=flipud(gray(256));
        he.Colormap(:,1)=1;
        he.GridVisible='off';
        he.XDisplayLabels=XTLE;
        he.YDisplayLabels=XTLE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nonexc=rcmv(:,1)>exc & rcmv(:,2)<chain_length-exc;
        sims=hic(rcmv(nonexc,1)+1+rcmv(nonexc,2)*chain_length);
        error_fun=@(x) norm(sims/x-rcmv(nonexc,3),4);
        [coff,all_errors(ss,ca)]=fminunc(error_fun,sims'*sims/(sims'*rcmv(nonexc,3)));
        norm_meas=min(rcmv(:,3)*coff,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fid,'%s\t%s\t%f\t%f\n',para{9},para{10},capture_ratios(ca),all_errors(ss,ca));
        
        if all_errors(ss,ca)<min_error
            min_error=all_errors(ss,ca);
            [stiffness,density,processivity,separation,D3_per_D1]=para{7:11};
            capture_ratio=capture_ratios(ca);
        end
        
        hic_exp=zeros(size(hic));
        hic_exp(rcmv(:,1)+1+rcmv(:,2)*chain_length)=norm_meas;
        hic_exp(rcmv(:,2)+1+rcmv(:,1)*chain_length)=norm_meas;
        
        for vv=1:length(vps)
            name=['4C_',sr_nows{ss},'_',num2str(capture_ratios(ca)),'_',num2str(vps(vv))];
            names=[names,{name}];
            handles=[handles,{figure('position',[0 100 1800 800],'NumberTitle','off', 'Name',name, 'Visible','off')}];
            subplot(3,1,1)
            bar(hic(vps(vv)+1,:),1);
            ax=gca;
            ax.XTickLabel=XTL;
            ax.XTick=XT;
            ylim([0,uplim(ca)]);
            subplot(3,1,2)
            bar(hic_exp(vps(vv)+1,:),1);
            ax=gca;
            ax.XTickLabel=XTL;
            ax.XTick=XT;
            ylim([0,uplim(ca)]);
            subplot(3,1,3)
            bar(log2(hic(vps(vv)+1,:)./hic_exp(vps(vv)+1,:)),1);
            ax=gca;
            ax.XTickLabel=XTL;
            ax.XTick=XT;
            ylim([-1,1]);
        end
    end
end
fclose(fid);

for hl=1:length(handles)
    saveas(handles{hl},fullfile(parent_dic,fit_dic,'figures',[names{hl},'.jpg']),'jpeg');
end
close('all');

set_pro_path='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/set_range_and_CTCF_ver3/build/set_range_and_ctcf_ver3';
sim_pro_path='/home/ljw/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/loop_extrusion_3d_sim_ver5/build/loop_extrusion_3d_sim_ver5';
expers={'HEC-1-B_WT','HEC-1-B_F','HEC-1-B_FF','HEC-1-B_FFRR','HEC-1-B_R1','HEC-1-B_R2','HEC-1-B_RF','HEC-1-B_RRFF'};
chips={'chip/HEC-1-B_WT_CTCF.bam_peaks.txt','chip/HEC-1-B_F_CTCF.bam_peaks.txt','chip/HEC-1-B_FF_CTCF.bam_peaks.txt',...
    'chip/HEC-1-B_FFRR_CTCF.bam_peaks.txt','chip/HEC-1-B_R1_CTCF.bam_peaks.txt','chip/HEC-1-B_R2_CTCF.bam_peaks.txt','chip/HEC-1-B_RF_CTCF.bam_peaks.txt','chip/HEC-1-B_RRFF_CTCF.bam_peaks.txt'};
c4_rep1={'4C/4C_HEC-1-B_a12_WT_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_WT_rep1.bedGraph','4C/4C_HEC-1-B_a12_F_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_F_rep1.bedGraph','4C/4C_HEC-1-B_a12_FF_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_FF_rep1.bedGraph',...
    '4C/4C_HEC-1-B_a12_FFRR_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_FFRR_rep1.bedGraph','4C/4C_HEC-1-B_a12_R1_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_R1_rep1.bedGraph','4C/4C_HEC-1-B_a12_R2_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_R2_rep1.bedGraph',...
    '4C/4C_HEC-1-B_a12_RF_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_RF_rep1.bedGraph','4C/4C_HEC-1-B_a12_RRFF_rep1.bedGraph,4C/4C_HEC-1-B_HS5-1_RRFF_rep1.bedGraph'};
c4_rep2={'4C/4C_HEC-1-B_a12_WT_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_WT_rep2.bedGraph','4C/4C_HEC-1-B_a12_F_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_F_rep2.bedGraph','4C/4C_HEC-1-B_a12_FF_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_FF_rep2.bedGraph',...
    '4C/4C_HEC-1-B_a12_FFRR_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_FFRR_rep2.bedGraph','4C/4C_HEC-1-B_a12_R1_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_R1_rep2.bedGraph','4C/4C_HEC-1-B_a12_R2_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_R2_rep2.bedGraph',...
    '4C/4C_HEC-1-B_a12_RF_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_RF_rep2.bedGraph','4C/4C_HEC-1-B_a12_RRFF_rep2.bedGraph,4C/4C_HEC-1-B_HS5-1_RRFF_rep2.bedGraph'};
nipbls='~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B/nipbl/Hec1bwt_nipbl_treat_afterfiting_chr5.bdg';

long_runs={'run1','run2','run3','run4','run5','run6','run7'};
long_record_1D=50000;
for ad=1:length(expers)
    fid=fopen(fullfile(parent_dic,expers{ad},'set_sim.sh'),'w');
    fprintf(fid,'cd %s\n',fullfile(parent_dic,expers{ad}));
    fprintf(fid,...
        '%s -start 140160700 -end 140920000 \\\n-C4s_rep1 %s -C4s_rep2 %s \\\n-CTCF_peaks %s -nipbls %s \\\n-mut_inform mut_inform -chromosome chr5 -stiffnesses %s -densities %s -processivities %s -separations %s -D3_per_D1s %s -capture_ratios %f -record_1D %d\n',...
        set_pro_path,c4_rep1{ad},c4_rep2{ad},chips{ad},nipbls,stiffness,density,processivity,separation,D3_per_D1,capture_ratio,long_record_1D);
    for lr=1:length(long_runs)
        fprintf(fid,'cd %s\n',fullfile(parent_dic,expers{ad},long_runs{lr}));
        fprintf(fid,'%s \\\n-pos_file %s -THR_MAX 1 &\n',sim_pro_path,fullfile(parent_dic,expers{ad},'pos_file'));
    end
    fprintf(fid,'wait\n');
    fclose(fid);
end