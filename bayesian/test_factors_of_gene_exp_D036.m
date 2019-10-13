clc
clear
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % load promoter loop
% %%%%%%%%%%%%%%%%%%%%%%%%%
load chr_lib
files={'GSE84660_CHiC_contacts_b2b.bed','GSE84660_CHiC_contacts_b2g.bed'};
bait_pos=[];
other_pos=[];
loop_count=[];
fi_ind=1:2; % change this to determine the file to use
spot_num=[213893144,196091990,195559022,169995874,155618089,159040099];
for fi=fi_ind
    file=files{fi};
    fid=fopen(file,'rb');
    readin=fread(fid,'*char')';
    fclose(fid);
    readin=splitlines(readin);
    readin([1,2,end])=[];
    readin=split(readin,char(9));
    if fi==1
        bait_pos=[join(readin(:,1:3),char(9),2);join(readin(:,4:6),char(9),2)];
        other_pos=bait_pos([end/2+1:end,1:end/2]);
        [~,IA,~]=unique(strcat(bait_pos,other_pos));
        loop_count=repmat(sum(reshape(str2double(readin(:,10:15))./spot_num,[],3,2),3),2,1);
        bait_pos=split(bait_pos(IA),char(9));
        other_pos=split(other_pos(IA),char(9));
        loop_count=loop_count(IA,:);
    else
        bait_pos=[bait_pos;readin(:,1:3)];
        other_pos=[other_pos;readin(:,4:6)];
        loop_count=[loop_count;sum(reshape(str2double(readin(:,10:15))./spot_num,[],3,2),3)];
    end
end
tape=~strcmp(bait_pos(:,1),other_pos(:,1));
bait_pos(tape,:)=[];
other_pos(tape,:)=[];
loop_count(tape,:)=[];
O2O=false; % change this to determine whether consider bait regulated by several enhancer
if O2O
    [~,~,IC]=unique(join(bait_pos,char(9),2));
    tape=accumarray(IC,loop_count(:,1)>0,[],[],[],false)<=1 & accumarray(IC,loop_count(:,2)>0,[],[],[],false)<=1 & accumarray(IC,loop_count(:,3)>0,[],[],[],false)<=1;
    tape=~tape(IC);
    bait_pos(tape,:)=[];
    other_pos(tape,:)=[];
    loop_count(tape,:)=[];
end
loop_thres=0; % change this to determine loop size threshold
tape=abs(mean(str2double(bait_pos(:,2:3)),2)-mean(str2double(other_pos(:,2:3)),2))<loop_thres;
bait_pos(tape,:)=[];
other_pos(tape,:)=[];
loop_count(tape,:)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % load enhancer
% %%%%%%%%%%%%%%%%%%%%%%%%%%
ind_chr=cellfun(@(x) reshape(find(strcmp(other_pos(:,1),x)),1,[]),chr_lib{1},'uniformoutput',false);
loop_enhancer=zeros(size(other_pos,1),6);
enhancer_mode='bedgraph'; % change this to determine use bedgraph or macs
enhancer_ind=1:6; % change this to determine which enhancer to use
switch enhancer_mode
    case 'bedgraph'
        for oe=enhancer_ind
            file=['GSM224724',num2str(oe),'_H3K27ac_ChIPseq_Donor',num2str(ceil(oe/3)),'_Day',num2str(mod(oe-1,3)*3),'.bedgraph'];
            [~,tape]=unix(['wc -l ',file]);
            tape=str2double(tape(1:find(tape==' ',1,'first')-1));
            chr=cell(tape,1);
            whole_pe=zeros(tape,3);
            fid=fopen(file,'rb');
            rownow=0;
            while ~feof(fid)
                readin=split(fgetl(fid),char(9));
                rownow=rownow+1;
                chr(rownow)=readin(1);
                whole_pe(rownow,:)=str2double(readin(2:4));
            end
            fclose(fid);
            for cc=1:length(chr_lib{1})
                tape_pe=whole_pe(strcmp(chr,chr_lib{1}{cc}),:);
                tape=zeros(1,chr_lib{2}(cc));
                for en=1:size(tape_pe,1)
                    tape(tape_pe(en,1)+1:tape_pe(en,2))=tape_pe(en,3);
                end
                for ii=ind_chr{cc}
                    loop_enhancer(ii,oe)=mean(tape(str2double(other_pos{ii,2}):str2double(other_pos{ii,3})));
                end
            end
        end
        enhancer_normal=false; % change this to determine whether normalize enhancer
        if enhancer_normal
            enhancer_fqsize=[74859282,36966224,37938741,40400178,31080296,52300443];
            loop_enhancer=loop_enhancer./enhancer_fqsize;
        end
    case 'macs'
        thres=100;
        for oe=enhancer_ind
            file=['GSM224724',num2str(oe),'_H3K27ac_ChIPseq_Donor',num2str(ceil(oe/3)),'_Day',num2str(mod(oe-1,3)*3),'.peak'];
            fid=fopen(file,'rb');
            readin=fread(fid,'*char')';
            fclose(fid);
            readin=splitlines(readin);
            readin([1,end])=[];
            readin=split(readin,char(9));
            for cc=1:length(chr_lib{1})
                tape=strcmp(readin(:,1),chr_lib{1}{cc});
                loop_enhancer(ind_chr{cc},oe)=attach_int_2_int(str2double(readin(tape,2:3)),str2double(readin(tape,5)),str2double(other_pos(ind_chr{cc},2:3)),thres);
            end
        end
end
loop_enhancer=loop_enhancer(:,1:3)+loop_enhancer(:,4:6);
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % load CTCF
% %%%%%%%%%%%%%%%%%%%%%%%%%%
inside=500;
loop_insulator=zeros(size(other_pos,1),3);
insulator_mode='bedgraph_mean'; % change this to determine use bedgraph_mean, bedgraph_sum, or macs
switch insulator_mode
    case {'bedgraph_mean','bedgraph_sum'}
        for li=1:3
            file=['GSM224724',num2str(6+li),'_CTCF_ChIPseq_Day',num2str(mod(li-1,3)*3),'.bedgraph'];
            [~,tape]=unix(['wc -l ',file]);
            tape=str2double(tape(1:find(tape==' ',1,'first')-1));
            chr=cell(tape,1);
            whole_pe=zeros(tape,3);
            fid=fopen(file,'rb');
            rownow=0;
            while ~feof(fid)
                readin=split(fgetl(fid),char(9));
                rownow=rownow+1;
                chr(rownow)=readin(1);
                whole_pe(rownow,:)=str2double(readin(2:4));
            end
            fclose(fid);
            for cc=1:length(chr_lib{1})
                tape_pe=whole_pe(strcmp(chr,chr_lib{1}{cc}),:);
                tape=zeros(1,chr_lib{2}(cc));
                for en=1:size(tape_pe,1)
                    tape(tape_pe(en,1)+1:tape_pe(en,2))=tape_pe(en,3);
                end
                for ii=ind_chr{cc}
                    switch insulator_mode
                        case 'bedgraph_mean'
                            loop_insulator(ii,li)=mean(tape(min(str2double(bait_pos{ii,3}),str2double(other_pos{ii,3}))+inside:max(str2double(bait_pos{ii,2}),str2double(other_pos{ii,2}))-inside));
                        case 'bedgraph_sum'
                            loop_insulator(ii,li)=sum(tape(min(str2double(bait_pos{ii,3}),str2double(other_pos{ii,3}))+inside:max(str2double(bait_pos{ii,2}),str2double(other_pos{ii,2}))-inside));
                    end
                end
            end
        end
        insulator_normal=false; % change this to determine whether normalize insulator
        if insulator_normal
            insulator_fqsize=[28214852,36121176,41994497];
            loop_insulator=loop_insulator./insulator_fqsize;
        end
    case 'macs'
        loop_int=sort([mean(str2double(bait_pos(:,2:3)),2),mean(str2double(other_pos(:,2:3)),2)],2);
        [header,seq]=fastaread('hg19.fa');
        extend=50;
        ext_peak_path='ext_peak.fa';
        output_path='output_path/';
        motif_path='MA0139.1.meme';
        for li=1:3
            file=['GSM224724',num2str(6+li),'_CTCF_ChIPseq_Day',num2str(mod(li-1,3)*3),'.peak'];
            fid=fopen(file,'rb');
            readin=fread(fid,'*char')';
            fclose(fid);
            readin=splitlines(readin);
            readin([1,end])=[];
            readin=split(readin,char(9));
            fid=fopen(ext_peak_path,'w');
            for ct=1:size(readin,1)
                ind=find(strcmp(header,readin{ct,1}));
                fprintf(fid,'>%s\n',num2str(ct));
                fprintf(fid,'%s\n',seq{ind}(max(str2double(readin{ct,2})-extend,1):min(str2double(readin{ct,3})+extend,length(seq{ind}))));
            end
            fclose(fid);
            [CBS_chr,CBS_se,CBS_ori,CBS_enrich]=my_call_CBS(readin(:,1),str2double(readin(:,2:3)),str2double(readin(:,5)),ext_peak_path,output_path,motif_path);
            for cc=1:length(chr_lib{1})
                tape=strcmp(CBS_chr,chr_lib{1}{cc});
                loop_insulator(ind_chr{cc},li)=insert_int_2_int(CBS_se(tape,:),CBS_enrich(tape),loop_int(ind_chr{cc},:),inside);
            end
        end
        clear('seq');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % get fpkm
% %%%%%%%%%%%%%%%%%%%%%%%%%
file='GSE84661_GV19_TSS_RNAseq_RSEM_counts.bed';
fid=fopen(file,'rb');
readin=fread(fid,'*char')';
fclose(fid);
readin=splitlines(readin);
readin([1,2,end])=[];
readin=split(readin,char(9));
fpkm=str2double(readin(:,8:13));
fpkm=fpkm(:,[1,3,5])+fpkm(:,[2,4,6]);
loop_fpkm=zeros(size(other_pos,1),3);
thres=0;
for cc=1:length(chr_lib{1})
    tape=strcmp(readin(:,1),chr_lib{1}{cc});
    loop_fpkm(ind_chr{cc},:)=attach_int_2_int(str2double(readin(tape,2)),fpkm(tape,:),str2double(bait_pos(ind_chr{cc},2:3)),thres);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % get loop size
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
loop_size=abs(mean(str2double(bait_pos(:,2:3)),2)-mean(str2double(other_pos(:,2:3)),2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculate factors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,IA,IC]=unique(join(bait_pos,char(9),2));
tape=length(IA);
sum_enhancer=zeros(tape,3);
mean_insulator=zeros(tape,3);
mean_size=zeros(tape,3);
mean_count=zeros(tape,3);
for ii=1:size(loop_count,2)
    sum_enhancer(:,ii)=accumarray(IC,loop_enhancer(:,ii),[],[],[],false);
    mean_insulator(:,ii)=accumarray(IC,loop_enhancer(:,ii).*loop_insulator(:,ii),[],[],[],false)./sum_enhancer(:,ii);
    mean_insulator(isnan(mean_insulator(:,ii)),ii)=mean(mean_insulator(:,ii),'omitnan');
    mean_size(:,ii)=accumarray(IC,loop_enhancer(:,ii).*loop_size,[],[],[],false)./sum_enhancer(:,ii);
    mean_size(isnan(mean_size(:,ii)),ii)=mean(mean_size(:,ii),'omitnan');
    mean_count(:,ii)=accumarray(IC,loop_enhancer(:,ii).*loop_count(:,ii),[],[],[],false)./sum_enhancer(:,ii);
    mean_count(isnan(mean_count(:,ii)),ii)=mean(mean_count(:,ii),'omitnan');
end
np_fpkm=loop_fpkm(IA,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save Data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=cat(3,np_fpkm,sum_enhancer,mean_insulator,mean_size,mean_count);
save(['Data_',num2str(fi_ind,'%d'),'_',num2str(O2O,'%d'),'_',num2str(loop_thres,'%d'),'_',enhancer_mode,'_',num2str(enhancer_ind,'%d'),'_',num2str(enhancer_normal,'%d'),'_',num2str(insulator_mode,'%d')],'Data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard nonexp bait
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tape=sum(Data(:,:,1),2)==0;
Data(tape,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard nonenhancer gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tape=sum(Data(:,:,2),2)==0;
Data(tape,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard noninsulator gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tape=sum(Data(:,:,3),2)==0;
Data(tape,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set names of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_names={'Fpkm','SE','MI','MS','MC'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % draw slope excel
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_excel=zeros(size(Data,1),length(var_names));
for se=1:size(slope_excel,1)
    for vn=1:length(var_names)
        tape=polyfit(1:3,log(Data(se,:,vn)),1);
        slope_excel(se,vn)=tape(1);
    end
end
slope_excel=sortrows(slope_excel,3);
for vn=1:size(slope_excel,2)
    figure('position',[100 100 1500 1000],'NumberTitle','off', 'Name',['slope_',var_names{vn}])
    plot(1:size(slope_excel,1),slope_excel(:,vn),'.k')
    print('-painters',gcf,['slope_',var_names{vn},'.eps'],'-depsc');
    saveas(gcf,['slope_',var_names{vn},'.jpg'],'jpg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretize
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,Data]=sort(Data,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape discretized data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=reshape(Data,[],length(var_names),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_ind=[1:3,5]; % change this to determine the variable to learn
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sorted heatmap
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data_ori=reshape(Data,[],3,length(var_names));
[Data_ori(:,:,1),dsi]=sort(Data_ori(:,:,1),2);
for ii=1:size(Data_ori,1)
    Data_ori(ii,:,2:end)=Data_ori(ii,dsi(ii,:),2:end);
end
xfour=repelem([0,1,2;1,2,3;1,2,3;0,1,2],1,size(Data_ori,1))+0.5;
yfour=[repmat(1:size(Data_ori,1),2,3);repmat(0:size(Data_ori,1)-1,2,3)];
def_colors=[255,0,0;135,206,250;154,205,50;218,112,214;238,118,33]/255;
for vn=var_ind
    figure('position',[100 100 1000 500],'NumberTitle','off', 'Name',var_names{vn})
    fill(xfour,yfour,reshape(Data_ori(:,:,vn),1,[]),'EdgeColor','none')
    set(gca,'xtick',1:3);
    CM=[linspace(1,def_colors(vn,1));linspace(1,def_colors(vn,2));linspace(1,def_colors(vn,3))]';
    colormap(CM);
    ax=gca;
    ax.CLim=[1,3];
    colorbar('eastoutside')
    axis([0.5,3.5,0,size(Data_ori,1)]);
    axis ij
    print('-painters',gcf,['sort_heapmap_',var_names{vn},'.eps'],'-depsc');
    saveas(gcf,['sort_heapmap_',var_names{vn},'.jpg'],'jpg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate AD-tree from data (no loop count)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
SizeUnit=50000;
dim=repmat(3,size(var_names(var_ind)));
tree=my_MakeADTree(dim,Data(:,var_ind),SizeUnit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% learn Bayesian Network by MMHC (no loop count)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_thres=0.05;
sample_thres=5;
subset_thres=5;

skeleton=my_MMPC(dim,tree,p_thres,sample_thres,subset_thres);
DAG=zeros(length(dim));
para=1;
list_length=100;
endur_fail=15;
DAG=my_discrete_Bayes_paraUni_strucUni_Greedy(DAG,skeleton,para,dim,tree,list_length,endur_fail);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw Bayesian Network and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edge_mat=corr(Data(:,var_ind),Data(:,var_ind));
% edge_mat=zeros(length(var_ind));
% for ii=1:length(var_ind)
%     for jj=ii+1:length(var_ind)
%         [~,~,tape]=crosstab(Data(:,var_ind(ii)),Data(:,var_ind(jj)));
%         edge_mat(ii,jj)=-log10(tape);
%         edge_mat(jj,ii)=edge_mat(ii,jj);
%     end
% end
figure('position',[100 100 1500 1000],'NumberTitle','off', 'Name','gene expression factor')
tape=DAG.*edge_mat;
tape(isnan(tape))=0;
tape=digraph(tape,var_names(var_ind));
plot(tape,'EdgeLabel',tape.Edges.Weight)
print('-painters',gcf,'gene_expression_factor_nc.eps','-depsc');
saveas(gcf,'gene_expression_factor_nc.jpg','jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% learn Bayesian Network by exact search (no loop count)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAG=my_discrete_Bayes_paraUni_strucUni_Exact(para,dim,tree);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw Bayesian Network and save (exact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[100 100 1500 1000],'NumberTitle','off', 'Name','gene expression factor (exact)')
tape=DAG.*edge_mat;
tape(isnan(tape))=0;
tape=digraph(tape,var_names(var_ind));
plot(tape,'EdgeLabel',tape.Edges.Weight)
print('-painters',gcf,'gene_expression_factor_exact_nc.eps','-depsc');
saveas(gcf,'gene_expression_factor_exact_nc.jpg','jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw crosstab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfour=repelem([0,1,2;1,2,3;1,2,3;0,1,2],1,3)+0.5;
yfour=[repmat(1:3,2,3);repmat(0:2,2,3)]+0.5;
def_colors=[255,0,0]/255;
for vn1=var_ind
    for vn2=var_ind
        if vn2>vn1
            figure('position',[100 100 1000 1000],'NumberTitle','off', 'Name',['contingency_table_',var_names{vn1},'_',var_names{vn2}])
            tape=crosstab(Data(:,vn1),Data(:,vn2));
            fill(xfour,yfour,tape(:)','EdgeColor','none')
            set(gca,'xtick',1:3,'ytick',1:3);
            CM=[linspace(1,def_colors(1));linspace(1,def_colors(2));linspace(1,def_colors(3))]';
            colormap(CM);
            ax=gca;
%             ax.CLim=[min(tape(:)),max(tape(:))];
            ax.CLim=[0,8000];
            colorbar('eastoutside')
            axis([0.5,3.5,0.5,3.5]);
            axis ij
            [II,JJ,VV]=find(tape);
            text(JJ,II,num2str(VV));
            print('-painters',gcf,['contab_',var_names{vn1},'_',var_names{vn2},'.eps'],'-depsc');
            saveas(gcf,['contab_',var_names{vn1},'_',var_names{vn2},'.jpg'],'jpg');
        end
    end
end



