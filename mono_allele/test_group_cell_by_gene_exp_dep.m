clc
clear
GP_fold='2019_10_12';
mkdir(GP_fold);
file_list={'mouse_cortex_23178.xlsx'};
file_list_sn={'23178'};
alpha_gene=strcat('pcdha',strtrim(cellstr(num2str((1:12)'))))';
beta_gene=strcat('pcdhb',strtrim(cellstr(num2str((1:22)'))))';
gamma_gene=[strcat('pcdhga',strtrim(cellstr(num2str((1:12)'))))',strcat('pcdhgb',strtrim(cellstr(num2str([1,2,4:8]'))))'];
% target_gene_list={alpha_gene,beta_gene,gamma_gene};
% target_gene_list_sn={'alpha','beta','gamma'};
% target_gene_list_ord=[2,3,3];
target_gene_list={beta_gene,gamma_gene};
target_gene_list_sn={'beta','gamma'};
target_gene_list_ord=[3,3];
% classes={'','Pri','Ant','GAB','Glu'};
% ref_num=[1,1,1,2,2];
classes={''};
ref_num=[1];
tol=1e-6;
uprow=inf;
[~,~,tape]=xlsread('GSE115746_complete_metadata_28467-cells.xlsx');
sample_name_ref=lower(tape(4:end,1));
sample_name_ref=cellfun(@(x) x(2:end-1),sample_name_ref,'uniformoutput',false);
class_ref=tape(4:end,2:3);
        
for fl=1:length(file_list)
    fl
    file=file_list{fl};
    [~,~,raw]=xlsread(file);
    gene=lower(raw(2,2:end));
    sample=lower(raw(3:end,1));
    raw=cell2mat(raw(3:end,2:end));
    raw(raw>1)=1;
    raw(raw<1)=0;
    for tl=1:length(target_gene_list)
        tl
        target_gene=lower(target_gene_list{tl});
        yfour=[0:length(target_gene)-1;0:length(target_gene)-1;1:length(target_gene);1:length(target_gene)];
        [LIA,LOCB]=ismember(target_gene,gene);
        if ~all(LIA)
            error('not all target genes are found');
        end
        Data_mix=raw(:,LOCB);
        for cl=1:length(classes)
            cl
            tic
            Data_class=Data_mix(ismember(sample,sample_name_ref(strncmpi(class_ref(:,ref_num(cl)),classes{cl},length(classes{cl})))),:);
            [Data_class,~,IC]=unique(Data_class,'row');
            Counts_class=accumarray(IC,ones(size(IC)),[],[],[],false);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Data_dep,Counts_dep,Data_merge,Counts_merge]=dep_data(Data_class,Counts_class,target_gene_list_ord(tl),tol,uprow);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meth={'class','merge'};
            for me=1:length(meth)
                eval(['Data=Data_',meth{me},';']);
                eval(['Counts=Counts_',meth{me},';']);
                expnum=sum(Data,2);
                Data(expnum==0,:)=[];
                Counts(expnum==0)=[];
                expnum(expnum==0)=[];
                [Data,ind]=sortrows([Data,expnum],[length(target_gene)+1,-(1:length(target_gene))]);
                Counts=Counts(ind);
                expnum=expnum(ind);
                Data(:,end)=[];
                tape=cumsum(Counts');
                tape=tape/tape(end);
                xfour=[0,tape(1:end-1);tape;tape;0,tape(1:end-1)];
                xtick=zeros(1,2*(target_gene_list_ord(tl)-1)+2);
                for do=1:2*(target_gene_list_ord(tl)-1)
                    xtick(do+1)=sum(Counts(expnum==do));
                end
                xtick(end)=sum(Counts(expnum>2*(target_gene_list_ord(tl)-1)));
                xtick=cumsum(xtick);
                xtick=xtick/xtick(end);
                
                [II,JJ]=find(Data);
                figure('position',[100 100 1000 1500],'NumberTitle','off', 'Name',[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'_bar'],'Visible', 'off')
                expnum_max=max(expnum);
                barnum=accumarray(expnum,Counts);
                total_barnum=sum(barnum);
                bar((1:expnum_max)',barnum,'r','EdgeColor','none')
                for bn=1:expnum_max
                    text(bn,barnum(bn),sprintf('%.2f%%',barnum(bn)*100/total_barnum),'HorizontalAlignment','center','VerticalAlignment','bottom');
                end
                print('-painters',gcf,fullfile(GP_fold,[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'_bar.eps']),'-depsc');
                saveas(gcf,fullfile(GP_fold,[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'_bar.jpg']),'jpg');
                figure('position',[100 100 1000 1500],'NumberTitle','off', 'Name',[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me}],'Visible', 'off')
                fill(xfour(:,II),yfour(:,JJ),'r','EdgeColor','none')
                set (gca,'Ydir','reverse');
                axis([0,1,0,length(target_gene)]);
                if xtick(end-1)>=xtick(end)
                    xtick(end)=[];
                end
                set(gca,'xtick',xtick,'xticklabel',strcat(strtrim(cellstr(num2str(100*xtick','%.2f'))),'%'),'ytick',(1:length(target_gene))-0.5,'yticklabel',target_gene)
%                 if strcmp(meth{me},'class')
%                     print('-painters',gcf,fullfile(fold,GP_fold,[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'.eps']),'-depsc');
%                     saveas(gcf,fullfile(fold,GP_fold,[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'.jpg']),'jpg');
%                 end
%                 print('-painters',gcf,fullfile(fold,GP_fold,[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'.eps']),'-depsc');
                saveas(gcf,fullfile(GP_fold,[file_list_sn{fl},'_',target_gene_list_sn{tl},'_',classes{cl},'_',meth{me},'.jpg']),'jpg');
                close('all');
            end
            toc
        end
    end
end