clc
clear

exc=10;
parent_dic='~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/human_HEC-1-B';
expers={'HEC-1-B_WT','HEC-1-B_F','HEC-1-B_FF','HEC-1-B_FFRR','HEC-1-B_R1','HEC-1-B_R2','HEC-1-B_RF','HEC-1-B_RRFF'};
shorten={'WT','F','FF','FFRR','R1','R2','RF','RRFF'};
data_name={'sim','entropy','exp'};
uplim=0.02;

vpss=cell(size(expers));
hics=cell(3,length(expers));
XTLEs=cell(size(expers));
XTLs=cell(size(expers));
XTs=cell(size(expers));
for ee=1:length(expers)
    [vpss{ee},hics{1,ee},hics{2,ee},hics{3,ee},XTLEs{ee},XTLs{ee},XTs{ee}]=my_draw(fullfile(parent_dic,expers{ee}),exc);
end

comp_pair=[1,2;1,3;1,4;1,5;1,6;1,7;1,8];
C4s=cell(1,2);
for cp=1:size(comp_pair,1)
    ai=comp_pair(cp,1);
    bi=comp_pair(cp,2);
    cvps=intersect(vpss{ai},vpss{bi});
    for vv=1:length(cvps)
        for dd=1:3
            figure('position',[100 100 1200 1000],'NumberTitle','off', 'Name',[data_name{dd},'_',shorten{ai},'_',shorten{bi},'_',num2str(cvps(vv))]);
            for ee=[ai,bi]
                tmp=find([ai,bi]==ee);
                nor(ee)=sum(hics{dd,ee}(cvps(vv)+1,[1:cvps(vv)+1-exc,cvps(vv)+1+exc:end]));
                subplot(3,1,tmp)
                bar(hics{dd,ee}(cvps(vv)+1,:)/nor(ee),1)
                ax=gca;
                ax.XTickLabel=XTLs{ee};
                ax.XTick=XTs{ee};
                ylim([0,uplim]);
            end
            
            subplot(3,1,3)
            bar(hics{dd,bi}(cvps(vv)+1,:)/nor(bi)-hics{dd,ai}(cvps(vv)+1,:)/nor(ai),1)
            ylim([-uplim/5,uplim/5]);
            
            print('-painters',fullfile(parent_dic,'figures',[data_name{dd},'_',shorten{ai},'_',shorten{bi},'_',num2str(cvps(vv))]),'-djpeg');
            print('-painters',fullfile(parent_dic,'figures',[data_name{dd},'_',shorten{ai},'_',shorten{bi},'_',num2str(cvps(vv))]),'-depsc');
        end
    end
end

exp_hic=[1];
for ee=1:length(exp_hic)
    hic=load(fullfile(parent_dic,expers{ee},'hic'));
    figure('position',[100 100 1200 1000],'NumberTitle','off', 'Name',['hic_sim_',shorten{ee}]);
    he=heatmap(hic);
    he.ColorLimits=[0,uplim];
    he.Colormap=flipud(gray(256));
    he.Colormap(:,1)=1;
    he.GridVisible='off';
    he.XDisplayLabels=XTLEs{ee};
    he.YDisplayLabels=XTLEs{ee};
    
    print('-painters',fullfile(parent_dic,'figures',['hic_sim_',shorten{ee}]),'-djpeg');
    print('-painters',fullfile(parent_dic,'figures',['hic_sim_',shorten{ee}]),'-depsc');
end