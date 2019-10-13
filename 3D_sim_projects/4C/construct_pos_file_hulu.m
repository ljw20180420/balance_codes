clc
clear
bin=600;
range=[37084400,38042900];
range(2)=ceil((range(2)-range(1))/bin)*bin+range(1);
chain_length=(range(2)-range(1))/bin;
Fpeaks=[];
Rpeaks=[];
chip='~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/mouse_cortex/p0-cortex_WT/cortex_ctcf_peaks.txt';
fid=fopen(chip,'r');
while(~feof(fid))
    tmp=split(fgetl(fid));
    if strcmp(tmp{6},'+')
        Fpeaks=[Fpeaks;str2double(tmp(2:3))'];
    else
        if strcmp(tmp{6},'-')
            Rpeaks=[Rpeaks;str2double(tmp(2:3))'];
        end
    end
end
fclose(fid);
Findex=ceil((mean(Fpeaks,2)-range(1))/bin);
Findex=unique(Findex);
Rindex=ceil((mean(Rpeaks,2)-range(1))/bin);
Rindex=unique(Rindex);
XTL=[repmat({char(8594)},length(Findex),1);repmat({char(8592)},length(Rindex),1)];
XT=[Findex;Rindex];
[XT,IT]=sort(XT);
XTL=XTL(IT);
figure
bar(zeros(chain_length,1),1);
ax=gca;
ax.XTickLabel=XTL;
ax.XTick=XT;

% nipbls={'~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/hulu/nipbl/NIPBL1_treat_afterfiting_chr18.bdg','~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/hulu/nipbl/NIPBL2_treat_afterfiting_chr18.bdg'};
% nip_cover=zeros(range(2)-range(1),1);
% for ni=1:length(nipbls)
%     fid=fopen(nipbls{ni},'r');
%     fgetl(fid);
%     while(~feof(fid))
%         tmp=split(fgetl(fid));
%         peak=str2double(tmp(2:3));
%         peak(1)=max(range(1),peak(1));
%         peak(2)=min(range(2),peak(2));
%         if(peak(1)<peak(2))
%             height=str2double(tmp{4});
%             peak=peak-range(1)+1;
%             nip_cover(peak(1):peak(2)-1)=nip_cover(peak(1):peak(2)-1)+height;
%         end
%     end
%     fclose(fid);
% end
% figure
% bar(nip_cover,1);


stall=0.97;
chain_length=1500;
VP=[321,381,813,836,1037,1114];
left_stall=zeros(chain_length,1);
left_stall(VP(1:4)+1)=stall;
right_stall=zeros(chain_length,1);
right_stall(VP(5:6)+1)=stall;
birth=ones(chain_length,1);
birth(VP(4)+1:VP(5))=20;
stiffness=2;
density=0.2;
processivity=400;
separation=400;
D3_per_D1=1250;
capture_ratio=2;
record_1D=50000;
fid=fopen("~/new_fold/old_desktop/wuqiang/wuyonghu/3D_sim_projects/JINGWEI/4C/hulu/pos_file",'w');
fprintf(fid,"VP\n");
for ii=1:length(VP)
    fprintf(fid,"%d\n",VP(ii));
end
fprintf(fid,"left_stall\n");
for ii=1:length(left_stall)
    fprintf(fid,"%f\n",left_stall(ii));
end
fprintf(fid,"right_stall\n");
for ii=1:length(right_stall)
    fprintf(fid,"%f\n",right_stall(ii));
end
fprintf(fid,"birth\n");
for ii=1:length(birth)
    fprintf(fid,"%f\n",birth(ii));
end
fprintf(fid,"stiffness\n");
for ii=1:length(stiffness)
    fprintf(fid,"%f\n",stiffness(ii));
end
fprintf(fid,"density\n");
for ii=1:length(density)
    fprintf(fid,"%f\n",density(ii));
end
fprintf(fid,"processivity\n");
for ii=1:length(processivity)
    fprintf(fid,"%d\n",processivity(ii));
end
fprintf(fid,"separation\n");
for ii=1:length(separation)
    fprintf(fid,"%d\n",separation(ii));
end
fprintf(fid,"D3_per_D1\n");
for ii=1:length(D3_per_D1)
    fprintf(fid,"%d\n",D3_per_D1(ii));
end
fprintf(fid,"capture_ratio\n");
for ii=1:length(capture_ratio)
    fprintf(fid,"%f\n",capture_ratio(ii));
end
fprintf(fid,"record_1D\n");
for ii=1:length(record_1D)
    fprintf(fid,"%d\n",record_1D(ii));
end
fclose(fid);