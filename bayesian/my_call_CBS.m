function [CBS_chr,CBS_se,CBS_ori,CBS_enrich]=my_call_CBS(CTCF_chr,CTCF_se,CTCF_enrich,ext_peak_path,output_path,motif_path)
try
    rmdir(output_path,'s');
end
system(['/home/ljw/meme/bin/fimo --o ',output_path,' --max-stored-scores 1000000000000 ',motif_path,' ',ext_peak_path]);
fid=fopen([output_path,'fimo.tsv'],'rb');
readin=fread(fid,'*char')';
fclose(fid);
readin=splitlines(readin);
if isempty(readin{end})
    readin(end)=[];
end
readin(end-3:end)=[];
readin=split(readin,char(9));

CBS_index=str2double(readin(2:end,3));
CBS_se=str2double(readin(2:end,4:5))+CTCF_se(CBS_index,1)-1;
CBS_ori=strcmp(readin(2:end,6),'+');
CBS_chr=CTCF_chr(CBS_index);
[uCBS_index,~,IC]=unique(CBS_index);
tape=accumarray(IC,1,[],[],[],false);
CBS_enrich=CTCF_enrich(CBS_index)./tape(IC);

fid=fopen(ext_peak_path,'rb');
readin=fread(fid,'*char')';
fclose(fid);
readin=splitlines(readin);
if isempty(readin{end})
    readin(end)=[];
end
readin=readin(2:2:end);
log_emp=true(size(readin));
log_emp(uCBS_index)=false;
ind_emp=reshape(find(log_emp),1,[]);
emp_CBS_se=zeros(length(ind_emp),2);
emp_CBS_ori=zeros(length(ind_emp),1);
load CTCF_Jaspar CTCF_Jaspar
for em=1:length(ind_emp)
    [emp_CBS_se(em,1),emp_CBS_ori(em)]=find_CBS_with_highest_score(readin{ind_emp(em)},CTCF_Jaspar);
    emp_CBS_se(em,2)=emp_CBS_se(em,1)+size(CTCF_Jaspar,2)-1;
    emp_CBS_se(em,:)=emp_CBS_se(em,:)+CTCF_se(ind_emp(em),1)-1;
end

CBS_chr=[CBS_chr;CTCF_chr(log_emp)];
CBS_se=[CBS_se;emp_CBS_se];
CBS_ori=[CBS_ori;emp_CBS_ori];
CBS_enrich=[CBS_enrich;CTCF_enrich(log_emp)];
[~,IA,IC]=unique(strcat(CBS_chr,cellstr(num2str(CBS_se))),'stable');
CBS_chr=CBS_chr(IA);
CBS_se=CBS_se(IA);
CBS_ori=CBS_ori(IA);
CBS_enrich=accumarray(IC,CBS_enrich,[],[],[],false);
end