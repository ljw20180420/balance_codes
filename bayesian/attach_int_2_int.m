function tar_val=attach_int_2_int(val_int,val,tar_int,thres)
if isempty(tar_int)
    tar_val=[];
    return
end
[tar_int,~,IC]=unique(tar_int,'rows');
[tar_int_mean,si]=sort(mean(tar_int,2),'ascend');
tar_int=tar_int(si,:);
val_int=mean(val_int,2);
ind=discretize(val_int,[0;(tar_int_mean(1:end-1)+tar_int_mean(2:end))/2;inf]);
log_tape=val_int<tar_int(ind,2)+thres & val_int>tar_int(ind,1)-thres;
tar_val=zeros(size(tar_int,1),size(val,2));
for ii=1:size(val,2)
    tar_val(si,ii)=accumarray(ind(log_tape),val(log_tape,ii),[size(tar_int,1),1],[],[],false);
end
tar_val=tar_val(IC,:);
end