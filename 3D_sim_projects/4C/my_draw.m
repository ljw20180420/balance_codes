function [vps,hic,hic_max_entropy,hic_exp,XTLE,XTL,XT]=my_draw(save_path,exc)
[rcmv,~,vps,XTLE,XTL,XT]=load_target_CTCF(fullfile(save_path,'pos_file'));
chain_length=length(XTLE);

hic=load(fullfile(save_path,'hic'));
hic=max(0,hic-median(hic(:)));

hic_max_entropy=load(fullfile(save_path,'hic_max_entropy'));
hic_max_entropy=max(0,hic_max_entropy-median(hic_max_entropy(:)));

nonexc=rcmv(:,1)>exc & rcmv(:,2)<chain_length-exc;
sims=hic(rcmv(nonexc,1)+1+rcmv(nonexc,2)*chain_length);
error_fun=@(x) norm(sims/x-rcmv(nonexc,3),4);
coff=fminunc(error_fun,sims'*sims/(sims'*rcmv(nonexc,3)));
norm_meas=min(rcmv(:,3)*coff,1);
hic_exp=zeros(size(hic));
hic_exp(rcmv(:,1)+1+rcmv(:,2)*chain_length)=norm_meas;
hic_exp(rcmv(:,2)+1+rcmv(:,1)*chain_length)=norm_meas;
end