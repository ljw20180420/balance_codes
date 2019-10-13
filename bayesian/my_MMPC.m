function skeleton=my_MMPC(dim,tree,p_thres,sample_thres,subset_thres)
skeleton=zeros(length(dim));
for T=1:length(dim)
    T
    CPC=reshape(find(skeleton(1:T-1,T)),1,[]);
    CX=T+1:length(dim);
    CPC=my_barMMPC(CX,T,CPC,dim,tree,p_thres,sample_thres,subset_thres);
    length(CPC)
    skeleton(T,CPC)=1;
end
skeleton=skeleton.*skeleton';
end