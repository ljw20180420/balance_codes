function CPC=my_barMMPC(CX,T,CPC,dim,tree,p_thres,sample_thres,subset_thres)
CXminassoc=inf(size(CX));
newPC=CPC;
[F,assocF,CX,CXminassoc]=my_MaxMinHeuristic(CX,T,CPC,newPC,CXminassoc,dim,tree,p_thres,sample_thres,subset_thres);
while assocF>0
    ins=find([CPC,inf]>F,1,'first');
    CPC=[CPC(1:ins-1),F,CPC(ins:end)];
    rmind=find(CX==F,1,'first');
    CX(rmind)=[];
    CXminassoc(rmind)=[];
    newPC=F;
    [F,assocF,CX,CXminassoc]=my_MaxMinHeuristic(CX,T,CPC,newPC,CXminassoc,dim,tree,p_thres,sample_thres,subset_thres);
end
kk=1;
while kk<=length(CPC)
    flag=0;
    for ss=0:min(subset_thres,length(CPC)-1)
        combos=my_nchoosek(CPC([1:kk-1,kk+1:end]),ss);
        for cc=1:size(combos,1)
            if my_Assoc(CPC(kk),T,combos(cc,:),dim,tree,p_thres,sample_thres)==0
                CPC(kk)=[];
                flag=1;
                break
            end
        end
        if flag==1
            break
        end
    end
    if flag==0
        kk=kk+1;
    end
end
end