function [F,assocF,CX,CXminassoc]=my_MaxMinHeuristic(CX,T,CPC,newPC,CXminassoc,dim,tree,p_thres,sample_thres,subset_thres)
for cx=1:length(CX)
    X=CX(cx);
    minassoc=CXminassoc(cx);
    for ss=0:min(subset_thres,length(CPC))
        combos=my_nchoosek(CPC,ss);
        for cc=1:size(combos,1)
            if length(CPC)==length(newPC) || any(ismember(newPC,combos(cc,:)))
                assoc=my_Assoc(X,T,combos(cc,:),dim,tree,p_thres,sample_thres);
                if assoc<minassoc
                    minassoc=assoc;
                end
                if minassoc==0
                    break
                end
            end
        end
        if minassoc==0
            break
        end
    end
    CXminassoc(cx)=minassoc;
end
[assocF,cx]=max(CXminassoc);
F=CX(cx);
CX(CXminassoc==0)=[];
CXminassoc(CXminassoc==0)=[];
end