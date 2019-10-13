function assoc=my_Assoc(X,T,S,dim,tree,p_thres,sample_thres)
if T==1 && X==5
    1
end
insX=find([S,inf]>X,1,'first');
insT=find([S,inf]>T,1,'first');
if insX<insT
    sub=[S(1:insX-1),X,S(insX:insT-1),T,S(insT:end)];
    ins1=insX;
    ins2=insT+1;
else
    if insX>insT
        sub=[S(1:insT-1),T,S(insT:insX-1),X,S(insX:end)];
    else
        if X<T
            sub=[S(1:insX-1),X,T,S(insT:end)];
        else
            sub=[S(1:insX-1),T,X,S(insT:end)];
        end
    end
    ins1=insT;
    ins2=insX+1;
end
tab_12S=my_MakeContab(sub,dim,tree);
if tree(1,1)/prod(dim(S))/(dim(X)-1)/(dim(T)-1)>=sample_thres
    tab_2S=sum(tab_12S,ins1);
    tab_1S=sum(tab_12S,ins2);
    tab_S=sum(tab_1S,ins1);
    G2=2*sum(reshape(tab_12S.*log(tab_12S.*tab_S./tab_1S./tab_2S),[],1),'omitnan');
    df=nnz(tab_12S)-nnz(tab_1S)-nnz(tab_2S)+nnz(tab_S);
    if df<=0
        p_value=0;
    else
        p_value=chi2cdf(G2,df,'upper');
    end
    assoc=max(0,p_thres-p_value);
else
    assoc=0;
end
end