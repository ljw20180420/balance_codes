function local_score=BDeU_local_score(var,pa,para,dim,tree)
pdimpa=prod(dim(pa));
para_pa=para/pdimpa;
para_pav=para_pa/dim(var);
ins=find([pa,inf]>var,1,'first');
sub=[pa(1:ins-1),var,pa(ins:end)];
tab_pav=my_MakeContab(sub,dim,tree);
tab_pa=sum(tab_pav,ins);
local_score=pdimpa*gammaln(para_pa)-sum(reshape(gammaln(para_pa+tab_pa),[],1))...
    +sum(reshape(gammaln(para_pav+tab_pav),[],1))-pdimpa*dim(var)*gammaln(para_pav);
end