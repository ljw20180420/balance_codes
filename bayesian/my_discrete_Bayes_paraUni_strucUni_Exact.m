function DAG_max=my_discrete_Bayes_paraUni_strucUni_Exact(para,dim,tree)
[II,JJ]=find(tril(true(length(dim)),-1));
tril_ind=II+(JJ-1)*length(dim);
triu_ind=JJ+(II-1)*length(dim);
DAG=zeros(length(dim));
var_score=zeros(1,length(dim));
for jj=1:length(dim)
    var_score(jj)=BDeU_local_score(jj,[],para,dim,tree);
end
now_score=sum(var_score);
DAG_max=DAG;
max_score=now_score;
iter_lim=3^length(triu_ind)-1
for iter=1:3^length(triu_ind)-1
    if mod(iter,1000)==1
        iter
    end
    DAG_old=DAG;
    tape=str2num(dec2base(iter,3)');
    tape=[zeros(length(tril_ind)-length(tape),1);tape];
    DAG=zeros(length(dim));
    DAG(tril_ind(tape==1))=1;
    DAG(triu_ind(tape==2))=1;
    if isdag(digraph(DAG))
        c_ind=reshape(find(any(DAG~=DAG_old,1)),1,[]);
        for jj=c_ind
            pa=reshape(find(DAG(:,jj)),1,[]);
            tape=BDeU_local_score(jj,pa,para,dim,tree);
            now_score=now_score+tape-var_score(jj);
            var_score(jj)=tape;
        end
        if now_score>max_score
            DAG_max=DAG;
            max_score=now_score;
        end
    else
        DAG=DAG_old;
    end
end
end