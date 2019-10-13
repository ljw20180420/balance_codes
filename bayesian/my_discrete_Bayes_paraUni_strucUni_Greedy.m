function DAG=my_discrete_Bayes_paraUni_strucUni_Greedy(DAG,skeleton,para,dim,tree,list_length,endur_fail)
delta_score=zeros(size(DAG));
cpa=cell(1,length(dim));
for jj=1:length(dim)
    cpa{jj}=find(skeleton(:,jj));
end
TABU_list=zeros([size(DAG),list_length]);
active_log=false(list_length,1);
tl=1;

now_score=0;
fail_num=0;
max_score=now_score;
max_tl=tl;
TABU_list(:,:,tl)=DAG;
active_log(tl)=true;
tl=mod(tl+1,list_length);

for jj=1:length(dim)
    pa_old=reshape(find(DAG(:,jj)),1,[]);
    local_score_old=BDeU_local_score(jj,pa_old,para,dim,tree);
    cpa_now=cpa{jj};
    for ii=cpa_now'
        if DAG(ii,jj)>0
            pa=setdiff(pa_old,ii);
            delta_score(ii,jj)=BDeU_local_score(jj,pa,para,dim,tree)-local_score_old;
        else
            ins=find([pa_old,inf]>ii,1,'first');
            pa=[pa_old(1:ins-1),ii,pa_old(ins:end)];
            delta_score(ii,jj)=BDeU_local_score(jj,pa,para,dim,tree)-local_score_old;
        end
    end
end

max_delta_score=-inf;
for jj=1:length(dim)
    cpa_now=cpa{jj};
    for ii=cpa_now'
        DAGij=DAG(ii,jj);
        DAGji=DAG(jj,ii);
        DAG(ii,jj)=1-DAG(ii,jj);
        DAG(jj,ii)=0;
        if (DAGij==1 || isdag(digraph(DAG))) && ~any(all(all(TABU_list(:,:,active_log)==DAG,1),2))
            delta_score_total=delta_score(ii,jj);
            if DAGji==1
                delta_score_total=delta_score_total+delta_score(jj,ii);
            end
            if delta_score_total>max_delta_score
                max_delta_score=delta_score_total;
                max_ii=ii;
                max_jj=jj;
            end
        end
        DAG(ii,jj)=DAGij;
        DAG(jj,ii)=DAGji;
    end
end
max_delta_score
if max_delta_score==-inf
    DAG=TABU_list(:,:,max_tl);
    return;
end

now_score=now_score+max_delta_score;
if now_score<=max_score
    fail_num=fail_num+1;
else
    fail_num=0;
    max_score=now_score;
    max_tl=tl;
end

while fail_num<endur_fail
    DAGji=DAG(max_jj,max_ii);
    DAG(max_ii,max_jj)=1-DAG(max_ii,max_jj);
    DAG(max_jj,max_ii)=0;
    TABU_list(:,:,tl)=DAG;
    active_log(tl)=true;
    tl=mod(tl+1,list_length);
    
    pa_old=reshape(find(DAG(:,max_jj)),1,[]);
    local_score_old=BDeU_local_score(max_jj,pa_old,para,dim,tree);
    cpa_now=cpa{max_jj};
    for ii=cpa_now'
        if ii==max_ii
            delta_score(ii,max_jj)=-delta_score(ii,max_jj);
        else
            if DAG(ii,max_jj)>0
                pa=setdiff(pa_old,ii);
                delta_score(ii,max_jj)=BDeU_local_score(max_jj,pa,para,dim,tree)-local_score_old;
            else
                ins=find([pa_old,inf]>ii,1,'first');
                pa=[pa_old(1:ins-1),ii,pa_old(ins:end)];
                delta_score(ii,max_jj)=BDeU_local_score(max_jj,pa,para,dim,tree)-local_score_old;
            end
        end
    end
    if DAGji==1
        pa_old=reshape(find(DAG(:,max_ii)),1,[]);
        local_score_old=BDeU_local_score(max_ii,pa_old,para,dim,tree);
        cpa_now=cpa{max_ii};
        for ii=cpa_now'
            if ii==max_jj
                delta_score(ii,max_ii)=-delta_score(ii,max_ii);
            else
                if DAG(ii,max_ii)>0
                    pa=setdiff(pa_old,ii);
                    delta_score(ii,max_ii)=BDeU_local_score(max_ii,pa,para,dim,tree)-local_score_old;
                else
                    ins=find([pa_old,inf]>ii,1,'first');
                    pa=[pa_old(1:ins-1),ii,pa_old(ins:end)];
                    delta_score(ii,max_ii)=BDeU_local_score(max_ii,pa,para,dim,tree)-local_score_old;
                end
            end
        end
    end
    
    max_delta_score=-inf;
    for jj=1:length(dim)
        cpa_now=cpa{jj};
        for ii=cpa_now'
            DAGij=DAG(ii,jj);
            DAGji=DAG(jj,ii);
            DAG(ii,jj)=1-DAG(ii,jj);
            DAG(jj,ii)=0;
            if (DAGij==1 || isdag(digraph(DAG))) && ~any(all(all(TABU_list(:,:,active_log)==DAG,1),2))
                delta_score_total=delta_score(ii,jj);
                if DAGji==1
                    delta_score_total=delta_score_total+delta_score(jj,ii);
                end
                if delta_score_total>max_delta_score
                    max_delta_score=delta_score_total;
                    max_ii=ii;
                    max_jj=jj;
                end
            end
            DAG(ii,jj)=DAGij;
            DAG(jj,ii)=DAGji;
        end
    end
    max_delta_score
    if max_delta_score==-inf
        DAG=TABU_list(:,:,max_tl);
        return;
    end
    
    now_score=now_score+max_delta_score;
    if now_score<=max_score
        fail_num=fail_num+1;
    else
        fail_num=0;
        max_score=now_score;
        max_tl=tl;
    end
end
DAG=TABU_list(:,:,max_tl);
end