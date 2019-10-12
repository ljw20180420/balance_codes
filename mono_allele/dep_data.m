function [Data_dep,Counts_dep,Data_merge,Counts_merge]=dep_data(Data,Counts,ord,tol,uprow)
v_num=size(Data,2);
Data_dep=cell(size(Data,1)+1,1);
Data_dep{end}=zeros(1,v_num);
for ii=1:size(Data,1)
    ind=find(Data(ii,:));
    if length(ind)<=2*ord
        Data_dep{ii}=zeros(sum(arrayfun(@(x) nchoosek(length(ind),x),1:min(ord,length(ind)))),v_num);
        kk=0;
        for oo=1:min(ord,length(ind))
            combos=nchoosek(ind,oo);
            for co=1:size(combos,1)
                kk=kk+1;
                Data_dep{ii}(kk,combos(co,:))=1;
            end
        end
    end
end
[Data_dep,~,IC]=unique(cell2mat(Data_dep),'row');
Data_multi=accumarray(IC,ones(size(IC)),[],[],[],false);
[~,sid]=sort(Data_multi,'descend');
Data_dep=Data_dep(sid,:);
Data_dep(uprow+1:end,:)=[];
c_num=size(Data_dep,1);
Counts_pair=zeros(c_num);
Data_merge=zeros(c_num^2,v_num);
for d1=1:c_num
    for d2=1:c_num
        Data_merge(d1+(d2-1)*c_num,:)=(Data_dep(d1,:)+Data_dep(d2,:)>0);
    end
end
[LIA,LOCB]=ismember(Data_merge,Data,'row');
LOCB=LOCB(LIA);
[~,~,IC]=unique(LOCB);
tape=accumarray(IC,1,[],[],[],false);
Counts_pair(LIA)=Counts(LOCB)./tape(IC);

toln=inf;
while toln>tol
    % M-step
    Counts_dep=sum(Counts_pair,2);
    Counts_dep=Counts_dep/sum(Counts_dep);
    % E-step
    Counts_pair_old=Counts_pair;
    tape=Counts_dep.*Counts_dep';
    wei=tape(LIA);
    tape=accumarray(IC,wei,[],[],[],false);
    Counts_pair(LIA)=Counts(LOCB).*wei./tape(IC);
    toln=max(max(abs(Counts_pair-Counts_pair_old)));
    toln
end

Counts_merge=reshape(Counts_dep.*Counts_dep',[],1);
Data_dep(Counts_dep==0,:)=[];
Counts_dep(Counts_dep==0)=[];
Data_merge(Counts_merge==0,:)=[];
Counts_merge(Counts_merge==0)=[];
[Data_merge,~,IC]=unique(Data_merge,'row');
Counts_merge=accumarray(IC,Counts_merge,[],[],[],false);
end