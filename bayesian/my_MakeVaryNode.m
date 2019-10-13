function tree=my_MakeVaryNode(dim,RecordNums)
if size(RecordNums,1)==0
    tree=[0,0];
else
    Childnums=cell(dim(1),1);
    for kk=1:dim(1)
        Childnums{kk}=RecordNums(RecordNums(:,1)==kk,2:end);
    end
    abs_Childnums=cellfun(@(x) size(x,1),Childnums);
    [~,MCV]=max(abs_Childnums);
    tree=cell(dim(1),1);
    for kk=1:dim(1)
        if kk~=MCV
            tree{kk}=my_MakeADTree(dim(2:end),Childnums{kk});
        end
    end
    tape=cellfun(@(x) size(x,1),tree);
    tape(MCV)=[];
    tape=cumsum(tape);
    tape=[[1;tape(1:end-1)+1],tape]+dim(1);
    tree=[MCV,dim(1)-1;tape;cell2mat(tree)];
end
end