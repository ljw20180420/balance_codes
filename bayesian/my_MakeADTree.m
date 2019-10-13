function tree=my_MakeADTree(dim,RecordNums,SizeUnit)
tree=zeros(SizeUnit,2);
tree_ptr=0;
Sub_ptr=0;
Sub=zeros(size(dim));
tree_Sub=zeros(size(dim));
ReRow=cell(size(dim));
Val_ptr=0;
Val=zeros(size(dim));
tree_Val=zeros(size(dim));

tree_ptr=tree_ptr+1;
if tree_ptr>size(tree,1)
    tree=[tree;zeros(ceil((tree_ptr-size(tree,1))/SizeUnit)*SizeUnit,size(tree,2))];
end
tree(tree_ptr,1)=size(RecordNums,1);
tree(tree_ptr,2)=length(dim);
tree_ptr=tree_ptr+tree(tree_ptr,2);

while true
    if Sub_ptr==Val_ptr
        if Val_ptr>0
            flag=tree(tree_Val(Val_ptr),2)==0;
        else
            flag=tree(1,2)==0;
        end
        if flag
            while Sub_ptr>0 && flag
                if Sub_ptr==Val_ptr
                    pa_ptr=tree_Sub(Sub_ptr)+Val(Val_ptr)-(Val(Val_ptr)>tree(tree_Sub(Sub_ptr),1));
                    tree(pa_ptr,2)=tree_ptr-tree_Sub(Sub_ptr)+1;
                    MCV=tree(tree_Sub(Sub_ptr),1);
                    if Val(Val_ptr)+(Val(Val_ptr)+1==MCV)<dim(Sub(Sub_ptr))
                        Val(Val_ptr)=Val(Val_ptr)+1+(Val(Val_ptr)+1==MCV);
                        flag=false;
                    else
                        Val_ptr=Val_ptr-1;
                    end
                else
                    if Val_ptr==0
                        ha_ptr=1;
                        pa_ptr=1+Sub(Sub_ptr);
                    else
                        ha_ptr=tree_Val(Val_ptr);
                        pa_ptr=tree_Val(Val_ptr)+Sub(Sub_ptr)-Sub(Sub_ptr-1);
                    end
                    tree(pa_ptr,2)=tree_ptr-ha_ptr+1;
                    if Sub(Sub_ptr)<length(dim)
                        Sub(Sub_ptr)=Sub(Sub_ptr)+1;
                        flag=false;
                    else
                        Sub_ptr=Sub_ptr-1;
                    end
                end
            end
            if Sub_ptr==0
                break;
            end
        else
            Sub_ptr=Sub_ptr+1;
            if Sub_ptr>1
                Sub(Sub_ptr)=Sub(Sub_ptr-1)+1;
            else
                Sub(Sub_ptr)=1;
            end
        end
    else
        MCV=tree(tree_Sub(Sub_ptr),1);
        Val_ptr=Val_ptr+1;
        if MCV>1
            Val(Val_ptr)=1;
        else
            Val(Val_ptr)=2;
        end
    end
    
    tree_ptr=tree_ptr+1;
    if tree_ptr>size(tree,1)
        tree=[tree;zeros(ceil((tree_ptr-size(tree,1))/SizeUnit)*SizeUnit,size(tree,2))];
    end
    if Sub_ptr>Val_ptr
        tree_Sub(Sub_ptr)=tree_ptr;
        if Val_ptr==0
            ha_ptr=1;
            pa_ptr=1+Sub(Sub_ptr);
        else
            ha_ptr=tree_Val(Val_ptr);
            pa_ptr=tree_Val(Val_ptr)+Sub(Sub_ptr)-Sub(Sub_ptr-1);
        end
        tree(pa_ptr,1)=tree_ptr-ha_ptr+1;
        
        if Val_ptr>0
            [~,tree(tree_ptr,1)]=max(histcounts(RecordNums(ReRow{Val_ptr},Sub(Sub_ptr)),(1:dim(Sub(Sub_ptr))+1)-0.5));
        else
            [~,tree(tree_ptr,1)]=max(histcounts(RecordNums(:,Sub(Sub_ptr)),(1:dim(Sub(Sub_ptr))+1)-0.5));
        end
        tree(tree_ptr,2)=dim(Sub(Sub_ptr))-1;
        tree_ptr=tree_ptr+tree(tree_ptr,2);
    else
        tree_Val(Val_ptr)=tree_ptr;
        pa_ptr=tree_Sub(Sub_ptr)+Val(Val_ptr)-(Val(Val_ptr)>tree(tree_Sub(Sub_ptr),1));
        tree(pa_ptr,1)=tree_ptr-tree_Sub(Sub_ptr)+1;
        if Val_ptr>1
            ReRow{Val_ptr}=ReRow{Val_ptr-1}(RecordNums(ReRow{Val_ptr-1},Sub(Sub_ptr))==Val(Val_ptr));
        else
            ReRow{Val_ptr}=find(RecordNums(:,Sub(Sub_ptr))==Val(Val_ptr));
        end
        
        tree(tree_ptr,1)=length(ReRow{Val_ptr});
        if tree(tree_ptr,1)>0
            tree(tree_ptr,2)=length(dim)-Sub(Sub_ptr);
        end
        tree_ptr=tree_ptr+tree(tree_ptr,2);
    end
end
tree(tree_ptr+1:end,:)=[];
end