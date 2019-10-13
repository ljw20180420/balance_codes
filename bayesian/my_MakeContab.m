function tab=my_MakeContab(sub,dim,tree)
cumdim=cumprod(dim(sub));
ptr=0;
str=0;
val=zeros(size(sub));
mcv=zeros(size(sub));
ro_ptr=zeros(size(sub));
tab=zeros([dim(sub),1]);
while true
%     ptr
%     str
%     val
%     mcv
%     ro_ptr
    if ptr>0
        flag=(tree(ro_ptr(ptr),1)==0 || ptr==length(sub));
    else
        flag=(tree(1,1)==0 || ptr==length(sub));
    end
    if flag
        if ptr==length(sub)
            tab(high_dim_index(cumdim,val))=tree(ro_ptr(ptr),1);
        end
        while ptr>0 && flag
            if val(ptr)<dim(sub(ptr))
                val(ptr)=val(ptr)+1;
                str=find(val(1:ptr-1)~=mcv(1:ptr-1),1,'last');
                if isempty(str)
                    str=0;
                end
                flag=false;
            else
                tape=reshape(tab(high_dim_index(cumdim,val(1:ptr-1))),dim(sub(ptr)),[]);
                tab(high_dim_index(cumdim,[val(1:ptr-1),mcv(ptr)]))=tape(mcv(ptr),:)-sum(tape([1:mcv(ptr)-1,mcv(ptr)+1:end],:),1);
                ptr=ptr-1;
            end
        end
        if ptr==0
            break;
        end
    else
        ptr=ptr+1;
        val(ptr)=1;
    end
    
    if str>0
        ha_ptr=ro_ptr(str);
        pa_ptr=ha_ptr+sub(ptr)-sub(str);
    else
        ha_ptr=1;
        pa_ptr=ha_ptr+sub(ptr);
    end
    tape=ha_ptr+tree(pa_ptr,1)-1;
    mcv(ptr)=tree(tape,1);
    
    if val(ptr)==mcv(ptr)
        if str>0
            ro_ptr(ptr)=ro_ptr(str);
        else
            ro_ptr(ptr)=1;
        end
    else
        ro_ptr(ptr)=tape+tree(tape+val(ptr)-(val(ptr)>mcv(ptr)),1)-1;
        str=ptr;
    end
end
end