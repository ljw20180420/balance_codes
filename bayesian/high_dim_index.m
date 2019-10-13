function ind=high_dim_index(cumdim,sub)
if isempty(sub)
    ind=1:cumdim(end);
else
    L=length(sub);
    ini=sub(1)+(sub(2:end)-1)*cumdim(1:L-1)';
    ind=ini:cumdim(L):cumdim(end);
end
end