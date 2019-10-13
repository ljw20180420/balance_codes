function combos=my_nchoosek(CPC,ss)
if length(CPC)>1
    combos=nchoosek(CPC,ss);
else
    if ss==0
        combos=zeros(1,0);
    else
        combos=CPC;
    end
end
end
    