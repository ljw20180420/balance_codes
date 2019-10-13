function [rcmv,capture_ratios,vps,XTLE,XTL,XT]=load_target_CTCF(pos_file)
[~,ln]=system(['wc -l ',pos_file]);
ln=split(ln);
ln=str2double(ln{1});
left_stall=zeros(ln,1);
right_stall=zeros(ln,1);
rcmv=zeros(ln,4);
vps=[];
capture_ratios=[];
fid=fopen(pos_file,'r');
line=fgetl(fid);
while ~feof(fid)
    if strcmp(line,'left_stall')
        ii=0;
        while ~feof(fid)
            line=fgetl(fid);
            words=split(line);
            if isnan(str2double(words))
                break;
            else
                ii=ii+1;
                left_stall(ii)=str2double(words);
            end
        end
    else
        if strcmp(line,'right_stall')
            jj=0;
            while ~feof(fid)
                line=fgetl(fid);
                words=split(line);
                if isnan(str2double(words))
                    break;
                else
                    jj=jj+1;
                    right_stall(jj)=str2double(words);
                end
            end
        else
            if strcmp(line,'target')
                kk=0;
                while ~feof(fid)
                    line=fgetl(fid);
                    words=split(line);
                    if length(words)~=4
                        break;
                    else
                        kk=kk+1;
                        rcmv(kk,:)=str2double(words)';
                    end
                end
            else
                if strcmp(line,'VP')
                    while ~feof(fid)
                        line=fgetl(fid);
                        tmp=str2double(line);
                        if(isnan(tmp))
                            break;
                        else
                            vps=[vps,tmp];
                        end
                    end
                else
                    if strcmp(line,'capture_ratio')
                        while ~feof(fid)
                            line=fgetl(fid);
                            tmp=str2double(line);
                            if(isnan(tmp))
                                break;
                            else
                                capture_ratios=[capture_ratios,tmp];
                            end
                        end
                    else
                        line=fgetl(fid);
                    end
                end
            end
        end
    end
end
fclose(fid);
left_stall(ii+1:end)=[];
right_stall(jj+1:end)=[];
rcmv(kk+1:end,:)=[];
chain_length=length(left_stall);

XTLE=repmat({''},chain_length,1);
XTLE(left_stall>0 & right_stall==0)={char(8595)};
XTLE(left_stall==0 & right_stall>0)={char(8593)};
XTLE(left_stall>0 & right_stall>0)={char(8597)};
for vv=1:length(vps)
    if(left_stall(vps(vv)+1)>0)
        XTLE(vps(vv)+1)={char(8650)};
    else
        if (right_stall(vps(vv)+1)>0)
            XTLE(vps(vv)+1)={char(8648)};
        end
    end
end
tmp=false(chain_length,1);
for aa=1:chain_length
    if(~isempty(XTLE{aa}))
        tmp(aa)=true;
    end
end
XTL=XTLE(tmp);
for aa=1:length(XTL)
    XTL{aa}=char(XTL{aa}-1);
end
XT=find(tmp);
end