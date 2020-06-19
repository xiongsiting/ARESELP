function [nums,inlist] = loadfilelist(txtname)

fid = fopen(txtname);
data=textscan(fid,'%s ');
fclose(fid);

data = data{1}; 
idx = cell({});
for i = 1:size(data,1)
   if strcmp(data{i},'<<')
       idx{i,1} = i;
   end
   if strcmp(data{i},'>>')
       idx{i,2} = i;
   end
end
idx = cell2mat(idx);
nf = size(idx,1);
nums = zeros(nf,1);
inlist = cell(nf,1);
%outlist = cell(nf,1);
for i = 1:nf
    idxin = idx(i,1);
    nums(i) = str2num(data{idxin + 1});
    for j = 1:nums(i)
        inlist{i,j} = data{idxin + 1 + j};
    end
    %idxout = idx(i,2);
    %outlist{i} = data{idxout + 1};
    
end

end


