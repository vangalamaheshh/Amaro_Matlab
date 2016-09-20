function d2 = copyD(d)
%copyD make a copy of a datastruct object
%
%   D2 = copyD(D)

midx = strmatch('memory',d.storagetype);
didx = strmatch('disk',d.storagetype);

if size(midx,1) > size(midx,2)
    midx = midx';
end

if size(didx,1) > size(didx,2)
    didx = didx';
end


dfielddata = [d.fielddata{didx}];

olddiskfiles = {dfielddata.Filename}; 
newdiskfiles = get_new_filenames(olddiskfiles);

[cpfiles,idx] = unique(olddiskfiles);

newfielddata = repmat({},1,length(didx));
for i = idx
    copyfile(olddiskfiles{i},newdiskfiles{i})
    newfielddata{i} = d.fielddata{didx(i)};
    newfielddata{i}.Filename = newdiskfiles{i};
end



d2 = datastruct('memfields',{d.fieldnames{midx}},{d.fielddata{midx}},...
    'diskfields',{d.fieldnames{didx}},newfielddata);

for k = didx
    idx = strmatch(d.fieldnames{k},d2.fieldnames,'exact');
    d2.rowmapping{idx} = d.rowmapping{k};
    d2.colmapping{idx} = d.colmapping{k};
end
