function [newSISstruct] = preproc_mergesis(Cin,mkuarray)
%
%           Revisions:
%               
%% Find the unique sis fields; verify fields exist and write them if not

if ~iscell(Cin)
    Ctemp = {Cin};
    Cin = Ctemp;
    clear Ctemp
 
end

sisfields = [];

for k = 1:length(Cin)
    sisfields = [sisfields; fieldnames(Cin{k}.sis)];
    SIS0{k} = Cin{k}.sis;
end

[s,i,j] = unique(sisfields);
sisfields = s;


for k = 1:length(SIS0)
    for fld = sisfields'
        if ~isfield(SIS0{k},fld)
            [SIS0{k}(char(fld))] = deal('Empty');
        end
    end
    SIS0{k} = orderfields(SIS0{k},sisfields);
end


%% Make sure sis samples are in the same order

emptystruct = SIS0{1}(1);

for fld = fieldnames(SIS0{1}(1))'
    emptystruct.(char(fld)) = 'Empty';
end


for k = 1 : length(SIS0)
    
    fillsis = setdiff([1:max(mkuarray{k})],[1:length(SIS0{k})]);
    SIS0{k}(fillsis) = emptystruct;
    SIS1{k} = SIS0{k}(mkuarray{k});
end


% 
% 
% for k = 1: length(SIS0)
%     tok = repmat({'_'},length(SIS0{k}),1);
%     [plat{k},sisnames{k}] = cellfun(@strtok,{SIS0{k}.array}',tok,'UniformOutput',0);
%     [match,mk,mu] = match_string_sets(sisnames{k},uniquenames);
%     missing = setdiff([1:length(uniquenames)],mu);
% 
%     for nm = uniquenames(missing)'
%         SIS0{k}(end+1) = emptystruct;
%         SIS0{k}(end).array = char(nm);
% 
%     end
%     SIS1{k}(mu) = SIS0{k}(mk);
%     SIS1{k}(missing) = SIS0{k}(length(mk)+1:end);
%     
% end



%% This is where we merge the fields

%To do: Array specific fields merge with '||'; sample specific fields are checked
%for same names


arrayspecific = {'array','dup','core','good','force','leaveout','rep','loh','platform','batch','histqc'} ;


cellSIS = cellfun(@struct2cell,SIS1,'UniformOutput',0);

colons = repmat({'###'},size(cellSIS{1}));

colSIS= cellfun(@strcat,cellSIS,repmat({colons},1,length(cellSIS)),'UniformOutput',0);


newSIScell = strcat(colSIS{:});

newSISstruct = cell2struct(newSIScell,sisfields,1)';
