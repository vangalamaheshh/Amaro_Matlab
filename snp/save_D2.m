function [D] = save_D2(fname,D)
%save_D save gistic data object D in file fname.
%
%   D = save_D(fname,D)
%
%    Revisions:
%
%           11 Dec 07:  Jen Dobson (jdobson@broad.mit.edu)
%           Added conversion of SIS field and file pointer redirection for
%           datastruct object.

%---
% $Id$
% $Date: 2008-08-12 12:14:38 -0400 (Tue, 12 Aug 2008) $
% $LastChangedBy: jdobson $
% $Rev$

if isa(D,'datastruct')
D = rewriteD(D);
end


if ~ischar(fname)
    error('FNAME must be a string');
end

if ~isstruct(D) && ~isa(D,'datastruct') && ~iscell(D)
    error('D must be a struct, datastruct, or cell of structs/datastructs');
end

if ~is_full_filename(fname)
    fname = [pwd filesep fname];
end

if isa(D,'cell')

    isstct  = cellfun(@isstruct,D);

    if isequal(isstct,ones(size(isstct)))  %if everything's a struct

        s = whos('D');

        if s.bytes < 2000000000
            save(fname,'D');
        else
            save(fname,'D','-v7.3')
        end


    else  %if at least one datastruct

        if ~isequal(isstct,zeros(size(isstct)))  %turn everything into datastruct
            for k = find(isstct)
                D{k} = datastruct(D{k});
            end
        end

        newloc = cell(2,1,1);
        oldloc = cell(2,1,1);
        for cc = 1:length(D)
            Ds{cc} = prepare_datastruct(D{cc});  %#ok

            [fn,ds,fieldnames] = gethdf5info_lite(Ds{cc});
            idx = 1:length(fn);
            if ~isempty(oldloc{1}) && cc~= 1
                idx = idx + size(oldloc,3);
            end

            oldloc(1,1,idx) = fn;  %#ok
            oldloc(2,1,idx) = ds;  %#ok
            newloc(2,1,idx) = ds;
            newloc(2,1,idx) = strcat(newloc(2,1,idx),'_');
            newloc(2,1,idx) = regexprep(newloc(2,1,idx),'_.*',['_' num2str(cc)]);
            newloc(1,1,idx) = cellstr(fname);  %#ok

            Ds{cc} = sethdf5info(Ds{cc},fieldnames,'Name',{newloc{2,1,idx}});  %#ok
            Ds{cc} = sethdf5info(Ds{cc},fieldnames,'Filename',{newloc{1,1,idx}});  %#ok

        end

        save(fname,'Ds','-v7.3')

        for k = 1:size(newloc,3)

            unixstr = ['~/CancerGenomeAnalysis/trunk/C/transh5toexisting18 ' oldloc{1,1,k} ' ' oldloc{2,1,k} ' ' newloc{1,1,k} ' ' newloc{2,1,k}];

            status = unix(unixstr);

        end

    end


elseif isa(D,'struct')

    s = whos('D');

    if s.bytes < 2000000000
        save(fname,'D');
    else
        save(fname,'D','-v7.3')
    end


    elseif isa(D,'datastruct')

        Ds = prepare_datastruct(D);
     
        
        
        [fn,ds,fieldnames] = gethdf5info_lite(Ds);
        
        if ~isempty(fn)  %If has diskfields

            idx = 1:length(fn);
            oldloc(1,1,idx) = fn;
            oldloc(2,1,idx) = ds;
            newloc(2,1,idx) = ds;
            newloc(1,1,idx) = cellstr(fname);  %#ok

            Ds = sethdf5info(Ds,fieldnames,'Name',newloc(2,1,idx));  %#ok
            Ds = sethdf5info(Ds,fieldnames,'Filename',newloc(1,1,idx));  %#ok
            
            save(fname,'Ds','-v7.3')

            for k = 1:size(newloc,3)

                unixstr = ['~/CancerGenomeAnalysis/trunk/C/transh5toexisting18 ' oldloc{1,1,k} ' ' oldloc{2,1,k} ' ' newloc{1,1,k} ' ' newloc{2,1,k}];

                status = unix(unixstr);

            end

        end

end

end

function Dps = prepare_struct(Dps)

Dps.convertedfields = [];

fields = fieldnames(Dps);

for fld = fields'

    if isnumeric(Dps.(char(fld))) || ischar(Dps.(char(fld)))

        continue;

    elseif iscell(Dps.(char(fld)))

        [Dps.(char(fld)),suc] = convert_cell_to_char(Dps.(char(fld)));
        if suc
            Dps.convertedfields = strvcat(Dps.convertedfields,char(fld));
        end


    elseif isstruct(Dps.(char(fld)))

        tmpstruct = strvcat(fieldnames(Dps.(char(fld))));

        tmpstruct = strvcat(tmpstruct,'>>>END FIELDS<<<');

        tmpstruct = strvcat(tmpstruct,strvcat(struct2cell(Dps.(char(fld)))));

        Dps.(char(fld)) = strvcat(['STRUCT:::' num2str(size(Dps.(char(fld)),1)) ' x '...
            num2str(size(Dps.(char(fld)),2))],tmpstruct);

        Dps.convertedfields = strvcat(Dps.convertedfields,char(fld));
    end


end

end


function Dpd = prepare_datastruct(Dpd)
Dpd.convertedfields = [];

fields = fieldnames(Dpd)';
dskflds = diskfieldnames(Dpd);

if size(fields,1)>size(fields,2)
    fields = fields';
end


for fld = setdiff(fields,dskflds)

    if isnumeric(Dpd.(char(fld))) || ischar(Dpd.(char(fld)))

        continue;

    elseif iscell(Dpd.(char(fld)))

        [Dpd.(char(fld)),suc] = convert_cell_to_char(Dpd.(char(fld)));
        if suc
            Dpd.convertedfields = strvcat(Dpd.convertedfields,char(fld));
        end

    elseif isstruct(Dpd.(char(fld)))

        tmpstruct = strvcat(fieldnames(Dpd.(char(fld))));

        tmpstruct = strvcat(tmpstruct,'>>>END FIELDS<<<');

        fielddata = convert_cell_to_char(struct2cell(Dpd.(char(fld))));

        tmpstruct = strvcat(tmpstruct,fielddata);

        Dpd.(char(fld)) = strvcat(['STRUCT:::' num2str(size(Dpd.(char(fld)),1)) ' x '...
            num2str(size(Dpd.(char(fld)),2))],tmpstruct);

        Dpd.convertedfields = strvcat(Dpd.convertedfields,char(fld));
    end


end

end



function [dfld,success] = convert_cell_to_char(dfld)

success = 0;
cellclasses = cellfun(@class,dfld,'UniformOutput',0);

if isempty(find(~strcmp('char',cellclasses)))  %are all cells 'char'?

    header = ['CELL of CHAR:::' num2str(size(dfld,1)) ' x ' ...
        num2str(size(dfld,2))];
    dfld = strvcat(header, strvcat( dfld));

    success = 1;

elseif isempty(find(~strcmp('cell',cellclasses)))

    header = ['CELL of CELL:::' num2str(size(dfld,1)) ' x ' ...
        num2str(size(dfld,2)) ];

    if isvector(cellclasses)

        tmpfld = header;

        for k = 1:length(cellclasses)

            tmpcell = dfld{k};

            subcellclasses = cellfun(@class,tmpcell,'UniformOutput',0);

            if isempty(find(~strcmp('char',subcellclasses)))

                header = ['CELL of CHAR:::' num2str(size(tmpcell,1)) ' x ' ...
                    num2str(size(tmpcell,2))];
                tmpcell = strvcat(header, strvcat( tmpcell));

                tmpfld = strvcat(tmpfld,tmpcell);

            elseif isempty(find(~strcmp('cell',subcellclasses)))

                header = ['CELL of CELL:::' num2str(size(tmpcell,1)) ' x ' ...
                    num2str(size(tmpcell,2))];
                
                tmpcell = convert_cell_to_char(tmpcell);
                
                tmpfld = strvcat(tmpfld,tmpcell);
                
            else
                
                return

            end
        end
        dfld = tmpfld;

        success = 1;
    end



end




end



