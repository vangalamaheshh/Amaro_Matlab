function D = convert_struct(D,opt)
%D = CONVERT_STRUCT(D)
%
%This function converts fields that were turned into char arrays for saving back into structs or cells

if isfield(D,'convertedfields') && ~isempty(D.convertedfields)

    for fld = cellstr(D.convertedfields)'

        if strmatch('STRUCT',D.(char(fld))(1,:))==1

            fldnamesend = strmatch('>>>END FIELDS<<<',D.(char(fld)));
            structfieldnames = cellstr(D.(char(fld))(2:fldnamesend-1,:));
            structfields = cellstr(D.(char(fld))(fldnamesend+1:end,:));

            if ~isempty(strmatch('CELL of CHAR',structfields{1}))
                structfields = {structfields{2:end}}';
                dims = regexp(D.(char(fld))(1,:),'\d+','match');
                dims = cellfun(@str2num,dims);
                tmpcell = reshape(structfields,length(structfieldnames),size(structfields,1)/length(structfieldnames));
                D.(char(fld)) = cell2struct(tmpcell,structfieldnames,1);

                if ~isequal(size(D.(char(fld))), dims)
                    D.(char(fld)) = D.(char(fld))';
                end



            elseif ~isempty(strmatch('CELL of CELL',structfields{1}))
             
                dims1 = {};


                while strmatch('CELL of CELL',structfields{1})
                    dims1 = [dims1 {cellfun(@str2num,regexp(structfields{1},'\d+','match'))}]; %#ok
                    structfields = structfields(2:end);
                end

                xx = strmatch('CELL of CELL',structfields);
                structfields = structfields(setdiff(1:length(structfields),xx));

                tmpcell = [];
                charstarts = strmatch('CELL of CHAR',structfields);
                for k = 1:length(charstarts)
                    if k ~= length(charstarts)
                        addflds = {structfields{charstarts(k)+1:charstarts(k+1)-1}}';
                    else
                        addflds = {structfields{charstarts(k)+1:end}}';
                    end

                    tmpcell = [tmpcell {addflds}];  %#ok
                end

                tmpcell = reshape(tmpcell,dims1{end-1}(1),dims1{end}(2));




                D.(char(fld)) = cell2struct(tmpcell,structfieldnames,1);

              
            else

                dims = regexp(D.(char(fld))(1,:),'\d+','match');
                dims = cellfun(@str2num,dims);
                tmpcell = reshape(structfields,length(structfieldnames),size(structfields,1)/length(structfieldnames));
                D.(char(fld)) = cell2struct(tmpcell,structfieldnames,1);

                if ~isequal(size(D.(char(fld))), dims)
                    D.(char(fld)) = D.(char(fld))';
                end




            end



        elseif ~isempty(strmatch('CELL of CHAR',D.(char(fld))(1,:)))

            dims = regexp(D.(char(fld))(1,:),'\d+','match');
            dims = cellfun(@str2num,dims);
            chardat = D.(char(fld))(2:end,:);
            D.(char(fld)) = cellstr(chardat);

            if ~isequal(size(D.(char(fld))), dims)
                D.(char(fld)) = D.(char(fld))';
            end


        elseif ~isempty(strmatch('CELL of CELL',D.(char(fld))(1,:))) && ~strcmp(opt,'skip-cell')
   disp('Converting cell array of cells -- this may take some time.  Use the ''skip-cell'' option to skip.')
            dims1 = regexp(D.(char(fld))(1,:),'\d+','match');
            dims1 = cellfun(@str2num,dims1);
            subcellstarts = strmatch('CELL of CHAR',D.(char(fld)));

            for k = 1:length(subcellstarts)
                dims2 = regexp(D.(char(fld))(subcellstarts(k),:),'\d+','match');
                dims2 = cellfun(@str2num,dims2);
                if k ~= length(subcellstarts)
                    subcell{k}=cellstr(D.(char(fld))((subcellstarts(k)+1):(subcellstarts(k+1)-2),:)); %#ok
                else
                    subcell{k}=cellstr(D.(char(fld))((subcellstarts(k)+1):end,:));  %#ok
                end
                if ~isequal(size(D.(char(fld))),dims2)
                    subcell{k} = subcell{k}';  %#ok
                end

            end
            D.(char(fld)) = subcell;

            if ~isequal(size(D.(char(fld))), dims1)
                D.(char(fld)) = D.(char(fld))';
            end
  disp('Conversion complete')
        end

    end


    D = rmfield(D,'convertedfields');
elseif isfield(D,'convertedfields') && isempty(D.convertedfields)
    D = rmfield(D,'convertedfields');
end

end  %Turns char arrays back into cells or structs


