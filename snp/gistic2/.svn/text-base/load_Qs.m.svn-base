function Qs = load_Qs( filename )
%LOAD_QS Load either old or new Qs segments into new Qs structure
%   Qs = load_Qs( filenname )
%   The data in the file FILENAME are loaded into Qs. If Qs is a cell
%   array (old form), it is translated into the new struct form.

    load(filename);
    if exist('Qs','var')
        if iscell(Qs)
            tempQs = Qs;
            qfields={'amp','del','aod','doa'};
            Qs = struct;
            for k = 1:min(length(tempQs),length(qfields))
                Qs(qfields(k)) = tempQs{k};
            end
        elseif ~isstruct(Qs)
            error('Qs in file ''%s'' is not a cell array or struct',filename);
        end
        Qs.header = {'chrn','pos_start (snp)','pos_end (snp)','scna_score',...
                     'sample_id','starting_cn_level','ending_cn_level','fract_chr_arm_length',...
                     'deconstruction_score','arm_score','EMPTY','amplitude'};
        
    else
        error('File ''%s'' does not contain a Qs.',filename);
    end
end

