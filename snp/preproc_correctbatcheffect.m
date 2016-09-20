function [C,P,batches,snps] = preproc_correctbatcheffect(C,CBEparams,output_dir)
%[C,P,batches,snps] = preproc_correctbatcheffect(C,CBEparams)
%
%       Revisions:
%           3 Dec 07:  Added reorder_D_cols before batch correction to
%           improve data write to HDF5 after batch correct
% 
% bcdir = [output_dir 'batchcorrection/'];
% mkdir(bcdir);

verbose('Correcting Batch Effect',20)




if ~iscell(C)
    Ctmp{1} = C;
    clear C;
    C = Ctmp;
    wascell = 0;
else
    wascell = 1;
end


for k = 1:length(C)


    %  ADD BATCH to SUPDAT (batch numbers should be assigned after plates merge)
    [ub,ui,uj]=unique(strvcat(get_sis(C{k},'batch')),'rows');
    C{k}=add_D_sup(C{k},'BATCH','Batch',uj','cols');
    % Order columns by batch to speed disk write after correction
    batrow = find(strmatch('BATCH',C{k}.supacc));
    [s,batidx] = sort(C{k}.supdat(batrow,:));
    C{k} = reorder_D_cols(C{k},batidx);
    [C{k},P,batches] = correct_batch_effect_new(C{k},CBEparams);
    C{k} = add_history(C{k},CBEparams);
%     save([bcdir 'P_D' num2str(k)],'P');

end

if ~wascell
    C = C{1};
end
