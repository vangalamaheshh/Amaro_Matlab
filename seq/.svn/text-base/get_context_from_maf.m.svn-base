function n = get_context_from_maf(fname,ncat,context_field_name)
% function n = get_context_from_maf(fname,ncat)
%
% Mike Lawrence 2010-01-27

if ~exist('context_field_name','var'), context_field_name='context'; end


M = load_struct(fname);

n = get_context_from_mafstruct(M,ncat,context_field_name);
