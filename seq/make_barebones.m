function [mut P] = make_barebones(mut,P)

if ~exist('P','var'), P=[]; end

if ~isfield(mut,'patient') && ~isfield(mut,'pat_idx')
  mut.patient = regexprep(mut.Tumor_Sample_Barcode,'-Tumor$','');
end

if ~isfield(mut,'Variant_Classification')
  mut.Variant_Classification = repmat({''},slength(mut),1);
end

if isfield(P,'build')
  pbuild = P.build;
  if isfield(mut,'NCBI_Build')
    mbuild = mut.NCBI_Build{1};
  elseif isfield(mut,'Build')
    mbuild = mut.Build{1};
  end
end
if ~exist('pbuild','var') && ~exist('mbuild','var')
  fprintf('Assuming hg19\n');
%  fprintf('Is this OK? (dbcont/dbquit)'); keyboard
  P.build = 'hg19';
elseif ~exist('pbuild','var') && exist('mbuild','var')
  fprintf('Inferring build %d from MAF');
  P.build = mbuild;
elseif exist('pbuild','var') && ~exist('mbuild','var')
  % (no problem)
else % both exist
  if ~strcmpi(pbuild,mbuild) && interpret_build(pbuild)~=interpret_build(mbuild)
    fprintf('WARNING: Possible build conflict:\tP=%s\n\tmaf=%s\n',pbuild,mbuild);
  end
end

b = interpret_build(P.build);
if b==18, P.build = 'hg18'; end
if b==19, P.build = 'hg19'; end 

try
  mut = add_and_convert_simple_fieldnames(mut,P);
catch me
  disp(me); disp(me.message);
  fprintf('WARNING: problem with add_and_convert_simple_fieldnames\n');
end
mut.newbase = find_newbase(mut);

if ~isfield(mut,'dataset')
  mut.dataset = repmat({''},slength(mut),1);
end

if ~isfield(mut,'i_tumor_f')
  mut.i_tumor_f = nan(slength(mut),1);
end
if ~isfield(mut,'Protein_Change')
  mut.Protein_Change = repmat({''},slength(mut),1);
end

flds = {'dataset','patient','gene','chr','start','end','type',...
        'classification','ref_allele','newbase','Protein_Change',...
        'i_tumor_f','t_alt_count','t_ref_count',...
        'context65'};

if isfield(mut,'pat_idx') && ~isfield(mut,'patient')
  flds = regexprep(flds,'^patient$','pat_idx');
end

mut = keep_fields_if_exist(mut,flds);


