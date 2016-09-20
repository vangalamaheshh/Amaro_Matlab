function convert_gaf_to_refseq_format(in,out)
% convert_gaf_to_refseq_format(in,out)
%
% in = GAF, e.g. transcript.genome.v3_0.gaf
% out = file with same format as R.mat, suitable for loading with load_refseq

n = get_numcols(in);
if n~=19, error('Was expecting GAF to have 19 columns... please check format'); end
f = fopen(in);
z = textscan(f,repmat('%s',1,n),'headerLines',0,'delimiter',char(9),'bufsize',5e6);
fclose(f);

if isempty(grep('^chr',z{15},1))
  error('please use transcript.genome.v*.gaf instead of transcript.gene.v*.gaf');
end

% relevant columns:
% col14 = cDNA coordinates
% col15 = DNA coordinates
% col16 = gene|ID
% col17 = chr:txstart-txend:strand
% col18 = NM_xxx
% col19 = Confidence,CDSstart,CDSstop

R=[];
R.col14 = z{14};
R.col15 = z{15};
R.col16 = z{16};
R.col17 = z{17};
R.col18 = z{18};
R.col19 = z{19};

% Need to convert to R.mat:
%             id: 778
%     transcript: NM_004985
%            chr: chr12
%         strand: -
%       tx_start: 25358180
%         tx_end: 25403854
%     code_start: 25362729
%       code_end: 25398318
%        n_exons: 5
%    exon_starts: 25358180 25378548 25380168 25398208 25403685
%      exon_ends: 25362845 25378707 25380346 25398329 25403854
%           gene: KRAS
%    exon_frames: 0 2 0 0 -1
%        version: n/a
%         tx_len: 5297
%       code_len: 567
%       n_codons: 189

R = parse_in(R,'col16','^([^\|]+)\|(.*)$',{'gene','id'}); R = rmfield(R,'col16');
R = rename_field(R,'col18','transcript');
R = parse_in(R,'col17','^(chr[^:]+):(\d+)\-(\d+):(\+|\-)',{'chr','tx_start','tx_end','strand'},2:3);
R.numloci = nan(slength(R),1); for i=1:slength(R), R.numloci(i) = sum(R.col17{i}==';')+1; end
R = parse_in(R,'col19','Confidence=(\d+)','confidence',1);
R = parse_in(R,'col19','CDSstart=(\d+)','code_start',1);
R = parse_in(R,'col19','CDSstop=(\d+)','code_end',1);
R = rmfield(R,{'col17','col19'});

% get rid of junk we are not interested in

R = reorder_struct_exclude(R,R.numloci>3);
R = reorder_struct_exclude(R,R.confidence==100);
R = reorder_struct_exclude(R,isnan(R.code_start)|isnan(R.code_end));















