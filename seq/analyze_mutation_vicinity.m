function analyze_mutation_vicinity(x,outname,P)
% analyze_mutation_vicinity(x,outname,P)
%
% for each of (A,C,G,T) -> (A,C,G,T)
%    make a sequence-logo-type diagram to show the local sequence context

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'radius',10);

x = require_fields_with_convert(x,{'chr','start','end'},...
   {'Chromosome','Start_position','End_position'});

if ~isnumeric(x.chr), x.chr = convert_chr(x.chr); end
x = make_numeric(x,{'start','end'});

