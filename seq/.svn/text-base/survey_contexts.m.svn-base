function c = survey_contexts(seq,context_type);
% returns vector c telling the context of each base along dna
% (edges are assigned "0" meaning "unknown"
%
% if context_type is ommited, default of 0 is imposed.
%
% Base categories:
% (all with respect to the transcribed strand)
%
% context_type 0
%
%     1. A
%     2. T
%     3. C in CpG
%     4. C in TpC but not TpCpG
%     5. other C
%     6. G in CpG
%     7. G in GpA but not CpGpA
%     8. other G
%
% context_type 1
%
%     A_A = 1     A_C = 2     A_G = 3     A_T = 4
%     C_A = 5     C_C = 6     C_G = 7     C_T = 8
%     G_A = 9     G_C = 10    G_G = 11    G_T = 12
%     T_A = 13    T_C = 14    T_G = 15    T_T = 16
%

if ~exist('context_type','var'), context_type=0; end

seq = upper(seq);
tlen = length(seq);

c = zeros(1,tlen);

if context_type==0

      for bp = 2:tlen-1
           if seq(bp) == 'A'
               category = 1;                        % A
           elseif seq(bp) == 'T'
               category = 2;                        % T
           elseif seq(bp) == 'C'
               if seq(bp+1) == 'G'
                    category = 3;                   % C in CpG
               elseif seq(bp-1) == 'T'
                    category = 4;                   % C in TpC but not TpCpG
               else
                    category = 5;                   % other C
               end
           elseif seq(bp) == 'G'
               if seq(bp-1) == 'C'
                    category = 6;                   % G in CpG
               elseif seq(bp+1) == 'A'
                    category = 7;                   % G in GpA but not CpGpA
               else
                    category = 8;                   % other G
               end
           else
               category = 9;                        % "N"?
           end
           c(bp)=category;
        end

elseif context_type==1

  for bp = 2:tlen-1
    x=0;
    if seq(bp-1)=='A', x=1;
    elseif seq(bp-1)=='C', x=5;
    elseif seq(bp-1)=='G', x=9;
    elseif seq(bp-1)=='T', x=13;
    else error('unknown base'); end
    if seq(bp+1)=='A', x=x+0;
    elseif seq(bp+1)=='C', x=x+1;
    elseif seq(bp+1)=='G', x=x+2;
    elseif seq(bp+1)=='T', x=x+3;
    else error('unknown base'); end
    c(bp)=x;
  end

else
  error('unknown context_type');
end
