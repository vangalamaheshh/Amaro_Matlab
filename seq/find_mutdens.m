function mutdens = find_mutdens(chr,pos,chr2,pos2)
% find_mutdens(chr,pos)
%
% for each mutation, finds number of mutations within 20bp
%
% find_mutdens(chr1,pos1,chr2,pos2)
%
% for each mutation (chr2,pos2), finds number of mutations (chr1,pos1) within 20bp
%
%

win = 20;

if length(chr)~=length(pos), error('length(chr)~=length(pos)'); end

if nargin==2
  chr2=chr;
  pos2=pos;
elseif nargin==4
  if length(chr2)~=length(pos2), error('length(chr2)~=length(pos2)'); end
else
  error('wrong argument structure');
end

if any(pos<1) || any(pos2<1), error('pos must be >=1'); end

n = length(chr2);
mutdens = nan(n,1);

[u ui uj] = unique(chr);
vj = listmap(chr2,u);
for i=1:length(u), fprintf('chr%d ',u(i));
  idx = find(uj==i);
  idx2 = find(vj==i);
  if isempty(idx) || isempty(idx2), continue; end
  mx = max(pos(idx));
  mx2 = max(pos2(idx2));
  mx = max(mx,mx2);
  h = histc(pos(idx),round(-win/2)+1:mx+round(win/2));
  ch = cumsum(h);
  d = ch(win+1:end)-ch(1:end-win);
  mutdens(idx2)=d(pos2(idx2));
end, fprintf('\n');


