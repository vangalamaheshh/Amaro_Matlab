function N = get_covered_territory_tally_from_FWB(cov_fwb,categ_fwb)

fprintf('Getting covered territory tally... ');

S = org.broadinstitute.cga.tools.seq.FixedWidthBinary(cov_fwb);
C = org.broadinstitute.cga.tools.seq.FixedWidthBinary(categ_fwb);
swidth = S.getWidth();

categ_list = regexprep(categ_fwb,'^(.*)/([^/]+)$','$1/categs.txt');
if exist(categ_list,'file')
  z = load_struct(categ_list);
  if isfield(z,'num')
    tmp = str2double(z.num);
    if any(tmp<1)
      error('This function will not work properly with tracks that have category numbers less than 1');
    end
  end
  ncat = slength(z);
else
  ncat = C.maxValForWidth();
end

N = zeros(ncat,1);

maxchunk = 1e6;
nr = S.getNumRegions;
step = round(nr/100);
for i=1:nr, if ~mod(i,step), fprintf('%d/%d ',i,nr); end
  chr = S.getRegionChr(i-1);   % java is 0-based
  if (chr<1 || chr>24), continue; end
  st = S.getRegionStart(i-1);
  en = S.getRegionEnd(i-1);
  for st2=st:maxchunk:en
    en2 = min(en,st2+maxchunk-1);
    cov = double(S.get(chr,st2,en2));
    cov(cov==-1) = 0;  % position with no data -> change to zero coverage
    cat = double(C.get(chr,st2,en2));
    if swidth==1 %faster
      cat = cat(cov==1);
      if ~isempty(cat)
        if length(cat)==1
          if cat>=1 && cat<=ncat
            N(cat)=N(cat)+1;
          end
        else
          N=N+histc(cat,1:ncat);
        end
      end
    else %slower
      for j=1:length(cat)
        if cat(j)>=1 & cat(j)<=ncat
          N(cat(j)) = N(cat(j)) + cov(j);
        end
      end
    end
  end
end

fprintf('\n');

S.close();
C.close();
