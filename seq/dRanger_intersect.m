function [x1 x2] = dRanger_intersect(R1,R2,tolerance)
if ~exist('tolerance','var'), tolerance = 100; end

flds = {'chr1','chr2','pos1','pos2','str1','str2'};
require_fields(R1,flds);
require_fields(R2,flds);
R1=make_numeric(R1,flds);
R2=make_numeric(R2,flds);

if ~isfield(R1,'individual')
  R1.individual = repmat({'this_individual'},slength(R1),1);
end
if ~isfield(R2,'individual')
  R2.individual = repmat({'this_individual'},slength(R2),1);
end

x1 = nan(slength(R1),1);
x2 = nan(slength(R2),1);

all_indivs = unique([R1.individual;R2.individual]);
R1.ino = listmap(R1.individual,all_indivs);
R2.ino = listmap(R2.individual,all_indivs);
for ii=1:length(all_indivs), fprintf('\n%s\n',all_indivs{ii});
 ii1 = find(R1.ino==ii);
 ii2 = find(R2.ino==ii);
 for c1=1:24, fprintf('chr%d ',c1);
  c1i1 = ii1(R1.chr1(ii1)==c1);
  c1i2 = ii2(R2.chr1(ii2)==c1);
  for c2=1:24
    c2i1 = c1i1(R1.chr2(c1i1)==c2);
    c2i2 = c1i2(R2.chr2(c1i2)==c2);
    for s1=0:1
      s1i1 = c2i1(R1.str1(c2i1)==s1);
      s1i2 = c2i2(R2.str1(c2i2)==s1);
      for s2=0:1
        s2i1 = s1i1(R1.str2(s1i1)==s2);
        s2i2 = s1i2(R2.str2(s1i2)==s2);
        for w=0:1
          win = w*round(tolerance/2);
          st1 = [round((R1.pos1(s2i1)+win)/tolerance) round((R1.pos2(s2i1)+win)/tolerance)];
          st2 = [round((R2.pos1(s2i2)+win)/tolerance) round((R2.pos2(s2i2)+win)/tolerance)];
          for ii=1:length(s2i1),i=s2i1(ii);
            if ~isnan(x1(i)), continue; end
            for jj=1:length(s2i2),j=s2i2(jj);
              if ~isnan(x2(j)), continue; end
              if st1(ii,1)==st2(jj,1) && st1(ii,2)==st2(jj,2)
                x1(i) = j;
                x2(j) = i;
end,end,end,end,end,end,end,end,end,fprintf('\n');

