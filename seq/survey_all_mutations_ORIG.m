function S = survey_all_mutations(dna)

%% examines all possible point mutations _in silico_
%% and counts how many would be synonymous, missense, and nonsense
%% 
%% assumes "dna" is already the transcribed strand, in +0 frame.
%%
%% for now, mutations that _remove_ a stop codon are called "missense"...
%% may want to change this in a future version.
%%
%% returns a matrix S where columns=positions and rows:(1=silent,2=missense,3=nonsense),
%% and sum of each column=3
%%
%% Mike Lawrence 2008/04/29

C = [656565;656567;656571;656584;656765;656767;656771;656784;657165;657167;657171;657184;658465;658467;658471;658484; ...
     676565;676567;676571;676584;676765;676767;676771;676784;677165;677167;677171;677184;678465;678467;678471;678484; ...
     716565;716567;716571;716584;716765;716767;716771;716784;717165;717167;717171;717184;718465;718467;718471;718484; ...
     846565;846567;846571;846584;846765;846767;846771;846784;847165;847167;847171;847184;848465;848467;848471;848484 ];

M = cat(3,[0,0,1;2,3,2;1,0,0],[0,0,1;3,3,2;0,0,0],[0,0,1;2,3,2;1,0,0],[0,0,1;3,...
3,2;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,...
3;3,3,0;0,0,0],[1,0,1;1,3,2;1,0,0],[0,0,1;3,3,2;0,0,0],[1,0,1;2,3,2;0,0,0],[0,...
0,1;3,3,2;0,0,0],[0,0,2;3,3,1;0,0,0],[0,0,2;3,3,1;0,0,0],[0,0,0;3,3,3;0,0,0],...
[0,0,2;3,3,1;0,0,0],[0,0,1;2,3,2;1,0,0],[0,0,1;3,3,2;0,0,0],[0,0,1;2,3,2;1,0,...
0],[0,0,1;3,3,2;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,...
0,0],[0,0,3;3,3,0;0,0,0],[1,0,3;1,3,0;1,0,0],[0,0,3;3,3,0;0,0,0],[1,0,3;2,3,...
0;0,0,0],[0,0,3;3,3,0;0,0,0],[1,0,3;2,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[1,0,3;2,...
3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,1;2,3,2;1,0,0],[0,0,1;3,3,2;0,0,0],[0,0,...
1;2,3,2;1,0,0],[0,0,1;3,3,2;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,...
0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;2,3,0;1,0,0],[0,0,3;3,3,0;0,0,0],...
[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,...
0],[0,0,3;3,3,0;0,0,0],[0,0,3;3,3,0;0,0,0],[0,1,1;3,2,2;0,0,0],[0,0,1;3,3,0;0,...
0,2],[0,0,1;3,3,2;0,0,0],[0,0,1;3,3,0;0,0,2],[0,0,3;3,1,0;0,2,0],[0,0,3;3,3,...
0;0,0,0],[0,0,3;3,2,0;0,1,0],[0,0,3;3,3,0;0,0,0],[0,1,0;3,2,3;0,0,0],[0,0,1;3,...
3,1;0,0,1],[0,0,0;3,2,2;0,1,1],[0,0,1;3,3,1;0,0,1],[1,0,1;2,1,2;0,2,0],[0,0,...
1;3,3,2;0,0,0],[1,0,1;2,2,2;0,1,0],[0,0,1;3,3,2;0,0,0]);

d = double(upper(dna));
d = d(1:3*floor(length(d)/3));
d = reshape(d,3,length(d)/3);
d = 10000*d(1,:) + 100*d(2,:) + d(3,:);

S = zeros(3,3*length(d));
s=1;
for i=1:length(d)
   pos = find(C==d(i));
   if isempty(pos), error('invalid DNA'); end
   S(:,s:s+2)=M(:,:,pos);
   s=s+3;
end

return



%%%%%%%%%%%%%%%%%%%%%%% this is the code that generated matrix M

SIL=1;
MIS=2;
NON=3;

S='M = cat(3,';
bases='ACGT';
for c=1:length(C)
  n=C(c);
  d='...';
  d(1)=char(floor(n/10000));
  n=mod(n,10000);
  d(2)=char(floor(n/100));
  n=mod(n,100);
  d(3)=char(n);
  oldaa=my_nt2aa(d);
  M=zeros(3,3);
  for i=1:3
    oldbase=d(i);
    for newbase=bases
      if newbase==oldbase, continue; end
      d(i)=newbase;
      newaa=my_nt2aa(d);
      if newaa==oldaa
         M(SIL,i)=M(SIL,i)+1;
      elseif newaa=='*'
         M(NON,i)=M(NON,i)+1;
      else
         M(MIS,i)=M(MIS,i)+1;
      end
    end
    d(i)=oldbase;
  end
  S=[S '['];
  for r=1:3
    for c=1:3
      S = [S num2str(M(r,c))];  
      if c<3, S = [S ',']; end
    end
    if r<3, S = [S ';']; end
  end
  S = [S ']'];
  if c<length(C), S = [S ',']; end
end
S = [S ');'];

ll=0;
for i=1:length(S)
  fprintf('%c', S(i));
  ll=ll+1;
  if ll>75 && S(i)==',', fprintf('...\n'); ll=0; end
end



