function Q = collapse_context_and_sense_categories(C,K)
% C is struct {num,name} from c65e29 track
% K is categ struct from category discovery
%
% returns Q,  rows = raw categories from C
%             columns = collapsed categories from K
%             pages = 4 newbases
%
%             each cell tells the sense of the mutation:
%                 0 = incompatible combination
%                 1 = synonymous
%                 2 = missense
%                 3 = nonsense / nonstop / splice

k = assign_65x4_to_categ_set(K);
c = map_categories_to_65(C);
w = nansub(k,c);

C = parse_in(C,'name','^([ACGT]).*:(.*)$',{'from','effect'});
C = parse_in(C,'effect','^(syn|mis|non)/(syn|mis|non)/(syn|mis|non)$',{'o1','o2','o3'});
s = [listmap(C.o1,{'syn','mis','non'}) listmap(C.o2,{'syn','mis','non'}) listmap(C.o3,{'syn','mis','non'})];
s(strcmp(C.effect,'splice-site'),:) = 3;

Q = zeros(size(w));
bases = 'ACGT';
for from=1:4
  rows = find(strcmp(C.from,bases(from)));
  scol = 1;
  for newbase=1:4
    if newbase==from, continue; end
    Q(rows,:,newbase) = bsxfun(@times,w(rows,:,newbase),s(rows,scol));
    scol=scol+1;
  end
end


% spot-checking

%z=[Q(:,:,1),nan(1885,1),Q(:,:,2),nan(1885,1),Q(:,:,3),nan(1885,1),Q(:,:,4)];
%pr(C.name,z,1:100:1885)

%A in A_A:syn/syn/syn 0    0    0    0  NaN    0    0    0    1  NaN    0    0    1    0  NaN    0    0    0    1
%A in A_T:mis/mis/mis 0    0    0    0  NaN    0    0    0    2  NaN    0    0    2    0  NaN    0    0    0    2
%A in C_G:non/non/non 0    0    0    0  NaN    0    0    0    3  NaN    0    0    3    0  NaN    0    0    0    3
%A in G_G:mis/syn/mis 0    0    0    0  NaN    0    0    0    2  NaN    0    0    1    0  NaN    0    0    0    2
%A in T_C:non/mis/non 0    0    0    0  NaN    0    0    0    3  NaN    0    0    2    0  NaN    0    0    0    3
%C in A_C:syn/non/mis 0    0    0    1  NaN    0    0    0    0  NaN    0    0    0    3  NaN    0    2    0    0
%C in C_A:non/syn/non 0    0    0    3  NaN    0    0    0    0  NaN    0    0    0    1  NaN    3    0    0    0
%C in G_A:syn/mis/mis 0    0    0    1  NaN    0    0    0    0  NaN    0    0    0    2  NaN    0    2    0    0
%C in G_T:mis/non/non 0    0    0    2  NaN    0    0    0    0  NaN    0    0    0    3  NaN    0    3    0    0
%C in T_T:syn/syn/mis 0    0    0    1  NaN    0    0    0    0  NaN    0    0    0    1  NaN    2    0    0    0
%G in A_G:mis/mis/non 2    0    0    0  NaN    0    0    0    2  NaN    0    0    0    0  NaN    0    0    0    3
%G in C_C:splice-site 0    3    0    0  NaN    0    0    0    3  NaN    0    0    0    0  NaN    0    0    0    3
%G in G_C:mis/syn/non 0    2    0    0  NaN    0    0    0    1  NaN    0    0    0    0  NaN    0    0    0    3
%G in T_A:non/non/syn 3    0    0    0  NaN    0    0    0    3  NaN    0    0    0    0  NaN    0    0    0    1
%T in A_A:syn/non/non 0    0    0    1  NaN    0    0    3    0  NaN    0    0    0    3  NaN    0    0    0    0
%T in A_T:non/mis/syn 0    0    0    3  NaN    0    0    2    0  NaN    0    0    0    1  NaN    0    0    0    0
%T in C_T:syn/mis/non 0    0    0    1  NaN    0    0    2    0  NaN    0    0    0    3  NaN    0    0    0    0
%T in G_G:non/syn/syn 0    0    0    3  NaN    0    0    1    0  NaN    0    0    0    1  NaN    0    0    0    0
%T in T_G:syn/syn/non 0    0    0    1  NaN    0    0    1    0  NaN    0    0    0    3  NaN    0    0    0    0
