function X = dRanger_find_densities(X,T,larger_window_size)

if ~exist('larger_window_size','var'), larger_window_size = 3000; end

% make indices
fprintf('Indexing...\n');
I2=cell(24,1);I6=cell(24,1);
for chr=1:24
  I2{chr} = find(T(:,2)==chr);
  I6{chr} = find(T(:,6)==chr);
end

nx = slength(X);
X.Lsupp1tum = zeros(nx,1);
X.tot1tum = zeros(nx,1);
X.tot2tum = zeros(nx,1);
X.S1 = zeros(nx,1);
X.S2 = zeros(nx,1);

% windowed positions (for S filter)
T8wind = round(T(:,8)/1000);
T4wind = round(T(:,4)/1000);

for i=1:nx
  if ~mod(i,100), fprintf('%d/%d ',i,nx); end

  % boundaries of the locality
  chr1 = X.chr1(i);
  str1 = X.str1(i);
  left1 = X.min1(i);
  right1 = X.max1(i);
  chr2 = X.chr2(i);
  str2 = X.str2(i);
  left2 = X.min2(i);
  right2 = X.max2(i);

  % LARGER-REGION DENSITIES
if 0
  mid1 = round((right1+left1)/2);
  Lleft1 = mid1 - round(larger_window_size/2);
  Lright1 = mid1 + round(larger_window_size/2);

  mid2 = round((right2+left2)/2);
  Lleft2 = mid2 - round(larger_window_size/2);
  Lright2 = mid2 + round(larger_window_size/2);
end

  % end1 in tumor

if 0
  idx1 = I2{chr1}(find(T(I2{chr1},4)<=Lright1 & T(I2{chr1},5)>=Lleft1));
  idx2 = I6{chr1}(find(T(I6{chr1},8)<=Lright1 & T(I6{chr1},9)>=Lleft1));
end

%  supp1 = find(T(idx1,6)==chr2 & T(idx1,8)<=right2 & T(idx1,9)>=left2 & T(idx1,3)==str1 & T(idx1,7)==str2);
%  supp2 = find(T(idx2,2)==chr2 & T(idx2,4)<=right2 & T(idx2,5)>=left2 & T(idx2,7)==str1 & T(idx2,3)==str2);

if 0
  X.Lsupp1tum(i)=sum(T(idx1,6)==chr2 & T(idx1,8)<=Lright2 & T(idx1,9)>=Lleft2 & T(idx1,3)==str1 & T(idx1,7)==str2) +...
                 sum(T(idx2,2)==chr2 & T(idx2,4)<=Lright2 & T(idx2,5)>=Lleft2 & T(idx2,7)==str1 & T(idx2,3)==str2);
end

  % end2 in tumor
  %idx1 = find(T(:,2)==chr2 & T(:,4)<=right2 & T(:,5)>=left2);
  %idx2 = find(T(:,6)==chr2 & T(:,8)<=right2 & T(:,9)>=left2);
  %supp1 = find(T(idx1,6)==chr1 & T(idx1,8)<=right1 & T(idx1,9)>=left1 & T(idx1,3)==str2 & T(idx1,7)==str1);
  %supp2 = find(T(idx2,2)==chr1 & T(idx2,4)<=right1 & T(idx2,5)>=left1 & T(idx2,7)==str2 & T(idx2,3)==str1);
  %X.Lsupp2tum(i,1)=length(supp1)+length(supp2); X.Ltot2tum(i,1)=length(idx1)+length(idx2);

% X.Lsupp2tum(i,1)=sum(T(idx1,6)==chr1 & T(idx1,8)<=Lright1 & T(idx1,9)>=Lleft1 & T(idx1,3)==str2 & T(idx1,7)==str1) +...
%                  sum(T(idx2,2)==chr1 & T(idx2,4)<=Lright1 & T(idx2,5)>=Lleft1 & T(idx2,7)==str2 & T(idx2,3)==str1);

  % LOCAL DENSITIES

  % end1 in tumor
%  idx1 = find(T(:,2)==chr1 & T(:,4)<=right1 & T(:,5)>=left1);
%  idx2 = find(T(:,6)==chr1 & T(:,8)<=right1 & T(:,9)>=left1);
%  supp1 = find(T(idx1,6)==chr2 & T(idx1,8)<=right2 & T(idx1,9)>=left2 & T(idx1,3)==str1 & T(idx1,7)==str2);
%  supp2 = find(T(idx2,2)==chr2 & T(idx2,4)<=right2 & T(idx2,5)>=left2 & T(idx2,7)==str1 & T(idx2,3)==str2);
%  X.supp1tum(i,1)=length(supp1)+length(supp2); X.tot1tum(i,1)=length(idx1)+length(idx2);

if 0
   X.tot1tum(i)=sum(T(idx1,4)<=right1 & T(idx1,5)>=left1) +... 
                sum(T(idx2,8)<=right1 & T(idx2,9)>=left1);
end

  % end2 in tumor
%  idx1 = find(T(:,2)==chr2 & T(:,4)<=right2 & T(:,5)>=left2);
%  idx2 = find(T(:,6)==chr2 & T(:,8)<=right2 & T(:,9)>=left2);
%  supp1 = find(T(idx1,6)==chr1 & T(idx1,8)<=right1 & T(idx1,9)>=left1 & T(idx1,3)==str2 & T(idx1,7)==str1);
%  supp2 = find(T(idx2,2)==chr1 & T(idx2,4)<=right1 & T(idx2,5)>=left1 & T(idx2,7)==str2 & T(idx2,3)==str1);
%  X.supp2tum(i,1)=length(supp1)+length(supp2); X.tot2tum(i,1)=length(idx1)+length(idx2);

%  supp1 = find(T(idx1,6)==chr2 & T(idx1,8)<=right2 & T(idx1,9)>=left2 & T(idx1,3)==str1 & T(idx1,7)==str2);
%  supp2 = find(T(idx2,2)==chr2 & T(idx2,4)<=right2 & T(idx2,5)>=left2 & T(idx2,7)==str1 & T(idx2,3)==str2);

if 0
   X.tot2tum(i)=sum(T(I2{chr2},4)<=right2 & T(I2{chr2},5)>=left2) +...
                sum(T(I6{chr2},8)<=right2 & T(I6{chr2},9)>=left2);
end % if 0

   % added 2009-06-23:
   %
   % "S" (Stripe) filter

   % end1 in tumor

   idx1 = I2{chr1}(T(I2{chr1},4)<=right1 & T(I2{chr1},5)>=left1);
   idx2 = I6{chr1}(T(I6{chr1},8)<=right1 & T(I6{chr1},9)>=left1);

%   nonsupp1 = idx1(T(idx1,6)~=chr2 | T(idx1,8)>right2 | T(idx1,9)<left2 | T(idx1,3)~=str1 | T(idx1,7)~=str2);
%   nonsupp2 = idx2(T(idx2,2)~=chr2 | T(idx2,4)>right2 | T(idx2,5)<left2 | T(idx2,7)~=str1 | T(idx2,3)~=str2);

%   offtarget = [T(nonsupp1,6) T(nonsupp1,8) T(nonsupp1,3) T(nonsupp1,7);...
%                T(nonsupp2,2) T(nonsupp2,4) T(nonsupp2,7) T(nonsupp2,3)];

   nonsupp1 = idx1(T(idx1,6)~=chr2 | T(idx1,8)>right2 | T(idx1,9)<left2);
   nonsupp2 = idx2(T(idx2,2)~=chr2 | T(idx2,4)>right2 | T(idx2,5)<left2);

   offtarget = [T(nonsupp1,6) T8wind(nonsupp1);...
                T(nonsupp2,2) T4wind(nonsupp2)];

   X.S1(i) = size(unique(offtarget,'rows'),1);

   % end2 in tumor

   idx1 = I2{chr2}(T(I2{chr2},4)<=right2 & T(I2{chr2},5)>=left2);
   idx2 = I6{chr2}(T(I6{chr2},8)<=right2 & T(I6{chr2},9)>=left2);

%   nonsupp1 = idx1(T(idx1,6)~=chr1 | T(idx1,8)>right1 | T(idx1,9)<left1 | T(idx1,3)~=str2 | T(idx1,7)~=str1);
%   nonsupp2 = idx2(T(idx2,2)~=chr1 | T(idx2,4)>right1 | T(idx2,5)<left1 | T(idx2,7)~=str2 | T(idx2,3)~=str1);

%   offtarget = [T(nonsupp1,6) T8wind(nonsupp1) T(nonsupp1,3) T(nonsupp1,7);...
%                T(nonsupp2,2) T4wind(nonsupp2) T(nonsupp2,7) T(nonsupp2,3)];

   nonsupp1 = idx1(T(idx1,6)~=chr1 | T(idx1,8)>right1 | T(idx1,9)<left1);
   nonsupp2 = idx2(T(idx2,2)~=chr1 | T(idx2,4)>right1 | T(idx2,5)<left1);

   offtarget = [T(nonsupp1,6) T8wind(nonsupp1);...
                T(nonsupp2,2) T4wind(nonsupp2)];

   X.S2(i) = size(unique(offtarget,'rows'),1);
end


