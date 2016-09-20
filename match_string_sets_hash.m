function [M,mi,mj,h,us1j]=match_string_sets_hash(set1,set2,h,us1j,no_sparse,use_verbose)

% fill hash with longer list
% go over short list and check if it is there

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


sets_unique=0;
if exist('h','var') && ~isempty(h)
  if iscell(set2)
    set2=fill_empty(set2,'EMPTY');
    set2=strvcat(set2);
  end
  [uset2,us2i,us2j]=unique(set2,'rows');
  sw=0;
else
  if iscell(set1)
    set1=fill_empty(set1,'EMPTY');
    set1=strvcat(set1);
  end
  
  if iscell(set2)
    set2=fill_empty(set2,'EMPTY');
    set2=strvcat(set2);
  end
  
  [uset1,us1i,us1j]=unique(set1,'rows');
  [uset2,us2i,us2j]=unique(set2,'rows');
  if (length(us1i)==length(us1j)) &&(length(us2i)==length(us2j))
    sets_unique=1;
  end
  
  if size(uset1,1)<size(uset2,1)
    [uset1,uset2]=deal(uset2,uset1);
    [us1j,us2j]=deal(us2j,us1j);
    sw=1;
  else
    sw=0;
  end

  h=java.util.Hashtable;
  if exist('use_verbose','var') && use_verbose
    for i=1:size(uset1,1)
      put(h,deblank(uset1(i,:)),i);
      if mod(i,use_verbose)==0
        verbose(num2str(i));
      end
    end
  else
    for i=1:size(uset1,1)
      put(h,deblank(uset1(i,:)),i);
    end    
  end
end

mi=zeros(size(uset2,1),1);
mj=1:size(uset2,1);

if exist('use_verbose','var') && use_verbose
  for i=1:size(uset2,1)
    t=get(h,deblank(uset2(i,:)));
    if ~isempty(t)
      mi(i)=t;
    end
    if mod(i,use_verbose)==0
      verbose(num2str(i));
    end
  end
else
  for i=1:size(uset2,1)
    t=get(h,deblank(uset2(i,:)));
    if ~isempty(t)
      mi(i)=t;
    end
  end  
end


nz=find(mi~=0);
mi=mi(nz);
mj=mj(nz);

if exist('no_sparse','var') && no_sparse
  M=[];
else
  M=sparse(mi,mj,ones(length(mi),1),length(us1j),size(uset2,1));
end

if ~sets_unique
  Mf=M(us1j,us2j);
  [mif,mjf,dum]=find(Mf);
  M=Mf;
  mi=mif;
  mj=mjf;
else
  [tmp,revus1j]=sort(us1j);
  [tmp,revus2j]=sort(us2j);
  mi=revus1j(mi);
  mj=revus2j(mj);
  [smj,smji]=sort(mj);
  mj=smj;
  mi=mi(smji);
end

if sw
  M=M';
  [mi,mj]=deal(mj,mi);
  [smj,smji]=sort(mj);
  mj=smj;
  mi=mi(smji);
end

%Mx=sparse(size(set1,1),size(set2,1)); % Is this OK?
%whos Mx
%Mx(1:size(M,1),1:size(M,2))=M; 
%M=Mx;
