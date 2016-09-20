function [Nn_new cat_new] = toggle_Nn_5col_2col(Nn_old)

X = generate_categ_context65_names();

oldsz = size(Nn_old);

if oldsz(1)==64 && oldsz(2)==5

  a=1:16; c=17:32; g=33:48; t=49:64;
  toA=1:48; toC=49:96; toG=97:144; toT=145:192;

  newsz = oldsz; newsz(1)=192; newsz(2)=2;
  Nn_new = zeros(newsz);

  Nn_new(:,1,:,:,:) = [Nn_old([c g t],1,:,:,:) ; Nn_old([a g t],1,:,:,:) ; Nn_old([a c t],1,:,:,:) ; Nn_old([a c g],1,:,:,:)];
  Nn_new(:,2,:,:,:) = [Nn_old([c g t],2,:,:,:) ; Nn_old([a g t],3,:,:,:) ; Nn_old([a c t],4,:,:,:) ; Nn_old([a c g],5,:,:,:)];

  cat_new=[];
  cat_new.context = X.name([c g t   a g t   a c t   a c g]);
  cat_new.newbase = [repmat({'->A'},48,1);repmat({'->C'},48,1);repmat({'->G'},48,1);repmat({'->T'},48,1)];
  cat_new.name = stringsplice([cat_new.context cat_new.newbase],' ');

elseif oldsz(1)==192 && oldsz(2)==2

  a=1:16; c=17:32; g=33:48; t=49:64;
  toA=1:48; toC=49:96; toG=97:144; toT=145:192;

  newsz = oldsz; newsz(1)=64; newsz(2)=5;
  Nn_new = zeros(newsz);

  Nn_new(:,1,:,:,:) = Nn_old([toT toA(33:end)],1,:,:,:);
  Nn_new([c g t],2,:,:,:) = Nn_old(toA,2,:,:,:);
  Nn_new([a g t],3,:,:,:) = Nn_old(toC,2,:,:,:);
  Nn_new([a c t],4,:,:,:) = Nn_old(toG,2,:,:,:);
  Nn_new([a c g],5,:,:,:) = Nn_old(toT,2,:,:,:);

  cat_new = reorder_struct(X,1:64);

elseif oldsz(1)==32 && oldsz(2)==5

  a=1:16; c=17:32;
  toA=1:16; toC=17:32; toG=33:64; toT=65:96;

  newsz = oldsz; newsz(1)=96; newsz(2)=2;
  Nn_new = zeros(newsz);

  Nn_new(:,1,:,:,:) = [Nn_old([c],1,:,:,:) ; Nn_old([a],1,:,:,:) ; Nn_old([a c],1,:,:,:) ; Nn_old([a c],1,:,:,:)];
  Nn_new(:,2,:,:,:) = [Nn_old([c],2,:,:,:) ; Nn_old([a],3,:,:,:) ; Nn_old([a c],4,:,:,:) ; Nn_old([a c],5,:,:,:)];

  cat_new=[];
  cat_new.context = X.name([c   a    a c    a c ]);
  cat_new.newbase = [repmat({'->A'},16,1);repmat({'->C'},16,1);repmat({'->G'},32,1);repmat({'->T'},32,1)];
  cat_new.name = stringsplice([cat_new.context cat_new.newbase],' ');

elseif oldsz(1)==96 && oldsz(2)==2

  a=1:16; c=17:32;
  toA=1:16; toC=17:32; toG=33:64; toT=65:96;

  newsz = oldsz; newsz(1)=32; newsz(2)=5;
  Nn_new = zeros(newsz);

  Nn_new(:,1,:,:,:) = Nn_old([toG],1,:,:,:);
  Nn_new([c],2,:,:,:) = Nn_old(toA,2,:,:,:);
  Nn_new([a],3,:,:,:) = Nn_old(toC,2,:,:,:);
  Nn_new([a c],4,:,:,:) = Nn_old(toG,2,:,:,:);
  Nn_new([a c],5,:,:,:) = Nn_old(toT,2,:,:,:);

  cat_new = reorder_struct(X,1:32);

else
  error('invalid dimensions of input Nn');
end

