function CL21c=arm_center(CL21_T,cyto)

if exist('cyto','var')
  CL21_T=add_cyto(CL21_T,cyto);
  CL21_T.chrarmn=CL21_T.armn+CL21_T.chrn*2;
end

collapse_type='median';
%collapse_type='mode';
Y_1=collapse_D(CL21_T,'chrarmn',collapse_type); % 'mode'
if isfield(Y_1,'orig')
  Y_1=rmfield(Y_1,'orig');
end

Y2_1=Y_1;
Y2_1.dat=2.^(Y_1.dat+1);

CL21a=CL21_T;

CL21a.dat=2.^(CL21a.dat+1);

CL21b=CL21a;
for i=1:size(Y_1.dat,1)
  in_arm=find(CL21a.chrarmn==Y_1.chrarmn(i));
  delta=-(repmat(Y2_1.dat(i,:),length(in_arm),1)-2);
%  delta=delta.*(sign(delta)~=sign(CL21_T.dat(in_arm,:))); % make sure delta w/ opposite sign compared to log2ratio
%                                                        % --> brings closer to 2
  delta2=delta;
  delta2=min(delta,(CL21a.dat(in_arm,:)<=2).*(2-CL21a.dat(in_arm,:)));
  delta2=max(delta,(CL21a.dat(in_arm,:)>2).*(2-CL21a.dat(in_arm,:)));
  
  CL21b.dat(in_arm,:)=CL21a.dat(in_arm,:)+delta2;
  [i length(in_arm)]
end
CL21c=CL21b;
CL21c.dat=log2(CL21c.dat)-1;
