function [QA QD num_bpts] = ziggurat(D,cyto)
% Ziggurat peel-off: takes D and cyto, and peels off segments in a
% ziggurat fashion to output QA (amps) and QD.  QA and QD fields are: 1)
% sample number, 2) SNP start, 3) SNP end, 4) amplitude (in copy-number)
%
%
% ---
% $Id$
% $Date: 2008-09-10 18:13:10 -0400 (Wed, 10 Sep 2008) $
% $LastChangedBy: cmermel $
% $Rev$

  %% get coordinates associated with chromosome start and end
  changeChrn = diff(D.chrn);
  chrnEnd = find(changeChrn == 1);
  chrnEnd(end+1) = size(D.chrn,1);
  chrnStart(1)=1;
  chrnStart(2:length(chrnEnd)) = chrnEnd(1:end-1)+1;

  %% Get coordinates associated with each chromosome arm
  %  load('/xchip/gistic/variables/cyto_rg_hg17_ucsc20070227.mat')
  armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9 X Y]+[p-q]+','match'),'UniformOutput',false));
%  armnames = unique(cellfun(@char,regexp({cyto.name},'^[0-9]+[p-q]+','match'),'UniformOutput',false));
  armnames = armnames(2:end);

  chrarms = {};
  for i=1:length(armnames)
    idx = strmatch(armnames(i),{cyto.name});
    band.name = char(armnames(i));
    band.start = cyto(idx(1)).start+1;
    band.end = cyto(idx(end)).end;
    band.chrn=cyto(idx(end)).chrn;
    band.length = band.end-band.start+1;
    chrarms{1}(i) = band;
    end
    
    chr=cat(1,chrarms{1}.chrn);
    [spos,sposi]=sort(chr);
    chrarms{1} = chrarms{1}(sposi);

    armlengths = [chrarms{1}.length];
    armstart = [chrarms{1}.start];

  %% find segments
  [breakpt_ids sample_ids] = find(diff(D.dat) ~= 0);
  num_bpts = zeros(1,size(D.dat,2));
  
  for j=1:length(unique(sample_ids))
    if mod(j,100)==0
      disp(j)
    end
    bpt = breakpt_ids(find(sample_ids == j));
    bpt = union(bpt,chrnEnd);
    num_bpts(j) = length(bpt);
    B = zeros(length(bpt),5);
    B(:,5) = repmat(j,length(bpt),1);
    B(1,2) = 1;
    B(2:end,2) = bpt(1:end-1)+1;
    B(:,3)=bpt;
    B(:,1)=D.chrn(B(:,2));
    %isq = D.pos(B(:,2))> armstart(2*B(:,1))'
    %B(:,5) = armlengths(2*B(:,1)-1+(D.pos(B(:,2)) > armstart(2*B(:,1))'));
    B(:,4)=D.dat(B(:,2),j);
    AmpB = B;
    DelB = B;
    AmpB(:,4) = B(:,4).*(B(:,4)>0);
    DelB(:,4) = 1*B(:,4).*(B(:,4)<0);
    AmpB(:,4)=(2.^(AmpB(:,4)+1))-2;
    Qamp{j} = zigg_deconstruction(AmpB);
    DelB(:,4)=-((2.^(DelB(:,4)+1))-2);
    Qdel{j} = zigg_deconstruction(DelB);
end

QA = cat(1,Qamp{:}); 
QD = cat(1,Qdel{:});
  
  
   
    
            
        
        
    
