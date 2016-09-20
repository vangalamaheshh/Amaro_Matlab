function twobytwo = check_broad_segs(D,Q,armlengths,armstart,fract_thresh,amp_thresh,isamp)
  
%% Function takes in a segmented data structure and its associated
%ziggurat matrix Q, two arrays (one containing
%the armlengths of each chr_arm [sorted by chr arm] and the other
%containing the start positions of each chr_arm), fraction and
%amplitude thresholds, and a amp_del toggle to look at amplifications
%or deletions.  It returns a 2x2 matrix: the number of
%chromosome arms with and without amps (or dels) greater than fract_thresh of the
%arm length amplified (or deleted) to greater than amp_thresh, vs the
%number of chromosome arms with (or without) segments longer than
%fract_thresh and higher than amp_thresh.
%
% ---
% $Id$
% $Date: 2007-11-23 14:09:40 -0500 (Fri, 23 Nov 2007) $
% $LastChangedBy: rameen $
% $Rev$

    
  %%% Convert D into 1-0 matrix for values exceeding amp_thresh
  D.dat=2.^(D.dat+1);
  if(isamp==0)
    D.dat=2-D.dat;
  else
    D.dat=D.dat-2;
  end
  D.dat(D.dat<amp_thresh)=0;
  D.dat(D.dat~=0)=1;
  
  %%% Find start, end SNPs for each chromosome arm
  snparmstart(1)=1;
  for i = 2:length(armlengths)
    chri=find(D.chrn==round(i/2));
    afterstart=find(D.pos>armstart(i));
    snparmstart(i)=min(intersect(chri,afterstart));
    snparmend(i-1)=snparmstart(i)-1;
    snparmlength(i-1)=snparmend(i-1)-snparmstart(i-1)+1;
  end
  snparmend(length(armlengths))=size(D.dat,1);
  snparmlength(length(armlengths))=snparmend(length(armlengths))-snparmstart(length(armlengths))+1;
 
  probedarms=find(snparmlength>0);  % identify arms with markers
  
 %%% Assay the fraction of each probed arm that surpasses amp_thresh 
  samp_chr=zeros(2,length(probedarms),size(D.dat,2)); % repository for
                                                      % fraction information
 
  % Calculate lengths from Q
  max_length=zeros(length(armlengths)/2,size(D.dat,2),2);
  for i = 1:size(Q,1)
    if(i==1)
      j=1;
    elseif(Q(i,1)==1&&Q(i-1,1)~=1)
      j=j+1;
    end
    if(Q(i,4)>amp_thresh)
     if(snparmstart(Q(i,1)*2-1)<=Q(i,2)&&snparmend(Q(i,1)*2-1)>=Q(i,2))
       temp_max_length=min(snparmend(Q(i,1)*2-1),Q(i,3))-Q(i,2)+1;
       if temp_max_length>max_length(Q(i,1),j,1)
         max_length(Q(i,1),j,1)=temp_max_length;
       end
     end
     if(snparmend(Q(i,1)*2)>=(Q(i,3)&&snparmstart(Q(i,1)*2)<=Q(i,3)))
       temp_max_length=Q(i,3)-max(snparmstart(Q(i,1)*2),Q(i,2))+1;
       if temp_max_length>max_length(Q(i,1),j,2)
         max_length(Q(i,1),j,2)=temp_max_length;
       end
     end
   end
  end
  
  % Calculate fractions from D and Q
  for i = 1:length(probedarms)   
    samp_chr(1,i,:)=sum(D.dat(snparmstart(probedarms(i)):snparmend(probedarms(i)),:))/ ...
        (snparmlength(probedarms(i)));
    if( probedarms(i)/2~=round(probedarms(i)/2))
      samp_chr(2,probedarms(i),:)=max_length(round(probedarms(i)/2),:,1)/ ...
          snparmlength(probedarms(i));
    else
      samp_chr(2,i,:)=max_length(probedarms(i)/2,:,2)/ ...
          snparmlength(probedarms(i));
    end
  end

  twobytwo=zeros(2,2);
  twobytwo(1,1)=length(intersect(find(samp_chr(1,:,:)>fract_thresh), ...
                                 find(samp_chr(2,:,:)>fract_thresh)));
  twobytwo(1,2)=length(setdiff(find(samp_chr(1,:,:)>fract_thresh), ...
                                 find(samp_chr(2,:,:)>fract_thresh))); 
  twobytwo(2,1)=length(setdiff(find(samp_chr(2,:,:)>fract_thresh), ...
                                 find(samp_chr(1,:,:)>fract_thresh)));
  twobytwo(2,2)=length(intersect(find(samp_chr(1,:,:)<fract_thresh), ...
                                 find(samp_chr(2,:,:)<fract_thresh)));
