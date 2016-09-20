function Creg = preproc_malechrn23(Craw,dataislog)
%PREPROC_MALECHRN23 take log of raw data and correct for male chromosome 23.
%
%CREG = PREPROC_MALECHRN23(CRAW,DATAISLOG) returns the data structure CRAW
%in CREG after correction for half-intensity on male chromosome 23. The
%difference between the median values for females and males is added (in
%linspace, quoitient is multiplied) to the intensity data for all males.
%DATAISLOG is 1 if the .dat in CRAW has been log transformed (default).
%DATAISLOG is 0 if the .dat in CRAW has not been log transformed.  CRAW can
%be a gistic data structure or a cell array of gistic data structures.  
%
%       History:
%           26 Sept 07:  Created by Jen Dobson (jdobson@broad.mit.edu)
%
%---
%  $Id$
%  $Date: 2007-11-19 17:28:20 -0500 (Mon, 19 Nov 2007) $
%  $LastChangedBy: jdobson $
%  $Rev$





verbose('Adding difference of ch.23 medians between males and females in log space (multiply in linear space) \nto the X chromosome data of males',30);
if ~exist('dataislog','var')
    dataislog = 1
end

if ~iscell(Craw)
    Ctemp = {Craw};
    Craw = Ctemp;
    clear Ctemp;
end


for k = 1:length(Craw)
    C = Craw{k};
    Xpos = find(C.chrn==23);  %increase intensity of male X chromosome
    gender = C.supdat(find_supid(C,'GENDER'),:);
    malesIDX = find(gender==1);
    femalesIDX = find(gender==2);

    if dataislog
        C.dat(Xpos,:)=C.dat(Xpos,:)+repmat(C.supdat(find_supid(C,'GENDER'),:)==1,length(Xpos),1); % Gender:
       
    else
        medquot = femalemed/malemed;
        C.dat(Xpos,malesIDX)=C.dat(Xpos,malesIDX).*repmat(medquot,length(Xpos),length(malesIDX));
        medchange = medquot;
    end

    C = add_history(C,mfilename);
    Craw{k} = C
end


Creg = Craw;