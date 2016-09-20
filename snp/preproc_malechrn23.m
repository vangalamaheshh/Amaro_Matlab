function C = preproc_malechrn23(C,dataislog)
%PREPROC_MALECHRN23 take log of raw data and correct for male chromosome 23.
%
%C = PREPROC_MALECHRN23(C,DATAISLOG) returns the data structure C
% after correction for half-intensity on male chromosome 23. The
%difference between the median values for females and males is added to the
%intensity data for all males. DATAISLOG is 1 if the .dat in C has been log
%transformed (default). DATAISLOG is 0 if the .dat in C has not been log
%transformed.  C can be a gistic data structure or a cell array of gistic
%data structures.
%
%       History:
%           26 Sept 07:  Created by Jen Dobson (jdobson@broad.mit.edu)
%
%---
%  $Id$
%  $Date: 2007-12-11 17:03:53 -0500 (Tue, 11 Dec 2007) $
%  $LastChangedBy: jdobson $
%  $Rev$

minfemales = 5;   %Function will double the median intensity of males if fewer than this number of females



verbose('Adding difference of ch.23 medians between males and females in log space\nto the X chromosome data of males',30);

if ~exist('dataislog','var')
    dataislog = 1
end

if iscell(C)

    for k = 1:length(C)

        Xpos = find(C{k}.chrn==23);  %increase intensity of male X chromosome
        gender = C{k}.supdat(find_supid(C{k},'GENDER'),:);
        malesIDX = find(gender==1);
        femalesIDX = find(gender==2);
        malemed = median(median(C{k}.dat(Xpos,malesIDX)));

        if length(femalesIDX)>=minfemales
            femalemed = median(median(C{k}.dat(Xpos,femalesIDX)));
            if dataislog
                meddiff = femalemed - malemed;
                C{k}.dat(Xpos,malesIDX)=C{k}.dat(Xpos,malesIDX)+repmat(meddiff,length(Xpos),length(malesIDX));
                medchange = meddiff;
            else
                medquot = femalemed/malemed;
                C{k}.dat(Xpos,malesIDX)=C{k}.dat(Xpos,malesIDX).*repmat(medquot,length(Xpos),length(malesIDX));
                medchange = medquot;
            end
        else
            if dataislog
                meddiff = malemed;
                C{k}.dat(Xpos,malesIDX)=C{k}.dat(Xpos,malesIDX)+repmat(meddiff,length(Xpos),length(malesIDX));
                medchange = meddiff;
            else
                medquot = 2;
                C{k}.dat(Xpos,malesIDX)=C{k}.dat(Xpos,malesIDX).*repmat(medquot,length(Xpos),length(malesIDX));
                medchange = medquot;
            end

            C{k} = add_history(C{k},mfilename,dataislog,medchange);

        end

    end

else
    Xpos = find(C.chrn==23);  %increase intensity of male X chromosome
    gender = C.supdat(find_supid(C,'GENDER'),:);
    malesIDX = find(gender==1);
    femalesIDX = find(gender==2);

    malemed = median(median(C.dat(Xpos,malesIDX)));

    if length(femalesIDX)>minfemales
        femalemed = median(median(C.dat(Xpos,femalesIDX)));

        if dataislog
            meddiff = femalemed - malemed;
            C.dat(Xpos,malesIDX)=C.dat(Xpos,malesIDX)+repmat(meddiff,length(Xpos),length(malesIDX));
            medchange = meddiff;
        else
            medquot = femalemed/malemed;
            C.dat(Xpos,malesIDX)=C.dat(Xpos,malesIDX).*repmat(medquot,length(Xpos),length(malesIDX));
            medchange = medquot;
        end
    else
        if dataislog
            meddiff = malemed;
            C.dat(Xpos,malesIDX)=C.dat(Xpos,malesIDX)+repmat(meddiff,length(Xpos),length(malesIDX));
            medchange = meddiff;
        else
            medquot = 2;
            C.dat(Xpos,malesIDX)=C.dat(Xpos,malesIDX).*repmat(medquot,length(Xpos),length(malesIDX));
            medchange = medquot;
        end


        C = add_history(C,mfilename,dataislog,medchange);
    end
    
end
