function D=unite_Ds(Ds,rc,check_names)
%UNITE_DS unite the data structures contained in a cell array of data
%structures.
%
%   D = UNITE_DS(DS,RC,CHECK_NAMES) unites the data structures in cell
%   array DS.  RC toggles between column-wise unite (default or set
%   RC='col') in which data structures are concatinated along dim 2 (i.e.
%   markers are matched) or row-wise unite in which structures are
%   concatinated along dim 1 (set RC='row'). (Warning: if used in 'row'
%   mode, sample names are NOT MATCHED.  See PREPROC_MERGEPLATFORMS
%   instead.  If the Ds have different numbers of fields, only the common
%   fields are included in the new D.
%
%   TO DO: 
%   1) modify to take union, rather than intersection, of fields
%   2) modify to ask whether to check for inconsistencies in union of
%   fields, and if yes, to check for those inconsistencies
%   
%
%   Updated 17 Sept 07 to include .ref field
%   Updated 02 Oct 07 to make .sdesc field a row cell array rather than column
%   (did not fix case for if ischar(D.sdesc)). Jen Dobson (jdobson@broad)
%   Updated 18 Oct 07 with help documentation.
%   Updated 24 Oct 07 to remove fields that are not common to all Ds.
%   Updated 12 Dec 07 to write new D as datastruct object if any of input Ds are
%   datastructs.
%
%---
% $Id$
% $Date: 2007-12-28 17:25:27 -0500 (Fri, 28 Dec 2007) $
% $LastChangedBy: jdobson $
% $Rev$


%  If any one of the input Ds is a datastruct, the output D will be a
%  datastruct
Ds_classes = cellfun(@class,Ds,'UniformOutput',0);
isdatastruct = strmatch('datastruct',Ds_classes);

if ~isempty(isdatastruct)
    make_datastruct = 1;
    dskflds = {};
    for k = isdatastruct'
        dskflds = [dskflds diskfieldnames(Ds{k})];
        hdf5dir = regexp(get_datafile(Ds{k},'dat'),'/.+/','match');  %will use dir of last D as hdf5dir
    end
    dskflds = unique(dskflds);    
else
    make_datastruct = 0;
end
%%%


if ~exist('rc','var')
  rc='cols';
end

if ~exist('check_names','var')
  check_names=0;
end

st=1;
while isempty(Ds{st})
  st=st+1;
end

% only unite the fields that are common to all data structures
fieldstounite = fieldnames(Ds{st});
unitesis = 1;
if ~isfield(Ds{st},'sis')
    unitesis = 0;
else
    sistounite = fieldnames(Ds{st}.sis);
end

for i=(st+1):length(Ds)
    fieldstounite = intersect(fieldstounite,fieldnames(Ds{i}));
    if (isfield(Ds{i},'sis') && unitesis)
        sistounite = intersect(sistounite,fieldnames(Ds{i}.sis));
    else
        unitesis = 0;
    end
end


D=Ds{st};
Ds{st}=[];
  
%remove fields that are not common to all
removefields = setdiff(fieldnames(D),fieldstounite);
if ~isempty(removefields)
    for fl = setdiff(fieldnames(D),fieldstounite)
        D = rmfield(D,char(fl));
    end
end

if unitesis
    removesis = setdiff(fieldnames(D.sis),sistounite);
    if ~isempty(removesis)
        for fl = removesis
            D.sis = rmfield(D.sis,char(fl));
        end
    end
end

for i=(st+1):length(Ds)
  verbose(num2str(i));

  if isempty(Ds{i})
    continue;
  end
  
  switch rc(1:min(length(rc),3))
   case {'col','sam','con'}
    if check_names
      if ~same(D.marker,Ds{i}.marker)
        error('Markers do not match');
      end
    end
    D.dat=[ D.dat Ds{i}.dat ];
    if isfield(D,'affy_call')
      D.affy_call=[ D.affy_call Ds{i}.affy_call ];
    end
    if isfield(D,'affy_calls')
      D.affy_calls=[ D.affy_calls Ds{i}.affy_calls ];
    end
    if ischar(D.sdesc)
      D.sdesc=strvcat(D.sdesc,Ds{i}.sdesc);
    else
      D.sdesc=[as_column(D.sdesc); as_column(Ds{i}.sdesc)]';
    end
    if isfield(D,'sscale')
      if ischar(D.sscale)
        D.sscale=strvcat(D.sscale,Ds{i}.sscale);
      else
        D.sscale=[D.sscale Ds{i}.sscale];
      end
    end
    if isfield(D,'scans')
      D.scans=[D.scans; Ds{i}.scans];
    end
    if isfield(D,'supdat');
      D.supdat=[ D.supdat Ds{i}.supdat ];
    end
    if isfield(D,'sup');
      D.sup=[ D.sup; Ds{i}.sup ];
    end
    if isfield(D,'scale2vec');
      D.scale2vec=[ D.scale2vec Ds{i}.scale2vec ];
    end
    if isfield(D,'prectrls');
      D.prectrls=[ D.prectrls Ds{i}.prectrls ];
    end
    if isfield(D,'gcm_name');
      D.gcm_name=[ D.gcm_name Ds{i}.gcm_name ];
    end
    if isfield(D,'origidx');
      D.origidx=[ D.origidx Ds{i}.origidx ];
    end
    if isfield(D,'residx');
      D.residx=[ D.residx Ds{i}.residx ];
    end
    if isfield(D,'cbs');
      D.cbs=[ D.cbs Ds{i}.cbs ];
    end
    if isfield(D,'hmm');
      D.hmm=[ D.hmm Ds{i}.hmm ];
    end
    if isfield(D,'flag');
      D.flag=[ D.flag Ds{i}.flag ];
    end
    if isfield(D,'cbs_rl');
      D.cbs_rl=[ D.cbs_rl Ds{i}.cbs_rl ];
    end
    if isfield(D,'si');
      D.si=[ D.si Ds{i}.si ];
    end
    if isfield(D,'sitab');
      D.sitab=[ D.sitab Ds{i}.sitab ];
    end
    if unitesis;
        for fl = setdiff(fieldnames(Ds{i}.sis),sistounite)
            if isempty(fl)
                continue
            else
            Ds{i}.sis = rmfield(Ds{i}.sis,char(fl));
            end
        end
      D.sis=[ as_row(D.sis) as_row(Ds{i}.sis) ];
    end
    if isfield(D,'orig');
      D=rmfield(D,'orig');
      disp('Removing orig field');
    end
    if isfield(D,'raw');
      D.raw=[ D.raw Ds{i}.raw ];
    end
    if isfield(D,'smooth');
      D.smooth=[ D.smooth Ds{i}.smooth ];
    end
    if isfield(D,'peaks');
      D.peaks=[ D.peaks Ds{i}.peaks ];
    end
    if isfield(D,'joins');
      D.joins=[ D.joins Ds{i}.joins ];
    end
    if isfield(D,'level');
      D.level=[ D.level Ds{i}.level ];
    end
    if isfield(D,'final_dist');
      D.final_dist=[ D.final_dist Ds{i}.final_dist ];
    end
    if isfield(D,'medians');
      D.medians=[ D.medians Ds{i}.medians ];
    end    
    if isfield(D,'used_normals');
      D.used_normals=[ D.used_normals Ds{i}.used_normals ];
    end    
   
   case {'row','gen','mir'}
    if check_names
      if ~same(D.sdesc,Ds{i}.sdesc)
        error('Sample names do not match');
      end
    end
    D.dat=[ D.dat; Ds{i}.dat ];
    Ds{i}.dat=[]; % save memory
    if isfield(D,'affy_call')
      D.affy_call=[ D.affy_call; Ds{i}.affy_call ];
    end
    if isfield(D,'affy_calls')
      D.affy_calls=[ D.affy_calls; Ds{i}.affy_calls ];
    end
    if isfield(D,'score');
      D.score=[ D.score; Ds{i}.score ];
    end
    if isfield(D,'marker');
      D.marker=[ D.marker; Ds{i}.marker ];
    end
    if isfield(D,'cM');
      D.cM=[ D.cM; Ds{i}.cM ];
    end
    if isfield(D,'pos');
      D.pos=[ D.pos; Ds{i}.pos ];
    end
    if isfield(D,'gorigidx');
      D=rmfield(D,'gorigidx');
      disp('Removing gorigidx field');
    end
    if isfield(D,'orig');
      D=rmfield(D,'orig');
      disp('Removing orig field');
    end
    if isfield(D,'chr');
      D.chr=[ D.chr; Ds{i}.chr ];
    end
    if isfield(D,'chrn');
      D.chrn=[ D.chrn; Ds{i}.chrn ];
    end
    if isfield(D,'gsymb')
      D.gsymb=[ as_column(D.gsymb); as_column(Ds{i}.gsymb)];
    end
    if isfield(D,'gdesc')
      if iscell(D.gdesc)
        if size(D.gdesc,2)==1
          D.gdesc=[ as_column(D.gdesc); as_column(Ds{i}.gdesc) ];
          D.gacc=[ as_column(D.gacc); as_column(Ds{i}.gacc) ];
        else
        D.gdesc=[ as_row(D.gdesc) as_row(Ds{i}.gdesc) ];
        D.gacc=[ as_row(D.gacc) as_row(Ds{i}.gacc) ];
        end
      else
        D.gdesc=strvcat(D.gdesc,Ds{i}.gdesc);
        D.gacc=strvcat(D.gacc,Ds{i}.gacc);
      end      
    end
    if isfield(D,'ref')
      D.ref=[D.ref; Ds{i}.ref];
    end
   otherwise
    error('no match');
  end

  Ds{i}=[];
end


if make_datastruct && ~isa(D,'datastruct')
    D = datastruct(D);
    dflds = intersect(dskflds,fieldnames(D));
    dfile = 'D_merged';
    D{jj} = convert_to_diskfield(D{jj},dflds,strcat(dfile,'_',dflds,'.h5'),...
        strcat(dflds,'_dataset'));
    try
        lst = ls([hdf5dir 'D_unite*.h5']);
        lstfnames = regexp(lst,'/[^/]+\.','match');
        fnums = regexp(lstfnames,'\d','match');
        fnums = fnums(~cellfun(@isempty,fnums));
        fnums = cellfun(@str2num,[fnums{:}]);
        newnum = max(fnums);
        dfile = [hdf5dir 'D_unite' num2str(newnum)];
    catch
        dfile = [hdf5dir 'D_unite'];
    end
    
    
    D{jj} = convert_to_diskfield(D{jj},dflds,strcat(dfile,'_',dflds,'.h5'),...
        strcat(dflds,'_dataset'));
    
end

