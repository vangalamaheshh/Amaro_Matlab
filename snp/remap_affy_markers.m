function [C markers] = remap_affy_markers(C,markers,annot_fname)
% REMAP_AFFY_MARKERS -- Map affy markers to new genomic locations
% according to the Affy CSV annotation file indicated by annot_fname, or
% by the structure, markers, that encapsulates the same info (and is
% output by this function).
%    
%
% ---
% $Id$
% $Date: 2008-07-21 13:53:38 -0400 (Mon, 21 Jul 2008) $
% $LastChangedBy: rameen $
% $Rev$

% Set defaults: delimiter used in affy annotation files and columns for
% marker ID, chromosome, and physical position. May need to add options
% here for non affy-CSV annotation files

if isempty(markers)&&~exist('annot_fname','var')
  error('Require input of either new marker mapping or CSV filename')
end

dlm = '","';
marker_id_col = 2; 
chr_col = 11;
pos_col = 14;
start_line = 3;

if ~exist('markers','var') || isempty(markers)

  annot = read_dlm_file(annot_fname,dlm);  % Read in annotation file

  markers.id=mat2cell(zeros(size(annot,2)-start_line+1,1),ones(size(annot,2)-start_line+1,1));
  markers.chr=markers.id;
  pos_cell = markers.id;
  markers.pos = zeros(size(annot,2)-start_line+1,1);


  % Extract markers, chr numbers, and positions
  for i = start_line:size(annot,2)
    if(mod(i,50000)==0)
      disp(i)
    end
    markers.id{i-start_line+1}=annot{i}{marker_id_col};
    markers.chr{i-start_line+1}=annot{i}{chr_col};
    pos_cell{i-start_line+1}=str2num(annot{i}{pos_col});
    if isempty(pos_cell{i-start_line+1})
      markers.chr{i-start_line+1}=[];
      markers.pos(i-start_line+1)=NaN;
    else
      markers.pos(i-start_line+1)=pos_cell{i-start_line+1};
    end
  end
end

% Match to markers in C
[M mi mj h us1j] = match_string_sets_hash(markers.id,C.marker);

clear M h us1j pos_cell

if length(mj) < size(C.marker,1)
  warning('Not all markers remapped, missing markers:')
  disp({C.marker{setdiff([1:size(C.marker,1)],mj)}})
end

if length(mi) < size(markers.id,1)
  warning('The following markers are missing:')
  disp({markers.id{setdiff([1:size(markers.id,1)],mi)}})
end

C.chr(mj)=markers.chr(mi);
C.pos(mj)=markers.pos(mi);

C.chr=cellfun(@char,C.chr,'UniformOutput',false);

if isfield(C,'chrn')
  C=rmfield(C,'chrn');
  C=add_chrn(C);
end
