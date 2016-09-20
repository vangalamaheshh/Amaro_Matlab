function SEG = read_glad_file( fname, min_seg_length )
%  read_glad_file.m
%  INPUT:   DNAcopy output file (ID, chr, start, end, size, mean)
%  OUTPUT:  SEGMENTS structure (same columns as above, except ID)
%
%  Derek Chiang
%  dchiang@broad.mit.edu
%


% Load .dnacopy file
fid = fopen(char(fname));
C = textscan(fid,'%s%f64%f64%f64%f64%f64','headerLines',1);
fclose(fid);

% Assign contents to variables 
chr = C{2};
left = C{3};
right = C{4};
lengths = C{5};
segmean = C{6};

SEG.chr=[];
SEG.left=[];
SEG.right=[];
SEG.size=[];
SEG.ratios=[];


% Loop across chromosomes
for c=1:23
    numShort = 0;
    idx=find(chr==c);
    chrData = [ left(idx) right(idx) lengths(idx) segmean(idx) ];
    short = find( chrData(:,3) < min_seg_length );
    
    % Merge short segments in left to right order along chromosome
    while ~isempty(short)
        i = 1;
        N = size(chrData,1);
        if ( short(i) == 1 )
            disp(['Merge FIRST seg in chr ' num2str(c) ': ' num2str(chrData(short(i),1)) '-' num2str(chrData(short(i),2)) ...
                  ' to ' num2str(chrData(short(i),1)) '-' num2str(chrData(short(i)+1,2))]);

            mergeSize = chrData(1,3) + chrData(2,3);
            weightedMean = chrData(1,3) * chrData(1,4) + ...
                           chrData(2,3) * chrData(2,4);
            weightedMean = weightedMean / mergeSize;
            mergeData = [ chrData(1,1) chrData(2,2) mergeSize weightedMean ];
            newData = [ mergeData; chrData(3:end,:) ];
            chrData = newData;
            numShort = numShort + 1;
        elseif ( short(i) > 1 && short(i) < N  )

            offset=1;
            % Find multiple short segments that are mutually adjacent
            if length(short) > 1 
                while( offset+1 <= length(short) && ...
                       short(i+offset) == short(i+offset-1) + 1 )
                    offset = offset+1;
                end
            end
            
            numShort = numShort + 1;
%            disp(['Merge chr ' num2str(c) ': ' num2str(chrData(short(i),1)) '-' num2str(chrData(short(i),2)) ...
%                  ' to ' num2str(chrData(short(i)-1,1)) '-' num2str(chrData(short(i)+offset,2))]);
            
            % Previous segment & current segment
            mergeSize = chrData(short(i)-1,3) + chrData(short(i),3);
            weightedMean = chrData(short(i)-1,3) * chrData(short(i)-1,4) + ...
                           chrData(short(i),3)   * chrData(short(i),4);
            for k=1:offset
                mergeSize = mergeSize + chrData(short(i)+k,3);
                weightedMean = weightedMean + chrData(short(i)+k,3) * chrData(short(i)+k,4);
            end
            weightedMean = weightedMean / mergeSize;
            mergeData = [ chrData(short(i)-1,1) chrData(short(i)+offset,2) mergeSize weightedMean ];

            % Edge effect on left
            if ( short(i) > 2 )
                newData = [ chrData(1:(short(i)-2),:); mergeData];
            else
                newData = [ mergeData ];   % Merged segment is first segment
            end
            
            % Edge effect on right
            if ( short(i) < N-2 )
                newData = [ newData; chrData((short(i)+offset+1):N,:)];
            end

            % Recalculate current chromosome data & short segments
            chrData = newData;
        end
        
        short = find( chrData(:,3) < min_seg_length );
        if length(short) == 1
            if ( short(1) == N )
                disp(['Reached end of chr ' num2str(c)]);
                short = [];
            end
        end
    end
    
%    disp(['chr ' num2str(c) ':  Merged ' num2str(numShort) ' segments']);
    
    SEG.chr   = [ SEG.chr; repmat(c,size(chrData,1),1) ];
    SEG.left = [ SEG.left; chrData(:,1) ];
    SEG.right   = [ SEG.right; chrData(:,2) ];
    SEG.size  = [ SEG.size; chrData(:,3) ];
    SEG.ratios  = [ SEG.ratios; 2.^chrData(:,4) ];
    
end
