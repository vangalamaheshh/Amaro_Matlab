function change = get_aa_change(build,chr,pos,old_nt,new_nt,protein_sequence)

old_nt = upper(old_nt);
new_nt = upper(new_nt);

silent = true;

change = 'no_match';

stretch=15;
while(stretch<300)
  if ~silent, fprintf('Trying stretch=%d\n', stretch); end
  match.count = 0;
  match.protpos = [];
  for direction=1:2
    if direction==1
      % try matching rightward
      target_end = -1;
      dna = upper(genome_region(chr,pos-2,pos+stretch+1,build));
      new_dna = dna;
      test_dna = dna;
      if new_dna(3)==old_nt
        new_dna(3)=new_nt;
        test_dna(3)='?';
      elseif new_dna(3)==my_seqrcomplement(old_nt)
        new_dna(3)=my_seqrcomplement(new_nt);
        test_dna(3)='?';
      else
        change='old_nt does not match';
        return;
      end
    else
      % try matching leftward
      target_end = 1;
      dna = upper(genome_region(chr,pos-stretch-1,pos+2,build));
      new_dna = dna;
      test_dna = dna;
      x=length(dna)-2;
      if new_dna(x)==old_nt
        new_dna(x)=new_nt;
        test_dna(x)='?';
      elseif new_dna(x)==my_seqrcomplement(old_nt)
        new_dna(x)=my_seqrcomplement(new_nt);
        test_dna(x)='?';
      else
        change='old_nt does not match';
        return;
      end
    end
    for strand=1:2
      for frame=0:2
        if strand==1
           msg = dna;
           new_msg = new_dna;
           test_msg = test_dna;
           % match top strand
        else
           msg = my_seqrcomplement(dna);
           new_msg = my_seqrcomplement(new_dna);
           test_msg = my_seqrcomplement(test_dna);
        end
        msg = msg(frame+1:end);
        new_msg = new_msg(frame+1:end);
        test_msg = test_msg(frame+1:end);

        frag = my_nt2aa(msg);
        new_frag = my_nt2aa(new_msg);
        test_frag = my_nt2aa(test_msg); 

        frag = regexprep(frag, '\*', 'X');      
        new_frag = regexprep(new_frag, '\*', 'X');
        test_frag = regexprep(test_frag, '\*', 'X');

        change_pos = find(test_frag=='?');

       
        % look for this protein fragment in protein_sequence
        matches = regexp(protein_sequence,frag);
        if ~isempty(matches)
          if ~silent, fprintf('%d %s strand%d frame%d\n', target_end, frag, strand, frame); end
          match.count = match.count + length(matches);
          protpos = matches + change_pos - 1;
          if match.count==1
            match.oldaa = protein_sequence(protpos);
            match.newaa = '?';
            for i=1:length(frag)
              if new_frag(i) ~= frag(i)
                match.newaa = new_frag(i);
                break;
              end
            end
            if isempty(match.newaa), match.newaa=match.oldaa; end
            match.protpos = [match.protpos protpos];
            if ~silent, fprintf('Matched: %s%d%s\n', match.oldaa, protpos, match.newaa); end
          end
        end
      end % next frame
    end % next strand
  end % next direction

  if match.count==0
    % no match
    return
  end
  if length(unique(match.protpos))==1
    protpos = unique(match.protpos);
    if ~silent, fprintf('Matched! %d\n', protpos); end
    oldaa = match.oldaa;
    newaa = match.newaa;
    change = [oldaa num2str(protpos) newaa];
    return;
  end
stretch = stretch + 3;
end % next stretch

end
