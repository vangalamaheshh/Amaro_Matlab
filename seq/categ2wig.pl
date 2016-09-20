# given an example wiggle file,
# reformats a set of per-chromosome category files in categdir
# into a wiggle file with the same intervals etc.
#
# procedure:
#   (1) read through wiggle file
#       -- save all "track" and "fixedStep" lines
#       -- give error if any "variableStep" lines are encountered: these are not supported yet
#       -- for each fixedStep line, save the number of positions listed after it
#   (2) allocate memory to store the category of every position in the wiggle file
#   (3) extract the positions from each chromosome's (binary) category file
#   (4) write the category wiggle file
#
# Mike Lawrence 2009-11-17

if ($#ARGV+1 != 3) {
  print "Usage: categ2wig wigfile categdir outfile\n";
  exit;
}

$WIGFILE = $ARGV[0];
$CATEGDIR = $ARGV[1];
$OUTFILE = $ARGV[2];

# open each chromosome binary category file (ensuring they exist)
my @chr_handle;
for ($chr=1;$chr<=24;$chr++) {
    local *FILE;
    open(FILE,"<$CATEGDIR/chr$chr.bin") || die;
    push(@chr_handle,*FILE);
}


if (!(-e $WIGFILE)) { print "$WIGFILE not found!\n"; exit; }
$nlines = `wc -l $WIGFILE | awk '{print \$1}'`;
chomp $nlines;

# allocate memory
@categ = ((0) x $nlines);
$categ[$nlines-1] = 0;

@sec = `grep -P -n '^\\D' $WIGFILE`;
$nsec = $#sec+1;
@sec_fs = ((0) x $nsec);
@sec_line = ((0) x $nsec);
@sec_chr = ((0) x $nsec);
@sec_start = ((0) x $nsec);
for ($i=0;$i<$nsec;$i++) {
  ($sec_line[$i]) = ($sec[$i] =~ m/(\d*):/);
  if ($sec[$i] =~ m/^\d*:track/) { next; }
  elsif ($sec[$i] =~ m/^\d*:variableStep/) { die "variableStep not supported"; }
  elsif ($sec[$i] =~ m/^\d*:fixedStep/) {
    $sec_fs[$i] = 1;
    $result = (($sec_chr[$i],$sec_start[$i]) = ($sec[$i] =~ m/\d*:fixedStep chrom=chr(\S*) start=(\d*) step=1/));
    die "Invalid fixedStep parameters: $sec[$i]\n" if $result!=2;
    $sec_chr[$i]=23 if $sec_chr[$i] eq "X";
    $sec_chr[$i]=24 if $sec_chr[$i] eq "Y";
    die "Invalid chromosome $sec_chr[$i]" if ($sec_chr[$i]<1 || $sec_chr[$i]>24);
    }
  else { die "Invalid line: $sec[$i]\n"; }
}

# compute length of each section

@sec_length = ((0) x $nsec);
for ($i=0;$i<$nsec;$i++) {
    next if !$sec_fs[$i];
    if ($i+1==$nsec) { $last_line = $nlines; }
    else { $last_line = $sec_line[$i+1]-1; }
    $sec_length[$i] = $last_line - $sec_line[$i] + 1;
}

# retrieve categories for each section

for ($i=0;$i<$nsec;$i++) {
    next if !$sec_fs[$i];
    seek($chr_handle[$sec_chr[$i]],4*($sec_start[$i]-1),0);
    for($j=0;$j<$sec_length[$i];$j++) {
	$int = read($chr_handle[$sec_chr[$i]],$buffer,4);
	print "j=$j int=$buffer[0] $buffer[1] $buffer[2] $buffer[3]\n";
    }
}
