#! /bin/csh -f

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
#   (3) load each chromosome's category file, extract the positions
#   (4) write the category wiggle file
#
# Mike Lawrence 2009-11-17

if ($#argv != 3) then
  echo "Usage: categ2wig wigfile categdir outfile"
  exit 1
endif

set wigfile = $1
set categdir = $2
set outfile = $3

