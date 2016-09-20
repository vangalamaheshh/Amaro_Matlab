#inputs: PLATE, SIF_path, bird_dir, by_allele_dir, out_dir
#collects and matches birdseed normals with byallele copy number files and outputs a txt file for use in scatter gather

args<-commandArgs()
args<-args[7:length(args)]
PLATE <- as.character(args[1])
SIF_path <- as.character(args[2])
bird_dir<-as.character(args[3])
by_allele_dir <- as.character(args[4])
out_dir<-as.characte r(args[5])

sif<-read.table(SIF_path, sep="\t", header=T)

#####match normals and tumor samples

normals<-sif[which(sif$TUMOR_NORMAL=="Normal"),]
tumors<-sif[which(sif$TUMOR_NORMAL=="Tumor"),]

match_check<-setdiff(normals$PARTICIPANT_ID, tumors$PARTICIPANT_ID)
if (length(match_check)>0)
  {
    print(c("warning ", match_check, " not matched!"))
  }

reorder<-match(normals$PARTICIPANT_ID, tumors$PARTICIPANT_ID)
tumors<-tumors[reorder,]

bird_filenames <- list.files(bird_dir, pattern=c(PLATE, "_.*"))
by_allele_filenames <- list.files(by_allele_dir, pattern=c(PLATE, "_.*"))

COMMANDS<-NULL
bird_paths<-NULL
by_allele_paths<-NULL

for(i in 1:length(normals))
  {
    bird_paths[i]<-paste(bird_dir, bird_filenames[grep(normals$ARRAY[i], bird_filenames)], sep="")
    by_allele_paths[i]<-paste(by_allele_dir, by_allele_filenames[grep(tumors$ARRAY[i], by_allele_filenames)], sep="")
    COMMANDS[i]<-paste("matlab -nodesktop -nosplash -r run_LOH_run_sample.sh ", bird_paths[i], " ", by_allele_paths[i], " ", out_dir)
  }

write.table(COMMANDS, file=paste(PLATE, "_scat_gat.txt", sep=""), sep="\n", quote=F, row.names=F)
