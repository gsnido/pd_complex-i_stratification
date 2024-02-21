if(!AnnoLoc %in% list.dirs(full.names = F, recursive = T)){
  dir.create(AnnoLoc)
}

AnnoLoc = paste0(AnnoLoc, "/")


URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz"
AssemblyFilename = "gencode.v35.annotation.gff3.gz"

if(!file.exists(paste0(AnnoLoc, AssemblyFilename))){
  download.file(URL, destfile = paste0(AnnoLoc, AssemblyFilename))
}

txdbFilename  = paste0(paste0(strsplit(AssemblyFilename, "\\.")[[1]][c(2,4)], collapse = "."), "_DB")
#AnnoFilename = paste0("annoFileCollapsed_", strsplit(txdbFilename, "\\.")[[1]][1], ".Rds")
GTF_file = paste0(AnnoLoc, AssemblyFilename)
geneNameFile = paste0(AnnoLoc, strsplit(AssemblyFilename, "\\.")[[1]][2], "_geneNames.Rds")

#Get genome annotations
if(!file.exists(paste0(AnnoLoc, txdbFilename))){
  txdb <- makeTxDbFromGFF(paste0(AnnoLoc, AssemblyFilename))
  saveDb(txdb, paste0(AnnoLoc, txdbFilename))
} else {
  txdb <- loadDb(paste0(AnnoLoc, txdbFilename))
}

# if(!file.exists(paste0(AnnoLoc, AnnoFilename))){
#   annoFileCollapsed <- GetGenomeAnno_GENCODE(GTF_file = paste0(AnnoLoc, AssemblyFilename), txdb = txdb)
#   saveRDS(annoFileCollapsed, paste0(AnnoLoc, AnnoFilename))
# } else {
#   annoFileCollapsed <- readRDS(paste0(AnnoLoc, AnnoFilename))
# }
