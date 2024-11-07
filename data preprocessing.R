"
PROJECT: METAGENOMICS ANALYSIS ON TEXTILE WASTEWATER TREATMENT PLANT 

DESCRIPTION: This is the first program script of entire data analysis, for MAGs data preprocessing.
"

# (1) Import packages
##########################################
# source("https://bioconductor.org/biocLite.R"); biocLite("Biostrings")
# install.packages("seqinr")
# install.packages("reshape")
# install.packages("tidyverse")
# install.packages("scales")
# install.packages('data.table')
library("Biostrings")
library("seqinr")
library("tidyverse")
library("scales")
library('data.table')
##########################################


# (2) Functions definition
##########################################
A_in_Codon <- function (seq, forceToLower = TRUE, exact = FALSE, NA.A = NA, oldGC = FALSE) 
{ # FUNCTION:: a function to return base-A ration in a coding sequence. Copied from and modified from seqinr::GC().
  if (length(seq) == 1 && is.na(seq)) 
    return(NA)
  if (nchar(seq[1]) > 1) 
    stop("sequence is not a vector of chars")
  if (forceToLower) 
    seq <- tolower(seq)
  nc <- sum(seq == "c")
  ng <- sum(seq == "g")
  na <- sum(seq == "a")
  nt <- sum(seq == "t")
  if (oldGC) {
    warning("argument oldGC is deprecated")
    return((nc + ng)/length(seq))
  }
  if (!exact) {
    if (na + nc + ng + nt == 0) {
      result <- NA.A
    }
    else {
      result <- (na)/(na + nc + ng + nt)
    }
  }
  else {
    ngc <- ng + nc
    nat <- na + nt
    ngc <- ngc + sum(seq == "s")
    nat <- nat + sum(seq == "w")
    if (na + nc != 0) {
      nm <- sum(seq == "m")
      ngc <- ngc + nm * nc/(na + nc)
      nat <- nat + nm * na/(na + nc)
    }
    if (ng + nt != 0) {
      nk <- sum(seq == "k")
      ngc <- ngc + nk * ng/(ng + nt)
      nat <- nat + nk * nt/(ng + nt)
    }
    if (ng + na != 0) {
      nr <- sum(seq == "r")
      ngc <- ngc + nr * ng/(ng + na)
      nat <- nat + nr * na/(ng + na)
    }
    if (nc + nt != 0) {
      ny <- sum(seq == "y")
      ngc <- ngc + ny * nc/(nc + nt)
      nat <- nat + ny * nt/(nc + nt)
    }
    if (na + nc + ng != 0) {
      nv <- sum(seq == "v")
      ngc <- ngc + nv * (nc + ng)/(na + nc + ng)
      nat <- nat + nv * na/(na + nc + ng)
    }
    if (na + nc + nt != 0) {
      nh <- sum(seq == "h")
      ngc <- ngc + nh * nc/(na + nc + nt)
      nat <- nat + nh * (na + nt)/(na + nc + nt)
    }
    if (na + ng + nt != 0) {
      nd <- sum(seq == "d")
      ngc <- ngc + nd * ng/(na + ng + nt)
      nat <- nat + nd * (na + nt)/(na + ng + nt)
    }
    if (nc + ng + nt != 0) {
      nb <- sum(seq == "b")
      ngc <- ngc + nb * (nc + ng)/(nc + ng + nt)
      nat <- nat + nb * nt/(nc + ng + nt)
    }
    if (ngc + nat == 0) {
      result <- NA.A
    }
    else {
      result <- na/(ngc + nat)
    }
  }
  return(result)
}

A3 <- function (seq, frame = 0, ...) 
{# FUNCTION:: a function to return base-A ratio in 3rd-position of a codon
  if (nchar(seq[1]) > 1) {
    warning("sequence is not a vector of chars, I'm trying to cast it into one")
    seq <- s2c(seq[1])
  }
  if (frame != 0) 
    seq <- seq[(1 + frame):length(seq)]
  A_in_Codon(seq[seq(3, length(seq), by = 3)], ...)
}

get_nucleo_comp <- function(seq, nucleotide='a')
{# FUNCTION: to count te number of specific nucleotide in a sequence `seq`.
  if(!nucleotide %in% c("a", "t", "g", "c", "u")){
    print("`nucleotide` should be either on of a/t/g/c/u.")
    return()
  }
  seqinr::count(seq = seq, 1)[nucleotide]
}

get_atoms_comp <- function(seq_df, atom="O")
{# FUNCTION: to calculate codon-average number of atoms (O/H/C/N) of each CDS in a sequence dataframe
  switch (atom,
          O = ((seq_df$a * 0) + (seq_df$g * 1) + (seq_df$c * 1) + (seq_df$t * 2))/(seq_df$Length/3),
          N = ((seq_df$a * 5) + (seq_df$g * 5) + (seq_df$c * 3) + (seq_df$t * 2))/(seq_df$Length/3),
          H = ((seq_df$a * 5) + (seq_df$g * 5) + (seq_df$c * 5) + (seq_df$t * 4))/(seq_df$Length/3),
          C = ((seq_df$a * 5) + (seq_df$g * 5) + (seq_df$c * 4) + (seq_df$t * 5))/(seq_df$Length/3),
          'Please give correct atom(N/O/C/H)')
}

check_EggNOG <- function(x)
{# FUNCTION: to check EggNOG category by the COG.cat, belonging to >=2 categories will be defined as "UNIDENTIFIED".
  Eggnog_cat <- list("INFORMATION STORAGE AND PROCESSING"=c("J", "A", "K", "L", "B"), 
                     "CELLULAR PROCESSES AND SIGNALING"=c("D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X"),
                     "METABOLISM"=c("C", "G", "E", "F", "H", "I", "P", "Q"),
                     "POORLY CHARACTERIZED"=c("R", "S"))
  ifelse(x%in%Eggnog_cat$`INFORMATION STORAGE AND PROCESSING`, names(Eggnog_cat[1]), 
         ifelse(x%in%Eggnog_cat$`CELLULAR PROCESSES AND SIGNALING`, names(Eggnog_cat[2]), 
                ifelse(x%in%Eggnog_cat$METABOLISM, names(Eggnog_cat[3]), 
                       ifelse(x%in%Eggnog_cat$`POORLY CHARACTERIZED`, names(Eggnog_cat[4]), 
                              "UNIDENTIFIED"))))
}

codon_counting_eff <- function(seq)
{# FUNCTION: to calculate effective number of codons in sequence `seq`.
  t(uco(seq, index = "eff"))
} 


isCDS <- function(seq)
{# FUNCTION: to check if the input sequence is a CDS (i.e. [length %% 3 == 0] AND [start with ATG] AND [stop with TAG/TAA/TGA])
  if (length(seq) %% 3 != 0){
    return(FALSE)
    }
  if (paste(seq[1:3], collapse="") != "atg") {
    return(FALSE)
    }
  if (!paste(seq[(length(seq)-2):length(seq)], collapse = "") %in% c("tag", "taa", "tga")){
    return(FALSE)
    }
  return(TRUE)
}
##########################################


# (3) Main Program
##########################################
# Step 0: Set working environment
this.dir <- getwd()
INTERMETIDATE_dir <- "./intermediate_RDS/"
dir.create(INTERMETIDATE_dir, showWarnings = FALSE)
FFN <- "./ffn"
FAA <- "./faa"
EGG <- "./emapper"
PHYLA <- "HQ_BINS_LIST.csv"

# Step 1: Create a dataframe summerizing all .ffn and .emapper.annotations.
print('STAGE 1:: Dataframe is being processed...')
seqs_df <- data.frame()
ffn_files <- list.files(FFN, pattern = 'ffn')
faa_files <- list.files(FAA, pattern = 'faa')
# egg_files <- list.files(EGG, pattern = 'annotations')
phyla_temp <- read_csv(PHYLA)
pb <- txtProgressBar(min = 0, max = length(ffn_files), char = "#", style = 3)
k <- 0               # progressbar counter
initializer <- TRUE  # dataframe-initializing indicator
for (ffn in ffn_files){
  k <- k + 1
  str.location <- gregexpr(pattern ='\\.', ffn)
  index <- substr(ffn, 1, str.location[[1]][2]-1)
  fasta_seqinr <- read.fasta(file = file.path(FFN, ffn))
  # extract data directly from .ffn file.
  new_df <- data.frame(seq_label = names(fasta_seqinr),
                       gc = sapply(X = fasta_seqinr, FUN = GC),
                       gc3 = sapply(X = fasta_seqinr, FUN = GC3),
                       a3 = sapply(X = fasta_seqinr, FUN = A3),
                       a = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='a'),
                       g = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='g'),
                       t = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='t'),
                       c = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='c'),
                       length = sapply(X = fasta_seqinr, FUN = length))
  new_df$Annotation <- sapply(getAnnot(fasta_seqinr), substr, 17, length(getAnnot(fasta_seqinr)))
  new_df$Carbon   <- ((new_df$a * 5) + (new_df$g * 5) + (new_df$c * 4) + (new_df$t * 5))/(new_df$length/3)   # U has 4 carbons and T has 5
  new_df$Hydrogen <- ((new_df$a * 5) + (new_df$g * 5) + (new_df$c * 5) + (new_df$t * 4))/(new_df$length/3)   # U has 4 Hydrogens and T has 6
  new_df$Nitrogen <- ((new_df$a * 5) + (new_df$g * 5) + (new_df$c * 3) + (new_df$t * 2))/(new_df$length/3)   # U and T have the same Nitrogen content = 2
  new_df$Oxygen   <- ((new_df$a * 0) + (new_df$g * 1) + (new_df$c * 1) + (new_df$t * 2))/(new_df$length/3)   # U and T have the same Oxygen content = 2
  new_df$Adenine  <- (new_df$a)/(new_df$length/3)
  new_df$Thymine  <- (new_df$t)/(new_df$length/3)
  new_df$Guanine  <- (new_df$g)/(new_df$length/3)
  new_df$Cytosine <- (new_df$c)/(new_df$length/3)
  new_df$gc3_to_gc <- new_df$gc3 / new_df$gc 
  new_df$sequence <- sapply(fasta_seqinr, getSequence)
  new_df$isCDS    <- sapply(new_df$sequence, isCDS) #ORFs must be started with ATG and ended with STOP-codon, as well as length%3!=0; CDSs is subset of ORFs.
  new_df$sample_id    <- index
  # No EggNOG file for AS/AD files, commented this section.
  # # extract data directly from .emapper.annotations file.
  # egg = read.delim(file = paste('emapper/', egg_files[grep(egg_files, pattern = index)], sep=""), header = TRUE, sep = "\t", skip = 3)
  # rownames(egg) <- egg$X.query_name
  # pre_df <- subset.data.frame(egg, subset = TRUE, select =  c(Preferred_name, eggNOG.free.text.desc., KEGG_ko, COG.Functional.cat.))
  # pre_df <- pre_df[-(nrow(pre_df)-2):-nrow(pre_df),] # Remove useless rows
  # new_df <- merge.data.frame(new_df, pre_df, by=0, all.x = TRUE)
  # combine all data together.
  if (initializer){
    seqs_df <<- new_df
    initializer <- FALSE
  }else {
    seqs_df <- rbind(seqs_df, new_df)
  }
  setTxtProgressBar(pb, k)
}
seqs_df$Row.names <- NULL
print("======STAGE 1 completed!======")

# # Step 2: Concatenate phyla dataframe with the dataframe above.
print("STAGE 2:: Loading phyla information for database concantenation..")
# phyla_temp$X1 <- NULL
phyla_temp$Sample_type <- phyla_temp$Type 
phyla_temp$Sample <- sapply(phyla_temp$BIN, function(x) gsub(".fasta", "", x))
phyla_temp <- phyla_temp %>%  drop_na()
# phyla_temp <- phyla_temp %>%
#   select(Sample, Sample_type, taxa, GTDBTK) %>%
#   rename(Phylophlan = taxa) %>%
#   drop_na()
# phyla_temp$phyla <-substr(phyla_temp$GTDBTK,
#                           gregexpr(pattern ="p_", phyla_temp$GTDBTK),
#                           regexpr(pattern =';c_', phyla_temp$GTDBTK)-1)

seqs_df <- merge(seqs_df, phyla_temp, by.x="sample_id", by.y="Sample", all.x=TRUE)
print('======STAGE 2 completed!======')

# Step 3: Generate MAG-level dataframe (averaging by MAGs).
  # Originally, group by (Sample_type, sample_id, phyla, GTDBTK, Phylophlan).
print("STAGE 3:: Computing bin-level data.frames...")
seqs_df <- seqs_df %>% filter(isCDS == TRUE)
bins_df <- seqs_df %>% 
  select(sample_id, Sample_type, Oxygen, Nitrogen, Hydrogen, Carbon, 
         Adenine, Thymine, Guanine, Cytosine, length, gc, gc3, gc3_to_gc) %>% 
  group_by(Sample_type, sample_id) %>% 
  summarise(Sum_length=sum(length), gc_mean=mean(gc), gc3_mean=mean(gc3), GC_avg_ratio=mean(gc3_to_gc), cds_count=n(),
            Adenine=mean(Adenine), Thymine=mean(Thymine), Guanine=mean(Guanine), Cytosine=mean(Cytosine),
            Oxygen=mean(Oxygen), Nitrogen=mean(Nitrogen), Hydrogen=mean(Hydrogen), Carbon=mean(Carbon))
#   # Set aerobic/anaerobic zone info (USAB = "Anaerobic")
# bins_df$Sample_type <- paste(bins_df$Sample_type, " ", sep = "")
# bins_df$zone <- ifelse(str_detect(bins_df$Sample_type, "UASB"), 
#                        "Anaerobic", 
#                        substr(x=bins_df$Sample_type, start=1, stop=str_locate(bins_df$Sample_type, " ")-1))
# bins_df$zone <- ifelse(str_detect(bins_df$zone, "Aerobic"), "Aerobic", "Anaerobic")

print("Saving all databases...")
saveRDS(seqs_df, file = paste(INTERMETIDATE_dir, "Database_CDSs.rds", sep=""))
saveRDS(bins_df, file = paste(INTERMETIDATE_dir, "Database_MAGs.rds", sep=""))
print("======STAGE 3 completed!======")
print('Data Preprocessing FINISHED!')
##########################################

