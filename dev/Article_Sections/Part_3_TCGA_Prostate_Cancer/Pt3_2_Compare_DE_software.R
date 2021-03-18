
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #  #




# Intention of Script

# Analysis of Prostate Cancer TCGA data by common DE tools 

# Comparing which genes are declared DE for edgeR, limma + voom, DESeq2, Bayseq, TABI 


# Section 1 Comparing which transcript is declared DE for different group numbers - edgeR 

# Section 2 Comparing which transcripts are declared DE for different software 

# Section 3 Comparing TABI DE transcripts, with common DE software





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #  #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Pt 1 Comparing which transcripts are decalred DE for different group numbers - edgeR

# Including purity.score as covariate

# 3 groups, 7 groups for each DE sofware

tidybulk_DE()


# Load TCGA Prostate cancer data 

ed_PC_TCGA <- readRDS("~/TABI/dev/ed_PC_TCGA.rds")

make_groups(ed_PC_TCGA, 
            .sample = Sample_ID,
            .transcript = transcript,
            .abundance = read_count_normalised, 
            n_group = 3,
            x_cord = seq(from = 0, to = 8, by = 1))
