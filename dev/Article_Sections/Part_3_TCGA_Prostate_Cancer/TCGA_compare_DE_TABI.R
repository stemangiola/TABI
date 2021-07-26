
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Intention of Script

# Analysis of Prostate Cancer TCGA data by common DE tools 

# Comparing which genes are declared DE for edgeR, limma + voom, DESeq2, Bayseq, TABI 


# Pt 1 Comparing which transcript is declared DE for different group numbers - edgeR 

# Pt 2 Comparing which transcripts are declared DE for different software 

# Pt 3 Comparing TABI DE transcripts, with common DE software

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Pt 1 Comparing which transcripts are decalred DE for different group numbers - edgeR

# Including purity.score as covariate

# 3 groups, 7 groups for each DE sofware


ed_PC_TCGA <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/ed_PC_TCGA.rds")


# Load functions from Functions_Simulate_RNA_seq_df using box 
options(box.path = setwd("~/TABI/dev/Article_Sections"))   


#' @export
box::use(./TABI_Article_Functions) 

# If you want details and example usage 
# of the functions used with module TABI_Article_Functions
# use box::help(TABI_Article_Functions$function_name)

# e.g. 

box::help(TABI_Article_Functions$sigmoidal_sim_df)

set.seed(150)

ed_PC_TCGA_sample = (
  ed_PC_TCGA %>%
  select(sample = Sample.ID, 
         CAPRA_S = 'CAPRA-S', 
         transcript,read_count_normalised = 'read count normalised', 
         purity.score) %>% 
  na.omit() %>% 
  split(.$transcript)
                      ) [sample.int(size = 1000,
                                    n = length(levels(ed_PC_TCGA$transcript)))]
  #group_by(transcript) %>% 
  #tidyr::nest() ) [sample.int(size = 100,
  #                            n = length(levels(ed_PC_TCGA$transcript))), ] 


# ed_PC_TCGA_sample_1 = ed_PC_TCGA_sample %>% 
#   tidyr::unnest() %>% 
#   ungroup() %>% 
#   #filter(transcript %in% sample(levels(transcript),100)) %>% 
#   split(.$transcript)


n<-15

plan(multisession, 
     workers = n)



TABI_TCGA_fit_slopes_df = ed_PC_TCGA_sample  %>% 
  future_map(
    ~TABI::TABI_glm(
      .data = .,
      ~CAPRA_S+ purity.score,
      .abundance = read_count_normalised, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = transcript, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 2000, # total number of iterations
      warmup = 1000
    ) %$% 
      fit) 

saveRDS(TABI_TCGA_fit_slopes_df, file = "TABI_TCGA_fit_slopes_df.rds")


# ed_PC_TCGA_sample 
  



