library(Biostrings)
library(sarlacc)
testthat::test_that("The error profile of Nanopore sequencing data should be found.",{
    CHECKFUN <- function(aln){
        #Deletion and mutation
        aln_df <- as.data.frame(as.matrix(aln))
        aln_df[] <-  lapply(aln_df, factor, levels=c('A','C', 'G', 'T',"-"))
        atcg_count <- do.call(rbind.data.frame,lapply(aln_df,function(x)(table(x))))
        atcg_count <- cbind.data.frame(atcg_count,NA)
        colnames(atcg_count) <- c('A','C', 'G', 'T',"deletion", "insertion")
        ref <- gsub("-","",alignedSubject(aln)[1])
        ref <- strsplit(ref, split="")[[1]]
        atcg_count <- cbind(base=ref, atcg_count)
        
        #Transition matrix
        transition_mat <- lapply(c("A","C","G","T"),function(x){
            trans_count <- colSums(atcg_count[atcg_count$base==x,][,c("A","C","G","T")])
        })
        transition_mat <- do.call(rbind,transition_mat)
        rownames(transition_mat) <- c("A","C","G","T")
        
        #Insertion
        subject_vec <- strsplit(as.character(alignedSubject(aln)), split="")
        ins_length <- lapply(subject_vec, function(x){
            rle_result<- rle(x)
            ins_pos <- ins_length <- c()
            true_pos <- 0
            for(i in seq_along(rle_result$lengths)){
                if(rle_result$values[i]=="-"){
                    ins_length <- c(ins_length,rle_result$lengths[i])
                    ins_pos <- c(ins_pos, true_pos+1)
                }else{
                    true_pos <- true_pos+rle_result$lengths[i]
                }
            }
            names(ins_length)<- ins_pos
            ins_length
        })
        ins_length <- unlist(ins_length)
        for_split <- factor(as.numeric(names(ins_length)), levels = 1:(length(ref)+1))
        pos_length <- split(ins_length,f=for_split)
        length_count <- lapply(pos_length, table)
        insertion_count <- as(length_count, "IntegerList")
        atcg_count <- DataFrame(rbind(atcg_count, NA))
        atcg_count$insertion <- insertion_count
        
        #output
        atcg_count$base<- as.character(atcg_count$base)
        list(full=atcg_count, transition=transition_mat)
    }
    aln <- pairwiseAlignment(subject=DNAString(c("GGAAACGATCAGCTACGAACACT")), 
                             pattern=DNAStringSet(c("GGAACGTCAGCGGTACGAAACACTAAAA",
                                                    "GGAACGTCAGCGGTACGAAACACT",
                                                    "GGGGGAAACGATCAGCTACGAACACT",
                                                    "ACGATCAGCTACGAACACT")))
    check_out <- CHECKFUN(aln)
    result <- errorFinder(aln)
    testthat::expect_equal(check_out, result)
    }
)




