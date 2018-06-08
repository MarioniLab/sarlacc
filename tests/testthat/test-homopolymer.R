#Check by finding consencutive bases.
test_that("homopolymer regions should be found", {
    CHECKFUN1 <- function(seq){
        result <- lapply(seq, function(q){
            q <- as.character(q)
            q <- gsub("-","", q)
            seq_letters <- unique(strsplit(q, "", perl = TRUE)[[1]])
            q <- DNAString(as.character(q))
            homo_range <- bplapply(seq_letters,function(x){
                match_loc <- matchPattern(x, q)
                repeat_group <- t(.consencutive_finder(match_loc@ranges@start))
                cbind.data.frame(repeat_group, x)
            })
            homo_range <- data.table(do.call(rbind.data.frame, homo_range))
            colnames(homo_range) <- c("start", "end", "base")
            homo_ir <- IRanges(start = homo_range$start, end = homo_range$end)
            elementMetadata(homo_ir)$base <- as.character(homo_range$base)
            homo_ir <- homo_ir[order(homo_ir@start)]
            homo_ir[homo_ir@width != 1]#& homo_ir@elementMetadata$base != "N"]
        })
        IRangesList(result)
    }
    .consencutive_finder <- function(vec){
        breaks <- c(0, which(diff(vec) != 1), length(vec)) 
        sapply(seq(length(breaks) - 1), 
               function(x)vec[c((breaks[x] + 1), breaks[x+1])])
    }
    #Insertions should not affect the result
    test_seq <- DNAStringSet(c("ATGG--GGCCGTTTAA",
                               "ATG-GGGCCGTTTAA",
                               "ATGGGGCCGT-TTAA",
                               "ATGGGGC-CGTTTAA",
                               "---ATGGGGCCGTTTAA"))
    check1_result <- CHECKFUN1(test_seq)
    fun_result <- homopolymerFinder(test_seq)
    expect_identical(check1_result, fun_result)
})

#Check homopolymer region using rle
test_that("homopolymer regions should be found", {
    CHECKFUN2 <- function(seq){
        seq <- as.character(seq)
        test_rle <- lapply(strsplit(seq, split = ""), rle)
        test_rle <- lapply(test_rle, function(x){
            x$end <- cumsum(x$lengths)
            x$start <- x$end-x$lengths + 1
            x})
        ir <- lapply(test_rle, function(x){
            ir <- IRanges(start = x$start, end = x$end)
            elementMetadata(ir)$base <- x$values
            ir})
        
        IRangesList(lapply(ir,function(x)x[x@width!=1]))#&elementMetadata(x)$base!="N"]))
    }
    test_seq <- DNAStringSet(c("GGCTTTNN","GGCTTTNNTTTT"))
    check2_result <- CHECKFUN2(test_seq)
    fun_result <- homopolymerFinder(DNAStringSet(test_seq))
    expect_identical(fun_result, check2_result)
})

#Check cxx_find_homopolymers by using naive loop to walk through the sequence
test_that("the details of homopoylmer in query sequence should be found", {
    CHECKFUN3 <- function(seq){
        seq <- strsplit(as.character(seq), split="")
        homo_start <- c()
        homo_len <- c()
        cur_pos <- 1L
        homo_base <- c()
        str_index <- c()
        for(i in as.integer(1:length(seq))){
            cur_pos <- 1L
            cur_seq <- seq[[i]]
            while(cur_pos < length(cur_seq)){
                last_base <- cur_seq[cur_pos]
                cur_pos <- cur_pos + 1L
                cur_base <- cur_seq[cur_pos]
                if(last_base == cur_base){#& last_base!="N"){
                    homo_start <- c(homo_start, cur_pos - 1L)
                    homo_base <- c(homo_base, last_base)
                    len_counter <- 1L
                    while(last_base == cur_base & cur_pos <= length(cur_seq)){
                        len_counter <- len_counter + 1L
                        last_base <- cur_seq[cur_pos]
                        cur_pos <- cur_pos + 1L
                        cur_base <- cur_seq[cur_pos]
                    }
                    homo_len <- c(homo_len, len_counter)
                    str_index <- c(str_index, i - 1L)
                }else{
                    len_counter <- c()
                }
            }
        }
        list(str_index, homo_start, homo_len, homo_base)
    }
    test_seq <- DNAStringSet(c("GGCTTTNN","GGCAATN"))
    fun_result <- .Call(sarlacc:::cxx_find_homopolymers, test_seq)
    check3_result <- CHECKFUN3(test_seq)
    expect_identical(fun_result, check3_result)
})

