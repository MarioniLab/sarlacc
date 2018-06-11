readsMerger <- function(reads, umi, minimap2_dir ="/Users/LawCheukTing/Desktop/tools/minimap2", 
                        minimap2_arg = "-t 6", min_match = 0.7){
    elementMetadata(reads)$UMI <- umi
    
    #Run the first round of minimap2
    writeQualityScaledXStringSet(reads, "cluster4minimap2.fastq")
    minimap2 <- file.path(minimap2_dir, "minimap2")
    minimap2_arg <- paste0(minimap2_arg, " -x ava-ont -c cluster4minimap2.fastq cluster4minimap2.fastq")
    cm <- paste(minimap2, minimap2_arg, sep = " ")
    paf_raw <- system(cm , intern = TRUE)
    paf_raw <- .read_paf(paf_raw)
    clean_paf <- .paf_filter(paf_raw, min_match = min_match)
    cluster_list <- .cluster_paf(clean_paf)
    
    #Recreate groups for sequences without pair
    is_missing <- !names(reads)%in%unlist(cluster_list)
    pseudo_cluster <- split(names(reads[is_missing]), f=1:length(reads[is_missing]))
    cluster_list <- c(cluster_list, pseudo_cluster)
    
    #Clustered by read similarity and then grouped by UMI
    clustered_reads <- lapply(cluster_list, function(x)reads[x])
    grouped_umi <- mclapply(clustered_reads, function(x){
        umiGroup2(elementMetadata(x)$UMI, max.lev1 = 2)+1},mc.cores = 6L)
    clustered_reads <- do.call(c, clustered_reads)
    
    #Initialization
    seq_group <- 0
    iteration <- 0
    finish <- FALSE
    while(!finish){
        
        #Combine and rename all small groups
        seq_group_before <- unique(seq_group)
        seq_group <- grouped_umi[[1]]
        if(length(grouped_umi)>1){
            for(i in 2:length(grouped_umi)){
                seq_group <- c(seq_group, grouped_umi[[i]] + max(seq_group))
            }  
        }
        seq_group_after <- unique(seq_group)
        
        #Check group length
        finish <- length(seq_group_before) == length(seq_group_after)
        
        #multiple sequence alignment and create consenus sequnece
        elementMetadata(clustered_reads)$grouped_by_UMI <- factor(seq_group)
        multi_aln <- multiReadAlign(clustered_reads, groups = clustered_reads@elementMetadata$grouped_by_UMI)
        con_seq <- consensusReadSeq(multi_aln)
        
        #Record the origins of the consenus sequence
        seq_history <- split(names(clustered_reads), f = clustered_reads@elementMetadata$grouped_by_UMI)
        elementMetadata(con_seq)$seq_hisotry <- CharacterList(seq_history)
        names(con_seq) <- names(multi_aln)
        
        #Run minimap2 again on consenus sequences
        writeQualityScaledXStringSet(con_seq, "cluster4minimap2.fastq")
        paf_raw <- system(cm , intern = TRUE)
        paf_raw <- .read_paf(paf_raw)
        clean_paf <- .paf_filter(paf_raw, min_match = min_match)
        cluster_list <- .cluster_paf(clean_paf)
        
        #Those missing consenus sequences were added again.
        #(This function can be speeded up by not merging again and 
        # only combine them after multiReadAlign but before running minimap2)
        
        still_in <- unlist(cluster_list)
        missing_group <- names(con_seq[!names(con_seq)%in%still_in])
        cluster_list <- c(cluster_list, split(missing_group, missing_group))
        
        
        #Harvest consenus sequences
        if(finish){
            names(con_seq) <- lapply(con_seq@elementMetadata$seq_hisotry, function(x)paste(x, collapse = " "))
        }else{
            #not finished: merge conseq --> combine seq history --> extract UMI
            cluster_list <- mclapply(cluster_list, function(x){
                unlist(con_seq[as.character(x)]@elementMetadata$seq_hisotry)
                },mc.cores = 6L)
            grouped_umi <- lapply(cluster_list, function(x){
                umiGroup2(clustered_reads@elementMetadata$UMI[x], max.lev1 = 2) + 1
                }) 
        }
        iteration <- iteration + 1
    }
    cat("iteration: ", iteration)
    con_seq
}

#cluster sequneces in paf format
.cluster_paf <- function(paf){
    igraph_df <- cbind(paf$qname, paf$tname)
    igraph_df <- graph.data.frame(igraph_df)
    igraph_com <- components(igraph_df)
    cluster_list <- lapply(1:igraph_com$no,function(x){
        names(igraph_com$membership[igraph_com$membership == x])
    })
    cluster_list
}

#filter paired alignment with low similarity
.paf_filter <- function(paf, min_match=0.7){
    mean_length <- (paf$qlength+paf$tlength)/2
    pw_id <- paf$matches/mean_length
    paf$pw_id <- pw_id
    matches_filter <- which(pw_id < min_match)
    discard_loc <- matches_filter
    if(length(discard_loc)!=0)paf <- paf[-discard_loc,]
    paf
}

#read paf format
.read_paf <- function(paf){
    paf <- sub("s2.*dv","dv", paf, perl = TRUE)
    paf <- sub("zd.*cg", "cg", paf, perl = TRUE)
    col_num <- length(unlist(gregexpr("\t",paf[1])))+1
    paf <- paste(paf,collapse = "\n")
    paf_colClasses <- c("character",rep("numeric",4),"character", 
                        rep("numeric",6), rep("character",col_num-12))
    paf <- fread(paf, colClasses = paf_colClasses, header= FALSE)
    colnames(paf)[1:12] <- c("qname","qlength","qstart","qend","qstrand",
                             "tname","tlength","tstart","tend","matches",
                             "aln_length","quality")
    paf
}


#Test
setwd("/Users/LawCheukTing/Dropbox (Cambridge University)/Univeristy of Cambridge/Academic/Internship/Nanopore_analysis/test_20180522/201805_ONT_ligation/fastq")
umi_seq <- readFasta("umi_seq")
names(umi_seq@sread) <- umi_seq@id
umi_seq <- umi_seq@sread
full_umi_loc <- which(nchar(umi_seq)==14)
umi_seq <- umi_seq[full_umi_loc]
reads <- readQualityScaledDNAStringSet("20180605_chopped.fastq")
reads <- reads[full_umi_loc]
reads <- reads[1:2000]
umi_seq <- umi_seq[1:2000]
system.time(merged_conseq <- readsMerger(reads, umi_seq))


