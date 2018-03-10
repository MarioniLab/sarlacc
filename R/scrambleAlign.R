#' @export
scrambleAlign <- function(adaptor1, adaptor2, reads, tolerance=100, gapOpening=1, gapExtension=5, match=1, mismatch=0)
# Scrambles the input sequence and performs the same 
# thing as adaptorAlign but with a scrambled input. 
#
# written by Aaron Lun
# created 10 March 2018
{
    pre.out <- .preprocess_input(adaptor1, adaptor2, reads)
    adaptor1 <- pre.out$adaptor1
    adaptor2 <- pre.out$adaptor2
    reads <- pre.out$reads

    # Scrambling the start and end of the read sequences.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back
    scrambled.start <- .scramble_input(reads.start)
    scrambled.end <- .scramble_input(reads.end)

    .internal_align(adaptor1, adaptor2, scrambled.start, scrambled.end,
                    gapOpening=gapOpening, gapExtension=gapExtension,
                    match=match, mismatch=mismatch)
}

.scramble_input <- function(seqs) 
# Scrambles the input sequences. Uses R loops,
# but it should be fast enough for our purposes.   
{
    for (i in seq_along(seqs)) {
        current <- seqs[[i]]
        seqs[[i]] <- current[sample(length(current))]
    }
    return(seqs)
}
