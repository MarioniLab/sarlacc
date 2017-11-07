#For initial quality assessment and adaptor chopping, the reads need to be supplied in a FASTA file.
#Additionally, the UMI sourrounding adaptor sequence has to be supplied as a string element with the UMI marked as N's.

UMI <- DNAString("AAGCAGTGGTATNNNNNNCAACGCAGAGT")
fastqseq <- readFastq("file.fastq", withIds=TRUE)
