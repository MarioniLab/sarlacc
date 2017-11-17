#bash commands to make alignment using minimap2
#takes a .fastq file as input and returns a .sam file
#for further pipeline analysis, .sam has to be transformed into a .txt file

bsub -e min1.err -o min1.out -M 32000 "/nfs/research2/marioni/florian/minimap2/minimap2 -ax map-ont Mus_musculus.GRCm38.dna.sm.toplevel.fa chopped_c4_reads.fastq > align_chop_c4.sam"

#to transform .sam into .txt the following needs to be executed

awk '{print $1,$2,$3,$4,$5,$6,$10}' alignment.sam > align_pos.txt

#it returns a .txt file containing the first six columns from the .sam file





