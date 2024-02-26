# put load library commands at the start of the script for easy loading
library(msa)
library(Biostrings)

# Sequences are being aligned and read using the ReadDNAStringSet function
# printed the variable to make sure all the DNA samples imported properly
# make sure to put sequences.fasta in your Bioinformatics folder
seq<- readDNAStringSet("/Users/arlinemartinez/Downloads/sequences.fasta", format="fasta")
seq

# loading the Package msa to run a multiple sequence alignment on the imported DNA samples
library(msa)
#Running a multiple sequence alignment on the DNA sequences imported
alignment<- msa(seq, method= "ClustalW")
# Prints the alignment
alignment
# Shows complete view of alignments
print(alignment, show= "complete")

# The DNAMulitpleAlignment function autocmatically adds color over the muliple sequence alignment,which helps with visualization and comaprison of sequences 
dnaseq<-DNAMultipleAlignment(alignment)
#I wanted to see the alignment with color to beter identify any mutations in the sequence
print(dnaseq)
# Using the detail dnaseq, I can see a full list of the DNA sequences and the DNA sequence for Homo_sapien_6 has a mutation, specifcally a deletion 
detail(dnaseq)

# The gene of the homo sapien with the mutation is hbb gene, the accession number with the best match was LC121775.1, I ran this on blastn. 
# I am assigning the sequence for Homo_sapien 6 to the variable to aa_seq
aa_seq <- seq$Homo_sapiens_6 # easy to extract from the original sequence variable
# aa_seq<- DNAString("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG",nchar=NA)
#Now the DNA sequence is being translated to an amino acid sequence
actualAA <- Biostrings ::translate(aa_seq, genetic.code = GENETIC_CODE, no.init.codon=FALSE, if.fuzzy.codon = "error")
# prints the amino acid sequence, if the previous code worked 
actualAA


# I am saving the amino acid sequence as a FASTA file, but first I am writing the amino acid as an Xstring set to be able to save the amino acid sequence as a fasta file.
# for (i in names(actualAA)) 
#   +   writeXStringSet(myAAStringSet[[i]], filepath = paste("/Users/arlinemartinez/Documents/GitHub/Bioinformatics/",i, ".fasta"), format = "fasta")
write.fasta(names="actualAA", sequences=actualAA, file.out="actualAA.fasta")

# After saving the amino acid sequence to the file, I ran the amino acid sequence in UniProt, and the accession number of the best match was A0A0J9YWK4.

# The disease associated with this gene is Beta Thalassemia

# Protein Image will be uploaded on Github


