mySequenceFile1 <- system.file("sequence_1", "/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_1.fasta", package="msa")
mySequences1 <- readDNAStringSet("/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_1.fasta", format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mySequences1

mySequenceFile2 <- system.file("sequence_2", "/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_2.fasta", package="msa")
mySequences2 <- readDNAStringSet("/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_2.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities =FALSE)
mySequences2

mySequenceFile3 <- system.file("sequence_3", "/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_3.txt", package="msa")
mySequences3 <- readDNAStringSet("/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_3.txt", format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mySequences3

mySequenceFile4 <- system.file("sequence_4", "/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_4.fasta", package="msa")
mySequences4 <- readDNAStringSet("/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_4.fasta", format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mySequences4

mySequenceFile5 <- system.file("sequence_5", "/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_5.fasta", package="msa")
mySequences5 <- readDNAStringSet("/Users/arlinemartinez/Documents/GitHub/Bioinformatics/sequence_5.fasta", format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE, with.qualities = FALSE)
mySequences5

#mysequences <-c(sequence_1, sequence_2, sequence_3, sequence_4, sequence_5)
# exists('sequence_1')
# exists('sequence_2'

mySequences <- c(mySequences1, mySequences2, mySequences3, mySequences4, mySequences5)
mySequences

names(mySequences) <- c("GS262331.1", "DN626083.1", "lcl|HQ003260.1", "DN626085.1", "DN626081.1")

library(msa)

#firstAlignment <- msa(mySequences)
#firstAlignment
#print(firstAlignment, show="complete")

myClustalWAlignment <- msa(mySequences, "ClustalW")
myClustalWAlignment

#print(firstAlignment, show="complete")
print(myClustalWAlignment, show="complete")

length(mySequences)
width(mySequences)
nchar(mySequences)

letterFrequency(mySequences[[1]], letters="CG", OR=0)
letterFrequency(mySequences[[2]],letters="CG", OR=0)
letterFrequency(mySequences[[3]],letters="CG", OR=0)
letterFrequency(mySequences[[4]],letters="CG", OR=0)
letterFrequency(mySequences[[5]],letters="CG", OR=0)

myAlignment2<- msaConvert(myClustalWAlignment, type="seqinr::alignment")
d<-dist.alignment(myAlignment2, matrix="identity")
d
rownames(myClustalWAlignment)



#phylo_Data <- msaConvert(myClustalWAlignment, type="phangorn::phyDat") 
#write.phyDat(phylo_Data, "myClustalWAlignment", format = "fasta") 

library("Biostrings")
myseq <- readDNAStringSet("DN626081.1", format="fasta")
head(myseq)

mySequences5 <- Biostrings ::translate(mySequences5, genetic.code = GENETIC_CODE, no.init.codon=FALSE, if.fuzzy.codon = "error")
#for(n )

Alignment_phyDat <- msaConvert(myClustalWAlignment, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")
