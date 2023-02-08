library(Biostrings)
library(AnnotationHub)
hub <- AnnotationHub(cache="./annotationhub")
query(hub, c("homo sapiens","EnsDb"))
query(hub, c("Homo sapiens", "release-94"))

ens94_human_dna <- hub[["AH65745"]]
# Ensembl release 94
txdbENS <- hub[["AH64923"]]

all.genes <- genes(txdbENS)
coding_genes <- all.genes[all.genes$gene_biotype == "protein_coding"]
coding_genes <- keepSeqlevels(coding_genes, c(seq(1:22), 'MT', 'X', 'Y'), pruning.mode="coarse")

# chr:start_flanking|start-end|end_flanking:width:geneid(s)
headers = apply(as.data.frame(coding_genes), 1, function(row) paste(
  paste('chr', row[1], sep=''),
  paste(
    paste(row[2], row[2] , sep='|'),
    paste(row[2], row[3], sep='|'),
    sep='-'),
  row[5],
  paste(row[6], collapse = ','),
  sep=':')
)

#extract the sequences from the genome 
seqs <- getSeq(ens94_human_dna, coding_genes)

# get sequences in 3 frame translation
seqs_3frames <- lapply(1:3, function(pos) subseq(seqs, start=pos))

protein_seqs = lapply(seqs_3frames, function(seq) translate(seq, no.init.codon = T, if.fuzzy.codon = 'solve'))
seqnumber = length(headers)

# writing the output
writefun <- function(header, orf, sequence, counter, output_file){
  line = paste('>seq', counter, ':', header, ':F', orf, '\n', sep='')
  line = paste(line, sequence, sep='')
  write(line, file=output_file, append=TRUE)
}

output_file <- "3FTgenes_coding.fasta"

mapply(function(x,y,z,w) writefun(x, y, z, w, output_file),
       list(headers,headers,headers),
       list(rep(1,seqnumber), rep(2,seqnumber), rep(3,seqnumber)),
       protein_seqs,
       list(seq(1,seqnumber), seq(seqnumber + 1 ,(seqnumber*2)), seq((seqnumber*2) + 1 ,(seqnumber*3))))