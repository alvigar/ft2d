library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

require(data.table)

genes_id <- as.data.frame(fread("gtexgenesinnetworks.txt"))

genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"external_gene_name",'chromosome_name','start_position','end_position', 'transcript_count', 'gene_biotype'), filters =c("ensembl_gene_id"), values=genes_id, mart = ensembl)

genes

write.table(genes,"busquedaHGNC_SYMBOL.txt", sep=" \t ")

genes <- read.table("busquedaHGNC_SYMBOL.txt", stringsAsFactors=FALSE)

genes_completos <- genes[-which(genes$hgnc_symbol == ""), ]

biotypes <- getBM(attributes=c("external_gene_name", "gene_biotype"), filters =c("start", "end", "chromosome_name"), values=list(genes_completos[1,'start_position']-5000, genes_completos[1,'end_position']+5000, genes_completos[1,'chromosome_name']), mart = ensembl)

biotypes$gene_name <- genes_completos[1,'external_gene_name']
biotypes$length_gene <- as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
biotypes$biotype_gene <- genes_completos[1, 'gene_biotype']


A_M = matrix ( nrow=nrow(genes_completos), ncol=22 )
dimnames(A_M) <- list(rownames(A_M, do.NULL = FALSE, prefix = "gen"), c("external_gene_name", "length", "biotype", "protein_coding", "antisense", "processed_pseudogene", "transcribed_processed_pseudogene", "3prime_overlapping_ncRNA", "lincRNA", "retained_intron", "sense_intronic", "sense_overlapping", "miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA", "tRNA", "vaultRNA", "macro_lncRNA", "bidlncRNA"))

biotypes_M <- getBM(attributes=c("external_gene_name", "gene_biotype"), filters =c("start", "end", "chromosome_name"), values=list(genes_completos[1,'start_position']-500000, genes_completos[1,'end_position']+500000, genes_completos[1,'chromosome_name']), mart = ensembl)

biotypes_M$gene_name <- genes_completos[1,'external_gene_name']
biotypes_M$length_gene <- as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
biotypes_M$biotype_gene <- genes_completos[1, 'gene_biotype']

A_M[index, "external_gene_name"] = genes_completos[1,'external_gene_name']
A_M[index, "length"] = as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
A_M[index, "biotype"] = genes_completos[1, 'gene_biotype']


for(index in 2:nrow(genes_completos)){

  biotype_march_M <- getBM(attributes=c("external_gene_name", "gene_biotype"), filters =c("start", "end", "chromosome_name"), values=list(genes_completos[index,'start_position']-500000, genes_completos[index,'end_position']+500000, genes_completos[index,'chromosome_name']), mart = ensembl)
  
  biotype_march_M$gene_name <- genes_completos[index,'external_gene_name']
  biotype_march_M$length_gene <- as.integer(genes_completos[index, 'end_position']) - as.integer(genes_completos[index, 'start_position'])
  biotype_march_M$biotype_gene <- genes_completos[index, 'gene_biotype']

  biotypes_M <- rbind(biotypes_M, biotype_march_M)
  
  print(paste(index, nrow(genes_completos), sep=" - "))

  
}

write.table(biotypes,"biotypes.txt", sep=" \t ")
write.table(biotypes_M,"biotypes_M.txt", sep=" \t ")
write.table(A_M,"n_biotypes_M.txt", sep=" \t ")

nrow(biotypes_M)

biotypes <- read.table("biotypes.txt",stringsAsFactors=FALSE)
biotypes_M <- read.table("biotypes_M.txt",stringsAsFactors=FALSE)

A = matrix ( nrow=nrow(genes_completos), ncol=22 )
dimnames(A) <- list(rownames(A, do.NULL = FALSE, prefix = "gen"), c("external_gene_name", "length", "biotype", "protein_coding", "antisense", "processed_pseudogene", "transcribed_processed_pseudogene", "3prime_overlapping_ncRNA", "lincRNA", "retained_intron", "sense_intronic", "sense_overlapping", "miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA", "tRNA", "vaultRNA", "macro_lncRNA", "bidlncRNA"))

for(gen in 1:nrow(biotypes)){

  index <- 1
  for(k in 1:nrow(genes_completos)){
    if (!is.na(A[k, "external_gene_name"])) {
      if (A[k, "external_gene_name"] == biotypes[gen, 'gene_name']){
        index <- k
        break
      }
    }
    else {
      A[k, "external_gene_name"] = ""
      A[k, "external_gene_name"] = biotypes[gen, 'gene_name']
      A[k, "length"] = biotypes[gen, 'length_gene']
      A[k, "biotype"] = biotypes[gen, 'biotype_gene']
      index <- k
      break
    }
  }

  if (biotypes[gen, 2] == "protein_coding"){
    if (is.na(A[index, "protein_coding"])){
      A[index, "protein_coding"] <- 0
    }
    A[index, "protein_coding"] <- as.integer(A[index, "protein_coding"]) + 1
  }
  if (biotypes[gen, 2] == "antisense"){
    if (is.na(A[index, "antisense"])){
      A[index, "antisense"] <- 0
    }
    A[index, "antisense"] <- as.integer(A[index, "antisense"]) + 1
  }
  if (biotypes[gen, 2] == "processed_pseudogene"){
    if (is.na(A[index, "processed_pseudogene"])){
      A[index, "processed_pseudogene"] <- 0
    }
    A[index, "processed_pseudogene"] <- as.integer(A[index, "processed_pseudogene"]) + 1
  }
  if (biotypes[gen, 2] == "transcribed_processed_pseudogene"){
    if (is.na(A[index, "transcribed_processed_pseudogene"])){
      A[index, "transcribed_processed_pseudogene"] <- 0
    }
    A[index, "transcribed_processed_pseudogene"] <- as.integer(A[index, "transcribed_processed_pseudogene"]) + 1
  }
  if (biotypes[gen, 2] == "3prime_overlapping_ncRNA"){
    if (is.na(A[index, "3prime_overlapping_ncRNA"])){
      A[index, "3prime_overlapping_ncRNA"] <- 0
    }
    A[index, "3prime_overlapping_ncRNA"] <- as.integer(A[index, "3prime_overlapping_ncRNA"]) + 1
  }
  if (biotypes[gen, 2] == "lincRNA"){
    if (is.na(A[index, "lincRNA"])){
      A[index, "lincRNA"] <- 0
    }
    A[index, "lincRNA"] <- as.integer(A[index, "lincRNA"]) + 1
  }
  if (biotypes[gen, 2] == "retained_intron"){
    if (is.na(A[index, "retained_intron"])){
      A[index, "retained_intron"] <- 0
    }
    A[index, "retained_intron"] <- as.integer(A[index, "retained_intron"]) + 1
  }
  if (biotypes[gen, 2] == "sense_intronic"){
    if (is.na(A[index, "sense_intronic"])){
      A[index, "sense_intronic"] <- 0
    }
    A[index, "sense_intronic"] <- as.integer(A[index, "sense_intronic"]) + 1
  }
  if (biotypes[gen, 2] == "sense_overlapping"){
    if (is.na(A[index, "sense_overlapping"])){
      A[index, "sense_overlapping"] <- 0
    }
    A[index, "sense_overlapping"] <- as.integer(A[index, "sense_overlapping"]) + 1
  }
  if (biotypes[gen, 2] == "miRNA"){
    if (is.na(A[index, "protein_coding"])){
      A[index, "miRNA"] <- 0
    }
    A[index, "miRNA"] <- as.integer(A[index, "miRNA"]) + 1
  }
  if (biotypes[gen, 2] == "piRNA"){
    if (is.na(A[index, "piRNA"])){
      A[index, "piRNA"] <- 0
    }
    A[index, "piRNA"] <- as.integer(A[index, "piRNA"]) + 1
  }
  if (biotypes[gen, 2] == "rRNA"){
    if (is.na(A[index, "rRNA"])){
      A[index, "rRNA"] <- 0
    }
    A[index, "rRNA"] <- as.integer(A[index, "rRNA"]) + 1
  }
  if (biotypes[gen, 2] == "siRNA"){
    if (is.na(A[index, "siRNA"])){
      A[index, "siRNA"] <- 0
    }
    A[index, "siRNA"] <- as.integer(A[index, "siRNA"]) + 1
  }
  if (biotypes[gen, 2] == "snRNA"){
    if (is.na(A[index, "snRNA"])){
      A[index, "snRNA"] <- 0
    }
    A[index, "snRNA"] <- as.integer(A[index, "snRNA"]) + 1
  }
  if (biotypes[gen, 2] == "snoRNA"){
    if (is.na(A[index, "snoRNA"])){
      A[index, "snoRNA"] <- 0
    }
    A[index, "snoRNA"] <- as.integer(A[index, "snoRNA"]) + 1
  }
  if (biotypes[gen, 2] == "tRNA"){
    if (is.na(A[index, "tRNA"])){
      A[index, "tRNA"] <- 0
    }
    A[index, "tRNA"] <- as.integer(A[index, "tRNA"]) + 1
  }
  if (biotypes[gen, 2] == "vaultRNA"){
    if (is.na(A[index, "vaultRNA"])){
      A[index, "vaultRNA"] <- 0
    }
    A[index, "vaultRNA"] <- as.integer(A[index, "vaultRNA"]) + 1
  }
  if (biotypes[gen, 2] == "macro_lncRNA"){
    if (is.na(A[index, "macro_lncRNA"])){
      A[index, "macro_lncRNA"] <- 0
    }
    A[index, "macro_lncRNA"] <- as.integer(A[index, "macro_lncRNA"]) + 1
  }
  if (biotypes[gen, 2] == "bidlncRNA"){
    if (is.na(A[index, "bidlncRNA"])){
      A[index, "bidlncRNA"] <- 0
    }
    A[index, "bidlncRNA"] <- as.integer(A[index, "bidlncRNA"]) + 1
  }

}

write.table(A, "n_biotypes.txt", sep=" \t ")

for(gen in 1:nrow(biotypes_M)){

  index <- 1
  for(k in 1:nrow(genes_completos)){
    if (!is.na(A_M[k, "external_gene_name"])) {
      if (A_M[k, "external_gene_name"] == biotypes_M[gen, 1]){
        index <- k
        break
      }
    }
    else {
      A_M[k, "external_gene_name"] = ""
      A_M[k, "external_gene_name"] = biotypes_M[gen, 1]
      A_M[k, "length"] = biotypes_M[gen, 'length_gene']
      A_M[k, "biotype"] = biotypes_M[gen, 'biotype_gene']
      for(o in 4:22){
        A_M[k,o] <- 0
      }
      index <- k
      break
    }
  }
  
  if (biotypes_M[gen, 2] == "protein_coding"){
    A_M[index, "protein_coding"] <- as.integer(A_M[index, "protein_coding"]) + 1
  }
  else if (biotypes_M[gen, 2] == "antisense"){
    A_M[index, "antisense"] <- as.integer(A_M[index, "antisense"]) + 1
  }
  else if (biotypes_M[gen, 2] == "processed_pseudogene"){
    A_M[index, "processed_pseudogene"] <- as.integer(A_M[index, "processed_pseudogene"]) + 1
  }
  else if (biotypes_M[gen, 2] == "transcribed_processed_pseudogene"){
    A_M[index, "transcribed_processed_pseudogene"] <- as.integer(A_M[index, "transcribed_processed_pseudogene"]) + 1
  }
  else if (biotypes_M[gen, 2] == "3prime_overlapping_ncRNA"){
    A_M[index, "3prime_overlapping_ncRNA"] <- as.integer(A_M[index, "3prime_overlapping_ncRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "lincRNA"){
    A_M[index, "lincRNA"] <- as.integer(A_M[index, "lincRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "retained_intron"){
    A_M[index, "retained_intron"] <- as.integer(A_M[index, "retained_intron"]) + 1
  }
  else if (biotypes_M[gen, 2] == "sense_intronic"){
    A_M[index, "sense_intronic"] <- as.integer(A_M[index, "sense_intronic"]) + 1
  }
  else if (biotypes_M[gen, 2] == "sense_overlapping"){
    A_M[index, "sense_overlapping"] <- as.integer(A_M[index, "sense_overlapping"]) + 1
  }
  else if (biotypes_M[gen, 2] == "miRNA"){
    A_M[index, "miRNA"] <- as.integer(A_M[index, "miRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "piRNA"){
    A_M[index, "piRNA"] <- as.integer(A_M[index, "piRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "rRNA"){
    A_M[index, "rRNA"] <- as.integer(A_M[index, "rRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "siRNA"){
    A_M[index, "siRNA"] <- as.integer(A_M[index, "siRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "snRNA"){
    A_M[index, "snRNA"] <- as.integer(A_M[index, "snRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "snoRNA"){
    A_M[index, "snoRNA"] <- as.integer(A_M[index, "snoRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "tRNA"){
    A_M[index, "tRNA"] <- as.integer(A_M[index, "tRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "vaultRNA"){
    A_M[index, "vaultRNA"] <- as.integer(A_M[index, "vaultRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "macro_lncRNA"){
    A_M[index, "macro_lncRNA"] <- as.integer(A_M[index, "macro_lncRNA"]) + 1
  }
  else if (biotypes_M[gen, 2] == "bidlncRNA"){
    A_M[index, "bidlncRNA"] <- as.integer(A_M[index, "bidlncRNA"]) + 1
  }

  print(paste(gen, nrow(biotypes_M), sep=" -"))
  
}

write.table(A_M, "n_biotypes_M.txt", sep=" \t ")

A_M <- read.table("n_biotypes_M.txt",stringsAsFactors=FALSE)

A <- read.table("n_biotypes.txt",stringsAsFactors=FALSE)


coeficientes_w = c(0.4, 0.2, 0.15, 0.05, 0.05, 0.05, 0.05, 0.05)

columnas = c("protein_coding", "antisense", "lincRNA", "retained_intron", "sense_intronic", "sense_overlapping", "macro_lncRNA", "bidlncRNA")

complex_march = matrix (nrow=nrow(A), ncol=4)
dimnames(complex_march) <- list(rownames(complex_march, do.NULL = FALSE, prefix = "gen"), c('external_gene_name', 'complexity', 'biotype', 'length'))

A <- A[-which(A$external_gene_name == '0'),]

A

for(gen_complex in 1:nrow(A)){

  complex_march[gen_complex, 1] = A[gen_complex, "external_gene_name"]
  complex_march[gen_complex, 3] = A[gen_complex, "biotype"]
  complex_march[gen_complex, 4] = A[gen_complex, "length"]
  aux = 0
  for(coef in 1:8){
    aux = aux + A[gen_complex, columnas[coef]]*coeficientes_w[coef]
  }
  complex_march[gen_complex, 2] <- aux
}

write.table(complex_march, "complex_A.txt", sep=" \t ")

library(corrplot)

cgenes = read.delim("controlForSonia.tsv",stringsAsFactors=F,header=F)$V1

comp_a = cbind(A,disease=rep(TRUE,nrow(A)))
comp_a$disease[match(cgenes,comp_a$external_gene_name)] = FALSE
comp_a_2 = comp_a[,unlist(apply(comp_a,2,function(x){ length(unique(x)) != 1}))]

corrplot(cor(comp_a_2[,-c(1,2,3,22)]))

comp_a

complex_march_M = matrix (nrow=nrow(genes_completos), ncol=2)

for(gen_complex in 1:nrow(genes_completos)){

  complex_march_M[gen_complex, 1] = A_M[gen_complex, "external_gene_name"]
  complex_march_M[gen_complex, 2] <- 0
  for(coef in 1:8){
    if (!is.na(A_M[gen_complex, columnas[coef]])){
      complex_march_M[gen_complex, 2] <- as.integer(complex_march_M[gen_complex, 2]) + as.integer(A_M[gen_complex, columnas[coef]])*coeficientes_w[coef]
    }
  }
}

fit = glm(formula = disease ~ ., data=comp_a[,c(-1,-2,-3)],family=binomial())
summary(fit)