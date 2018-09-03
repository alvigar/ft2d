#------------------------------------------------------------------------#

##Author: Alfonso Vidal García
##Co-Author: Juan Antonio Botía Blaya

#------------------------------------------------------------------------#

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

require(data.table)

genes_id <- as.data.frame(fread("gtexgenesinnetworks.txt"))

genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"external_gene_name",
                            'chromosome_name','start_position','end_position', 'transcript_count', 
                            'gene_biotype'), filters =c("ensembl_gene_id"), values=genes_id, mart = ensembl)

genes_completos <- genes[-which(genes$hgnc_symbol == ""), ]

A_m = matrix ( nrow=nrow(genes_completos), ncol=22 ) #Radius 5000
dimnames(A_m) <- list(rownames(A_m, do.NULL = FALSE, prefix = "gen"), 
                      c("external_gene_name", "length", "biotype", "protein_coding", "antisense", 
                        "processed_pseudogene", "transcribed_processed_pseudogene", 
                        "3prime_overlapping_ncRNA", "lincRNA", "retained_intron", "sense_intronic", 
                        "sense_overlapping", "miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA", 
                        "tRNA", "vaultRNA", "macro_lncRNA", "bidlncRNA"))

A_M = matrix ( nrow=nrow(genes_completos), ncol=22 ) #Radius 500000
dimnames(A_M) <- list(rownames(A_M, do.NULL = FALSE, prefix = "gen"), 
                      c("external_gene_name", "length", "biotype", "protein_coding", "antisense", 
                        "processed_pseudogene", "transcribed_processed_pseudogene", 
                        "3prime_overlapping_ncRNA", "lincRNA", "retained_intron", "sense_intronic", 
                        "sense_overlapping", "miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA",
                        "tRNA", "vaultRNA", "macro_lncRNA", "bidlncRNA"))


#Radius 5000
biotypes_m <- getBM(attributes=c("external_gene_name", "gene_biotype"), filters =c("start", "end", "chromosome_name"), 
                    values=list(genes_completos[1,'start_position']-5000, genes_completos[1,'end_position']+5000, genes_completos[1,'chromosome_name']), 
                    mart = ensembl)

#Radius 500000
biotypes_M <- getBM(attributes=c("external_gene_name", "gene_biotype"), filters =c("start", "end", "chromosome_name"),
                    alues=list(genes_completos[1,'start_position']-500000, genes_completos[1,'end_position']+500000, genes_completos[1,'chromosome_name']),
                    mart = ensembl)


biotypes_m$gene_name <- genes_completos[1,'external_gene_name']
biotypes_m$length_gene <- as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
biotypes_m$biotype_gene <- genes_completos[1, 'gene_biotype']

biotypes_M$gene_name <- genes_completos[1,'external_gene_name']
biotypes_M$length_gene <- as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
biotypes_M$biotype_gene <- genes_completos[1, 'gene_biotype']

A_m[index, "external_gene_name"] = genes_completos[1,'external_gene_name']
A_m[index, "length"] = as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
A_m[index, "biotype"] = genes_completos[1, 'gene_biotype']

A_M[index, "external_gene_name"] = genes_completos[1,'external_gene_name']
A_M[index, "length"] = as.integer(genes_completos[1, 'end_position']) - as.integer(genes_completos[1, 'start_position'])
A_M[index, "biotype"] = genes_completos[1, 'gene_biotype']

for(index in 2:nrow(genes_completos)){
  biotype_march_m <- getBM(attributes=c("external_gene_name", "gene_biotype"), 
                           filters=c("start", "end", "chromosome_name"), 
                           values=list(genes_completos[index,'start_position']-5000, 
                                       genes_completos[index,'end_position']+5000, 
                                       genes_completos[index,'chromosome_name']),
                           mart = ensembl)
  
  biotype_march_m$gene_name <- genes_completos[index,'external_gene_name']
  biotype_march_m$length_gene <- as.integer(genes_completos[index, 'end_position']) - as.integer(genes_completos[index, 'start_position'])
  biotype_march_m$biotype_gene <- genes_completos[index, 'gene_biotype']
  
  biotypes_m <- rbind(biotypes_m, biotype_march_m)
  
  biotype_march_M <- getBM(attributes=c("external_gene_name", "gene_biotype"), 
                           filters =c("start", "end", "chromosome_name"), 
                           values=list(genes_completos[index,'start_position']-500000, 
                                       genes_completos[index,'end_position']+500000, 
                                       genes_completos[index,'chromosome_name']), 
                           mart = ensembl)
  
  biotype_march_M$gene_name <- genes_completos[index,'external_gene_name']
  biotype_march_M$length_gene <- as.integer(genes_completos[index, 'end_position']) - as.integer(genes_completos[index, 'start_position'])
  biotype_march_M$biotype_gene <- genes_completos[index, 'gene_biotype']
  
  biotypes_M <- rbind(biotypes_M, biotype_march_M)
}

## To this point we have all the transcription for every gen.

for(gen in 1:nrow(biotypes_m)){
  
  index <- 1
  for(k in 1:nrow(genes_completos)){
    if (!is.na(A_m[k, "external_gene_name"])) {
      if (A_m[k, "external_gene_name"] == biotypes_m[gen, 1]){
        index <- k
        break
      }
    }
    else {
      A_m[k, "external_gene_name"] = ""
      A_m[k, "external_gene_name"] = biotypes_m[gen, 1]
      A_m[k, "length"] = biotypes_m[gen, 'length_gene']
      A_m[k, "biotype"] = biotypes_m[gen, 'biotype_gene']
      for(o in 4:22){
        A_m[k,o] <- 0
      }
      index <- k
      break
    }
  }
  
  if (biotypes_m[gen, 2] == "protein_coding"){
    A_m[index, "protein_coding"] <- as.integer(A_m[index, "protein_coding"]) + 1
  }
  else if (biotypes_m[gen, 2] == "antisense"){
    A_m[index, "antisense"] <- as.integer(A_m[index, "antisense"]) + 1
  }
  else if (biotypes_m[gen, 2] == "processed_pseudogene"){
    A_m[index, "processed_pseudogene"] <- as.integer(A_m[index, "processed_pseudogene"]) + 1
  }
  else if (biotypes_m[gen, 2] == "transcribed_processed_pseudogene"){
    A_m[index, "transcribed_processed_pseudogene"] <- as.integer(A_m[index, "transcribed_processed_pseudogene"]) + 1
  }
  else if (biotypes_m[gen, 2] == "3prime_overlapping_ncRNA"){
    A_m[index, "3prime_overlapping_ncRNA"] <- as.integer(A_m[index, "3prime_overlapping_ncRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "lincRNA"){
    A_m[index, "lincRNA"] <- as.integer(A_m[index, "lincRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "retained_intron"){
    A_m[index, "retained_intron"] <- as.integer(A_m[index, "retained_intron"]) + 1
  }
  else if (biotypes_m[gen, 2] == "sense_intronic"){
    A_m[index, "sense_intronic"] <- as.integer(A_m[index, "sense_intronic"]) + 1
  }
  else if (biotypes_m[gen, 2] == "sense_overlapping"){
    A_m[index, "sense_overlapping"] <- as.integer(A_m[index, "sense_overlapping"]) + 1
  }
  else if (biotypes_m[gen, 2] == "miRNA"){
    A_m[index, "miRNA"] <- as.integer(A_m[index, "miRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "piRNA"){
    A_m[index, "piRNA"] <- as.integer(A_m[index, "piRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "rRNA"){
    A_m[index, "rRNA"] <- as.integer(A_m[index, "rRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "siRNA"){
    A_m[index, "siRNA"] <- as.integer(A_m[index, "siRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "snRNA"){
    A_m[index, "snRNA"] <- as.integer(A_m[index, "snRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "snoRNA"){
    A_m[index, "snoRNA"] <- as.integer(A_m[index, "snoRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "tRNA"){
    A_m[index, "tRNA"] <- as.integer(A_m[index, "tRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "vaultRNA"){
    A_m[index, "vaultRNA"] <- as.integer(A_m[index, "vaultRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "macro_lncRNA"){
    A_m[index, "macro_lncRNA"] <- as.integer(A_m[index, "macro_lncRNA"]) + 1
  }
  else if (biotypes_m[gen, 2] == "bidlncRNA"){
    A_m[index, "bidlncRNA"] <- as.integer(A_m[index, "bidlncRNA"]) + 1
  }
}

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
}

write.table(A_m, "n_biotype_mod.txt", sep="\t")
write.table(A_M, "n_biotype_mod_M.txt", sep="\t")