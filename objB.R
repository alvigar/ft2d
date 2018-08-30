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

## To this point we have the two tables of the species of RNA. Now we select the elements that have protein_coding

comp_a = cbind(A_m,disease=rep(0,nrow(A_m)))
comp_a_M = cbind(A_M,disease=rep(0,nrow(A_M)))

comp_a$disease[match(discapacidad,trimws(comp_a$external_gene_name))] = 1
comp_a$disease[match(control,trimws(comp_a$external_gene_name))] = 2

comp_a_M$disease[match(discapacidad,trimws(comp_a_M$external_gene_name))] = 1
comp_a_M$disease[match(control,trimws(comp_a_M$external_gene_name))] = 2

comp_a <- comp_a[-which(comp_a$disease == 0), ]
comp_a$disease[which(comp_a$disease == 2)] = 0

comp_a_M <- comp_a_M[-which(comp_a_M$disease == 0), ]
comp_a_M$disease[which(comp_a_M$disease == 2)] = 0

comp_a[c(2,4:ncol(A_m))] <- scale(comp_a[c(2,4:ncol(A_m))])
comp_a_M[c(2,4:ncol(A_M))] <- scale(comp_a_M[c(2,4:ncol(A_M))])

discapacidad <- read.delim("DI_Oct2016.csv", stringsAsFactors = FALSE, header = FALSE)$V1
control <- read.delim("controlForSonia.tsv", stringsAsFactors = FALSE, header = FALSE)$V1

protein <- subset(comp_a, biotype == 'protein_coding' )

protein_disease <- subset(protein, disease == 1)

vector <- c(protein_disease['external_gene_name'])$external_gene_name

for(i in 1:length(vector)){
  vector[i] <- trimws(vector[i])
}

genes_disease <- intersect(discapacidad, vector)

protein_M <- subset(comp_a_M, biotype == 'protein_coding' )

protein_disease_M <- subset(protein_M, disease == 1)

vector_M <- c(protein_disease['external_gene_name'])$external_gene_name

for(i in 1:length(vector_M)){
  vector_M[i] <- trimws(vector_M[i])
}

genes_disease_M <- intersect(discapacidad, vector_M)

## Beyond this point start the split of genes in tissues

A_m <- read.table("n_biotypes.txt",stringsAsFactors=FALSE)
A_M <- read.table("n_biotypes_M.txt",stringsAsFactors=FALSE)

A_m <- A_m[which(A_m$length > 0),]
sumas_m <- colSums (A_m[4:17], na.rm = FALSE, dims = 1)
A_m <- A_m[c('external_gene_name', 'length', 'biotype',names(sumas_m[which(sumas_m > 0)]))]

A_M <- A_M[which(A_M$length > 0),]
sumas_M <- colSums (A_M[4:17], na.rm = FALSE, dims = 1)
A_M <- A_M[c('external_gene_name', 'length', 'biotype',names(sumas_M[which(sumas_M > 0)]))]

tejidos <- list.files(path="data/data/.")

for(tej in 1:length(tejidos)){
    
    fichero <- paste('data/data/', tejidos[tej], sep='')
    
    tabla_individuos <- readRDS(fichero)
    
    if(ncol(tabla_individuos) > 0){
      
      expressed.genes = rowSums(tabla_individuos > 0.1)  > (0.2 * ncol(tabla_individuos))
      
      genes = tabla_individuos[expressed.genes,]
      
      genes <- unlist(strsplit(c(rownames(genes)), "\\."))[c(TRUE, FALSE)]
      
      genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"external_gene_name",'chromosome_name','start_position','end_position', 'transcript_count', 'gene_biotype'), filters =c("ensembl_gene_id"), values=genes, mart = ensembl)
      
      A <- A_m #Con radio (-5000, +5000)
      A <- A[which(A$length > 0),]
      
      for(i in 1:nrow(A)){
        if (length(A[i, 'external_gene_name']) > 0){
          if (trimws(A[i, 'external_gene_name']) %in% c(genes['external_gene_name'])$external_gene_name & trimws(A[i, 'external_gene_name']) %in% genes_disease){
            genes_area_gen <- biotype[which(biotype$gene_name == trimws(A[i, 'external_gene_name'])), ]
            for(j in 1:nrow(genes_area_gen)){
              if (!(trimws(genes_area_gen[j, 'external_gene_name']) %in% c(genes['external_gene_name'])$external_gene_name)){
                cantidad <- as.integer(A[i, trimws(genes_area_gen[j, 'gene_biotype'])]) - 1
                if (length(cantidad) > 0){
                  A[i, trimws(genes_area_gen[j, 'gene_biotype'])] <- cantidad
                }
                else {
                  A[i, trimws(genes_area_gen[j, 'gene_biotype'])] <- 0
                }
              }
            }
          }
          else{
            A <- A[-i,]
          }
        }
      }
      
      write.table(A, paste('biotypes_', gsub('.rds', '', tejidos[tej]), '_m.txt', sep=''), sep="\t")
      
      A <- A_M #Con radio (-500000, +500000)
      A <- A[which(A$length > 0),]
      
      for(i in 1:nrow(A)){
        if (length(A[i, 'external_gene_name']) > 0){
          if (trimws(A[i, 'external_gene_name']) %in% c(genes['external_gene_name'])$external_gene_name & trimws(A[i, 'external_gene_name']) %in% genes_disease_M){
            genes_area_gen <- biotype[which(biotype$gene_name == trimws(A[i, 'external_gene_name'])), ]
            for(j in 1:nrow(genes_area_gen)){
              if (!(trimws(genes_area_gen[j, 'external_gene_name']) %in% c(genes['external_gene_name'])$external_gene_name)){
                cantidad <- as.integer(A[i, trimws(genes_area_gen[j, 'gene_biotype'])]) - 1
                if (length(cantidad) > 0){
                  A[i, trimws(genes_area_gen[j, 'gene_biotype'])] <- cantidad
                }
                else {
                  A[i, trimws(genes_area_gen[j, 'gene_biotype'])] <- 0
                }
              }
            }
          }
          else{
            A <- A[-i,]
          }
        }
      }
      
      write.table(A, paste('biotypes_', gsub('.rds', '', tejidos[tej]), '_M.txt', sep=''), sep="\t")
    }
}

types_m <- names(A_m)[-c(1,2,3)]
types_M <- names(A_M)[-c(1,2,3)]

tejidos_m = matrix ( nrow=length(tejidos), ncol=(length(types_m)+1) )
dimnames(tejidos_m) <- list(tejidos, c("genes_count", types_m))
tejidos_M = matrix ( nrow=length(tejidos), ncol=length(types_M)+1 )
dimnames(tejidos_M) <- list(tejidos, c("genes_count", types_M))

for(tej in 1:length(tejidos)){
  
  tejido_matrix_m <- read.table(paste('biotypes_', gsub('.rds', '', tejidos[tej]), '_m.txt', sep=''), stringsAsFactors = FALSE)
  sumas_m <- colSums (tejido_matrix_m[4:ncol(tejido_matrix_m)], na.rm = FALSE, dims = 1)
  tejido_matrix_M <- read.table(paste('biotypes_', gsub('.rds', '', tejidos[tej]), '_M.txt', sep=''), stringsAsFactors = FALSE)
  sumas_M <- colSums (tejido_matrix_M[4:ncol(tejido_matrix_M)], na.rm = FALSE, dims = 1)
  
  tejidos_m[tejidos[tej], 'genes_count'] <- nrow(tejido_matrix_m)
  tejidos_M[tejidos[tej], 'genes_count'] <- nrow(tejido_matrix_M)
  
  for(i in 1:length(types_m)){
    tejidos_m[tejidos[tej], types_m[i]] <- unname(sumas_m[types_m[i]])
  }
  
  for(i in 1:length(types_M)){
    tejidos_M[tejidos[tej], types_M[i]] <- unname(sumas_M[types_M[i]])
  }
  
}

library(xtable)

print(xtable(tejidos_m, type="latex"), file="tejidos_m.tex")
print(xtable(tejidos_M, type="latex"), file="tejidos_MM.tex")

## To this point we have all the diferents all the RNA species splited by tissues of Project GTEx.

