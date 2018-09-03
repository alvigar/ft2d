#------------------------------------------------------------------------#

##Author: Alfonso Vidal García
##Co-Author: Juan Antonio Botía Blaya

#------------------------------------------------------------------------#

A_m <- read.delim("n_biotype_mod.txt", stringsAsFactors = FALSE) #Con radio (-5000, +5000)
A_M <- read.delim("n_biotype_mod_M.txt", stringsAsFactors = FALSE) #Con radio (-500000, +500000)

A_m <- A_m[which(A_m$length > 0),]
sumas_m <- colSums (A_m[4:17], na.rm = FALSE, dims = 1)
A_m <- A_m[c('external_gene_name', 'length', 'biotype',names(sumas_m[which(sumas_m > 0)]))]

A_M <- A_M[which(A_M$length > 0),]
sumas_M <- colSums (A_M[4:17], na.rm = FALSE, dims = 1)
A_M <- A_M[c('external_gene_name', 'length', 'biotype',names(sumas_M[which(sumas_M > 0)]))]

discapacidad <- read.delim("DI_Oct2016.csv", stringsAsFactors = FALSE, header = FALSE)$V1
control <- read.delim("controlForSonia.tsv", stringsAsFactors = FALSE, header = FALSE)$V1

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
    
    write.table(A, paste('biotypes_tissues/biotypes_', gsub('.rds', '', tejidos[tej]), '_m.txt', sep=''), sep="\t")
    
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
    
    write.table(A, paste('biotypes_tissues/biotypes_', gsub('.rds', '', tejidos[tej]), '_M.txt', sep=''), sep="\t")
  }
}
