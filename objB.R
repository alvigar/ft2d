#------------------------------------------------------------------------#

##Author: Alfonso Vidal García
##Co-Author: Juan Antonio Botía Blaya

#------------------------------------------------------------------------#


library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

require(data.table)

A_m <- read.delim("n_biotype_mod.txt", stringsAsFactors = FALSE) #Con radio (-5000, +5000)
A_M <- read.delim("n_biotype_mod_M.txt", stringsAsFactors = FALSE) #Con radio (-500000, +500000)

## To this point we have the two tables of the species of RNA. Now we select the elements that have protein_coding

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

## Beyond this point start the split of genes in tissues

tejidos <- list.files(path="data/data/.")

types_m <- names(A_m)[-c(1,2,3)]
types_M <- names(A_M)[-c(1,2,3)]

tejidos_m = matrix ( nrow=length(tejidos), ncol=(length(types_m)+1) )
dimnames(tejidos_m) <- list(tejidos, c("genes_count", types_m))
tejidos_M = matrix ( nrow=length(tejidos), ncol=length(types_M)+1 )
dimnames(tejidos_M) <- list(tejidos, c("genes_count", types_M))

for(tej in 1:length(tejidos)){
  
  tejido_matrix_m <- read.table(paste('biotypes_tissues/biotypes_', gsub('.rds', '', tejidos[tej]), '_m.txt', sep=''), stringsAsFactors = FALSE)
  sumas_m <- colSums (tejido_matrix_m[4:ncol(tejido_matrix_m)], na.rm = FALSE, dims = 1)
  tejido_matrix_M <- read.table(paste('biotypes_tissues/biotypes_', gsub('.rds', '', tejidos[tej]), '_M.txt', sep=''), stringsAsFactors = FALSE)
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

## To this point we have all the diferents all the RNA species splited by tissues of Project GTEx. Now we make de analysis of the data for each tissue

library(InformationValue)

for(tej in 1:length(tejidos)){
  
  
  A_m <- read.table(paste('biotypes_', gsub('.rds', '', tejidos[tej]), '.txt', sep=''), stringsAsFactors = FALSE)
  A_M <- read.table(paste('biotypes_', gsub('.rds', '', tejidos[tej]), '_M.txt', sep=''), stringsAsFactors = FALSE)
  
  sumas <- colSums (A_m[4:ncol(A_m)], na.rm = FALSE, dims = 1)
  A_m <- A_m[c('external_gene_name', 'length', 'biotype',names(sumas[which(sumas > 0)]))]
  
  comp_A_tejido_m = cbind(A_m,disease=rep(0,nrow(A_m)))
  
  comp_A_tejido_m$disease[match(discapacidad,trimws(comp_A_tejido_m$external_gene_name))] = 1
  comp_A_tejido_m$disease[match(control,trimws(comp_A_tejido_m$external_gene_name))] = 2
  
  
  comp_A_tejido_m <- comp_A_tejido_m[-which(comp_A_tejido_m$disease == 0), ]
  comp_A_tejido_m$disease[which(comp_A_tejido_m$disease == 2)] = 0
  
  sumas <- colSums (A_M[4:ncol(A_M)], na.rm = FALSE, dims = 1)
  A_M <- A_M[c('external_gene_name', 'length', 'biotype',names(sumas[which(sumas > 0)]))]
  A_M <- A_M[-ncol(A_M)]
  
  comp_A_tejido_M = cbind(A_M,disease=rep(0,nrow(A_M)))
  
  comp_A_tejido_M$disease[match(discapacidad,trimws(comp_A_tejido_M$external_gene_name))] = 1
  comp_A_tejido_M$disease[match(control,trimws(comp_A_tejido_M$external_gene_name))] = 2
  
  
  comp_A_tejido_M <- comp_A_tejido_M[-which(comp_A_tejido_M$disease == 0), ]
  comp_A_tejido_M$disease[which(comp_A_tejido_M$disease == 2)] = 0
  
  comp_A_tejido_m[c(2,4:ncol(A_m))] <- scale(comp_A_tejido_m[c(2,4:ncol(A_m))])
  comp_A_tejido_M[c(2,4:ncol(A_M))] <- scale(comp_A_tejido_M[c(2,4:ncol(A_M))])
  
  comp_a_disease_m <- comp_A_tejido_m[which(comp_A_tejido_m[,'disease'] == 1),]
  comp_a_no_disease_m <- comp_A_tejido_m[which(comp_A_tejido_m[,'disease'] == 0),]
  set.seed(100)
  
  comp_a_disease_training_rows <- sample(1:nrow(comp_a_disease_m), 0.7*nrow(comp_a_disease_m))
  comp_a_no_disease_training_rows <- sample(1:nrow(comp_a_no_disease_m), 0.7*nrow(comp_a_no_disease_m))
  
  training_ones <- comp_a_disease_m[comp_a_disease_training_rows, ]
  training_zeros <- comp_a_no_disease_m[comp_a_no_disease_training_rows, ]
  
  trainingData <- rbind(training_ones, training_zeros)
  
  test_ones <- comp_a_disease_m[-comp_a_disease_training_rows, ]
  test_zeros <- comp_a_no_disease_m[-comp_a_no_disease_training_rows, ]
  
  testData_m <- rbind(test_ones, test_zeros)
  
  formula <- as.formula(disease ~ length)
  
  formula <- as.formula(
    paste(
      paste(deparse(formula), collapse=""), 
      paste(colnames(comp_A_tejido_m[4:(ncol(A_m)-1)]), collapse="+"),
      sep="+"
    )
  )
  
  logitMod_m <- glm(formula, data=trainingData, family=binomial(link="logit"))
  
  predicted_m <- plogis(predict(logitMod_m, testData_m)) 
  
  optCutOff_m <- optimalCutoff(testData_m$disease, predicted_m)[1] 
  
  comp_a_disease_M <- comp_A_tejido_M[which(comp_A_tejido_M[,'disease'] == 1),]
  comp_a_no_disease_M <- comp_A_tejido_M[which(comp_A_tejido_M[,'disease'] == 0),]
  set.seed(100)
  
  comp_a_disease_training_rows <- sample(1:nrow(comp_a_disease_M), 0.7*nrow(comp_a_disease_M))
  comp_a_no_disease_training_rows <- sample(1:nrow(comp_a_no_disease_M), 0.7*nrow(comp_a_no_disease_M))
  
  training_ones <- comp_a_disease_M[comp_a_disease_training_rows, ]
  training_zeros <- comp_a_no_disease_M[comp_a_no_disease_training_rows, ]
  
  trainingData_M <- rbind(training_ones, training_zeros)
  
  
  test_ones <- comp_a_disease_M[-comp_a_disease_training_rows, ]
  test_zeros <- comp_a_no_disease_M[-comp_a_no_disease_training_rows, ]
  
  testData_M <- rbind(test_ones, test_zeros)
  
  formula <- as.formula(disease ~ length)
  
  formula <- as.formula(
    paste(
      paste(deparse(formula), collapse=""), 
      paste(colnames(comp_A_tejido_M[4:(ncol(A_M))]), collapse="+"),
      sep="+"
    )
  )
  
  for(k in 1:nrow(trainingData_M)){
    for(t in 1:ncol(trainingData_M)){
      if (is.nan(trainingData_M[k, t])){
        trainingData_M[k, t] <- 0.0
      }
    }
  }
  
  logitMod_M <- glm(formula, data=trainingData_M, family=binomial(link="logit"))
  
  predicted_M <- plogis(predict(logitMod_M, testData_M)) 
  
  optCutOff_M <- optimalCutoff(testData_M$disease, predicted_M)[1]
  
  
  if (tej == 1){
    var_m <- names(summary(logitMod_m)$coefficients[,4])[-1]
    var_M <- names(summary(logitMod_M)$coefficients[,4])[-1]
    
    values_m <- matrix(nrow=length(var_m), ncol=length(c(tejidos)))
    dimnames(values_m) <- list(var_m, c(tejidos))
    
    values_M <- matrix(nrow=length(var_M), ncol=length(c(tejidos)))
    dimnames(values_M) <- list(var_M, c(tejidos))
  }
  
  pvalue_m <- summary(logitMod_m)$coefficients[,4]
  pvalue_m <- pvalue_m[-1]
  
  for(p in 1:length(var_m)){
    values_m[p, tejidos[tej]] <- unname(pvalue_m[p])
  }
  
  pvalue_M <- summary(logitMod_M)$coefficients[,4]
  pvalue_M <- pvalue_M[-1]
  
  for(p in 1:length(var_M)){
    values_M[p, tejidos[tej]] <- unname(pvalue_M[p])
  }
  
  print(tejidos[tej])
  
  plotROC(testData_m$disease, predicted_m)
  
  plotROC(testData_M$disease, predicted_M)
}

