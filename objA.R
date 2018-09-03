#------------------------------------------------------------------------#

##Author: Alfonso Vidal García
##Co-Author: Juan Antonio Botía Blaya

#------------------------------------------------------------------------#


library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

require(data.table)

A_m <- read.delim("n_biotype_mod.txt", stringsAsFactors = FALSE) #Con radio (-5000, +5000)
A_M <- read.delim("n_biotype_mod_M.txt", stringsAsFactors = FALSE) #Con radio (-500000, +500000)

## To this point we have the two tables of the species of RNA

discapacidad <- read.delim("DI_Oct2016.csv", stringsAsFactors = FALSE, header = FALSE)$V1
control <- read.delim("controlForSonia.tsv", stringsAsFactors = FALSE, header = FALSE)$V1

A_m <- A_m[which(A_m$length > 0),]
sumas_m <- colSums (A_m[4:17], na.rm = FALSE, dims = 1)
A_m <- A_m[c('external_gene_name', 'length', 'biotype',names(sumas_m[which(sumas_m > 0)]))]

A_M <- A_M[which(A_M$length > 0),]
sumas_M <- colSums (A_M[4:17], na.rm = FALSE, dims = 1)
A_M <- A_M[c('external_gene_name', 'length', 'biotype',names(sumas_M[which(sumas_M > 0)]))]

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

comp_a_disease <- comp_a[which(comp_a[,'disease'] == 1),]
comp_a_no_disease <- comp_a[which(comp_a[,'disease'] == 0),]
set.seed(100)

comp_a_disease_training_rows <- sample(1:nrow(comp_a_disease), 0.7*nrow(comp_a_disease))
comp_a_no_disease_training_rows <- sample(1:nrow(comp_a_no_disease), 0.7*nrow(comp_a_no_disease))

training_ones <- comp_a_disease[comp_a_disease_training_rows, ]
training_zeros <- comp_a_no_disease[comp_a_no_disease_training_rows, ]

trainingData <- rbind(training_ones, training_zeros)

test_ones <- comp_a_disease[-comp_a_disease_training_rows, ]
test_zeros <- comp_a_no_disease[-comp_a_no_disease_training_rows, ]

testData <- rbind(test_ones, test_zeros)

comp_a_M_disease <- comp_a_M[which(comp_a_M[,'disease'] == 1),]
comp_a_M_no_disease <- comp_a_M[which(comp_a_M[,'disease'] == 0),]
set.seed(100)

comp_a_M_disease_training_rows <- sample(1:nrow(comp_a_M_disease), 0.7*nrow(comp_a_M_disease))
comp_a_M_no_disease_training_rows <- sample(1:nrow(comp_a_M_no_disease), 0.7*nrow(comp_a_M_no_disease))

training_ones_M <- comp_a_M_disease[comp_a_M_disease_training_rows, ]
training_zeros_M <- comp_a_M_no_disease[comp_a_M_no_disease_training_rows, ]

trainingData_M <- rbind(training_ones_M, training_zeros_M)

test_ones_M <- comp_a_M_disease[-comp_a_M_disease_training_rows, ]
test_zeros_M <- comp_a_M_no_disease[-comp_a_M_no_disease_training_rows, ]

testData_M <- rbind(test_ones_M, test_zeros_M)

formula <- as.formula(disease ~ length)

formula <- as.formula(
  paste(
    paste(deparse(formula), collapse=""), 
    paste(colnames(comp_a[4:(ncol(A)-1)]), collapse="+"),
    sep="+"
  )
)

logitMod <- glm(formula, data=trainingData, family=binomial(link="logit"))

predicted <- plogis(predict(logitMod, testData)) 

formula_M <- as.formula(disease ~ length)

formula_M <- as.formula(
  paste(
    paste(deparse(formula_M), collapse=""), 
    paste(colnames(comp_a_M[4:(ncol(A_M)-1)]), collapse="+"),
    sep="+"
  )
)

logitMod_M <- glm(formula_M, data=trainingData_M, family=binomial(link="logit"))

predicted_M <- plogis(predict(logitMod_M, testData_M)) 


library(InformationValue)
optCutOff <- optimalCutoff(testData$disease, predicted)[1] 
optCutOff_M <- optimalCutoff(testData_M$disease, predicted_M)[1]

summary(logitMod)
summary(logitMod_M)

predicted <- plogis(predict(logitMod, testData))  

sensitivity(testData$disease, predicted, threshold = optCutOff)
specificity(testData$disease, predicted, threshold = optCutOff)

plotROC(testData$disease, predicted)

predicted_M <- plogis(predict(logitMod_M, testData_M))  

sensitivity(testData_M$disease, predicted_M, threshold = optCutOff_M)
specificity(testData_M$disease, predicted_M, threshold = optCutOff_M)

plotROC(testData_M$disease, predicted_M)