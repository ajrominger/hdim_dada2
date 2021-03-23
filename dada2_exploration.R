library(dada2)

dataPath <- '~/Dropbox/HDIM_beating'

rawForward <- list.files(dataPath, pattern="_R1_001.fastq", full.names = TRUE)[1:10]
rawReverse <- list.files(dataPath, pattern="_R2_001.fastq", full.names = TRUE)[1:10]
sampNames <- gsub('.*/|_R1_001.*', '', rawForward)

plotQualityProfile(rawForward)
plotQualityProfile(rawReverse)



filterPath <- file.path(dataPath, "filtered")

filterForward <- file.path(filterPath,
                           paste0(sampNames, "_R1_trimmed.fastq.gz"))

filterReverse <- file.path(filterPath,
                           paste0(sampNames, "_R2_trimmed.fastq.gz"))


out <- filterAndTrim(fwd = rawForward, filt = filterForward, 
                     rev = rawReverse, filt.rev = filterReverse, 
                     truncLen = c(230, 230),
                     maxN = 0, truncQ = 2, maxEE = c(1, 2), 
                     compress = TRUE, multithread = TRUE)
head(out)



errorsForward <- learnErrors(filterForward, multithread = TRUE)
errorsReverse <- learnErrors(filterReverse, multithread = TRUE)

plotErrors(errorsForward, nominalQ = TRUE)
plotErrors(errorsReverse, nominalQ = TRUE)


derepForward <- derepFastq(filterForward, verbose = TRUE)
derepReverse <- derepFastq(filterReverse, verbose = TRUE)

names(derepForward) <- sampNames
names(derepReverse) <- sampNames

dadaForward <- dada(derepForward, err = errorsForward, multithread = TRUE)
dadaReverse <- dada(derepReverse, err = errorsReverse, multithread = TRUE)


merged <- mergePairs(dadaForward, derepForward, 
                     dadaReverse, derepReverse, 
                     verbose = TRUE)

# inspect the merger data.frame from the first sample
head(merged[[1]])



seqTab <- makeSequenceTable(merged)
dim(seqTab)

# inspect distribution of sequence lengths
plot(density(nchar(getSequences(seqTab))))

seqTabClean <- removeBimeraDenovo(seqTab, method = 'consensus',
                                  multithread = TRUE, verbose = TRUE)
dim(seqTabClean)

plot(density(nchar(getSequences(seqTabClean))))

sum(seqTabClean)
sum(seqTab)
