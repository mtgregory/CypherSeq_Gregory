## Targeted CypherSeq script
## Written by Mark Gregory
## 11/19/2014

#######
##initialize the program
#####
library(BSgenome)
library(Biostrings)
library(ShortRead)
library(Rsamtools)
library(latticeExtra)

###########
#parameters
#################

#Adaptor trimming and Fastq filtering settings
inputFile1="R1.fastq.gz"  # Fastq file for reads1 (.fastq required)
inputFile2="R2.fastq.gz"  # Fastq file for reads2 (.fastq required)
adaptor="adaptor.txt"     # adaptor sequence in fasta format. Remember to change adaptor in line 70 also if you change it here!
barcodeLength=7        # Miniumum value is 5
tagLength=3      		  # e.g. CCC
uncorrectLength=13		# The # 5prime leading bases of the template (after tag) where no error is tolerated. Also used to determine molecule orientations.
mismatch.rate=0.05      # rate of allowed mismatches for adaptor trimming
minReadLength=50      # minimum read length of a paired end read

#Barcode distribution settings
BCDistPreAlgn=TRUE          #If true, saves pre-alignment barcode distributions
BCDistPostAlgn=TRUE     #If true, saves post-alignment barcode distributions
BCDistQC=TRUE           #If true, saves post quality control barcode distributions

#Bowtie Alignment settings
reference="genome_rCRS" # Base name for reference genome
bwa="bwa" # Specify the aligner directory path
alnCores = 3  # Number of cores to parallelize on during BWA alignment


#BAM QC filtering settings
chromosome='chr17'  #chromosome location of reference sequence. Must be in UCSC genome format, ie 'chr#'
gene_start=7578371  #start and ending chromosome coordinates of reference sequence 
gene_end=7578554
phredEncoding=33           # Enter 64 for HiSeq
qualType="FastqQuality"    # Do NOT change
qualityThreshold=30        # Valid range is [2,40]
maxLQFraction=0.2         # Fraction of bases in a read allowed to have Qscore < qualityThreshold. This setting filters out entire reads that have a large fraction of low quality basecalls.
qualityControl=TRUE      # Enter TRUE low quality bases should be masked. This setting filters out low quality bases from reads that survive the maxLQfraction.
maxReadN=0                #Maximum number of N's allowed per read
leftEndMask <- 3          #Number of bases to be masked on left end.  Used to eliminate end repair substitutions and end proximal indels incorrectly called as mismatches.
RightEndMask <- 3         #Number of bases to be masked on right end.  Used to eliminate end repair substitutions and end proximal indels incorrectly called as mismatches.

#Barcode family settings
barcodeFamilySize=10  	# Miniumum number of reads in a barcode  family. for example, Family size =10 is equivalent to 5 readpairs, or 10 reads.
consensusLevel=0.9  		# Valid range is [0.5,1]
minFamNum=1           #Min number of families a mutation must be found in to be validated

# Processing settings
blocks = 4    # number of blocks to break the .sam file into for processing
chunks=30 #number of pieces to break the data into for parallel processing during QC filtering
consensusCores = 2  # Number of cores to parallelize on during barcode family consensus building

##file structure settings
if(!file.exists('results')) dir.create('results/')


#################
#functions
#################

#Function used to trim adaptor sequences from barcodes
FastqTrimmer <- function (fastq, barcodeLength, 
                          tagLength, uncorrectLength, 
                          mismatch.rate) {
  #Isolate reads
  fqreads<-sread(fastq)
  
  ## Find barcodes
  barcode<-subseq(fqreads, start=1, end=barcodeLength+tagLength+uncorrectLength)
  
  #Trim barcodes and adaptor sequences
  adaptor<-readDNAStringSet('adaptor.txt')[[1]]
  adaptorLength <- nchar(adaptor)
  max.Rmismatch<-mismatch.rate*1:adaptorLength
  adaptor<-c(adaptor, DNAString(paste(rep('N',(151-barcodeLength-tagLength-adaptorLength)), collapse='')))
  max.Rmismatch<-append(max.Rmismatch, rep(max(max.Rmismatch), (151-barcodeLength-tagLength-adaptorLength)))
  max.Rmismatch[1:(barcodeLength+tagLength)]<-0
  fqreads.trim<-subseq(fqreads, start=barcodeLength+tagLength+1)
  trim<-trimLRPatterns(Rpattern=adaptor, subject=fqreads.trim, max.Rmismatch=max.Rmismatch, Rfixed='subject', ranges=TRUE, with.Rindels=TRUE)
  
  fqreads.trim<-subseq(fqreads.trim, start=start(trim), end=end(trim))
  qual.trim<-subseq(quality(quality(fastq)), start=barcodeLength+tagLength+start(trim), end=barcodeLength+tagLength+end(trim))
  
  
  dfFastq<-DataFrame(sread=fqreads.trim, quality=qual.trim, id=id(fastq), barcode=barcode)
  return (dfFastq)
}

#Function used to filter bam files
myFilter2<-function(df){
  idx<-rep(TRUE, nrow(df))
  
  ##extract read info
  qual <- df$qual
  qseq<-df$seq
  rwidth<-width(qseq)
  
  ##filter based on quality
  idx[rowSums(as(FastqQuality(qual),'matrix')<qualityThreshold, na.rm=TRUE)/rwidth>maxLQFraction]<-FALSE
  
  ##filter based on number of Ns
  numberOfN<-as.vector(letterFrequency(qseq,"N"))
  idx[which(numberOfN>maxReadN)]<-FALSE
  
  message("Records read: ", nrow(df), "; Records accepted: ", sum(idx))  
  return (idx)
}


#Function used to mark low quality bases as missing data
replaceLowQualityBases<-function(reads,quals,min.qual,letter) {
  if (!is.character(min.qual) || length(min.qual) != 1L || is.na(min.qual) || nchar(min.qual) != 1L)
    stop("'min.qual' must be a single letter")
  if (!is.character(letter) || length(letter) != 1L || is.na(letter) || nchar(letter) != 1L)
    stop("'letter' must be a single letter")
  raw1<-charToRaw(min.qual)
  
  
  unlisted_reads<-unlist(reads)
  unlisted_quals<-unlist(quals)
  at<-which(as.raw(unlisted_quals) < raw1)
  unlisted_ans<-replaceLetterAt(unlisted_reads,at,rep.int(letter,length(at)))
  ans<-as(successiveViews(unlisted_ans,width(reads)),"DNAStringSet")
  return (ans)
}






###########
## Process fastq files and generate new, trimmed files for alignment ## 0:07:00 for 21mn paired reads
############
#Test to see if trimmed fastq files have already been generated. 
test<-file.exists('results/R1.trim.fastq') & file.exists('results/R2.trim.fastq')

if (!test){
  fq1<-FastqStreamer(inputFile1)
  fq2<-FastqStreamer(inputFile2)
  dfBarcode<-DataFrame()
  
  #Read in fastq files one line at a time and perform filtering and trimming, then append to new file
  while (length(fastq1<-yield(fq1)) & length(fastq2<-yield(fq2))) {
    message("fastq1: ", length(fastq1), "; fastq2: ", length(fastq2))
    
    #Trim barcodes and adaptor sequences
    dfFastq<-mclapply(c(fastq1, fastq2), FastqTrimmer, barcodeLength, tagLength, uncorrectLength, mismatch.rate)
    
    #filter barcodes with Ns
    numberOfN<-lapply(dfFastq, function (fq) as.vector(letterFrequency(fq$barcode,"N")))
    idx<-numberOfN[[1]]==0 & numberOfN[[2]]==0
    
    #filter reads shorter than minReadLength
    idx[which(width(dfFastq[[1]]$sread)<minReadLength | width(dfFastq[[2]]$sread)<minReadLength)]<-FALSE
    dfFastq<-lapply(dfFastq, function (fq) fq[idx,])
    
    #append trimmed reads to new fastq files
    FastqList<-lapply(dfFastq, function(fq)
      ShortReadQ(sread=fq$sread, quality=fq$quality, id=fq$id))
    
    fqout<-c('results/R1.trim.fastq', 'results/R2.trim.fastq')
    
    mclapply(c(1,2), function(i) 
      writeFastq(FastqList[[i]], fqout[i], mode='a', full=TRUE, compress=FALSE))
    
    
    #add barcodes to dataframe
    dfBarcode<<-rbind(dfBarcode, DataFrame(ID1=dfFastq[[1]]$id, ID2=dfFastq[[2]]$id, 
                                           barcode=xscat(dfFastq[[1]]$barcode, dfFastq[[2]]$barcode)))
    
  }
  
  #save barcode file
  save(dfBarcode, file='results/dfBarcode.rda')}

############
## Alignment ## 1.5 hours
###############
#Test to see if alignament has been completed
test2<-file.exists('results/precorrection.sorted.bam.bai')

if (!test2){
  system2(bwa,c("mem", "-t", alnCores ,paste0('BWA_Genome_index/',reference,'.fa'),"results/R1.trim.fastq",'results/R2.trim.fastq',">",'results/precorrection.sam'))
  system2("samtools",c("view","-bS",'results/precorrection.sam',">",'results/precorrection.bam'))
  system2("samtools",c("sort",'results/precorrection.bam','results/precorrection.sorted'))
  system2("samtools",c("index",'results/precorrection.sorted.bam'))
}


########################
#load the reference gene sequence from fasta file and format
#######################
refseq<-getSeq(FaFile(dir(file.path(getwd()), paste0(reference,'.fa'), full=TRUE)))
refseq<-subseq(refseq[chromosome], start=gene_start, end=gene_end)
refBase<-alphabetByCycle(refseq)
refseqlength<-width(refseq)

refRange<-GRanges(seqnames=chromosome,ranges=IRanges(start=gene_start, end=gene_end))
mcols(refRange)$qseq<-refseq

######################
#Check for barcode file and index barcodes. (About 5 min?)
#Optional: calculate post-alignment barcode distribution
####################
if (!exists('dfBarcode')) {
  if(file.exists('results/dfBarcode.rda')) { 
    load ('results/dfBarcode.rda')} else {
      stop("Barcode table does not exist yet")   
    }
}

family<-with(dfBarcode, match(barcode, barcode))
clusterID<-BStringSet(gsub('(.+?)\\s.+$', '\\1', dfBarcode$ID1))

#Option to save pre alignment distribution
if (BCDistPreAlgn){
  famCount<-table(family)
  famTally<-table(famCount)
  PreAlgnFamTally<-data.frame(FamilySize=as.integer(names(famTally)), 
                              Freq=as.integer(famTally))
  
  out.file<-file(paste0('results/BCDistPreAlgn_',format(Sys.time(), "%b_%d_%Y"),'.csv'),'w')
  writeLines('',out.file)
  write.table(famTally, file=out.file, sep=',', row.names=FALSE)
  close(out.file)
}

#Option to save post alignment barcode distribution
if (BCDistPostAlgn){
  param<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isProperPair=TRUE, hasUnmappedMate=FALSE),
                      what=c('qname'), which=refRange)
  gal<-readGAlignments('results/precorrection.sorted.bam', param=param)  
  qfam<-family[match(mcols(gal)$qname, clusterID)]  
  famSize<-table(qfam)/2  
  famTally<-table(famSize)
  PostAlgnFamTally<-data.frame(FamilySize=as.integer(names(famTally)), 
                               Freq=as.integer(famTally))
  
  out.file<-file(paste0('results/BCDistPostAlgn_',format(Sys.time(), "%b_%d_%Y"),'.csv'),'w')
  writeLines('',out.file)
  write.table(famTally, file=out.file, sep=',', row.names=FALSE)
  close(out.file)
}


####################
# Filter the results of the bam alignment based on alignment to target gene, quality
# scores and number of N's and save results to new filtered.bam file
####################
fl<-dir(file.path(getwd(), 'results'), 'precorrection.sorted.bam$', full=TRUE)
filtered.bam<-'results/precorrection.sorted.filtered.bam'

if (!file.exists(filtered.bam)){
  
  bf<-BamFile(fl, yieldSize=1500000)
  
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), 
                        what=c("seq", "qual"),
                        which=refRange)
  filterBam(bf, filtered.bam, param=param, filter=FilterRules(myFilter2))
}

###################
#read filtered.bam file into a GAlignments object
##################
param<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isProperPair=TRUE, hasUnmappedMate=FALSE),
                    what=c('seq', 'qual', 'qname', 'flag'))
gal<-readGAlignments(filtered.bam, param=param)

##################
## Failsafe against crashes due to too many reads at a single genomic location,
## which can occur in some versions of R.
## If the program gives an error about 'long vectors' increase 
## the 'blocks' parameter at the top of the file.
##################
if (length(gal)<1){
  
  system2("samtools",c("view","-HS",'results/precorrection.sam', '>', 'results/header.txt'))
  system2("samtools",c("view","-S",'results/precorrection.sam', '>', 'results/samBody.txt'))
  headerLength <- as.numeric(countLines('results/header.txt'))
  samTotal <- as.numeric(countLines('results/precorrection.sam'))
  samBlockLength <- ceiling(as.integer((samTotal - headerLength)/blocks))+1
  system2("split",c("-l",samBlockLength,'results/samBody.txt', 'results/samBlock'))
  samBlockList <- list.files('results')[grepl('samBlock..$', list.files('results'))]
  for(i in samBlockList) {
    system2("cat",c('results/header.txt', paste0('results/', i), '>', paste0('results/',i, '.sam')))
  }
  
  gal <- GAlignments()
  for(i in samBlockList) {
    system2("samtools",c("view","-bS",paste0('results/', i, '.sam'),">",'results/precorrection.bam'))
    system2("samtools",c("sort",'results/precorrection.bam','results/precorrection.sorted'))
    system2("samtools",c("index",'results/precorrection.sorted.bam'))
    
    fl<-dir(file.path(getwd(), 'results'), 'precorrection.sorted.bam$', full=TRUE)
    filtered.bam<-'results/precorrection.sorted.filtered.bam'
    
    bf<-BamFile(fl, yieldSize=1500000)
    
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), 
                          what=c("seq", "qual"),
                          which=refRange)
    filterBam(bf, filtered.bam, param=param, filter=FilterRules(myFilter2))
    
    param<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isProperPair=TRUE, hasUnmappedMate=FALSE),
                        what=c('seq', 'qual', 'qname', 'flag'))
    galTemp<-readGAlignments(filtered.bam, param=param)
    gal <- c(gal, galTemp)
  }
}

# mask read ends to eliminate false mismatch calls due to lower fidelity
# end repair and end-adjacent indels
gal <- qnarrow(gal, start= leftEndMask + 1, end=-(RightEndMask + 1))
lengthIDX <- nchar(mcols(gal)$seq)
mcols(gal)$seq <- DNAStringSet(substr(mcols(gal)$seq, leftEndMask + 1, lengthIDX-RightEndMask))
mcols(gal)$qual <- PhredQuality(substr(mcols(gal)$qual, leftEndMask + 1, lengthIDX-RightEndMask))

#################
#Optional: Use quality filter to mark low quality bases as missing data 
#(i.e.replace bases with low quality with +'s). (About 1 minute)
#################
qseq<-mcols(gal)$seq
qual<-mcols(gal)$qual


if(qualityControl) {
  qualityChar<-rawToChar(as.raw(phredEncoding+qualityThreshold))
  
  #If number of reads is large (>20 million?) the DNAStringSet has to be broken 
  #into smaller chunks for replaceLowQualityBases to work. This should also
  #speed things up by allowing for parallelization
  
  #define cut points
  nseqs<-length(qseq)
  seqBreaks<-cut(seq(1:nseqs), chunks)
  
  #cut sequence and quality scores into chunks
  tmp<-split(qseq, seqBreaks)
  tmpqual<-split(qual,seqBreaks)
  
  #replace low quality bases and recombine into single DNAStringSet
  qseqlist<-lapply(seq(chunks),function(x) 
    replaceLowQualityBases(tmp[[x]],tmpqual[[x]],qualityChar,"+" )) #,
  #     mc.cores=chunks/2)   
  qseq<-do.call(c,qseqlist) 
}

#################
# Use sequenceLayer() to "lay" the query sequences on the reference. 
# This allows family members to pile up properly. It also trims away sequences that were
#'soft trimmed' during the alignment. About 1.5 minutes
#################

#extract cigar and split into chunks
cigar<-cigar(gal)
cigarlist<-split(cigar,seqBreaks)

#overlay sequence on reference. mclapply did not work on very large data sets,
#so had to implement serially
qseq_on_ref_list<-lapply(seq(chunks), function(i) 
  sequenceLayer(qseqlist[[i]],cigarlist[[i]],from='query',to='reference'))

mcols(gal)$seq<-do.call(c, qseq_on_ref_list)


###########
#filter short reads
###########
gal<-gal[width(gal)>minReadLength]

###########
#filter unpaired reads
##########
qname<-mcols(gal)$qname
readID<-match(qname, qname)
tableReadID<-table(readID)
freq<-as.vector(tableReadID[match(readID, names(tableReadID))])

gal<-gal[freq>1]


#############3 
#Assign reads to barcode families, find family sizes, ID good
#families, filter gal on family size
###########
mcols(gal)$qfam<-family[match(mcols(gal)$qname, clusterID)]

famSize<-table(mcols(gal)$qfam)
goodFamilies<-as.numeric(names(famSize)[which(famSize>=barcodeFamilySize)])

gal<-gal[mcols(gal)$qfam %in% goodFamilies]

#Option to save barcode distribution after quality control filters
if(BCDistQC) {
  famCount<-table(mcols(gal)$qfam)/2
  famTally<-table(famCount)
  QCFamTally<-data.frame(FamilySize=as.integer(names(famTally)), 
                         Freq=as.integer(famTally))
  
  out.file<-file(paste0('results/BCDistQC_',format(Sys.time(), "%b_%d_%Y"),'.csv'),'w')
  writeLines('',out.file)
  write.table(famTally, file=out.file, sep=',', row.names=FALSE)
  close(out.file)
  save(famTally, file=paste0('results/BCDistQC_',format(Sys.time(), "%b_%d_%Y"),'.rda'))
  
  
  print(xyplot(Freq~FamilySize, QCFamTally,
               scales=list(y=list(log=TRUE)),
               yscale.components=yscale.components.log10ticks,
               main='Family Distribution'))
}


###########
#split gal into families and find family consensus sequences
##########
qfam<-mcols(gal)$qfam
gal_by_fam<-split(gal, qfam)
famSizes<-elementLengths(gal_by_fam)
galFam<-unlist(gal_by_fam)

index<-rep(seq(length(famSizes)), famSizes)
cutindex<-cut(index, chunks)
tmp<-split(galFam, cutindex)

#If num of families is big, may have to implement serially. With 37 million
#aligned reads, I did not have enough memory to parallelize on 10 cores, but was
#able to do it on 5 (chunks/2)
famConsensusList<-mclapply(tmp, function(galFam){
  
  qseq<-mcols(galFam)$seq
  qstart<-start(galFam)-gene_start+1
  famSizes<-elementLengths(split(galFam, mcols(galFam)$qfam))
  
  shift00<-seq(0,length(famSizes)-1)*refseqlength
  shift0<-rep(shift00, famSizes)
  shift<-shift0+qstart-1
  mat<-consensusMatrix(qseq, shift=shift, width=length(famSizes)*refseqlength,
                       as.prob=FALSE)
  letters<-uniqueLetters(qseq)
  mat<-mat[letters,]
  mat["+",]<-0
  
  #Don't consider positions with readDepth less than barcodeFamilySize
  readDepth<-colSums(mat)
  mat[,readDepth<barcodeFamilySize]<-0
  
  #'Fix' positions with no depth
  readDepth<-colSums(mat)
  mat["+",readDepth==0]<-1
  
  readDepth<-colSums(mat)
  cm<-sweep(mat,2,readDepth,"/")
  
  famConsensus<-consensusString(cm, ambiguityMap='+', threshold=consensusLevel)
  
  famInterval<-rep(refseqlength,length(famSizes))
  v<-Views(famConsensus,successiveIRanges(famInterval))
  
  famConsSet<-as(v, "DNAStringSet")
  return(famConsSet)
  
}, mc.cores=2)

famConsSet<-do.call('c', unlist(famConsensusList, use.names=FALSE))

#Save the famConsSet object. This saves it to the results folder and appends a time
#stamp.
save(famConsSet, file=paste0('results/famConsSet_',format(Sys.time(), "%b_%d_%Y"),'.rda'))


###############
##Generate consensus matrix of family consensus sequences and 
##calculate total site specific mutational frequency.
##############
#Letters used in variant calling
seqLetters<-uniqueLetters(famConsSet)

#Generate consensus matrix across all families. Set all missing bases ('+') to 0.
cm.fam<-consensusMatrix(famConsSet, as.prob=FALSE)
cm.fam["+",]<-0
cm.fam.raw<-cm.fam

#Find family depth per position
readDepth<-colSums(cm.fam)

#Do not consider mutations that arent supported by minFamNum families 
cm.fam[cm.fam<minFamNum]<-0
#Calculate basecall frequencies
cm.fam<-prop.table(cm.fam,2)  
#Do not consider non-mutational basecalls
cm.fam[refBase==1]<-0
#Calculate total site specific mutational frequency
snpFreq<-colSums(cm.fam)


##################
##Create and output summary data frames of results
#################

#Create a character vector containing all reference bases
refPosition<-seq(1,width(refRange),1)
RefSeq <- unlist(strsplit(as.character(refseq), ""))

#Vector of chromosome coordinates
coordinate=seq(start(refRange), end(refRange),1)

#Make a data frame summarizing the variant analysis for each position
dfErrFreq<-data.frame(Position=refPosition,ErrFreq=snpFreq, Position=refPosition, 
                      Coordinate= coordinate,
                      Coverage=readDepth, 
                      RefSeq=RefSeq)
dfErrFreq<-cbind(dfErrFreq, t(cm.fam[seqLetters,]))

#Make a data frame summarizing the number of families supporting each raw variant call
dfFamFreq<-data.frame(Position=refPosition, 
                      Coordinate= coordinate,
                      Coverage=readDepth, 
                      RefSeq=RefSeq)
dfFamFreq<-cbind(dfFamFreq, t(cm.fam.raw[seqLetters,]))

#Save the results
label<-"ConsensusFamilies"
dffile<-paste0('results/',label,'_',format(Sys.time(), "%b_%d_%Y"),'.csv')
out.file<-file(dffile,'w')
writeLines(as.character(Sys.time()), out.file)
writeLines(paste('Analysis stage:', label, sep=","), out.file) 
writeLines(paste('Alignment:',reference, sep=","), out.file)
writeLines(paste('Quality Threshold:', qualityThreshold, sep=","), out.file)
writeLines(paste('Max Low Quality Threshold:', maxLQFraction, sep=","), out.file)
writeLines(paste('Max allowed Ns:', maxReadN, sep=","), out.file)
writeLines(paste('Min read length:', minReadLength, sep=","), out.file)
writeLines(paste('Barcode Family Size:', barcodeFamilySize, sep=","), out.file)
writeLines(paste('Min Num of Families per mutation:', minFamNum, sep=","), out.file)
writeLines(paste('Consensus level:', consensusLevel, sep=","), out.file)
writeLines('',out.file)
write.table(dfErrFreq, file=out.file, sep=',', row.names=FALSE)
writeLines('',out.file)
writeLines("Number of families supporting each raw variant call", out.file)
write.table(dfFamFreq, file=out.file, sep=',', row.names=FALSE)
close(out.file) 
