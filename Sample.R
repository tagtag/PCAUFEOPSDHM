#This is  a sample R script that performs analyses in the paper
# Tensor decomposition and principal component analysis-based unsupervised feature extraction outperforms state-of-the-art methods when applied to histone modification profiles
#Sanjiban Sekhar Roy, Y-h. Taguchi
#bioRxiv 2022.04.29.490081; doi: https://doi.org/10.1101/2022.04.29.490081

#GSE24850
#Current directory should include the following eleven files
#GSM611267_sa.rep1.bed.gz     GSM612980_input.rep2.bed.gz
#GSM611268_sa.rep2.bed.gz     GSM612981_input.rep3.bed.gz
#GSM611269_sa.rep3.bed.gz     GSM612982_input.rep4.bed.gz
#GSM611270_co.rep1.bed.gz     GSM612983_input.rep5.bed.gz
#GSM611271_co.rep2.bed.gz     GSM612984_input.rep6.bed.gz
#GSM612979_input.rep1.bed.gz


require(rtracklayer)
L<-25000
files <- list.files("./",pattern="bed.gz")
i<-1
x <- import(files[i])
seqlengths(x) <- getChromInfoFromUCSC("mm9")[match(names(seqlengths(x)),getChromInfoFromUCSC("mm9")[,1]),2]
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
gr0.data.binnedAvg <- binnedAverage(gr.windows, coverage(x), "value")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("mm9")[match(names(seqlengths(x)),getChromInfoFromUCSC("mm9")[,1]),2]
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    gr.data.binnedAvg <- binnedAverage(gr.windows, coverage(x), "value")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}

Z <- array(NA,c(dim(x_all)[1],5,2))
Z[,,1] <- x_all[,1:5]
Z[,,2] <- x_all[,6:10]

PCA <- prcomp(scale(Z[,,1]))


th <- function(sd){
    P2 <- pchisq((PCA$x[,1]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(1e-3,th)$par
P <- pchisq((PCA$x[,1]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,1])))
#[1] 0.04264199 0.04128722

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()
table(p.adjust(P,"BH")<0.01)
#FALSE   TRUE 
#104902   1302 

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene #shorthand (for convenience) txdb
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
write.table(file="gene.csv",mcols(GR),col.names=F,row.names=F,quote=F,sep="\t")

#GSE159075
# Current directory should include the following ten files
# "GSM5026266_1_input.bigWig"           "GSM5026267_2_time0_H3K4me3.bigwig"  
# "GSM5026268_3_time15_H3K4me3.bigWig"  "GSM5026269_4_time60_H3K4me3.bigWig" 
# "GSM5026270_5_time0_H3K27me3.bigWig"  "GSM5026271_6_time15_H3K27me3.bigWig"
# "GSM5026272_7_time60_H3K27me3.bigWig" "GSM5026273_8_time0_H3K27ac.bigwig"  
# "GSM5026274_9_time15_H3K27ac.bigwig"  "GSM5026275_10_time60_H3K27ac.bigWig"

L<-1000 
files <- list.files("./",pattern=".bigwig")
files <- c(files,list.files("./",pattern=".bigWig"))
files<-sort(files)
files <- files[-c(2:4)] #H3K27me3
files <- files[c(1,8:10)]#H3K27ac
i<-1
x <- import(files[i])
seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
x <- trim(x)
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
score <- score[1:25]
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
    x <- trim(x)
    #seqlengths(x)[26] <- 5e+6 #not for #H3K27ac
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


X <- x_all[rowSums(x_all[,2:4]==0)==0,2:4]
PCA <- prcomp(scale(X))

th <- function(sd){
    P2 <- pchisq((PCA$x[,1]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.1,th)$par
P <- pchisq((PCA$x[,1]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,1])))
#[1] 0.19218492 0.07708293 #H3K4me3
#[1] 0.8650455 0.7584596 #H3K27me3
#[1] 0.3911769 0.2348521 #H3K27ac


pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()
table(p.adjust(P,"BH")<0.01)
#FALSE    TRUE 
#1999143   34538 #H3K4me3

#FALSE    TRUE 
#2527728   62141 #H3K27me3

#FALSE    TRUE 
#1119779   61306 #H3K27ac

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[rowSums(x_all[,2:4]==0)==0][p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
genes <- unique(mcols(GR))
write.table(file="gene.csv",genes,col.names=F,row.names=F,quote=F,sep="\t")

#GSE74055
# Current directory should include the following 32 files
# "GSM1908888_ChIPseq-H3K4me1_E14_0831.H3K4me1_E14.bw"
# "GSM1908889_ChIPseq-H3K4me1_DKO_0831.H3K4me1_DKO.bw"
# "GSM1908890_ChIPseq-Input_E14_0831.Input_E14.bw"    
# "GSM1908891_ChIPseq-Input_DKO_0831.Input_DKO.bw"    
# "GSM1908894_ChIPseq-H3K4me1_E14_0803.AY2151.bw"     
# "GSM1908895_ChIPseq-H3K4me1_DKO_0803.AY2152.bw"     
# "GSM1908896_ChIPseq-Input_E14_0803.AY2145.bw"       
# "GSM1908897_ChIPseq-Input_DKO_0803.AY2146.bw"       
# "GSM2099808_ChIP_D0_E14_H3K4me1.bw"                 
# "GSM2099810_ChIP_D0_E14_input.bw"                   
# "GSM2099812_ChIP_D1_E14_H3K4me1.bw"                 
# "GSM2099814_ChIP_D1_E14_input.bw"                   
# "GSM2099816_ChIP_D1_5_E14_H3K4me1.bw"               
# "GSM2099818_ChIP_D1_5_E14_input.bw"                 
# "GSM2099820_ChIP_D2_0_E14_H3K4me1.bw"               
# "GSM2099822_ChIP_D2_0_E14_input.bw"                 
# "GSM2099824_ChIP_D2_5_E14_H3K4me1.bw"               
# "GSM2099826_ChIP_D2_5_E14_input.bw"                 
# "GSM2099828_ChIP_D3_0_E14_H3K4me1.bw"               
# "GSM2099830_ChIP_D3_0_E14_input.bw"                 
# "GSM2099832_ChIP_D0_DKO_H3K4me1.bw"                 
# "GSM2099834_ChIP_D0_DKO_input.bw"                   
# "GSM2099836_ChIP_D1_DKO_H3K4me1.bw"                 
# "GSM2099838_ChIP_D1_DKO_input.bw"                   
# "GSM2099840_ChIP_D1_5_DKO_H3K4me1.bw"               
# "GSM2099842_ChIP_D1_5_DKO_input.bw"                 
# "GSM2099844_ChIP_D2_0_DKO_H3K4me1.bw"               
# "GSM2099846_ChIP_D2_0_DKO_input.bw"                 
# "GSM2099848_ChIP_D2_5_DKO_H3K4me1.bw"               
# "GSM2099850_ChIP_D2_5_DKO_input.bw"                 
# "GSM2099852_ChIP_D3_0_DKO_H3K4me1.bw"               
# "GSM2099854_ChIP_D3_0_DKO_input.bw"    

files <- list.files("./",pattern=".bw")
L<-1000
i<-1
x <- import(files[i])
seqlengths(x) <- getChromInfoFromUCSC("mm9")[match(names(seqlengths(x)),getChromInfoFromUCSC("mm9")[,1]),2]
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("mm9")[match(names(seqlengths(x)),getChromInfoFromUCSC("mm9")[,1]),2]
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}

X <- x_all[,-grep("input",files,ignore.case=T)]/x_all[,grep("input",files,ignore.case=T)]
index0 <- rowSums(is.infinite(X))==0
X <- X[index0,]
X[is.na(X)]<-0
PCA<-prcomp(scale(X))

th <- function(sd){
    P2 <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.1,th)$par
P <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,2])))
#[1] 0.1873585 0.1577631 

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

table(p.adjust(P,"BH")<0.01)

#FALSE    TRUE 
#1466911   61329 

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene #shorthand (for convenience) txdb
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[index0][p.adjust(P,"BH")<0.01]) 
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t")

#GSE124690
#Current directory should include the following 14 files
#GSM3536497_H1_K27ac_Rep1.bed.gz    GSM3536514_K562_K27ac_Rep2.bed.gz
#GSM3536497_H1_K27ac_Rep2.bed.gz    GSM3536516_K562_K4me1_Rep1.bed.gz
#GSM3536499_H1_K4me1_Rep1.bed.gz    GSM3536516_K562_K4me1_Rep2.bed.gz
#GSM3536499_H1_K4me1_Rep2.bed.gz    GSM3536518_K562_K4me3_Rep1.bed.gz
#GSM3536501_H1_K4me3_Rep1.bed.gz    GSM3536518_K562_K4me3_Rep2.bed.gz
#GSM3536501_H1_K4me3_Rep2.bed.gz    GSM3680223_K562_H3K4me1_Abcam_8895.bed.gz
#GSM3536514_K562_K27ac_Rep1.bed.gz  GSM3680224_K562_H3K4me1_ActMot_39113.bed.gz

require(rtracklayer)
#files <- list.files("./",pattern="K4me1") #H3K4me1
#files <- list.files("./",pattern="K4me3") #H3K4me3
files <- list.files("./",pattern="K27ac") #H3K27ac
L<-1000
i<-1
x <- import(files[i])
seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
gr0.data.binnedAvg <- binnedAverage(gr.windows, coverage(x), "value")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    gr.data.binnedAvg <- binnedAverage(gr.windows, coverage(x), "value")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


PCA<-prcomp(scale(x_all))

k<-1 #H3K4me3, H3K27ac
k<-2 #H3K4me1
th <- function(sd){
    P2 <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th)$par 
P <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,k])))
#[1] 0.5468109 0.4257019 # H3K4me1 
#[1] 0.5600861 0.1208255 # H3K4me3 
#[1] 0.4509331 0.3174211 # H3K27ac 


table(p.adjust(P,"BH")<0.01)
#FALSE    TRUE 
#2931240  164466 # H3K4me1

#FALSE    TRUE 
#3058172   37534 # H3K4me3 

#FALSE    TRUE 
#3014457   81249  # H3K27ac

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t") #k=1


#GSE188173
#Current directory should include the following 16 files
#GSM5671438_K27ac105CR.bw    GSM5671444_K27ac167CR.bw    GSM5671450_K27ac77CR.bw
#GSM5671439_K27ac105CR_T.bw  GSM5671445_K27ac167CR_T.bw  GSM5671451_K27ac77CR_T.bw
#GSM5671440_K27ac136CR.bw    GSM5671446_K27ac35CR.bw     GSM5671452_K27ac96CR.bw
#GSM5671441_K27ac136CR_T.bw  GSM5671447_K27ac35CR_T.bw   GSM5671453_K27ac96CR_T.bw
#GSM5671442_K27ac147CR.bw    GSM5671448_K27ac73CR.bw
#GSM5671443_K27ac147CR_T.bw  GSM5671449_K27ac73CR_T.bw

require(rtracklayer)
L<-1000 
files <- list.files("./",pattern=".bw")
i<-1
x <- import(files[i])
seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
x <- trim(x)
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
score <- score[1:25]
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
    x <- trim(x)
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


PCA <- prcomp(scale(x_all))

k<-1
th <- function(sd){
    P2 <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th)$par 
P <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,k])))
#[1] 0.5319398 0.2854939

table(p.adjust(P,"BH")<0.01)
#FALSE    TRUE 
#2990268  105438 

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t") #k=1

#GSE159022
#Current directory should include the following four files
#GSM4817504_Het-1_H3K27me3_rpkm.bw  GSM4817506_KO-1_H3K27me3_rpkm.bw
#GSM4817505_Het-2_H3K27me3_rpkm.bw  GSM4817507_KO-2_H3K27me3_rpkm.bw

require(rtracklayer)
L<-1000 
files <- list.files("./",pattern=".bw")
i<-1
x <- import(files[i])
seqlengths(x) <- getChromInfoFromUCSC("mm9")[match(names(seqlengths(x)),getChromInfoFromUCSC("mm9")[,1]),2]
x <- trim(x)
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
score <- score[1:25]
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("mm9")[match(names(seqlengths(x)),getChromInfoFromUCSC("mm9")[,1]),2]
    x <- trim(x)
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


PCA <- prcomp(scale(x_all))

k<-1
th <- function(sd){
    P2 <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th)$par 
P <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,k])))
#[1] 0.06544008 0.05732543

table(p.adjust(P,"BH")<0.01)

#FALSE    TRUE 
#2599001   55923 

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene #shorthand (for convenience) txdb
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t")

#GSE168971
#Current directory should include the following eight files
#GSM5173256_3065_9.ALL.nomulti.bw  GSM5173259_3068_9.ALL.nomulti.bw  GSM5173286_3066_Inp.ALL.nomulti.bw
#GSM5173257_3066_9.ALL.nomulti.bw  GSM5173260_3078_9.ALL.nomulti.bw  GSM5173287_3068_Inp.ALL.nomulti.bw
#GSM5173258_A165_9.ALL.nomulti.bw  GSM5173261_A170_9.ALL.nomulti.bw

require(rtracklayer)
L<-1000 
files <- list.files("./",pattern=".bw")
i<-1
x <- import(files[i])
x <- keepStandardChromosomes(x,pruning.mode="coarse")
seqlengths(x) <- getChromInfoFromUCSC("mm10")[1:length(seqlengths(x)),2]
x <- trim(x)
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    x <- keepStandardChromosomes(x,pruning.mode="coarse")
    seqlengths(x) <- getChromInfoFromUCSC("mm10")[1:length(seqlengths(x)),2]
    x <- trim(x)
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


PCA <- prcomp(scale(x_all))

k<-1
th <- function(sd){
    P2 <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th)$par 
P <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,k])))
#[1] 0.3697877 0.2528360 

table(p.adjust(P,"BH")<0.01)

#FALSE    TRUE 
#2667059   58490

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene #shorthand (for convenience) txdb
GR <- genes(txdb)
seqlevels(gr0.data.binnedAvg) <- seqlevels(GR)[1:22]
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t")

#GSE159411
#Current directory should include the following four files
#GSM4829300_h3k36me3_cpc.1_hg19.chr_tagDensity.bw  GSM4829307_h3k36me3_ipsc.1_hg19.chr_tagDensity.bw
#GSM4829301_h3k36me3_cpc.2_hg19.chr_tagDensity.bw  GSM4829308_h3k36me3_ipsc.2_hg19.chr_tagDensity.bw

require(rtracklayer)
L<-1000 
files <- list.files("./",pattern=".bw")
i<-1
x <- import(files[i])
x <- keepStandardChromosomes(x,pruning.mode="coarse")
seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
x <- trim(x)
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    seqlengths(x) <- getChromInfoFromUCSC("hg19")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg19")[,1]),2]
    x <- trim(x)
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


PCA <- prcomp(scale(x_all))

k<-1
th <- function(sd){
    P2 <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th)$par 
P <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,k])))
# [1] 0.6005592 0.4890187 

table(p.adjust(P,"BH")<0.01)

#FALSE    TRUE 
#2842363  253326 

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t") #k=1


#GSE181596
#Current directory should include the following four files
#GSM5507791_ss1463_1_rm_dups.bw  GSM5507793_ss1463_3_rm_dups.bw
#GSM5507792_ss1463_2_rm_dups.bw  GSM5507794_ss1463_4_rm_dups.bw

require(rtracklayer)
L<-1000 
files <- list.files("./",pattern=".bw")
i<-1
x <- import(files[i])
x <- keepStandardChromosomes(x,pruning.mode="coarse")
seqlengths(x) <- getChromInfoFromUCSC("hg38")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg38")[,1]),2]
x <- trim(x)
gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
score <-  mcolAsRleList(x, "score")
score[is.na(score)]<-0
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
x_all <- data.matrix( mcols(gr0.data.binnedAvg))
for (i in c(2:length(files)))
{
    cat(i," ")
    x <- import(files[i])
    x <- keepStandardChromosomes(x,pruning.mode="coarse")
    seqlengths(x) <- getChromInfoFromUCSC("hg38")[match(names(seqlengths(x)),getChromInfoFromUCSC("hg38")[,1]),2]
    x <- trim(x)
    gr.windows <- tileGenome(seqinfo(x), tilewidth=L,cut.last.tile.in.chrom=TRUE)
    score <-  mcolAsRleList(x, "score")
    score[is.na(score)]<-0
    gr.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
    index <- match(granges(gr0.data.binnedAvg),granges(gr.data.binnedAvg))
    x_all <- cbind(x_all,mcols(gr.data.binnedAvg)[index,])
}


PCA<- prcomp(scale(x_all))

k<-1
th <- function(sd){
    P2 <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th)$par 
P <- pchisq((PCA$x[,k]/sd)^2,1,lower.tail=F)
print(c(sd,sd(PCA$x[p.adjust(P,"BH")>0.01,k])))
#[1] 1.855494 1.654783


table(p.adjust(P,"BH")<0.01)
#FALSE    TRUE 
#3051309   36972 

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #shorthand (for convenience) txdb
GR <- genes(txdb)
seqlevels(gr0.data.binnedAvg) <- seqlevels(GR)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01])
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t")

#GSE74055 with TD
#till get x_all, see GSE74055

require(rTensor)
Z <- array(NA,c(dim(x_all)[1],16,2))
Z[,,1] <- data.matrix(x_all[,c(1,2,5,6,seq(9,32,by=2))])
Z[,,2] <- data.matrix(x_all[,c(3,4,7,8,seq(10,32,by=2))])
Z <- apply(Z,2:3,scale)
HOSVD <- hosvd(as.tensor(Z),c(10,16,2))
plot(HOSVD$U[[2]][,1],type="h",ylim=c(-0.3,0))
plot(HOSVD$U[[3]][,2],type="h")
HOSVD$Z@data[,1,2]
#[1]  154.775233 3948.725962   87.294546  182.883755   18.789742  -67.358369
#[7]   42.548520    1.880876  -12.800926   65.154321


th <- function(sd){
    P2 <- pchisq((HOSVD$U[[1]][,2]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=100,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.1,th)$par
P <- pchisq((HOSVD$U[[1]][,2]/sd)^2,1,lower.tail=F)
print(c(sd,sd(HOSVD$U[[1]][p.adjust(P,"BH")>0.01,2])))
#[1] 0.0003451553 0.0003089417

pdf(file="hist.pdf")
hist(1-P,breaks=100)
dev.off()

table(p.adjust(P,"BH")<0.01)
#FALSE    TRUE 
#2584737   70187 

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene #shorthand (for convenience) txdb
GR <- genes(txdb)
ov <- findOverlaps(GR,granges(gr0.data.binnedAvg)[p.adjust(P,"BH")<0.01]) 
GR <- GR[queryHits(ov)]
gene <- unique(mcols(GR))
write.table(file="gene.csv",gene,col.names=F,row.names=F,quote=F,sep="\t")


