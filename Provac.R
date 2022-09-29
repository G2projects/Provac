# --------------------------------------------------------------------#

# Title: Provac
# Developer: Fernando Palluzzi
# PI: Daniela Gallo
# Centre: Agostino Gemelli
# Since: 09/07/2021

# --------------------------------------------------------------------#



#### edgeR analysis


# R

library(limma)
#library(affy)
#library(affyPLM)
library(edgeR)
library(gplots)
library(ggplot2)
library(reshape2)
library(preprocessCore)
source("~/euNet/euNet_core/rutils.R")

samples <- read.table("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_sample.txt", header = TRUE, stringsAsFactors = FALSE)
grp <- c("sensitive", "sensitive", "sensitive", "sensitive", "sensitive", "sensitive", "sensitive", "resistant", "resistant", "resistant", "resistant", "resistant", "resistant", "resistant")
colors <- c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "red", "red", "red", "red", "red", "red", "red")

counts <- readDGE(samples$filename)$counts
colnames(counts) <- samples$samplename
head(counts)
counts <- counts[-20803:-20812,]

d <- DGEList(counts=counts, group=grp)
d <- calcNormFactors(d)
nrow(d)
d2 <- d[rowSums(cpm(d) > 1) > 2,]
nrow(d2)

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/MDS_plot.jpg", width = 15, height = 10, units = 'in', res = 300)
plotMDS(d2, col=colors)
abline(v = 0, lty = 2, col = "gray")
dev.off()

library(sva)
batch <- c("A", "D", "A", "C", "D", "B", "B", "A", "A", "A", "C", "B", "B", "B")
counts <- ComBat(d2$counts, batch)
counts <- counts[apply(counts, 1, min) > 0,]
head(counts)
dim(counts)

d <- DGEList(counts=counts, group=grp)
d <- calcNormFactors(d)
nrow(d)
d2 <- d[rowSums(cpm(d) > 1) > 2,]
nrow(d2)

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/MDS_plot_batchCorrection.jpg", width = 15, height = 10, units = 'in', res = 300)
plotMDS(d2, col=colors)
abline(v = 0, lty = 2, col = "gray")
dev.off()

counts <- counts[, c(-2, -8)]
d2$sample <- d2$sample[c(-2, -8),]

d2 <- estimateCommonDisp(d2)
d2 <- estimateTrendedDisp(d2)
d2 <- estimateTagwiseDisp(d2)

d2$counts <- d2$counts[, c(-2, -8)]
d2$pseudo.counts <- d2$pseudo.counts[, c(-2, -8)]
d2$samples <- d2$samples[c(-2, -8),]

#max(table(d2$tagwise.dispersion))
#blacklisted <- which(round(d2$tagwise.dispersion, 5) == 0.00173)
#d2$counts <- d2$counts[-blacklisted,]
#d2$pseudo.counts <- d2$pseudo.counts[-blacklisted,]
#d2$AveLogCPM <- d2$AveLogCPM[-blacklisted]
#d2$trended.dispersion <- d2$trended.dispersion[-blacklisted]
#d2$tagwise.dispersion <- d2$tagwise.dispersion[-blacklisted]

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/MeanVariance_plot.jpg", width = 15, height = 10, units = 'in', res = 300)
plotMeanVar(d2, show.raw.vars=TRUE, show.tagwise.vars=TRUE, xlab = 'Mean expression level (log10 scale)', ylab = 'Pooled gene-level variance (log10 scale)')
dev.off()

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/DispersionMean_plot.jpg", width = 15, height = 10, units = 'in', res = 300)
plotBCV(d2, log = "y")
dev.off()

DEG <- exactTest(d2, pair = c("sensitive", "resistant"))

rstt <- topTags(DEG, n=nrow(d2))
dim(rstt$table[rstt$table$FDR <= 0.2 & rstt$table$logFC > 0,])
dim(rstt$table[rstt$table$FDR <= 0.2 & rstt$table$logFC < 0,])
dim(rstt$table[rstt$table$FDR <= 0.05 & rstt$table$logFC > 1,])
dim(rstt$table[rstt$table$FDR <= 0.05 & rstt$table$logFC < -1,])

den <- rstt$table
den$symbol <- rownames(den)
den <- data.frame(symbol = den[, 5], logFC = den[, 1], logCPM = den[, 2], pvalue = den[, 3], FDR = den[, 4])

write.delim(den, "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_ResistantVsSensitive_DEGs.txt", rownames = TRUE)

den <- rstt$table[rstt$table$FDR <= 0.2,]

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/MAplot_FDR5.jpg", width = 15, height = 10, units = 'in', res = 300)
plotSmear(DEG, de.tags = rownames(den), cex = 0.8, col="gold", deCol = "red3")
abline(h = -0.585, col = "blue", lty = 2, lwd = 1.5)
abline(h = 0.585, col = "blue", lty = 2, lwd = 1.5)
abline(h = 0, lty = 3)
dev.off()


# R

deg <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/edgeR/Provac_ResistantVsSensitive_DEGs.txt", stringsAsFactors = FALSE)

old.up <- read.delim("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Archive/old/Provac_oldSelection_up.list", stringsAsFactors = FALSE, header = FALSE)
old.up <- old.up$V1
old.down <- read.delim("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Archive/old/Provac_oldSelection_down.list", stringsAsFactors = FALSE, header = FALSE)
old.down <- old.down$V1

sig.level <- 0.05
magnitude <- 0.585

up <- deg$symbol[deg$FDR <= sig.level & deg$logFC > magnitude]
down <- deg$symbol[deg$FDR <= sig.level & deg$logFC < -1*magnitude]
length(up)
length(down)

up[up %in% old.up]
100*length(up[up %in% old.up])/length(old.up)

down[down %in% old.down]
100*length(down[down %in% old.down])/length(old.down)



#### DESeq2 analysis


library(DESeq2)

meta <- samples[c(-2, -8),]

counts.ds2 <- as.data.frame(counts)
counts.ds2 <- data.frame(symbol = rownames(counts.ds2), counts.ds2)
rownames(counts.ds2) <- NULL
head(counts.ds2)

ds2 <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~condition,
                              tidy = TRUE)
ds2 <- DESeq(ds2)
deg <- results(ds2, tidy = TRUE)


# R

#library(DESeq2)

#load("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/DESeq2/DESeq2_results.RData")

#ds2 <- DESeq(ds2.object)
#deg <- results(ds2)

#deg <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/DESeq2/DESeq2_results_7R7S.txt", stringsAsFactors = FALSE)
deg <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/DESeq2/DESeq2_results_6R6S.txt", stringsAsFactors = FALSE)
colnames(deg)[c(1, 3)] <- c("symbol", "logFC")

old.up <- read.delim("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Archive/old/Provac_oldSelection_up.list", stringsAsFactors = FALSE, header = FALSE)
old.up <- old.up$V1
old.down <- read.delim("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Archive/old/Provac_oldSelection_down.list", stringsAsFactors = FALSE, header = FALSE)
old.down <- old.down$V1

# Remove not-expressed genes
deg <- na.omit(deg)
dim(deg)

# Resistant vs. Sensitive
deg$logFC <- -1*deg$logFC
deg$stat <- -1*deg$stat
head(deg)


# Descriptives


dim(deg[deg$padj <= 0.2 & deg$logFC > 0,])
dim(deg[deg$padj <= 0.2 & deg$logFC < 0,])
dim(deg[deg$padj <= 0.05 & deg$logFC > 1,])
dim(deg[deg$padj <= 0.05 & deg$logFC < -1,])
dim(deg[deg$padj <= 0.01 & deg$logFC > 1,])
dim(deg[deg$padj <= 0.01 & deg$logFC < -1,])
dim(deg[deg$padj <= 0.005 & deg$logFC > 1,])
dim(deg[deg$padj <= 0.005 & deg$logFC < -1,])


sig.level <- 0.05
magnitude <- 0.585

up <- deg$symbol[deg$padj <= sig.level & deg$logFC > magnitude]
down <- deg$symbol[deg$padj <= sig.level & deg$logFC < -1*magnitude]

length(up)
length(down)

up[up %in% old.up]
100*length(up[up %in% old.up])/length(old.up)

down[down %in% old.down]
100*length(down[down %in% old.down])/length(old.down)


top <- 5

deg <- deg[order(deg$padj),]

top.up <- head(deg[deg$logFC > 0,], top)$symbol
top.down <- head(deg[deg$logFC < 0,], top)$symbol

top.up
top.up[top.up %in% old.up]

top.down
top.down[top.down %in% old.down]



## Predictive models

# R

library(randomForest)
library(e1071)
library(CMA)
library(OptimalCutpoints)
#library(boot)
#library(mosaic)
#library(parallel)
#library(robustbase)
library(ggplot2)
library(reshape2)
#source("~/euNet/Likelihood_functions.R")
source("~/euNet/euNet_core/rutils.R")

randomString <- function(L, N = 1, alphabet = "generic", seed = "None", sep = "")
{
	if (seed != "None") {
		set.seed(as.numeric(seed))
	}
	if (alphabet == "generic") {
		alphabet <- c(0:9, letters, LETTERS)
	} else {
		alphabet <- strsplit(alphabet, sep)[[1]]
	}
	
	K <- c(1:N)
	for (i in 1:N) {
		K[i] <- paste(sample(alphabet, L, replace = TRUE),
		              collapse = "")
	}
	print(table(K))
	cat("\n")
	return(K)
}

#tset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_extendedDataset_N43.txt", stringsAsFactors = FALSE)
tset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_extendedDataset_N43_NOV2021.txt", stringsAsFactors = FALSE)
vset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/CellLines12_RTqPCR_42DEGs.txt", stringsAsFactors = FALSE)

names(tset)[c(19, 23)] <- c("LINC01816", "NOL4L")

data <- tset[, -1]
data$y <- c(rep(0, 25), rep(1, 19))

#library(ComplexHeatmap)
#jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Subjects_42DEGs.jpg", width = 15, height = 10, units = 'in', res = 300)
#colors <- list(Phenotype = c("0" = "lightblue", "1" = "gold"))
#hann <- HeatmapAnnotation(Phenotype = data$y, col = colors)
#Heatmap(t(as.matrix(data[, 1:42])), name = "RTqPCR",
#        top_annotation = hann)
#dev.off()

L <- GenerateLearningsets(n = nrow(data), data$y,
	                      method = "CV",
	                      fold = 4,
	                      strat = TRUE)

C <- classification(as.matrix(data[, -43]), data$y,
	                learningsets = L,
	                classifier = "rfCMA")

E <- evaluation(C, scheme = "observationwise", measure = "brier score")

data$E <- E@score
data$sample <- tset$sample
data$sample[data$E > 1]

data <- x
data$K <- randomString(1, nrow(data), alphabet = 'ABCD')
#write.delim(data$K, "/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/randomKey.txt")
#keys <- read.delim("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/randomKey.txt", stringsAsFactors = FALSE)
#data$K <- keys$x

data <- data[data$E <= 1,]

T1 <- data[data$K == "A" | data$K == "B" | data$K == "C",]
T2 <- data[data$K == "A" | data$K == "B" | data$K == "D",]
T3 <- data[data$K == "A" | data$K == "C" | data$K == "D",]
T4 <- data[data$K == "B" | data$K == "C" | data$K == "A",]

set.seed(0)

n <- 5000
#n <- 10000

F1 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T1, ntree = n, importance = TRUE, mtry = 3)
F2 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T2, ntree = n, importance = TRUE, mtry = 3)
F3 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T3, ntree = n, importance = TRUE, mtry = 3)
F4 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T4, ntree = n, importance = TRUE, mtry = 3)

#F1 <- randomForest(as.factor(y) ~ RNF24 + SLC15A3 + UQCC1 + PLCG2 + GNG11 + TACC1 + CTNNBL1 + CKB + ABDHD5,
#                   data = T1, ntree = n, importance = TRUE, mtry = 4)
#F2 <- randomForest(as.factor(y) ~ RNF24 + SLC15A3 + UQCC1 + PLCG2 + GNG11 + TACC1 + CTNNBL1 + CKB + ABDHD5,
#                   data = T2, ntree = n, importance = TRUE, mtry = 4)
#F3 <- randomForest(as.factor(y) ~ RNF24 + SLC15A3 + UQCC1 + PLCG2 + GNG11 + TACC1 + CTNNBL1 + CKB + ABDHD5,
#                   data = T3, ntree = n, importance = TRUE, mtry = 4)
#F4 <- randomForest(as.factor(y) ~ RNF24 + SLC15A3 + UQCC1 + PLCG2 + GNG11 + TACC1 + CTNNBL1 + CKB + ABDHD5,
#                   data = T4, ntree = n, importance = TRUE, mtry = 4)


V1 <- data[data$K == "D",]
V2 <- data[data$K == "C",]
V3 <- data[data$K == "B",]
V4 <- data[data$K == "A",]


R1 <- data.frame(importance(F1))[, 3:4]
R2 <- data.frame(importance(F2))[, 3:4]
R3 <- data.frame(importance(F3))[, 3:4]
R4 <- data.frame(importance(F4))[, 3:4]
R <- cbind(R1, R2, R3, R4)
R$MDA <- apply(R[, c(1, 3, 5, 7)], 1, mean)
R$MDG <- apply(R[, c(2, 4, 6, 8)], 1, mean)
R <- R[, c(9, 10)]
#Fx <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
#                   data = data, ntree = 10000, importance = TRUE)
#R <- data.frame(importance(Fx))[, 3:4]
#names(R) <- c("MDA", "MDG")

R$fA <- 100*(R$MDA - min(R$MDA))/(max(R$MDA) - min(R$MDA))
R$fG <- 100*(R$MDG - min(R$MDG))/(max(R$MDG) - min(R$MDG))
R$f <- apply(R[, c(3, 4)], 1, mean)
#R[order(R$MDA, decreasing = TRUE),]
#R[order(R$MDG, decreasing = TRUE),]
R <- R[order(R$f, decreasing = TRUE),]
R


C1 <- predict(F1, V1)
C2 <- predict(F2, V2)
C3 <- predict(F3, V3)
C4 <- predict(F4, V4)


TP <- rbind(V1[C1 == V1$y & V1$y == 1,], V2[C2 == V2$y & V2$y == 1,], V3[C3 == V3$y & V3$y == 1,], V4[C4 == V4$y & V4$y == 1,])
TN <- rbind(V1[C1 == V1$y & V1$y == 0,], V2[C2 == V2$y & V2$y == 0,], V3[C3 == V3$y & V3$y == 0,], V4[C4 == V4$y & V4$y == 0,])
FP <- rbind(V1[C1 == 1 & V1$y == 0,], V2[C2 == 1 & V2$y == 0,], V3[C3 == 1 & V3$y == 0,], V4[C4 == 1 & V4$y == 0,])
FN <- rbind(V1[C1 == 0 & V1$y == 1,], V2[C2 == 0 & V2$y == 1,], V3[C3 == 0 & V3$y == 1,], V4[C4 == 0 & V4$y == 1,])

confusion <- matrix(c(nrow(TP), nrow(FP), nrow(FN), nrow(TN)), nrow = 2, byrow = TRUE)
rownames(confusion) <- c("Pred1", "Pred0")
colnames(confusion) <- c("Real1", "Real0")
confusion
P <- performance(confusion)
P


if (FALSE) {
j <- as.numeric(rownames(tset)[!(rownames(tset) %in% rownames(data))])
tset$sample[j]
data <- cbind(sample = tset$sample[-j], data)

rfc.provac <- list(data = data,
                   models = list(RFC1 = F1, RFC2 = F2, RFC3 = F3, RFC4 = F4),
                   rank = R,
                   confusion = confusion,
                   performance = P)

save("rfc.provac", file = "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC.RData")
}


# R

library(randomForest)
library(e1071)
library(CMA)
library(OptimalCutpoints)
#library(boot)
#library(mosaic)
#library(parallel)
#library(robustbase)
library(ggplot2)
library(reshape2)
#source("~/euNet/Likelihood_functions.R")
source("~/euNet/euNet_core/rutils.R")

seedTest <- function(x0, x1, mu = 0.5, alpha = 0.05) {
	y <- vector()
	for (j in 1:ncol(x0)) {
		m0 <- median(x0[, j])
		m0 <- ifelse(m0 >= 0, 1, -1)
		m1 <- median(x1[, j])
		m1 <- ifelse(m1 >= 0, 1, -1)
		#One subject was excluded, given their Brier score > 1.
		if (sign(m0) != sign(m1)) {
			y <- c(y, FALSE)
		} else {
			w0 <- wilcox.test(x0[, j], mu = m0*mu, conf.int = TRUE)
			w0$P <- ifelse(w0$p.value <= alpha, 0, 1)
			w1 <- wilcox.test(x1[, j], mu = m0*mu, conf.int = TRUE)
			w1$P <- ifelse(w1$p.value <= alpha, 0, 1)
			if (sign(w0$estimate) != sign(w1$estimate)) {
				y <- c(y, FALSE)
			} else if (w0$P != w1$P) {
				y <- c(y, FALSE)
			} else {
				y <- c(y, TRUE)
			}
		}
	}
	return(y)
}

tset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_extendedDataset_N43.txt", stringsAsFactors = FALSE)
names(tset)[c(19, 23)] <- c("LINC01816", "NOL4L")

vset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/CellLines12_RTqPCR_42DEGs.txt", stringsAsFactors = FALSE)
vset$cellLine <- c("HOSE", "FTI194", "COV318", "HEY", "OV90", "PEO1",
                   "OVCAR3", "OVGEM2", "OVGEM6", "OVGEM8", "OVGEM9",
                   "OVGEM11", "OVGEM16", "OVGEM20")

if (FALSE) {
vset$E <- c(0.5, 0.5, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0)
vset$y <- c(0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0)

vset$K <- randomString(1, nrow(vset), alphabet = 'ABCD')

T1 <- vset[vset$K == "A" | vset$K == "B" | vset$K == "C",]
T2 <- vset[vset$K == "A" | vset$K == "B" | vset$K == "D",]
T3 <- vset[vset$K == "A" | vset$K == "C" | vset$K == "D",]
T4 <- vset[vset$K == "B" | vset$K == "C" | vset$K == "A",]

#set.seed(0)

n <- 10000
m <- 3

F1 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T1, ntree = n, importance = TRUE, mtry = m)
F2 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T2, ntree = n, importance = TRUE, mtry = m)
F3 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T3, ntree = n, importance = TRUE, mtry = m)
F4 <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
                   data = T4, ntree = n, importance = TRUE, mtry = m)

R1 <- data.frame(importance(F1))[, 3:4]
R2 <- data.frame(importance(F2))[, 3:4]
R3 <- data.frame(importance(F3))[, 3:4]
R4 <- data.frame(importance(F4))[, 3:4]
R <- cbind(R1, R2, R3, R4)
R$MDA <- apply(R[, c(1, 3, 5, 7)], 1, mean)
R$MDG <- apply(R[, c(2, 4, 6, 8)], 1, mean)
R <- R[, c(9, 10)]

R$fA <- 100*(R$MDA - min(R$MDA))/(max(R$MDA) - min(R$MDA))
R$fG <- 100*(R$MDG - min(R$MDG))/(max(R$MDG) - min(R$MDG))
R$f <- apply(R[, c(3, 4)], 1, mean)
R <- R[order(R$f, decreasing = TRUE),]
R

V1 <- vset[vset$K == "D",]
V2 <- vset[vset$K == "C",]
V3 <- vset[vset$K == "B",]
V4 <- vset[vset$K == "A",]

C1 <- predict(F1, V1)
C2 <- predict(F2, V2)
C3 <- predict(F3, V3)
C4 <- predict(F4, V4)

TP <- rbind(V1[C1 == V1$y & V1$y == 1,], V2[C2 == V2$y & V2$y == 1,], V3[C3 == V3$y & V3$y == 1,], V4[C4 == V4$y & V4$y == 1,])
TN <- rbind(V1[C1 == V1$y & V1$y == 0,], V2[C2 == V2$y & V2$y == 0,], V3[C3 == V3$y & V3$y == 0,], V4[C4 == V4$y & V4$y == 0,])
FP <- rbind(V1[C1 == 1 & V1$y == 0,], V2[C2 == 1 & V2$y == 0,], V3[C3 == 1 & V3$y == 0,], V4[C4 == 1 & V4$y == 0,])
FN <- rbind(V1[C1 == 0 & V1$y == 1,], V2[C2 == 0 & V2$y == 1,], V3[C3 == 0 & V3$y == 1,], V4[C4 == 0 & V4$y == 1,])

confusion <- matrix(c(nrow(TP), nrow(FP), nrow(FN), nrow(TN)), nrow = 2, byrow = TRUE)
rownames(confusion) <- c("Pred1", "Pred0")
colnames(confusion) <- c("Real1", "Real0")
confusion
P <- performance(confusion)
P
}

load("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC.RData")
#data <- rfc.provac$data
data <- tset[, -1]

R <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RFC_rank.txt", stringsAsFactors = FALSE)
R$P <- 100*(R$padj - min(R$padj))/(max(R$padj) - min(R$padj))
R$fP <- apply(R[, c(6, 8)], 1, mean)
R <- R[order(R$fP, decreasing = TRUE),]
head(R)

V <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RFC_vrank.txt", stringsAsFactors = FALSE)
head(V)

mda <- 0
#mda <- median(rfc.provac$rank$MDA)
#seeds <- rownames(rfc.provac$rank)[rfc.provac$rank$MDA > mda]
#weeds <- V$symbol[V$MDA > mda]

fs <- quantile(rfc.provac$rank$fG)[3]
fv <- quantile(V$fG)[3]
seeds <- rownames(rfc.provac$rank)[rfc.provac$rank$fG > fs & rfc.provac$rank$MDA > mda]   # 19
weeds <- V$symbol[V$fG > fv & V$MDA > mda]   # 11

weeds[weeds %in% seeds]            # 7 (MYO5A, C3AR1, WDPCP, LDLR, TACC1, TNFSF13B, RNF24)
seeds <- weeds[weeds %in% seeds]   # 5 (MYO5A, C3AR1, WDPCP, TNFSF13B, RNF24)
#seeds <- weeds

seeds.old <- c("CTNNBL1", "IGSF9", "NOL4L", "RNF24", "UQCC1",   # UP
               "GNG11", "IGFBP7", "LDLR", "PLCG2", "SLC15A3")   # DOWN
seeds
seeds.old[seeds.old %in% seeds]

#Y <- vset[, c(-1, -2)]
Y <- vset[, colnames(vset) %in% seeds]
rownames(Y) <- vset$cellLine

#x0 <- apply(data, 2, mean)
#x1 <- apply(Y, 2, mean)
#x0 <- x0[names(x0) %in% names(x1)]
#x0 <- x0[order(names(x0))]
#x1 <- x1[order(names(x1))]

data <- data[, colnames(data) %in% colnames(Y)]
data <- data[, order(colnames(data))]
Y <- Y[, order(colnames(Y))]

include <- seedTest(data, Y)
data <- data[, include]
Y <- Y[, include]
seeds <- colnames(Y)
seeds

D <- dist(Y)
h <- hclust(D)
#h <- hclust(D, "ave")
#jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/CellLines_hclust_avg.jpg", width = 20, height = 10, units = 'in', res = 300)
plot(h)
#dev.off()

data$y <- c(rep(0, 25), rep(1, 18))

Fy <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = data, ntree = 10000, importance = TRUE)
predict(Fy, Y)

seeds <- colnames(Y)
Y$y <- c(0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0)
Y$K <- randomString(1, nrow(Y), alphabet = 'ABCD')

T1 <- Y[Y$K == "A" | Y$K == "B" | Y$K == "C",]
T2 <- Y[Y$K == "A" | Y$K == "B" | Y$K == "D",]
T3 <- Y[Y$K == "A" | Y$K == "C" | Y$K == "D",]
T4 <- Y[Y$K == "B" | Y$K == "C" | Y$K == "A",]

set.seed(0)

n <- 10000

F1 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T1, ntree = n, importance = TRUE)
F2 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T2, ntree = n, importance = TRUE)
F3 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T3, ntree = n, importance = TRUE)
F4 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T4, ntree = n, importance = TRUE)

Rv1 <- data.frame(importance(F1))[, 3:4]
Rv2 <- data.frame(importance(F2))[, 3:4]
Rv3 <- data.frame(importance(F3))[, 3:4]
Rv4 <- data.frame(importance(F4))[, 3:4]
Rv <- cbind(Rv1, Rv2, Rv3, Rv4)
Rv$MDA <- apply(Rv[, c(1, 3, 5, 7)], 1, mean)
Rv$MDG <- apply(Rv[, c(2, 4, 6, 8)], 1, mean)
Rv <- Rv[, c(9, 10)]

Rv$fA <- 100*(Rv$MDA - min(Rv$MDA))/(max(Rv$MDA) - min(Rv$MDA))
Rv$fG <- 100*(Rv$MDG - min(Rv$MDG))/(max(Rv$MDG) - min(Rv$MDG))
Rv$f <- apply(Rv[, c(3, 4)], 1, mean)
Rv <- Rv[order(Rv$f, decreasing = TRUE),]
Rv

V1 <- Y[Y$K == "D",]
V2 <- Y[Y$K == "C",]
V3 <- Y[Y$K == "B",]
V4 <- Y[Y$K == "A",]

C1 <- predict(F1, V1)
C2 <- predict(F2, V2)
C3 <- predict(F3, V3)
C4 <- predict(F4, V4)

TP <- rbind(V1[C1 == V1$y & V1$y == 1,], V2[C2 == V2$y & V2$y == 1,], V3[C3 == V3$y & V3$y == 1,], V4[C4 == V4$y & V4$y == 1,])
TN <- rbind(V1[C1 == V1$y & V1$y == 0,], V2[C2 == V2$y & V2$y == 0,], V3[C3 == V3$y & V3$y == 0,], V4[C4 == V4$y & V4$y == 0,])
FP <- rbind(V1[C1 == 1 & V1$y == 0,], V2[C2 == 1 & V2$y == 0,], V3[C3 == 1 & V3$y == 0,], V4[C4 == 1 & V4$y == 0,])
FN <- rbind(V1[C1 == 0 & V1$y == 1,], V2[C2 == 0 & V2$y == 1,], V3[C3 == 0 & V3$y == 1,], V4[C4 == 0 & V4$y == 1,])

confusion <- matrix(c(nrow(TP), nrow(FP), nrow(FN), nrow(TN)), nrow = 2, byrow = TRUE)
rownames(confusion) <- c("Pred1", "Pred0")
colnames(confusion) <- c("Real1", "Real0")
confusion
P <- performance(confusion)
P
#        Real.R Real.S
# Pred.R      6      1
# Pred.S      0      7
#  Se: 1.000
#  Sp: 0.875
# PPV: 0.857
# FPR: 0.125
# FDR: 0.143
#  F1: 0.923
#   A: 0.929

predicted <- c(C1, C2, C3, C4)
predicted <- predicted[order(names(predicted))]
seeds
Y <- Y[order(rownames(Y)),]
Y <- data.frame(Y[, 1:6], yhat = predicted)
Y

library(ComplexHeatmap)

minmax <- function(x, k = 100) {
	v <- k*(x - min(x))/(max(x) - min(x))
	return(v)
}

x <- Y[, 1:5]
#x <- apply(x, 1, minmax)
#rownames(x) <- colnames(Y[, 1:5])

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RTqPCR_markers.jpg", width = 15, height = 10, units = 'in', res = 300)
colors <- list(Experimental = c("0" = "green3", "1" = "orange"),
               Predicted = c("0" = "lightblue", "1" = "gold"))
hann <- HeatmapAnnotation(Predicted = Y$yhat, Experimental = Y$y, col = colors)
Heatmap(t(as.matrix(x)), name = "RTqPCR",
        top_annotation = hann,
        clustering_distance_columns = "maximum",
        clustering_distance_rows = "maximum")
dev.off()




## DEGs heatmaps

# R

library(ggplot2)
library(reshape2)
library(ggrepel)
library(circlize)
library(ComplexHeatmap)
library(randomForest)

source("~/euNet/euNet_core/rutils.R")

load("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC.RData")

exprs <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/DEBrowser/Provac_ResistantVsSensitive_DEGs.txt", stringsAsFactors = FALSE)
pcr <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_extendedDataset_N43.txt", stringsAsFactors = FALSE)

exprs$Regulation <- "NR"
exprs$Regulation[exprs$log2FoldChange > 1 & exprs$padj < 0.05] <- "UP"
exprs$Regulation[exprs$log2FoldChange < -1 & exprs$padj < 0.05] <- "DOWN"
cols <- c("UP" = "red2", "DOWN" = "blue", "NR" = "grey60")

top <- exprs$symbol[exprs$Tier == "Top42"]
top.degs <- exprs[exprs$Tier == "Top42" & exprs$padj < 0.001,]
up <- exprs[exprs$log2FoldChange > 1 & exprs$padj < 0.001,]
down <- exprs[exprs$log2FoldChange < -1 & exprs$padj < 0.001,]

if (FALSE) {
#png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_volcano.png", width = 20, height = 10, units = 'in', res = 450)
pdf("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_volcano.pdf", width = 20, height = 10)
ggplot(data = exprs,
       aes(x = log2FoldChange,
           y = log10padj)) +
  theme_bw() +
  theme(panel.border = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 28, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 28)) +
  geom_point(aes(colour = Regulation),
             alpha = 0.8,
             shape = 16,
             size = 3) +
  scale_colour_manual(values = cols) + 
#  geom_point(data = up,
#             shape = 21,
#             size = 3,
#             fill = "firebrick",
#             colour = "black") +
#  geom_point(data = down,
#             shape = 21,
#             size = 3,
#             fill = "steelblue",
#             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
#  geom_hline(yintercept = 3,
#             linetype = "dotted") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
#  geom_label_repel(data = top.degs,
#                   aes(label = symbol),
#                   force = 2,
#                   nudge_y = 1) +
  annotate("text", x = -3.8, y = 1.42, size = 8,
           label = "P.adj = 0.05",
           parse = FALSE) +
#  annotate("text", x = -4, y = 3.1, size = 5,
#           label = "P.adj = 0.001",
#           parse = FALSE) +
  labs(x = "log2(Resistant/Sensitive)", y = "-log10(BH-adjusted P-value)") +
  scale_x_continuous(breaks = c(seq(-10, 10, 1)),     
                     limits = c(-4, 4))
dev.off()
}

markers <- c("C3AR1", "TNFSF13B", "WDPCP", "MYO5A", "RNF24")

x <- log2(exprs[, 2:15] + 1)
#x <- x[exprs$symbol %in% top, -13]
x <- x[exprs$symbol %in% top,]
x <- t(apply(x, 1, scale))
rownames(x) <- exprs$symbol[exprs$symbol %in% top]
#colnames(x) <- colnames(exprs[, c(2:13, 15)])
colnames(x) <- colnames(exprs[, c(2:15)])
x <- x[, c(1:3, 7, 4:6, 8:10, 14, 11:13)]
colnames(x) <- c("1", "2", "3", "4", "5", "6", "7",
                 "8", "9", "10", "11", "12", "13", "14")
head(x)

if (FALSE) {
#png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_Top42.png", width = 15, height = 10, units = 'in', res = 450)
pdf("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_Top42.pdf", width = 15, height = 10)
colors <- list(Phenotype = c("Sensitive" = "lightblue", "Resistant" = "gold"))
hann <- HeatmapAnnotation(Phenotype = c(rep("Sensitive", 7), rep("Resistant", 7)),
                          col = colors,
                          show_legend = c(FALSE, FALSE))
#col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("red", "gold", "white", "deepskyblue", "dodgerblue"))
col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("blue", "blue", "white", "red", "red"))
Heatmap(as.matrix(x), name = "Standardized\nlog2(counts)",
        top_annotation = hann,
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 24),
        column_names_rot = 0,
        heatmap_legend_param = list(title_gp = gpar(col = "black", fontsize = 16),
                                    labels_gp = gpar(col = "black", fontsize = 16),
                                    direction = "horizontal"))
pheno <- Legend(labels = c("Sensitive", "Resistant"),
                legend_gp = gpar(fill = c("lightblue", "gold")),
                title = "Phenotype",
                title_gp = gpar(col = "black", fontsize = 16),
                labels_gp = gpar(col = "black", fontsize = 16),
                direction = "horizontal")
draw(pheno, x = unit(36, "cm"), y = unit(15, "cm"))
dev.off()
}

x1 <- x[rownames(x) %in% markers,]

if (FALSE) {
jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_markers.jpg", width = 15, height = 10, units = 'in', res = 300)
colors <- list(Phenotype = c("Sensitive" = "lightblue", "Resistant" = "gold"))
hann <- HeatmapAnnotation(Phenotype = c(rep("Sensitive", 7), rep("Resistant", 6)), col = colors)
Heatmap(as.matrix(x1), name = "Standardized\nlog2(counts)",
        top_annotation = hann,
        clustering_distance_columns = "maximum",
        clustering_distance_rows = "maximum")
dev.off()

# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(cyl = c("4" = "green", "6" = "gray", "8" = "darkred"),
            am = c("0" = "yellow", "1" = "orange"),
            mpg = circlize::colorRamp2(c(17, 25), 
                                       c("lightblue", "purple")) )
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  cyl = mtcars$cyl, am = mtcars$am, mpg = mtcars$mpg,
  col = col
)

# Combine the heatmap and the annotation
Heatmap(df, name = "mtcars",
        top_annotation = ha)

q <- as.matrix(pcr[, 2:43])
rownames(q) <- pcr$sample
colnames(q)[c(18, 22)] <- c("LINC01816", "NOL4L")
z <- t(apply(pcr[, 2:43], 1, scale))
rownames(z) <- pcr$sample
colnames(z) <- colnames(pcr[, 2:43])
head(z)
colnames(z)[c(18, 22)] <- c("LINC01816", "NOL4L")

#clustering_distance_rows = "maximum",
#clustering_distance_columns = "maximum",

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RTqPCR_human.jpg", width = 15, height = 10, units = 'in', res = 300)
colors <- list(Phenotype = c("Sensitive" = "lightblue", "Resistant" = "gold"))
hann <- HeatmapAnnotation(Phenotype = c(rep("Sensitive", 25), rep("Resistant", 18)), col = colors)
Heatmap(t(as.matrix(z)), name = "Standardized\nRTqPCR",
        top_annotation = hann)
dev.off()

top10 <- c(up$symbol, down$symbol)

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RTqPCR_human_top10_DEGs.jpg", width = 15, height = 10, units = 'in', res = 300)
colors <- list(Phenotype = c("Sensitive" = "lightblue", "Resistant" = "gold"))
hann <- HeatmapAnnotation(Phenotype = c(rep("Sensitive", 25), rep("Resistant", 18)), col = colors)
Heatmap(t(as.matrix(z[, colnames(z) %in% top10])), name = "Standardized\nRTqPCR",
        top_annotation = hann)
dev.off()
}

Q <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_avgExprs.txt", stringsAsFactors = FALSE)
Z <- apply(Q[, 2:5], 2, scale)
rownames(Z) <- Q$symbol
Z[, 3] <- -1*Z[, 3]
Z[, 4] <- -1*Z[, 4]
cor(Z[, 1], Z[, 3])
cor(Z[, 2], Z[, 4])
Z <- data.frame(symbol = rownames(Z),
                RNAseq = c(Z[, 1], Z[, 2]),
                RTqPCR = c(Z[, 3], Z[, 4]),
                Condition = c(rep("sensitive", 42), rep("resistant", 42)))
head(Z)

sen <- lm(RTqPCR ~ RNAseq, data = Z[Z$Condition == "SENSITIVE",])
summary(sen)

res <- lm(RTqPCR ~ RNAseq, data = Z[Z$Condition == "RESISTANT",])
summary(res)

if (FALSE) {
jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_RTqPCR_Top42_corr_overlayed.jpg", width = 15, height = 10, units = 'in', res = 300)
ggplot(Z, aes(x = RNAseq, y = RTqPCR, color = Condition, shape = Condition)) +
geom_point(size = 3) +
geom_smooth(method = lm, se = FALSE, fullrange = TRUE, lty = 1, lwd = 0.8) +
scale_colour_manual(values = c("red3", "royalblue")) +
scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 16)) +
annotate("text", x = 1.5, y = -1.6, size = 6,
         label = "Pearson correlation\nSensitive: 0.519\nResistant: 0.532",
         parse = FALSE)
dev.off()

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_RTqPCR_Top42_corr.jpg", width = 15, height = 10, units = 'in', res = 300)
txt <- data.frame(Condition = c("SENSITIVE", "RESISTANT"),
                  RNAseq = c(1, 1),
                  RTqPCR = c(-2.5, -2.5))
lbl <- ifelse(Z$symbol %in% c("IGSF9", "DSN1", "WDPCP", "CTNNBL1",
                              "UQCC1", "TMEM140", "CPZ", "TPP1",
                              "IGFBP7", "PLCG2"), Z$symbol, "")
Z$Condition[Z$Condition == "sensitive"] <- "SENSITIVE"
Z$Condition[Z$Condition == "resistant"] <- "RESISTANT"
ggplot(Z, aes(x = RNAseq, y = RTqPCR, color = Condition, shape = Condition)) +
geom_point(size = 4) +
geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
            lty = 1, lwd = 0.8) +
scale_colour_manual(values = c("red3", "royalblue")) +
scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16)) +
#annotate("text", x = 1.5, y = -2.1, size = 5,
#         label = "Pearson correlation\nSensitive: 0.519\nResistant: 0.532",
#         parse = FALSE) +
facet_wrap(facets =  vars(Condition)) +
geom_text(data = txt,
          label = c("Pearson corr. (beta) = 0.519\n95% C.I. = (0.255, 0.784)\nP-value = 4.25e-04",
                    "Pearson corr. (beta) = 0.529\n95% C.I. = (0.266, 0.792)\nP-value = 3.18e-04"),
          size = 5,
          show.legend = FALSE) +
geom_text(label = lbl, nudge_y = 0.1, show.legend = FALSE)
dev.off()

png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_RTqPCR_Top42_corr_noLabels.png", width = 15, height = 10, units = 'in', res = 450)
txt <- data.frame(Condition = c("SENSITIVE", "RESISTANT"),
                  RNAseq = c(1, 1),
                  RTqPCR = c(-2.5, -2.5))
lbl <- ifelse(Z$symbol %in% c("IGSF9", "DSN1", "WDPCP", "CTNNBL1",
                              "UQCC1", "TMEM140", "CPZ", "TPP1",
                              "IGFBP7", "PLCG2"), Z$symbol, "")
Z$Condition[Z$Condition == "sensitive"] <- "SENSITIVE"
Z$Condition[Z$Condition == "resistant"] <- "RESISTANT"
ggplot(Z, aes(x = RNAseq, y = RTqPCR, color = Condition, shape = Condition)) +
geom_point(size = 4) +
geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
            lty = 1, lwd = 0.8) +
scale_colour_manual(values = c("darkblue", "darkblue")) +
scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 20),
  axis.title = element_text(size = 22, face = "bold"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 18)) +
#annotate("text", x = 1.5, y = -2.1, size = 5,
#         label = "Pearson correlation\nSensitive: 0.519\nResistant: 0.532",
#         parse = FALSE) +
facet_wrap(facets =  vars(Condition)) +
geom_text(data = txt,
          label = c("Pearson corr. (beta) = 0.519\n95% C.I. = (0.255, 0.784)\nP-value = 4.25e-04",
                    "Pearson corr. (beta) = 0.529\n95% C.I. = (0.266, 0.792)\nP-value = 3.18e-04"),
          size = 6.5,
          show.legend = FALSE) +
ylab("HT RT-qPCR")
#+
#geom_text(label = lbl, nudge_y = 0.1, show.legend = FALSE)
dev.off()
}



# WILCOXON VALIDATION

# R
library(randomForest)
library(ggplot2)
source("~/euNet/euNet_core/rutils.R")

# save(rfc.provac, file = "/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC2.RData")
load("/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC2.RData")

#top.rfc <- rownames(rfc.provac$rank)[rfc.provac$rank$MDA > 0 & rfc.provac$rank$fG > median(rfc.provac$rank$fG)]
#top.rfc <- rownames(rfc.provac$rank)[rfc.provac$rank$MDA > 0 & rfc.provac$rank$fG > 34]

data <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_extendedDataset_N43_NOV2021.txt", stringsAsFactors = FALSE)
#data <- rfc.provac$data

#write.delim(Q, "/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_avgExprs.txt")

W <- data.frame(symbol = vector(), W = vector(), estimate = vector(), ci95 = vector(), pvalue = vector())
for (j in 2:43) {
	w <- wilcox.test(data[1:25, j], data[26:44, j], conf.int = TRUE)
	W <- rbind(W, data.frame(symbol = colnames(data)[j],
	                         W = w$statistic,
	                         estimate = w$estimate,
	                         ci95 = paste0(w$conf.int, collapse = ", "),
	                         pvalue = w$p.value))
}

#write.delim(W, "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RTqPCR_wilcoxon.txt")

seeds <- W$symbol[W$pvalue < 0.05]
#seeds <- seeds[seeds %in% rownames(rfc.provac$rank)[rfc.provac$rank$f > 40]]
seeds

#save("rfc10.provac", file = "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC10.RData")
load("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC.RData")
load("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC10.RData")

if (FALSE) {
model <- formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + ")))

data <- rfc.provac$data
data$K <- randomString(1, nrow(data), alphabet = 'ABCD')

T1 <- data[data$K == "A" | data$K == "B" | data$K == "C",]
T2 <- data[data$K == "A" | data$K == "B" | data$K == "D",]
T3 <- data[data$K == "A" | data$K == "C" | data$K == "D",]
T4 <- data[data$K == "B" | data$K == "C" | data$K == "A",]

set.seed(12345)

n <- 5000
#n <- 10000

F1 <- randomForest(model, data = T1, ntree = n, importance = TRUE, mtry = 3)
F2 <- randomForest(model, data = T2, ntree = n, importance = TRUE, mtry = 3)
F3 <- randomForest(model, data = T3, ntree = n, importance = TRUE, mtry = 3)
F4 <- randomForest(model, data = T4, ntree = n, importance = TRUE, mtry = 3)

V1 <- data[data$K == "D",]
V2 <- data[data$K == "C",]
V3 <- data[data$K == "B",]
V4 <- data[data$K == "A",]

R1 <- data.frame(importance(F1))[, 3:4]
R2 <- data.frame(importance(F2))[, 3:4]
R3 <- data.frame(importance(F3))[, 3:4]
R4 <- data.frame(importance(F4))[, 3:4]
R <- cbind(R1, R2, R3, R4)
R$MDA <- apply(R[, c(1, 3, 5, 7)], 1, mean)
R$MDG <- apply(R[, c(2, 4, 6, 8)], 1, mean)
R <- R[, c(9, 10)]
#Fx <- randomForest(as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + UQCC1 + WDPCP + ZNF738,
#                   data = data, ntree = 10000, importance = TRUE)
#R <- data.frame(importance(Fx))[, 3:4]
#names(R) <- c("MDA", "MDG")

R$fA <- 100*(R$MDA - min(R$MDA))/(max(R$MDA) - min(R$MDA))
R$fG <- 100*(R$MDG - min(R$MDG))/(max(R$MDG) - min(R$MDG))
R$f <- apply(R[, c(3, 4)], 1, mean)
#R[order(R$MDA, decreasing = TRUE),]
#R[order(R$MDG, decreasing = TRUE),]
R <- R[order(R$f, decreasing = TRUE),]
R

C1 <- predict(F1, V1)
C2 <- predict(F2, V2)
C3 <- predict(F3, V3)
C4 <- predict(F4, V4)

TP <- rbind(V1[C1 == V1$y & V1$y == 1,], V2[C2 == V2$y & V2$y == 1,], V3[C3 == V3$y & V3$y == 1,], V4[C4 == V4$y & V4$y == 1,])
TN <- rbind(V1[C1 == V1$y & V1$y == 0,], V2[C2 == V2$y & V2$y == 0,], V3[C3 == V3$y & V3$y == 0,], V4[C4 == V4$y & V4$y == 0,])
FP <- rbind(V1[C1 == 1 & V1$y == 0,], V2[C2 == 1 & V2$y == 0,], V3[C3 == 1 & V3$y == 0,], V4[C4 == 1 & V4$y == 0,])
FN <- rbind(V1[C1 == 0 & V1$y == 1,], V2[C2 == 0 & V2$y == 1,], V3[C3 == 0 & V3$y == 1,], V4[C4 == 0 & V4$y == 1,])

confusion <- matrix(c(nrow(TP), nrow(FP), nrow(FN), nrow(TN)), nrow = 2, byrow = TRUE)
rownames(confusion) <- c("Pred1", "Pred0")
colnames(confusion) <- c("Real1", "Real0")
confusion
P <- performance(confusion)
P

rfc10.provac <- list(data = data,
                     models = list(RFC1 = F1, RFC2 = F2, RFC3 = F3, RFC4 = F4),
                     rank = R,
                     confusion = confusion,
                     performance = P)
}


# R
library(randomForest)
library(ggplot2)
source("~/euNet/euNet_core/rutils.R")

W <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RNAseq_RTqPCR_summary.txt", stringsAsFactors = FALSE)
W$Regulation <- "NR"
W$Regulation[W$wPvalue <= 0.05 & W$shift > 0] <- "UP"
W$Regulation[W$wPvalue <= 0.05 & W$shift < 0] <- "DOWN"
head(W)

cols <- c("UP" = "red2", "DOWN" = "blue", "NR" = "grey60")

if (FALSE) {
png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RTqPCR_volcano.png", width = 20, height = 10, units = 'in', res = 450)
lbl <- ifelse(W$symbol %in% c("GNG11", "RNF24", "UQCC1", "CTNNBL1",
                              "PLCG2", "TTI1", "IGFBP7", "SLC15A3",
                              "TSPAN31", "CKB"), W$symbol, "")
ggplot(data = W,
       aes(x = shift,
           y = -log10(wPvalue))) +
  theme_bw() +
  theme(panel.border = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26)) +
  geom_point(aes(colour = Regulation),
             alpha = 0.8,
             shape = 16,
             size = 5) +
  scale_colour_manual(values = cols) + 
#  geom_point(data = up,
#             shape = 21,
#             size = 3,
#             fill = "firebrick",
#             colour = "black") +
#  geom_point(data = down,
#             shape = 21,
#             size = 3,
#             fill = "steelblue",
#             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
#  geom_hline(yintercept = 3,
#             linetype = "dotted") +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
#  geom_label_repel(data = top.degs,
#                   aes(label = symbol),
#                   force = 2,
#                   nudge_y = 1) +
  annotate("text", x = -1.2, y = 1.35, size = 7,
           label = "P-value = 0.05",
           parse = FALSE) +
#  annotate("text", x = -4, y = 3.1, size = 5,
#           label = "P.adj = 0.001",
#           parse = FALSE) +
  geom_text(label = lbl, cex = 6, nudge_x = -0.04, nudge_y = 0.06, show.legend = FALSE) +
  labs(x = "Wilcoxon estimated shift", y = "-log10(P-value)") +
  scale_x_continuous(breaks = c(seq(-2, 2, 0.25)),
                     limits = c(-1.25, 1.25))
dev.off()
}


# R

library(randomForest)
library(e1071)
library(CMA)
library(OptimalCutpoints)
#library(boot)
#library(mosaic)
#library(parallel)
#library(robustbase)
library(ggplot2)
library(reshape2)
library(SEMgraph)
library(SEMdata)
library(huge)
library(org.Hs.eg.db)
#source("~/euNet/Likelihood_functions.R")
source("~/euNet/euNet_core/rutils.R")

vset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/CellLines12_RTqPCR_42DEGs.txt", stringsAsFactors = FALSE)
vset$cellLine <- c("HOSE", "FTI194", "COV318", "HEY", "OV90", "PEO1",
                   "OVCAR3", "OVGEM2", "OVGEM6", "OVGEM8", "OVGEM9",
                   "OVGEM11", "OVGEM16", "OVGEM20")
v <- vset[8:14,]
v$Paclitaxel <- c(0, 0, 1, 0, 0, 1, 0)
v$CisplatinIC150uM <- c(1, 1, 1, 0, 1, 1, 0)
v$Carboplatin <- c(0, 1, 1, 0, 1, 0, 0)
v$cispty <- c(0, 1, 1, 0, 1, 1, 0)

vp0 <- v[v$Paclitaxel == 0, 3:44]
vp1 <- v[v$Paclitaxel == 1, 3:44]
vi0 <- v[v$CisplatinIC150uM == 0, 3:44]
vi1 <- v[v$CisplatinIC150uM == 1, 3:44]
vc0 <- v[v$Carboplatin == 0, 3:44]
vc1 <- v[v$Carboplatin == 1, 3:44]

V42 <- colnames(v)[3:44]
V42[1] <- "ABHD5"
V42 <- mapIds(org.Hs.eg.db, V42, column = "ENTREZID", keytype = "SYMBOL")
#pcr <- v[3:44]
pcr <- huge.npn(v[3:44])
colnames(pcr) <- V42
#G42k <- activeModule(kegg, type = "kou", seed = V42)
#G42r <- activeModule(reactome, type = "kou", seed = V42)
#G42s <- activeModule(string, type = "kou", seed = V42)
#G42 <- graph.union(G42k, G42r, as.directed(G42s))
r <- cor(pcr)
G42 <- corr2graph(r, n = nrow(pcr), type = "tmfg")

pactx.fit <- SEMrun(G42, pcr, v$Paclitaxel, algo = "ricf")
pactx.deg <- pactx.fit$gest
pactx.deg <- data.frame(Symbol = mapIds(org.Hs.eg.db, rownames(pactx.deg), column = "SYMBOL", keytype = "ENTREZID"), pactx.deg)
write.delim(pactx.deg, "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_primaryCells_paclitaxel_gest.txt")
Paclitaxel <- pactx.deg$Symbol[pactx.deg$pvalue < 0.05 & !(rownames(pactx.deg) %in% c("D", "A", "E"))]
Paclitaxel

#cispt.fit <- SEMrun(G42, pcr, v$CisplatinIC150uM, algo = "ricf")
cispt.fit <- SEMrun(G42, pcr, v$cispty, algo = "ricf")
cispt.deg <- cispt.fit$gest
cispt.deg <- data.frame(Symbol = mapIds(org.Hs.eg.db, rownames(cispt.deg), column = "SYMBOL", keytype = "ENTREZID"), cispt.deg)
write.delim(cispt.deg, "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_primaryCells_cisplatin_gest.txt")
Cisplatin <- cispt.deg$Symbol[cispt.deg$pvalue < 0.05 & !(rownames(cispt.deg) %in% c("D", "A", "E"))]
Cisplatin

carpt.fit <- SEMrun(G42, pcr, v$Carboplatin, algo = "ricf")
carpt.deg <- carpt.fit$gest
carpt.deg <- data.frame(Symbol = mapIds(org.Hs.eg.db, rownames(carpt.deg), column = "SYMBOL", keytype = "ENTREZID"), carpt.deg)
write.delim(carpt.deg, "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_primaryCells_carboplatin_gest.txt")
Carboplatin <- carpt.deg$Symbol[carpt.deg$pvalue < 0.05 & !(rownames(carpt.deg) %in% c("D", "A", "E"))]
Carboplatin


# R

library(randomForest)
library(e1071)
library(CMA)
library(OptimalCutpoints)
#library(boot)
#library(mosaic)
#library(parallel)
#library(robustbase)
library(ggplot2)
library(reshape2)
library(SEMgraph)
library(SEMdata)
library(huge)
library(org.Hs.eg.db)
library(VennDiagram)
#source("~/euNet/Likelihood_functions.R")
source("~/euNet/euNet_core/rutils.R")
source("~/UniPV_StatMedGen/Lectures/LezioniMaster2020/Statistica genomica 2020/SG.function.R")

vset <- read.delim("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/CellLines12_RTqPCR_42DEGs.txt", stringsAsFactors = FALSE)
vset$cellLine <- c("HOSE", "FTI194", "COV318", "HEY", "OV90", "PEO1",
                   "OVCAR3", "OVGEM2", "OVGEM6", "OVGEM8", "OVGEM9",
                   "OVGEM11", "OVGEM16", "OVGEM20")

v <- vset[3:14, c(1, 3:44)]
v$y <- c(2, 2, 2, 2, 2, 0, 0, 1, 0, 0, 1, 0)
names(v)[2] <- "ABHD5"

#pcr <- v[, 2:43]
pcr <- huge.npn(v[, 2:43])
rownames(pcr) <- v$cellLine
pcr <- as.data.frame(pcr)
pcr$y <- v$y

x <- pcr[pcr$y == 2, colnames(pcr) %in% c("GNG11", "SLC15A3", "UQCC1",
                                          "PLCG2", "RNF24", "LDLR",
                                          "IGFBP7", "CTNNBL1", "TTI1",
                                          "y")]
r <- cor(x[, 1:9])
#x <- pcr[pcr$y == 2,]
#r <- cor(x[, 1:42])
G <- corr2graph(r, n = nrow(x[, 1:9]), type = "tmfg")
dat <- pcr[pcr$y != 2,]
#dat <- pcr
#dat$y[dat$y == 0] <- 1
#dat$y[dat$y == 2] <- 0
dat$CTNNBL1 <- -1*dat$CTNNBL1
G.fit <- SEMrun(G, dat, dat$y, algo = "ricf")
#x1 <- pcr[pcr$y != 0,]
#x1$y[x1$y == 2] <- 0
#G.fit <- SEMrun(G, x1, x1$y, algo = "ricf")
G.fit$gest
#write.delim(G.fit$gest, "~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_stabilizedCells_gest.txt", rownames = TRUE)

#jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_stabilizedCells_covariance.jpg", width = 20, height = 10, units = 'in', res = 300)
gplot(G.fit$graph, cex.main = 2.5, fontsize = 22, main = "Stabilized cell RTqPCR covariance structure")
#dev.off()


if (FALSE) {
x0 <- pcr[pcr$y != 1, colnames(pcr) %in% c("GNG11", "SLC15A3", "UQCC1",
                                           "PLCG2", "RNF24", "LDLR",
                                           "IGFBP7", "CTNNBL1", "TTI1",
                                           "y")]
x0$y[x0$y == 2] <- 1
r0 <- cor(x0[, 1:9])
G0 <- corr2graph(r0, n = nrow(x0), type = "marg")
G0.fit <- SEMrun(G0, x0, x0$y, algo = "ricf")
gplot(G0.fit$graph)

x1 <- pcr[pcr$y != 0, colnames(pcr) %in% c("GNG11", "SLC15A3", "UQCC1",
                                           "PLCG2", "RNF24", "LDLR",
                                           "IGFBP7", "CTNNBL1", "TTI1",
                                           "y")]
x1$y[x1$y == 2] <- 0
r1 <- cor(x1[, 1:9])
G1 <- corr2graph(r1, n = nrow(x1), type = "marg")
G1.fit <- SEMrun(G1, x1, x1$y, algo = "ricf")
gplot(G1.fit$graph)

x2 <- pcr[, colnames(pcr) %in% c("GNG11", "SLC15A3", "UQCC1",
                                 "PLCG2", "RNF24", "LDLR",
                                 "IGFBP7", "CTNNBL1", "TTI1",
                                 "y")]
x2$y[x2$y == 1] <- 0
x2$y[x2$y == 2] <- 1
r2 <- cor(x2[, 1:9])
G2 <- corr2graph(r2, n = nrow(x2), type = "marg")
G2.fit <- SEMrun(G2, x2, x2$y, algo = "ricf")
gplot(G2.fit$graph)
}


pcr <- huge.npn(v[, 2:43])
#pcr <- v[, 2:43]
pcr <- pcr[, colnames(pcr) %in% c("GNG11", "SLC15A3", "UQCC1",
                                  "PLCG2", "RNF24", "LDLR",
                                  "IGFBP7", "CTNNBL1", "TTI1")]
rownames(pcr) <- v$cellLine
pcr

pca <- princomp(pcr)
summary(pca, loadings = TRUE, cutoff = 0.0, digits = 3) 

#screeplot(pca)
#screeplot(pca, type = "lines")

k <- 9
dimensions <- c(1:k)
eigenvalues <- pca$sdev[1:k]^2
#plot(dimensions, eigenvalues, type = "b", main = "PCA scree plot")
#abline(h = 1, lty = 2)

if (FALSE) {
x <- pca$scores[, 1]
y <- pca$scores[, 2]
group <- v$y
colors <- c("0" = "blue", "1" = "red", "2" = "green3")

plot(x, y, pch = 16, col = colors[as.character(group)],
     xlab = "PC1",
     ylab = "PC2",
     type = "p",
     asp = 1,
     xlim = c(-5, 5))
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
legend("topleft",
       c("Sensitive", "Resistant", "Stabilized"),
       fill = colors)

C <- pca$scores
corr.circle(x = cor(pcr, C)[, 1], y = cor(pcr, C)[, 2], dim = c(1, 2))
}

#jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_cellLines_biplot.jpg", width = 20, height = 10, units = 'in', res = 300)
biplot(pca, choices = c(1, 2), scale = 1, pc.biplot = FALSE)
#dev.off()

jpeg("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_cellLines_venn3.jpg", width = 20, height = 20, units = 'in', res = 300)
grid.newpage()
draw.triple.venn(area1 = 15,
                 area2 = 16,
                 area3 = 9,
                 n12 = 5,
                 n23 = 2,
                 n13 = 2,
                 n123 = 2,
          		 col = c("deepskyblue2", "gold", "#21908dff"),
          		 fill = c("deepskyblue2", "gold", "#21908dff"),
          		 alpha = c(0.2, 0.2, 0.2),
          		 fontfamily = "sans",
          		 cex = 3,
				 cat.fontfamily = "sans",
          		 cat.col = c("darkblue", "firebrick", "darkgreen"),
          		 cat.cex = 3,
				 rotation = 1,
                 category = c("Human\nsamples", "Primary\ncell lines", "Stabilized\n cell lines"),
                 cat.dist = c(0.05, 0.05, 0.05))
dev.off()

png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_cellLines_venn2.png", width = 20, height = 20, units = 'in', res = 450)
grid.newpage()
draw.pairwise.venn(area1 = 15,
                   area2 = 16,
                   cross.area = 5,
          		   #col = c("#440154ff", "gold"),
          		   #fill = c("#440154ff", "gold"),
          		   label.col = c("darkblue", "darkgreen", "firebrick"),
          		   col = c("deepskyblue2", "gold"),
          		   fill = c("deepskyblue2", "gold"),
          		   alpha = c(0.2, 0.2),
          		   fontfamily = "sans",
          		   cex = 7,
				   cat.fontfamily = "sans",
          		   cat.col = c("darkblue", "firebrick"),
          		   cat.cex = 3,
                   category = c("Human\nsamples", "Primary\ncell lines"))
dev.off()



# R

library(dplyr)
library(ggplot2)
library(stringr)
#library(UpSetR)

if (FALSE) {
# Old version
x <- data.frame(Gene = c("GNG11", "TACC1", "ABDHD5", "SLC15A3", "UQCC1",
                         "PLCG2", "TSPAN31", "RNF24", "LDLR", "IGFBP7",
                         "CKB", "CTNNBL1", "MYO5A", "TTI1", "WDPCP"),
                Estimate = c(-1.13, -0.65, -0.52, -0.98, 0.50, -0.92,
                             0.48, 0.51, -0.75, -0.39, 0.89, 0.43,
                             -0.22, 0.47, 0.21),
                P = c(2.17E-03, 0.0631, 0.0871, 0.0375, 0.0181, 0.0322,
                      0.0648, 0.0091, 0.0461, 0.0410, 0.0518, 0.0312,
                      0.4601, 0.0474, 0.2137),
                CI95 = c("-1.79, -0.43", "-1.33, 0.09", "-1.10, 0.07",
                         "-1.79, -0.08", "0.07, 0.89", "-1.73, -0.07",
                         "-0.05, 0.94", "0.16, 0.99", "-1.51, -0.02",
                         "-0.92, -0.03", "-0.01, 1.76", "0.03, 0.96",
                         "-0.89, 0.38", "5.96E-06, 1.12", "-0.15, 0.70"),
                Marker = c("RCFW", "RFC", "RFC", "RCFW", "RCFW", "RCFW",
                           "RCF", "RCFW", "RCFW", "RCFW", "RCF", "RCFW",
                           "RCF", "RCFW", "RCF"),
                Human = c(1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0),
                Primary = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0),
                Stabilized = c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0))
}

x <- data.frame(Gene = c("GNG11", "RNF24", "UQCC1", "CTNNBL1", "PLCG2",
                         "TTI1", "IGFBP7", "SLC15A3", "TSPAN31", "CKB"),
                Estimate = c(-1.09, 0.59, 0.54, 0.53, -0.91, 0.59, -0.42,
                             -0.98, 0.53, 0.94),
                P = c(0.00475, 0.00498, 0.01013, 0.01781, 0.02592,
                      0.02751, 0.02837, 0.02926, 0.04038, 0.04040),
                CI95 = c("-1.71, -0.34", "0.20, 1.12", "0.15, 0.93",
                         "0.05, 1.03", "-1.70, -0.13", "0.06, 1.21",
                         "-0.89, -0.06", "-1.75, -0.17", "0.03, 0.98",
                         "0.05, 1.79"),
                Marker = c("W", "WP", "W", "WP", "W", "WP", "W", "W", "W", "W"),
                Human = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                Primary = c(0, 1, 0, 1, 0, 1, 0, 0, 0, 0),
                Stabilized = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0))

#x$Human <- 1
#x$N <- apply(x[, 6:8], 1, sum)
x <- x[x$Human == 1,]
x$N <- apply(x[, 6:7], 1, sum)

#png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_polar_barplot.png", width = 20, height = 18, units = 'in', res = 450)
pdf("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_polar_barplot.pdf", width = 20, height = 18)
plt <- ggplot(x) +
  # Make custom panel grid
  theme_bw() +
  theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  #axis.text = element_text(size = 30),
  axis.text = element_text(size = c(29, 29, 29, 29, 29, 29, 29, 26, 29, 29)),
  #axis.text.x = element_text(angle = 0, hjust = 2, vjust = 2),
  axis.text.y = element_blank(),
  axis.title = element_text(size = 36, face = "bold"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 26),
  legend.title = element_text(size = 30)) +
    
  geom_hline(
    aes(yintercept = y), 
    #data.frame(y = c(0, 1, 1.3, 2, 3)),
    #color = c("grey", "grey", "red3", "grey", "grey")
    data.frame(y = c(0, 1, 1.3, 2)),
    color = c("grey", "grey", "red3", "grey")
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(str_wrap(Gene, 5), N),
      y = N,
      fill = Estimate
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = 0.9
  ) +
  
  # Add dots to represent the mean gain
  geom_point(
    aes(
      x = reorder(str_wrap(Gene, 7), -log10(P)),
      y = -log10(P)
    ),
    size = 4,
    color = "gray12"
  ) +
  
  # Lollipop shaft for mean gain per region
  geom_segment(
    aes(
      x = reorder(str_wrap(Gene, 5), -log10(P)),
      y = -1,
      xend = reorder(str_wrap(Gene, 5), -log10(P)),
      #yend = 3
      yend = 2
    ),
    linetype = "dashed",
    color = "gray12"
  ) + 
  
  scale_fill_gradientn(
    "Shift",
     colours = c("deepskyblue", "white", "darkorange"),
     #colours = c("dodgerblue", "deepskyblue", "white", "orange")
     limits = c(-1.1, 1.1)
  ) +
  
  #labs(x = "Wilcoxon estimated shift", y = "1: Patients,   2: Patients and OV.GEM,   3: Patients, OV.GEM, and Stabilized") +
  labs(x = "Selected markers", y = "1: Patients,   2: Patients and OV.GEM") +
  
  # Make it circular!
  coord_polar()
plt
dev.off()


p <- data.frame(Pathway = c("Fc epsilon RI signaling pathway ",
                            "Inositol phosphate metabolism ",
                            "Chemokine signaling pathway ",
                            "Phosphatidylinositol signaling system "),
                observed = c(3, 3, 4, 3),
                n = c(68, 73, 192, 97),
                pvalue = c(0.000392452875552, 0.000483513128126,
                           0.000692703575717, 0.001106829069923),
                FDR = c(0.027246340644872, 0.027246340644872,
                        0.027246340644872, 0.032651457562737),
                OR = c(23.5420118343195, 21.8549450549451,
                       11.0694288913774, 16.2553191489362),
                score = c(184.642213890043, 166.850093878057,
                          80.5290811249001, 110.63786422992),
                genes = c("PLCG2, INPP5D, ALOX5AP",
                          "PLCG2, PLCB1, INPP5D",
                          "GNG11, PLCG2, PLCB1, CXCL16",
                          "PLCG2, PLCB1, INPP5D"))
#p$P <- -log10(p$FDR)


#png("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_kegg_barplot.png", width = 20, height = 20, units = 'in', res = 450)
pdf("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_kegg_barplot.pdf", width = 20, height = 20)
p %>%
  arrange(score) %>%
  mutate(Pathway = factor(Pathway, levels = c("Chemokine signaling pathway ",
                                              "Phosphatidylinositol signaling system ",
                                              "Inositol phosphate metabolism ",
                                              "Fc epsilon RI signaling pathway "))) %>%
  ggplot(aes(x = Pathway, y = score)) +
    geom_segment(aes(xend = Pathway, yend = 0), color = "deepskyblue2", lwd = 2) +
    geom_point(size = 16, color = "deepskyblue2") +
    coord_flip() +
    
    theme_bw() +
    theme(panel.border = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 2, colour = "black"),
    axis.text = element_text(size = 38),
    axis.title = element_text(size = 38, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 38),
    legend.title = element_text(size = 38)) +
    
    annotate("text", x = 4.13, y = 48, size = 10,
             label = "PLCG2, INPP5D, ALOX5AP",
             parse = FALSE) +
    
    annotate("text", x = 3.13, y = 43, size = 10,
             label = "PLCG2, PLCB1, INPP5D",
             parse = FALSE) +
    
    annotate("text", x = 2.13, y = 43, size = 10,
             label = "PLCG2, PLCB1, INPP5D",
             parse = FALSE) +
    
    annotate("text", x = 1.13, y = 59, size = 10,
             label = "GNG11, PLCG2, PLCB1, CXCL16",
             parse = FALSE) +
    
    xlab("")
dev.off()


# R

library(randomForest)
library(e1071)
library(CMA)
library(OptimalCutpoints)
#library(boot)
#library(mosaic)
#library(parallel)
#library(robustbase)
library(ggplot2)
library(reshape2)
library(SEMgraph)
library(SEMdata)
library(huge)
library(org.Hs.eg.db)
library(VennDiagram)
#source("~/euNet/Likelihood_functions.R")
source("~/euNet/euNet_core/rutils.R")

load("~/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/RFC/Provac_RFC.RData")

seeds <- c("GNG11", "SLC15A3", "UQCC1", "PLCG2", "RNF24", "LDLR", "IGFBP7", "CTNNBL1", "TTI1")
#seeds <- c("CTNNBL1", "RNF24", "TTI1", "UQCC1", "GNG11", "PLCG2")
#seeds <- c("CTNNBL1", "RNF24", "TTI1", "UQCC1", "GNG11")
#seeds <- c("CTNNBL1", "RNF24", "TTI1", "UQCC1", "PLCG2")

data <- rfc.provac$data
Y <- data[, colnames(data) %in% seeds]
Y$y <- data$y
Y$K <- data$K
head(Y)

T1 <- Y[Y$K == "A" | Y$K == "B" | Y$K == "C",]
T2 <- Y[Y$K == "A" | Y$K == "B" | Y$K == "D",]
T3 <- Y[Y$K == "A" | Y$K == "C" | Y$K == "D",]
T4 <- Y[Y$K == "B" | Y$K == "C" | Y$K == "A",]

set.seed(0)

n <- 10000

F1 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T1, ntree = n, importance = TRUE)
F2 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T2, ntree = n, importance = TRUE)
F3 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T3, ntree = n, importance = TRUE)
F4 <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
                   data = T4, ntree = n, importance = TRUE)

#Fx <- randomForest(formula(paste0("as.factor(y) ~ ", paste(seeds, collapse = " + "))),
#                   data = Y, ntree = n, importance = TRUE)
#R <- data.frame(importance(Fx))[, 3:4]
#names(R) <- c("MDA", "MDG")
#R$fA <- 100*(R$MDA - min(R$MDA))/(max(R$MDA) - min(R$MDA))
#R$fG <- 100*(R$MDG - min(R$MDG))/(max(R$MDG) - min(R$MDG))
#R$f <- apply(R[, c(3, 4)], 1, mean)

Rv1 <- data.frame(importance(F1))[, 3:4]
Rv2 <- data.frame(importance(F2))[, 3:4]
Rv3 <- data.frame(importance(F3))[, 3:4]
Rv4 <- data.frame(importance(F4))[, 3:4]
Rv <- cbind(Rv1, Rv2, Rv3, Rv4)
Rv$MDA <- apply(Rv[, c(1, 3, 5, 7)], 1, mean)
Rv$MDG <- apply(Rv[, c(2, 4, 6, 8)], 1, mean)
Rv <- Rv[, c(9, 10)]

Rv$fA <- 100*(Rv$MDA - min(Rv$MDA))/(max(Rv$MDA) - min(Rv$MDA))
Rv$fG <- 100*(Rv$MDG - min(Rv$MDG))/(max(Rv$MDG) - min(Rv$MDG))
Rv$f <- apply(Rv[, c(3, 4)], 1, mean)
Rv <- Rv[order(Rv$f, decreasing = TRUE),]
Rv

#write.delim(Rv, "/home/fernando/GemelliBioinfoUnit/z_archive/PUBLISHED/Ovarian_carcinoma_Gallo/Provac_RFC9_ranking.txt", rownames = TRUE)

V1 <- Y[Y$K == "D",]
V2 <- Y[Y$K == "C",]
V3 <- Y[Y$K == "B",]
V4 <- Y[Y$K == "A",]

C1 <- predict(F1, V1)
C2 <- predict(F2, V2)
C3 <- predict(F3, V3)
C4 <- predict(F4, V4)

TP <- rbind(V1[C1 == V1$y & V1$y == 1,], V2[C2 == V2$y & V2$y == 1,], V3[C3 == V3$y & V3$y == 1,], V4[C4 == V4$y & V4$y == 1,])
TN <- rbind(V1[C1 == V1$y & V1$y == 0,], V2[C2 == V2$y & V2$y == 0,], V3[C3 == V3$y & V3$y == 0,], V4[C4 == V4$y & V4$y == 0,])
FP <- rbind(V1[C1 == 1 & V1$y == 0,], V2[C2 == 1 & V2$y == 0,], V3[C3 == 1 & V3$y == 0,], V4[C4 == 1 & V4$y == 0,])
FN <- rbind(V1[C1 == 0 & V1$y == 1,], V2[C2 == 0 & V2$y == 1,], V3[C3 == 0 & V3$y == 1,], V4[C4 == 0 & V4$y == 1,])

confusion <- matrix(c(nrow(TP), nrow(FP), nrow(FN), nrow(TN)), nrow = 2, byrow = TRUE)
rownames(confusion) <- c("Pred1", "Pred0")
colnames(confusion) <- c("Real1", "Real0")
confusion
P <- performance(confusion)
P




