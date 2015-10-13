#### This is one file with all previous R scripts combined into one. LW 4-2-15. Will comment on.


###CODE FROM bootstrap.R (IAN)

#not sure what this code is for.
bootstrap <- function(damdata, samdata, n = 100, bins = seq(0,1000,10)) {
	hd <- hist(damdata, breaks = bins, plot = F)
	results <- lapply(1:n, function(i) {
		hs <- hist(sample(samdata, replace = T), breaks = bins, plot = F)
		val <- sum(hs$counts/sum(hs$counts) - (hs$counts/sum(hs$counts))[1] * (hd$counts/sum(hd$counts)/max(hd$counts/sum(hd$counts))))
	})
	results <- unlist(results, use.names = F)
}



#read in data
dam <- read.table("damnew_integrated.txt", sep="\t", header=TRUE, as.is=TRUE)
d8 <- read.table("8hrnew_integrated.txt", sep="\t", header=TRUE, as.is=TRUE)
d16 <- read.table("16hrnew_integrated.txt", sep="\t", header=TRUE, as.is=TRUE)
d24 <- read.table("24hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d48 <- read.table("48hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d72 <- read.table("72hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d96 <- read.table("96hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d120 <- read.table("120hrnew_integrated.txt", sep="\t", header=T, as.is=T)

#read in coverages
all_cov <- read.csv("~/Dropbox/SMRT data (1)/coverageALL.csv", header = T, sep = ",") # sites seem to be missing in this

#spit out only rows where IPD is below 50,000
dam <- dam[(dam$IPD < 50000),]
d8 <- d8[(d8$IPD < 50000),]
d16 <- d16[(d16$IPD < 50000),]
d24 <- d24[(d24$IPD < 50000),]
d48 <- d48[(d48$IPD < 50000),]
d72 <- d72[(d72$IPD < 50000),]
d96 <- d96[(d96$IPD < 50000),]
d120 <- d120[(d120$IPD < 50000),]

#assign 1000 value to any IPD greater than 1000
dam[dam$IPD > 1000,]$IPD <- 1000
d8[d8$IPD > 1000,]$IPD <- 1000
d16[d16$IPD > 1000,]$IPD <- 1000
d24[d24$IPD > 1000,]$IPD <- 1000
d48[d48$IPD > 1000,]$IPD <- 1000
d72[d72$IPD > 1000,]$IPD <- 1000
d96[d96$IPD > 1000,]$IPD <- 1000
d120[d120$IPD > 1000,]$IPD <- 1000

bins <- seq(0,1000,1)
layout(matrix(1:8, ncol = 2))
hdam <- hist(dam$IPD, breaks = bins)
h8 <- hist(d8$IPD, breaks = bins)
h16 <- hist(d16$IPD, breaks = bins)
h24 <- hist(d24$IPD, breaks = bins)
h48 <- hist(d48$IPD, breaks = bins)
h72 <- hist(d72$IPD, breaks = bins)
h96 <- hist(d96$IPD, breaks = bins)
h120 <- hist(d120$IPD, breaks = bins)

layout(1)
plot(hdam$mids, hdam$counts/sum(hdam$counts), type = "b", pch = ".", xlim = c(0,500), ylim = c(0,1))
points(h8$mids, h8$counts/sum(h8$counts), type = "b", pch = ".", col = 2)
points(h16$mids, h16$counts/sum(h16$counts), type = "b", pch = ".", col = 2)
points(h24$mids, h24$counts/sum(h24$counts), type = "b", pch = ".", col = 4)
points(h48$mids, h48$counts/sum(h48$counts), type = "b", pch = ".", col = 5)
points(h72$mids, h72$counts/sum(h72$counts), type = "b", pch = ".", col = 6)
points(h96$mids, h96$counts/sum(h96$counts), type = "b", pch = ".", col = 7)
points(h120$mids, h120$counts/sum(h120$counts), type = "b", pch = ".", col = 8)
legend('topright', legend = c('dam', 8, 16, 24, 48, 72, 96, 120), fill = 1:8)

plot(hdam$mids, h8$counts/sum(h8$counts) * (1 - (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts)))), type = "b", pch = ".", col = 3)

plot(hdam$mids, hdam$counts/sum(hdam$counts), type = "b", pch = ".", xlim = c(0,500), ylim = c(0,1))

plot(h8$counts/sum(h8$counts), ylim = c(0,1))
points((hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))), col = "red")
points((h8$counts/sum(h8$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))), col = "green")
points(h8$counts/sum(h8$counts) -(h8$counts/sum(h8$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))), col = "orange")
legend('topright', legend = c('t8', 'dam-', 'rescaled dam-', 'residual area'), fill = c("black", "red", "green", "orange"))
sum(h8$counts/sum(h8$counts) -(h8$counts/sum(h8$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))


h8m <- sum(h8$counts/sum(h8$counts) - (h8$counts/sum(h8$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h8m_bs <- bootstrap(dam$IPD, d8$IPD, bins = bins)
h16m <- sum(h16$counts/sum(h16$counts) - (h16$counts/sum(h16$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h16m_bs <- bootstrap(dam$IPD, d16$IPD, bins = bins)
h24m <- sum(h24$counts/sum(h24$counts) - (h24$counts/sum(h24$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h24m_bs <- bootstrap(dam$IPD, d24$IPD, bins = bins)
h48m <- sum(h48$counts/sum(h48$counts) - (h48$counts/sum(h48$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h48m_bs <- bootstrap(dam$IPD, d48$IPD, bins = bins)
h72m <- sum(h72$counts/sum(h72$counts) - (h72$counts/sum(h72$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h72m_bs <- bootstrap(dam$IPD, d72$IPD, bins = bins)
h96m <- sum(h96$counts/sum(h96$counts) - (h96$counts/sum(h96$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h96m_bs <- bootstrap(dam$IPD, d96$IPD, bins = bins)
h120m <- sum(h120$counts/sum(h120$counts) - (h120$counts/sum(h120$counts))[1] * (hdam$counts/sum(hdam$counts)/max(hdam$counts/sum(hdam$counts))))
h120m_bs <- bootstrap(dam$IPD, d120$IPD, bins = bins)
plot(c(h8m,h16m,h24m,h48m,h72m,h96m,h120m), type = "b", pch = 20, xaxt = "n", ylab = "percent of dam sites methylated", xlab = "time", ylim = c(0.85,0.9))
axis(1, 1:7, c(8,16,24,48,72,96,120))
segments(1,min(h8m_bs),1,max(h8m_bs))
segments(2,min(h16m_bs),2,max(h16m_bs))
segments(3,min(h24m_bs),3,max(h24m_bs))
segments(4,min(h48m_bs),4,max(h48m_bs))
segments(5,min(h72m_bs),5,max(h72m_bs))
segments(6,min(h96m_bs),6,max(h96m_bs))
segments(7,min(h120m_bs),7,max(h120m_bs))



###### Code from Ian_coverage.R
#Filtered out low coverage sites, excluded pair-wise or complete exclusion.


uni_positions <- union(names(dam_2), union(names(d8_2), union(names(d16_2), union(names(d24_2), union(names(d48_2), union(names(d72_2), union(names(d96_2), names(d120_2))))))))
int_positions <- intersect(names(dam_2), intersect(names(d8_2), intersect(names(d16_2), intersect(names(d24_2), intersect(names(d48_2), intersect(names(d72_2), intersect(names(d96_2), names(d120_2))))))))

#stopped editing here

#positions <- intersect(names(dam_2), names(d96_2))
#make a list of position by coverage (#reads at position)
dam_cov <-lapply(dam_2[positions], nrow)
d8_cov <-lapply(d8_2[positions], nrow)
d16_cov <-lapply(d16_2[positions], nrow)
d24_cov <-lapply(d24_2[positions], nrow)
d48_cov <-lapply(d48_2[positions], nrow)
d72_cov <-lapply(d72_2[positions], nrow)
d96_cov <-lapply(d96_2[positions], nrow)
d120_cov <-lapply(d120_2[positions], nrow)

#unlist d_cov so that it can be data.framed
dam_cov <- unlist(dam_cov)
d8_cov <- unlist(d8_cov)
d16_cov <- unlist(d16_cov)
d24_cov <- unlist(d24_cov)
d48_cov <- unlist(d48_cov)
d72_cov <- unlist(d72_cov)
d96_cov <- unlist(d96_cov)
d120_cov <-unlist(d120_cov)

#restarted editing here

#combine all coverage data into one dataset (already combined into one file through "merge")
#find rows where coverage is below 16 across more than 3 timepoints
all_cov$X <- NULL
bad_cov <- all_cov[rowSums(all_cov[,2:8] < 16) > 0,]

#take problem rows (low coverage across multiple timepoints) out of timepoint data:

dam_3 <- dam[(dam$RefPos %in% bad_cov$position) == FALSE,]
d8_3 <- d8[(d8$RefPos %in% bad_cov$position) == FALSE,]
d16_3 <- d16[(d16$RefPos %in% bad_cov$position) == FALSE,]
d24_3 <- d24[(d24$RefPos %in% bad_cov$position) == FALSE,]
d48_3 <- d48[(d48$RefPos %in% bad_cov$position) == FALSE,]
d72_3 <- d72[(d72$RefPos %in% bad_cov$position) == FALSE,]
d96_3 <- d96[(d96$RefPos %in% bad_cov$position) == FALSE,]
d120_3 <- d120[(d120$RefPos %in% bad_cov$position) == FALSE,]

#split datasets without problem area
dam_4 <- split(dam_3, dam_3$RefPos)
d8_4 <- split(d8_3, d8_3$RefPos)
d16_4 <- split(d16_3, d16_3$RefPos)
d24_4 <- split(d24_3, d24_3$RefPos)
d48_4 <- split(d48_3, d48_3$RefPos)
d72_4 <- split(d72_3, d72_3$RefPos)
d96_4 <- split(d96_3, d96_3$RefPos)
d120_4 <- split(d120_3, d120_3$RefPos)

#filter out coverage list for coverage > 16
dam_cov16 <- subset(dam_cov, dam_cov > 16)
d8_cov16 <- subset(d8_cov, d8_cov > 16)
d16_cov16 <- subset(d16_cov, d16_cov > 16)
d24_cov16 <- subset(d24_cov, d24_cov > 16)
d48_cov16 <- subset(d48_cov, d48_cov > 16)
d72_cov16 <- subset(d72_cov, d72_cov > 16)
d96_cov16 <- subset(d96_cov, d96_cov > 16)
d120_cov16 <- subset(d120_cov, d120_cov > 16)

#only get rows where coverage is greater than 16 from datasets
dam_6 <- dam_3[(dam_3$RefPos %in% names(dam_cov16)) == TRUE,]
d16_5 <- d16_3[(d16_3$RefPos %in% names(d16_cov16)) == TRUE,]
d8_5 <- d8_3[(d8_3$RefPos %in% names(d8_cov16)) == TRUE,]
d24_5 <- d24_3[(d24_3$RefPos %in% names(d24_cov16)) == TRUE,]
d48_5 <- d48_3[(d48_3$RefPos %in% names(d48_cov16)) == TRUE,]
d72_5 <- d72_3[(d72_3$RefPos %in% names(d72_cov16)) == TRUE,]
d96_5 <- d96_3[(d96_3$RefPos %in% names(d96_cov16)) == TRUE,]
d120_5 <- d120_3[(d120_3$RefPos %in% names(d120_cov16)) == TRUE,]
#split datasets for wilcoxon

d8_6 <- split(d8_5, d8_5$RefPos)
dam_7 <- split(dam_6, dam_6$RefPos)
d16_6 <- split(d16_5, d16_5$RefPos)
d24_6 <- split(d24_5, d24_5$RefPos)
d48_6 <- split(d48_5, d48_5$RefPos)
d72_6 <- split(d72_5, d72_5$RefPos)
d96_6 <- split(d96_5, d96_5$RefPos)
d120_6 <- split(d120_5, d120_5$RefPos)

length(d16_cov16)
[1] 22102
str(d16_6)
List of 22102


#t.test on sites with cov >16
positions <- intersect(names(d8_6), names(dam_7))
pvals8 <- sapply(positions, function(x) {wilcox.test(d8_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d16_6), names(dam_7))
pvals16 <- sapply(positions, function(x) {wilcox.test(d16_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d24_6), names(dam_7))
pvals24 <- sapply(positions, function(x) {wilcox.test(d24_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d48_6), names(dam_7))
pvals48 <- sapply(positions, function(x) {wilcox.test(d48_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d72_6), names(dam_7))
pvals72 <- sapply(positions, function(x) {wilcox.test(d72_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d96_6), names(dam_7))
pvals96 <- sapply(positions, function(x) {wilcox.test(d96_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d120_6), names(dam_7))
pvals120 <- sapply(positions, function(x) {wilcox.test(d120_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})

#plot
> layout(matrix(1:8, 2, byrow=TRUE))
> plot(-log10(pvals8), main="wilcoxon test, 8hr v. dam (cov>16)", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals8)), col='red')
> 
> plot(-log10(pvals16), main="wilcoxon test, 16hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals16)), col='red')
> 
> plot(-log10(pvals24), main="wilcoxon test, 24hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals24)), col='red')
> 
> plot(-log10(pvals48), main="wilcoxon test, 48hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals48)), col='red')
> 
> plot(-log10(pvals72), main="wilcoxon test, 72hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals72)), col='red')
> 
> plot(-log10(pvals96), main="wilcoxon test, 96hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals96)), col='red')
> 
> plot(-log10(pvals120), main="wilcoxon test, 120hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals120)), col='red')


## FOR COMPLETE EXCLUSION OF BAD COVERAGE (above is pair-wise)

#make sure that all sites from every dataset with cov <16 are excluded in 16hr timepoint
list <- intersect(names(dam_cov16), names(d16_cov16))
list2 <- intersect(names(d8_cov16), list)
list3 <- intersect(list2,names(d24_cov16))
list4 <- intersect(list3, names(d48_cov16))
list5 <- intersect(list4, names(d72_cov16))
list6 <- intersect(list5, names(d96_cov16))
list7 <- intersect(list6, names(d120_cov16))

#only get rows where coverage is greater than 16 at ALL time points 
dam_8 <- dam_3[(dam_3$RefPos %in% list7) == TRUE,]
d16_8 <- d16_3[(d16_3$RefPos %in% list7) == TRUE,]
d8_8 <- d8_3[(d8_3$RefPos %in% list7) == TRUE,]
d24_8 <- d24_3[(d24_3$RefPos %in% list7) == TRUE,]
d48_8 <- d48_3[(d48_3$RefPos %in% list7) == TRUE,]
d72_8 <- d72_3[(d72_3$RefPos %in% list7) == TRUE,]
d96_8 <- d96_3[(d96_3$RefPos %in% list7) == TRUE,]
d120_8 <- d120_3[(d120_3$RefPos %in% list7) == TRUE,]

#split the datasets by refpos

dam_9 <- split(dam_8, dam_8$RefPos)
d8_9 <- split(d8_8, d8_8$RefPos)
d16_9 <- split(d16_8, d16_8$RefPos)
d24_9 <- split(d24_8, d24_8$RefPos)
d48_9 <- split(d48_8, d48_8$RefPos)
d72_9 <- split(d72_8, d72_8$RefPos)
d96_9 <- split(d96_8, d96_8$RefPos)
d120_9 <- split(d120_8, d120_8$RefPos)

#wilcoxon test on all pairs
pvals_e8 <- sapply(list7, function(x) {wilcox.test(d8_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e16 <- sapply(list7, function(x) {wilcox.test(d16_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e24 <- sapply(list7, function(x) {wilcox.test(d24_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e48 <- sapply(list7, function(x) {wilcox.test(d48_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e72 <- sapply(list7, function(x) {wilcox.test(d72_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e96 <- sapply(list7, function(x) {wilcox.test(d96_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e120 <- sapply(list7, function(x) {wilcox.test(d120_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})

#plot that shit
layout(matrix(1:8, 2, byrow=TRUE))
plot(-log10(pvals_e8), main="wilcoxon test, 8hr v. dam, exclusive", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e8)), col='red')

plot(-log10(pvals_e16), main="wilcoxon test, 16hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e16)), col='red')

plot(-log10(pvals_e24), main="wilcoxon test, 24hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e24)), col='red')

plot(-log10(pvals_e48), main="wilcoxon test, 48hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e48)), col='red')

plot(-log10(pvals_e72), main="wilcoxon test, 72hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e72)), col='red')

plot(-log10(pvals_e96), main="wilcoxon test, 96hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e96)), col='red')

plot(-log10(pvals_e120), main="wilcoxon test, 120hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e120)), col='red')


##### Code from ian_coverage2.R

# no threshold

pvals8 <- sapply(int_positions, function(x) {wilcox.test(d8_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals16 <- sapply(int_positions, function(x) {wilcox.test(d16_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals24 <- sapply(int_positions, function(x) {wilcox.test(d24_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals48 <- sapply(int_positions, function(x) {wilcox.test(d48_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals72 <- sapply(int_positions, function(x) {wilcox.test(d72_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals96 <- sapply(int_positions, function(x) {wilcox.test(d96_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals120 <- sapply(int_positions, function(x) {wilcox.test(d120_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})

threshold <- 0.05/37000

qvals8 <- qvalue(pvals8, pi0.method = 'bootstrap')
qvals16 <- qvalue(pvals16, pi0.method = 'bootstrap')


plot(1:7, c(sum(pvals8 <= threshold)/length(pvals8),sum(pvals16 <= threshold)/length(pvals16),sum(pvals24 <= threshold)/length(pvals24),sum(pvals48 <= threshold)/length(pvals48),sum(pvals72 <= threshold)/length(pvals72),sum(pvals96 <= threshold)/length(pvals96),sum(pvals120 <= threshold)/length(pvals120)), ylab = 'number of methylated sites', xlab = 'time point')


layout(matrix(1:8, ncol = 2))
hist(pvals8, breaks = 100)
hist(pvals16, breaks = 100)
hist(pvals24, breaks = 100)
hist(pvals48, breaks = 100)
hist(pvals72, breaks = 100)
hist(pvals96, breaks = 100)
hist(pvals120, breaks = 100)


#combine all coverage data into one dataset (already combined into one file through "merge")
#find rows where coverage is below 16 across more than 3 timepoints

all_cov$X <- NULL
good_cov <- all_cov[(rowSums(all_cov[,2:8] <= 20) > 0) == F,]

#t.test on sites with cov >16

pvals8 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d8_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals16 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d16_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals24 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d24_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals48 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d48_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals72 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d72_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals96 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d96_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})
pvals120 <- sapply(as.character(good_cov$position), function(x) {wilcox.test(d120_2[[x]]$IPD, dam_2[[x]]$IPD)$p.value})

threshold <- 0.05/30000
plot(1:7, c(sum(pvals8 <= threshold)/length(pvals8),sum(pvals16 <= threshold)/length(pvals16),sum(pvals24 <= threshold)/length(pvals24),sum(pvals48 <= threshold)/length(pvals48),sum(pvals72 <= threshold)/length(pvals72),sum(pvals96 <= threshold)/length(pvals96),sum(pvals120 <= threshold)/length(pvals120)), ylab = 'number of methylated sites', xlab = 'time point')


#take problem rows (low coverage across multiple timepoints) out of timepoint data:

dam_3 <- dam[(dam$RefPos %in% bad_cov$position) == FALSE,]
d8_3 <- d8[(d8$RefPos %in% bad_cov$position) == FALSE,]
d16_3 <- d16[(d16$RefPos %in% bad_cov$position) == FALSE,]
d24_3 <- d24[(d24$RefPos %in% bad_cov$position) == FALSE,]
d48_3 <- d48[(d48$RefPos %in% bad_cov$position) == FALSE,]
d72_3 <- d72[(d72$RefPos %in% bad_cov$position) == FALSE,]
d96_3 <- d96[(d96$RefPos %in% bad_cov$position) == FALSE,]
d120_3 <- d120[(d120$RefPos %in% bad_cov$position) == FALSE,]

#split datasets without problem area
dam_4 <- split(dam_3, dam_3$RefPos)
d8_4 <- split(d8_3, d8_3$RefPos)
d16_4 <- split(d16_3, d16_3$RefPos)
d24_4 <- split(d24_3, d24_3$RefPos)
d48_4 <- split(d48_3, d48_3$RefPos)
d72_4 <- split(d72_3, d72_3$RefPos)
d96_4 <- split(d96_3, d96_3$RefPos)
d120_4 <- split(d120_3, d120_3$RefPos)

#filter out coverage list for coverage > 16
dam_cov16 <- subset(dam_cov, dam_cov > 16)
d8_cov16 <- subset(d8_cov, d8_cov > 16)
d16_cov16 <- subset(d16_cov, d16_cov > 16)
d24_cov16 <- subset(d24_cov, d24_cov > 16)
d48_cov16 <- subset(d48_cov, d48_cov > 16)
d72_cov16 <- subset(d72_cov, d72_cov > 16)
d96_cov16 <- subset(d96_cov, d96_cov > 16)
d120_cov16 <- subset(d120_cov, d120_cov > 16)

#only get rows where coverage is greater than 16 from datasets
dam_6 <- dam_3[(dam_3$RefPos %in% names(dam_cov16)) == TRUE,]
d16_5 <- d16_3[(d16_3$RefPos %in% names(d16_cov16)) == TRUE,]
d8_5 <- d8_3[(d8_3$RefPos %in% names(d8_cov16)) == TRUE,]
d24_5 <- d24_3[(d24_3$RefPos %in% names(d24_cov16)) == TRUE,]
d48_5 <- d48_3[(d48_3$RefPos %in% names(d48_cov16)) == TRUE,]
d72_5 <- d72_3[(d72_3$RefPos %in% names(d72_cov16)) == TRUE,]
d96_5 <- d96_3[(d96_3$RefPos %in% names(d96_cov16)) == TRUE,]
d120_5 <- d120_3[(d120_3$RefPos %in% names(d120_cov16)) == TRUE,]
#split datasets for wilcoxon

d8_6 <- split(d8_5, d8_5$RefPos)
dam_7 <- split(dam_6, dam_6$RefPos)
d16_6 <- split(d16_5, d16_5$RefPos)
d24_6 <- split(d24_5, d24_5$RefPos)
d48_6 <- split(d48_5, d48_5$RefPos)
d72_6 <- split(d72_5, d72_5$RefPos)
d96_6 <- split(d96_5, d96_5$RefPos)
d120_6 <- split(d120_5, d120_5$RefPos)

length(d16_cov16)
[1] 22102
str(d16_6)
List of 22102


#t.test on sites with cov >16
positions <- intersect(names(d8_6), names(dam_7))
pvals8 <- sapply(positions, function(x) {wilcox.test(d8_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d16_6), names(dam_7))
pvals16 <- sapply(positions, function(x) {wilcox.test(d16_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d24_6), names(dam_7))
pvals24 <- sapply(positions, function(x) {wilcox.test(d24_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d48_6), names(dam_7))
pvals48 <- sapply(positions, function(x) {wilcox.test(d48_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d72_6), names(dam_7))
pvals72 <- sapply(positions, function(x) {wilcox.test(d72_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d96_6), names(dam_7))
pvals96 <- sapply(positions, function(x) {wilcox.test(d96_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d120_6), names(dam_7))
pvals120 <- sapply(positions, function(x) {wilcox.test(d120_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})

#plot
> layout(matrix(1:8, 2, byrow=TRUE))
> plot(-log10(pvals8), main="wilcoxon test, 8hr v. dam (cov>16)", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals8)), col='red')
> 
> plot(-log10(pvals16), main="wilcoxon test, 16hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals16)), col='red')
> 
> plot(-log10(pvals24), main="wilcoxon test, 24hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals24)), col='red')
> 
> plot(-log10(pvals48), main="wilcoxon test, 48hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals48)), col='red')
> 
> plot(-log10(pvals72), main="wilcoxon test, 72hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals72)), col='red')
> 
> plot(-log10(pvals96), main="wilcoxon test, 96hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals96)), col='red')
> 
> plot(-log10(pvals120), main="wilcoxon test, 120hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals120)), col='red')

#Bonferroni correction/% methylation 
> sum(pvals16 <0.05)
[1] 22032
> sum(pvals16)
[1] 15.80195
> length(pvals16)
[1] 22076
> 22032/22076
[1] 0.9980069
> length(pvals8)
[1] 37020
> sum(pvals8 < 0.05)
[1] 36996
> 36996/37020
[1] 0.9993517
> sum(pvals8 < 0.05/length(pvals8))
[1] 27125
> 27125/37020
[1] 0.732712
> sum(pvals16 < 0.05/length(pvals16))
[1] 9967
> length(pvals16)
[1] 22076
> 9967/22076
[1] 0.4514858
> sum(pvals24 < 0.05/length(pvals24))
[1] 27425
> length(pvals24)
[1] 36449
> 27425/36449
[1] 0.7524212
> sum(pvals48 < 0.05/length(pvals48))
[1] 24446
> length(pvals48)
[1] 35682
> 24446/35682
[1] 0.6851073
> sum(pvals72 < 0.05/length(pvals72))
[1] 31134
> length(pvals72)
[1] 37218
> 31134/37218
[1] 0.8365307
> sum(pvals96 < 0.05/length(pvals96))
[1] 30424
> length(pvals96)
[1] 37165
> 30424/37165
[1] 0.8186197
> sum(pvals120 < 0.05/length(pvals120))
[1] 28448
> length(pvals120)
[1] 36939
> 28448/36939
[1] 0.7701345
> 

## FOR COMPLETE EXCLUSION OF BAD COVERAGE (above is pair-wise)

#make sure that all sites from every dataset with cov <16 are excluded in 16hr timepoint
list <- intersect(names(dam_cov16), names(d16_cov16))
list2 <- intersect(names(d8_cov16), list)
list3 <- intersect(list2,names(d24_cov16))
list4 <- intersect(list3, names(d48_cov16))
list5 <- intersect(list4, names(d72_cov16))
list6 <- intersect(list5, names(d96_cov16))
list7 <- intersect(list6, names(d120_cov16))



#only get rows where coverage is greater than 16 at ALL time points 
dam_8 <- dam_3[(dam_3$RefPos %in% list7) == TRUE,]
d16_8 <- d16_3[(d16_3$RefPos %in% list7) == TRUE,]
d8_8 <- d8_3[(d8_3$RefPos %in% list7) == TRUE,]
d24_8 <- d24_3[(d24_3$RefPos %in% list7) == TRUE,]
d48_8 <- d48_3[(d48_3$RefPos %in% list7) == TRUE,]
d72_8 <- d72_3[(d72_3$RefPos %in% list7) == TRUE,]
d96_8 <- d96_3[(d96_3$RefPos %in% list7) == TRUE,]
d120_8 <- d120_3[(d120_3$RefPos %in% list7) == TRUE,]

#split the datasets by refpos

dam_9 <- split(dam_8, dam_8$RefPos)
d8_9 <- split(d8_8, d8_8$RefPos)
d16_9 <- split(d16_8, d16_8$RefPos)
d24_9 <- split(d24_8, d24_8$RefPos)
d48_9 <- split(d48_8, d48_8$RefPos)
d72_9 <- split(d72_8, d72_8$RefPos)
d96_9 <- split(d96_8, d96_8$RefPos)
d120_9 <- split(d120_8, d120_8$RefPos)

#wilcoxon test on all pairs
pvals_e8 <- sapply(list7, function(x) {wilcox.test(d8_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e16 <- sapply(list7, function(x) {wilcox.test(d16_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e24 <- sapply(list7, function(x) {wilcox.test(d24_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e48 <- sapply(list7, function(x) {wilcox.test(d48_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e72 <- sapply(list7, function(x) {wilcox.test(d72_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e96 <- sapply(list7, function(x) {wilcox.test(d96_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e120 <- sapply(list7, function(x) {wilcox.test(d120_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})

#plot that shit
layout(matrix(1:8, 2, byrow=TRUE))
plot(-log10(pvals_e8), main="wilcoxon test, 8hr v. dam, exclusive", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e8)), col='red')

plot(-log10(pvals_e16), main="wilcoxon test, 16hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e16)), col='red')

plot(-log10(pvals_e24), main="wilcoxon test, 24hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e24)), col='red')

plot(-log10(pvals_e48), main="wilcoxon test, 48hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e48)), col='red')

plot(-log10(pvals_e72), main="wilcoxon test, 72hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e72)), col='red')

plot(-log10(pvals_e96), main="wilcoxon test, 96hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e96)), col='red')

plot(-log10(pvals_e120), main="wilcoxon test, 120hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e120)), col='red')



##### CODE from coverage.R

#set working directory
setwd("~/Desktop/ Lacey/Bioinformatics")

#read in data
dam <- read.table("damnew_integrated.txt", sep="\t", header=TRUE, as.is=TRUE)
d8 <- read.table("8hrnew_integrated.txt", sep="\t", header=TRUE, as.is=TRUE)
d16 <- read.table("16hrnew_integrated.txt", sep="\t", header=TRUE, as.is=TRUE)
d24 <- read.table("24hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d48 <- read.table("48hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d72 <- read.table("72hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d96 <- read.table("96hrnew_integrated.txt", sep="\t", header=T, as.is=T)
d120 <- read.table("120hrnew_integrated.txt", sep="\t", header=T, as.is=T)

#spit out only rows where IPD is below 50,000
dam <- dam[(dam$IPD < 50000),]
d8 <- d8[(d8$IPD < 50000),]
d16 <- d16[(d16$IPD < 50000),]
d24 <- d24[(d24$IPD < 50000),]
d48 <- d48[(d48$IPD < 50000),]
d72 <- d72[(d72$IPD < 50000),]
d96 <- d96[(d96$IPD < 50000),]
d120 <- d120[(d120$IPD < 50000),]

#Add new column named "Strand" where if in dataset, RefPos = 0, strand= -1
dam[dam$Strand==0,]$Strand <- -1
d8[d8$Strand==0,]$Strand <- -1
d16[d16$Strand==0,]$Strand <- -1
d24[d24$Strand==0,]$Strand <- -1
d48[d48$Strand==0,]$Strand <- -1
d72[d72$Strand==0,]$Strand <- -1
d96[d96$Strand==0,]$Strand <- -1
d120[d120$Strand==0,]$Strand <- -1
#assign 1000 value to any IPD greater than 1000
dam[dam$IPD > 1000,]$IPD <- 1000
d8[d8$IPD > 1000,]$IPD <- 1000
d16[d16$IPD > 1000,]$IPD <- 1000
d24[d24$IPD > 1000,]$IPD <- 1000
d48[d48$IPD > 1000,]$IPD <- 1000
d72[d72$IPD > 1000,]$IPD <- 1000
d96[d96$IPD > 1000,]$IPD <- 1000
d120[d120$IPD > 1000,]$IPD <- 1000


#split data by position - this makes a list of data.frames
#dam_2 <- split(dam, dam$RefPos)
dam_2 <- split(dam, dam$RefPos)
d8_2 <- split(d8, d8$RefPos)
d16_2 <- split(d16, d16$RefPos)
d24_2 <- split(d24, d24$RefPos)
d48_2 <- split(d48, d48$RefPos)
d72_2 <- split(d72, d72$RefPos)
d96_2 <- split(d96, d96$RefPos)
d120_2 <- split(d120, d120$RefPos)

positions <- intersect(names(dam_2), names(d96_2))
#make a list of position by coverage (#reads at position)
dam_cov <-lapply(dam_2[positions], nrow)
d8_cov <-lapply(d8_2[positions], nrow)
d16_cov <-lapply(d16_2[positions], nrow)
d24_cov <-lapply(d24_2[positions], nrow)
d48_cov <-lapply(d48_2[positions], nrow)
d72_cov <-lapply(d72_2[positions], nrow)
d96_cov <-lapply(d96_2[positions], nrow)
d120_cov <-lapply(d120_2[positions], nrow)

#unlist d_cov so that it can be data.framed
dam_cov <- unlist(dam_cov)
d8_cov <- unlist(d8_cov)
d16_cov <- unlist(d16_cov)
d24_cov <- unlist(d24_cov)
d48_cov <- unlist(d48_cov)
d72_cov <- unlist(d72_cov)
d96_cov <- unlist(d96_cov)
d120_cov <-unlist(d120_cov)

#combine all coverage data into one dataset (already combined into one file through "merge")
#find rows where coverage is below 16 across more than 3 timepoints
bad_cov <- all_cov[rowSums(all_cov[,2:8] < 16) > 3,]

#this is the fixed coverage dataset:
all_cov <- read.csv("coverage_fixed.csv", header=TRUE, as.is=TRUE, sep=',')
all_cov$X <_ NULL

#take problem rows (low coverage across multiple timepoints) out of timepoint data:

dam_3 <- dam[(dam$RefPos %in% bad_cov$position) == FALSE,]
d8_3 <- d8[(d8$RefPos %in% bad_cov$position) == FALSE,]
d16_3 <- d16[(d16$RefPos %in% bad_cov$position) == FALSE,]
d24_3 <- d24[(d24$RefPos %in% bad_cov$position) == FALSE,]
d48_3 <- d48[(d48$RefPos %in% bad_cov$position) == FALSE,]
d72_3 <- d72[(d72$RefPos %in% bad_cov$position) == FALSE,]
d96_3 <- d96[(d96$RefPos %in% bad_cov$position) == FALSE,]
d120_3 <- d120[(d120$RefPos %in% bad_cov$position) == FALSE,]

#split datasets without problem area
dam_4 <- split(dam_3, dam_3$RefPos)
d8_4 <- split(d8_3, d8_3$RefPos)
d16_4 <- split(d16_3, d16_3$RefPos)
d24_4 <- split(d24_3, d24_3$RefPos)
d48_4 <- split(d48_3, d48_3$RefPos)
d72_4 <- split(d72_3, d72_3$RefPos)
d96_4 <- split(d96_3, d96_3$RefPos)
d120_4 <- split(d120_3, d120_3$RefPos)

#filter out coverage list for coverage > 16
dam_cov16 <- subset(dam_cov, dam_cov > 16)
d8_cov16 <- subset(d8_cov, d8_cov > 16)
d16_cov16 <- subset(d16_cov, d16_cov > 16)
d24_cov16 <- subset(d24_cov, d24_cov > 16)
d48_cov16 <- subset(d48_cov, d48_cov > 16)
d72_cov16 <- subset(d72_cov, d72_cov > 16)
d96_cov16 <- subset(d96_cov, d96_cov > 16)
d120_cov16 <- subset(d120_cov, d120_cov > 16)

#only get rows where coverage is greater than 16 from datasets
dam_6 <- dam_3[(dam_3$RefPos %in% names(dam_cov16)) == TRUE,]
d16_5 <- d16_3[(d16_3$RefPos %in% names(d16_cov16)) == TRUE,]
d8_5 <- d8_3[(d8_3$RefPos %in% names(d8_cov16)) == TRUE,]
d24_5 <- d24_3[(d24_3$RefPos %in% names(d24_cov16)) == TRUE,]
d48_5 <- d48_3[(d48_3$RefPos %in% names(d48_cov16)) == TRUE,]
d72_5 <- d72_3[(d72_3$RefPos %in% names(d72_cov16)) == TRUE,]
d96_5 <- d96_3[(d96_3$RefPos %in% names(d96_cov16)) == TRUE,]
d120_5 <- d120_3[(d120_3$RefPos %in% names(d120_cov16)) == TRUE,]
#split datasets for wilcoxon

d8_6 <- split(d8_5, d8_5$RefPos)
dam_7 <- split(dam_6, dam_6$RefPos)
d16_6 <- split(d16_5, d16_5$RefPos)
d24_6 <- split(d24_5, d24_5$RefPos)
d48_6 <- split(d48_5, d48_5$RefPos)
d72_6 <- split(d72_5, d72_5$RefPos)
d96_6 <- split(d96_5, d96_5$RefPos)
d120_6 <- split(d120_5, d120_5$RefPos)

length(d16_cov16)
[1] 22102
str(d16_6)
List of 22102


#t.test on sites with cov >16
positions <- intersect(names(d8_6), names(dam_7))
pvals8 <- sapply(positions, function(x) {wilcox.test(d8_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d16_6), names(dam_7))
pvals16 <- sapply(positions, function(x) {wilcox.test(d16_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d24_6), names(dam_7))
pvals24 <- sapply(positions, function(x) {wilcox.test(d24_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d48_6), names(dam_7))
pvals48 <- sapply(positions, function(x) {wilcox.test(d48_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d72_6), names(dam_7))
pvals72 <- sapply(positions, function(x) {wilcox.test(d72_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d96_6), names(dam_7))
pvals96 <- sapply(positions, function(x) {wilcox.test(d96_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})
positions <- intersect(names(d120_6), names(dam_7))
pvals120 <- sapply(positions, function(x) {wilcox.test(d120_6[[x]]$IPD, dam_7[[x]]$IPD)$p.value})

#plot
> layout(matrix(1:8, 2, byrow=TRUE))
> plot(-log10(pvals8), main="wilcoxon test, 8hr v. dam (cov>16)", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals8)), col='red')
> 
> plot(-log10(pvals16), main="wilcoxon test, 16hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals16)), col='red')
> 
> plot(-log10(pvals24), main="wilcoxon test, 24hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals24)), col='red')
> 
> plot(-log10(pvals48), main="wilcoxon test, 48hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals48)), col='red')
> 
> plot(-log10(pvals72), main="wilcoxon test, 72hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals72)), col='red')
> 
> plot(-log10(pvals96), main="wilcoxon test, 96hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals96)), col='red')
> 
> plot(-log10(pvals120), main="wilcoxon test, 120hr v. dam", ylab="-log10(p-values)", xlim=c(0,37500), ylim=c(0,15))
> abline(h=-log10(0.05/length(pvals120)), col='red')

#Bonferroni correction/% methylation 
> sum(pvals16 <0.05)
[1] 22032
> sum(pvals16)
[1] 15.80195
> length(pvals16)
[1] 22076
> 22032/22076
[1] 0.9980069
> length(pvals8)
[1] 37020
> sum(pvals8 < 0.05)
[1] 36996
> 36996/37020
[1] 0.9993517
> sum(pvals8 < 0.05/length(pvals8))
[1] 27125
> 27125/37020
[1] 0.732712
> sum(pvals16 < 0.05/length(pvals16))
[1] 9967
> length(pvals16)
[1] 22076
> 9967/22076
[1] 0.4514858
> sum(pvals24 < 0.05/length(pvals24))
[1] 27425
> length(pvals24)
[1] 36449
> 27425/36449
[1] 0.7524212
> sum(pvals48 < 0.05/length(pvals48))
[1] 24446
> length(pvals48)
[1] 35682
> 24446/35682
[1] 0.6851073
> sum(pvals72 < 0.05/length(pvals72))
[1] 31134
> length(pvals72)
[1] 37218
> 31134/37218
[1] 0.8365307
> sum(pvals96 < 0.05/length(pvals96))
[1] 30424
> length(pvals96)
[1] 37165
> 30424/37165
[1] 0.8186197
> sum(pvals120 < 0.05/length(pvals120))
[1] 28448
> length(pvals120)
[1] 36939
> 28448/36939
[1] 0.7701345
> 

## FOR COMPLETE EXCLUSION OF BAD COVERAGE (above is pair-wise)

#make sure that all sites from every dataset with cov <16 are excluded in 16hr timepoint
list <- intersect(names(dam_cov16), names(d16_cov16))
list2 <- intersect(names(d8_cov16), list)
list3 <- intersect(list2,names(d24_cov16))
list4 <- intersect(list3, names(d48_cov16))
list5 <- intersect(list4, names(d72_cov16))
list6 <- intersect(list5, names(d96_cov16))
list7 <- intersect(list6, names(d120_cov16))



#only get rows where coverage is greater than 16 at ALL time points 
dam_8 <- dam_3[(dam_3$RefPos %in% list7) == TRUE,]
d16_8 <- d16_3[(d16_3$RefPos %in% list7) == TRUE,]
d8_8 <- d8_3[(d8_3$RefPos %in% list7) == TRUE,]
d24_8 <- d24_3[(d24_3$RefPos %in% list7) == TRUE,]
d48_8 <- d48_3[(d48_3$RefPos %in% list7) == TRUE,]
d72_8 <- d72_3[(d72_3$RefPos %in% list7) == TRUE,]
d96_8 <- d96_3[(d96_3$RefPos %in% list7) == TRUE,]
d120_8 <- d120_3[(d120_3$RefPos %in% list7) == TRUE,]

#split the datasets by refpos

dam_9 <- split(dam_8, dam_8$RefPos)
d8_9 <- split(d8_8, d8_8$RefPos)
d16_9 <- split(d16_8, d16_8$RefPos)
d24_9 <- split(d24_8, d24_8$RefPos)
d48_9 <- split(d48_8, d48_8$RefPos)
d72_9 <- split(d72_8, d72_8$RefPos)
d96_9 <- split(d96_8, d96_8$RefPos)
d120_9 <- split(d120_8, d120_8$RefPos)

#wilcoxon test on all pairs
pvals_e8 <- sapply(list7, function(x) {wilcox.test(d8_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e16 <- sapply(list7, function(x) {wilcox.test(d16_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e24 <- sapply(list7, function(x) {wilcox.test(d24_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e48 <- sapply(list7, function(x) {wilcox.test(d48_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e72 <- sapply(list7, function(x) {wilcox.test(d72_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e96 <- sapply(list7, function(x) {wilcox.test(d96_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})
pvals_e120 <- sapply(list7, function(x) {wilcox.test(d120_9[[x]]$IPD, dam_9[[x]]$IPD)$p.value})

#plot that shit
layout(matrix(1:8, 2, byrow=TRUE))
plot(-log10(pvals_e8), main="wilcoxon test, 8hr v. dam, exclusive", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e8)), col='red')

plot(-log10(pvals_e16), main="wilcoxon test, 16hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e16)), col='red')

plot(-log10(pvals_e24), main="wilcoxon test, 24hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e24)), col='red')

plot(-log10(pvals_e48), main="wilcoxon test, 48hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e48)), col='red')

plot(-log10(pvals_e72), main="wilcoxon test, 72hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e72)), col='red')

plot(-log10(pvals_e96), main="wilcoxon test, 96hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e96)), col='red')

plot(-log10(pvals_e120), main="wilcoxon test, 120hr v. dam", ylab="-log10(p-values)", xlim=c(0,20500), ylim=c(0,15))
abline(h=-log10(0.05/length(pvals_e120)), col='red')

######CODE from heatpmap.R

d_m <- data.frame(pos = as.numeric(names(dam_2)), means = unlist(lapply(dam_2, function(i) {mean(i$IPD)})))
d8_m <- data.frame(pos = as.numeric(names(d8_2)), means = unlist(lapply(d8_2, function(i) {mean(i$IPD)})))
d16_m <- data.frame(pos = as.numeric(names(d16_2)), means = unlist(lapply(d16_2, function(i) {mean(i$IPD)})))
d24_m <- data.frame(pos = as.numeric(names(d24_2)), means = unlist(lapply(d24_2, function(i) {mean(i$IPD)})))
d48_m <- data.frame(pos = as.numeric(names(d48_2)), means = unlist(lapply(d48_2, function(i) {mean(i$IPD)})))
d72_m <- data.frame(pos = as.numeric(names(d72_2)), means = unlist(lapply(d72_2, function(i) {mean(i$IPD)})))
d96_m <- data.frame(pos = as.numeric(names(d96_2)), means = unlist(lapply(d96_2, function(i) {mean(i$IPD)})))
d120_m <- data.frame(pos = as.numeric(names(d120_2)), means = unlist(lapply(d120_2, function(i) {mean(i$IPD)})))

by <- c('pos', 'pos')
dall_m <- merge(merge(merge(merge(merge(merge(merge(d_m, d8_m, by = by), d16_m, by = by), d24_m, by = by), d48_m, by = by), d72_m, by = by), d96_m, by = by), d120_m, by = by)
names(dall_m) <- c('pos', 'ct', 'tp8', 'tp16', 'tp24', 'tp48', 'tp72', 'tp96', 'tp120')

heatmap(as.matrix(dall_m[,2:9]), NA)
heatmap(as.matrix(dall_m[,3:9]/dall_m[,2]), NA)

pdf('timecourses.pdf')
plots <- apply(dall_m, 1, function(rd) {
	plot(rd[3:9]/rd[2], type = 'b', xlab = 'time point', ylab = 'IPD ratio', main = rd[1], ylim = c(0,100), xaxt ='n')
	axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
})
dev.off()

layout(rbind(c(1:2)))
hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', main = '')
hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', ylim = c(0,50), main = '')

d_m <- data.frame(pos = as.numeric(names(dam_2)), means = unlist(lapply(dam_2, function(i) {median(i$IPD)})))
d8_m <- data.frame(pos = as.numeric(names(d8_2)), means = unlist(lapply(d8_2, function(i) {median(i$IPD)})))
d16_m <- data.frame(pos = as.numeric(names(d16_2)), means = unlist(lapply(d16_2, function(i) {median(i$IPD)})))
d24_m <- data.frame(pos = as.numeric(names(d24_2)), means = unlist(lapply(d24_2, function(i) {median(i$IPD)})))
d48_m <- data.frame(pos = as.numeric(names(d48_2)), means = unlist(lapply(d48_2, function(i) {median(i$IPD)})))
d72_m <- data.frame(pos = as.numeric(names(d72_2)), means = unlist(lapply(d72_2, function(i) {median(i$IPD)})))
d96_m <- data.frame(pos = as.numeric(names(d96_2)), means = unlist(lapply(d96_2, function(i) {median(i$IPD)})))
d120_m <- data.frame(pos = as.numeric(names(d120_2)), means = unlist(lapply(d120_2, function(i) {median(i$IPD)})))

by <- c('pos', 'pos')
dall_m <- merge(merge(merge(merge(merge(merge(merge(d_m, d8_m, by = by), d16_m, by = by), d24_m, by = by), d48_m, by = by), d72_m, by = by), d96_m, by = by), d120_m, by = by)
names(dall_m) <- c('pos', 'ct', 'tp8', 'tp16', 'tp24', 'tp48', 'tp72', 'tp96', 'tp120')

heatmap(as.matrix(dall_m[,2:9]), NA)
heatmap(as.matrix(dall_m[,3:9]/dall_m[,2]), NA)

pdf('timecourses_median.pdf')
plots <- apply(dall_m, 1, function(rd) {
	plot(rd[3:9]/rd[2], type = 'b', xlab = 'time point', ylab = 'IPD ratio', main = rd[1], ylim = c(0,50), xaxt ='n')
	axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
})
dev.off()

layout(rbind(c(1:2)))
hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', main = '')
hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', ylim = c(0,50), main = '')


###CODE FROM Histograms for ALL IPD.R

#create histogram of all IPDS for data.frame but don't plot it
dam_hist <- hist(dam$IPD, plot=F, breaks=seq(0,1000,25))
d8_hist <- hist(d8$IPD, plot=F, breaks=seq(0,1000,25))
d16_hist <- hist(d16$IPD, plot=F, breaks=seq(0,1000,25))
d24_hist <- hist(d24$IPD, plot=F, breaks=seq(0,1000,25))
d48_hist <- hist(d48$IPD, plot=F, breaks=seq(0,1000,25))
d72_hist <- hist(d72$IPD, plot=F, breaks=seq(0,1000,25))
d96_hist <- hist(d96$IPD, plot=F, breaks=seq(0,1000,25))
d120_hist <- hist(d120$IPD, plot=F, breaks=seq(0,1000,25))

#plot histogram output in plot() 
#plot(dam_hist$mids, dam_hist$counts/length(dam$IPD), type = "b", ylim = c(0,1),xlim=c(0,1000), main = "dam All IPD hist", xlab = 'IPD bin', ylab = 'relative frequency')
#open pdf, plot both dam and 8hr plots on one page

#pdf("8hr_dam_IPD_hist.pdf")
#split.screen(c(2,1))
#screen(1)
#plot(dam_hist$mids, dam_hist$counts/length(dam$IPD), type = "b", ylim = c(0,1),xlim=c(0,1000), main = "dam All IPD hist", xlab = 'IPD bin', ylab = 'relative frequency')
#screen(2)
#plot(d8_hist$mids, d8_hist$counts/length(d8$IPD), type = "b", ylim = c(0,1),xlim=c(25,200), main = "8hr All IPD hist", xlab = 'IPD bin', ylab = 'relative frequency')
#close.screen(all=T)
#dev.off()

#plot one on top of the other
pdf("up to 24 and dam hist.pdf")
plot(dam1_hist$mids, dam1_hist$counts/length(dam_1$IPD), type = "b", ylim = c(0,1),xlim=c(0,1000), main = "8-24 and dam (+/-) All IPD hist", xlab = 'IPD bin', ylab = 'relative frequency')
points(dam2_hist$mids, dam2_hist$counts/length(dam_2$IPD), type = "b", col = "grey")
points(d82_hist$mids, d82_hist$counts/length(d8_2$IPD), type = "b", col = "pink")
points(d81_hist$mids, d81_hist$counts/length(d8_1$IPD), type = "b", col = "red")
points(d161_hist$mids, d161_hist$counts/length(d16_1$IPD), type = "b", col = "cyan")
points(d162_hist$mids, d162_hist$counts/length(d16_2$IPD), type = "b", col = "blue")
points(d242_hist$mids, d242_hist$counts/length(d24_2$IPD), type = "b", col = "green3")
points(d241_hist$mids, d241_hist$counts/length(d24_1$IPD), type = "b", col = "green")
points(d482_hist$mids, d482_hist$counts/length(d48_2$IPD), type = "b", col = "cornflowerblue")
points(d481_hist$mids, d481_hist$counts/length(d48_1$IPD), type = "b", col = "cadetblue")
points(d722_hist$mids, d722_hist$counts/length(d72_2$IPD), type = "b", col = "coral")
points(d721_hist$mids, d721_hist$counts/length(d72_1$IPD), type = "b", col = "coral3")
points(d962_hist$mids, d962_hist$counts/length(d96_2$IPD), type = "b", col = "darkgoldenrod")
points(d961_hist$mids, d961_hist$counts/length(d96_1$IPD), type = "b", col = "darkgoldenrod1")
points(d1202_hist$mids, d1202_hist$counts/length(d120_2$IPD), type = "b", col = "darkorchid")
points(d1201_hist$mids, d1201_hist$counts/length(d120_1$IPD), type = "b", col = "darkorchid4")
legends <- c("dam+", "dam-", "8hr+", "8hr-", "16hr+", "16hr-", "24hr+", "24hr-", "48hr+", "48hr-", "72hr+", "72hr-", "96hr+", "96hr-","120hr+", "120hr-" )
legend('topright', legends, lty=1, col=c('black', 'grey', 'pink', 'red', 'cyan', 'blue', 'green3', 'green', 'cornflowerblue', 'cadetblue', 'coral', 'coral3', 'darkgoldenrod', 'darkgoldenrod1', 'darkorchid', 'darkorchid4'), bty='n', cex=.75)
dev.off()

#t-test between dam and other timepoints
# first, get all unique positions (use split datasets)

positions8 <- intersect(names(dam_2), names(d8_2))
pvals8 <- sapply(positions, function(x) {wilcox.test(d8_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
positions16 <- intersect(names(dam_2), names(d16_2))
pvals16 <- sapply(positions, function(x) {wilcox.test(d16_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
positions24 <- intersect(names(dam_2), names(d24_2))
pvals24 <- sapply(positions, function(x) {wilcox.test(d24_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
positions48 <- intersect(names(dam_2), names(d48_2))
pvals48 <- sapply(positions, function(x) {wilcox.test(d48_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
positions72 <- intersect(names(dam_2), names(d72_2))
pvals72 <- sapply(positions, function(x) {wilcox.test(d72_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
positions96 <- intersect(names(dam_2), names(d96_2))
pvals96 <- sapply(positions, function(x) {wilcox.test(d96_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
positions120 <- intersect(names(dam_2), names(d120_2))
pvals120 <- sapply(positions, function(x) {wilcox.test(d120_2$IPD[[x]], dam_2$IPD[[x]])$p.value})

#plot pvals
#plot(-log10(pvals))
#hist(pvals)
#sum(pvals < 0.05/length(pvals))
#[1] 18631
#length(pvals)
#[1] 18726
#> 18631/18726
#[1] 0.9949268
#> abline(h = -log10(0.05/length(pvals)))

#make sure x range is covered by binning parameters
#d8 d8_3_hist <- hist(d8_3$IPD, plot=F, breaks=seq(0,1000,25))

pdf("8newHIST.pdf")
jnk<- lapply(positions, function(x) {
	i <- new8_3[[x]]
	j <- new8_4[[x]]
	i_p <- hist(i$IPD, plot = FALSE, breaks = seq(0,1000,10))
	j_p <- hist(j$IPD, plot = FALSE, breaks = seq(0,1000,10))
	plot(i_p$mids, i_p$counts/length(i$IPD), type = "b", xlim = c(0,1000), ylim = c(0,1), main = x, xlab = 'IPD bin', ylab = 'relative frequency')
	points(j_p$mids, j_p$counts/length(j$IPD), type = "b", col = "red")
})
dev.off()

pdf("test.pdf")
jnk<- lapply(five_split, function(i) {hist(i$IPD, xlim=c(0,300), ylim=c(0,50), main= i$RefPos[1])})
dev.off()


#### CODE from PCA analysis.R

heatmap(as.matrix(all_cov[,2:8]))
image(as.matrix(all_cov[,2:8]))
pca <- prcomp(t(all_cov[,2:8]))
plot(pca$rotation)
plot(pca$rotation[,2:3])
plot(pca$rotation[,1:3])
plot(pca$rotation[,1:4])
plot(pca$rotation[,1:5])
all_cov2 <- all_cov
all_cov2[,2:8] <- (all_cov2[,2:8] < 10)

image(as.matrix(all_cov[,2:8]))
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 15)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 20)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 11)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 12)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 15)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 16)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
all_cov2[,2:8] <- (all_cov[,2:8] < 17)
image(as.matrix(all_cov2[,2:8]), col = c("white", "black"))
plot(rowSums(all_cov[,2:8] < 16))
plot(rowSums(all_cov[,2:8] < 17))
plot(rowSums(all_cov[,2:8] < 17))
all_cov[,2:8] < 17
plot(rowSums(all_cov[,2:8] < 18))
plot(rowSums(all_cov[,2:8] < 19))
plot(rowSums(all_cov[,2:8] < 20))
plot(rowSums(all_cov[,2:8] < 21))
plot(rowSums(all_cov[,2:8] < 15))


sum(sapply(0:30, function(x) {rowSums(all_cov[,2:8] < x)}) > 3)
[1] 50800
plot(0:30, sapply(0:30, function(x) {sum(rowSums(all_cov[,2:8] < x) > 3)}))
plot(0:30, sapply(0:30, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(0:20, sapply(0:20, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(0:20, sapply(0:20, function(x) {sum(rowSums(all_cov[,2:8] < x) > 5)}))
plot(0:20, sapply(0:20, function(x) {sum(rowSums(all_cov[,2:8] < x) > 7)}))
plot(0:20, sapply(0:20, function(x) {sum(rowSums(all_cov[,2:8] < x) > 6)}))
plot(0:20, sapply(0:20, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(0:30, sapply(0:40, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))

plot(0:30, sapply(0:30, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 3)}))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 5)}))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 6)}))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 4)}))
plot(rowSums(all_cov[,2:8] < 17))
plot(rowSums(all_cov[,2:8] < 16))
plot(0:25, sapply(0:25, function(x) {sum(rowSums(all_cov[,2:8] < x) > 6)}))
sum(rowSums(all_cov[,2:8] < 16) > 6
+ )
[1] 39
sum(rowSums(all_cov[,2:8] < 16) > 3)
[1] 120
all_cov[rowSums(all_cov[,2:8] < 16) > 3,]
     position dam X8hr X16hr X24hr X48hr X96hr X120hr
9259  1109622  25   15     5     8    18    27     15
9269  1111101  33   26     7    12     9    13     13
9271  1111539  33   28     6    10    12    11     19
9274  1111565  34   20    12    12    15    13     24
9278  1113232  33   20    10     7    10    15     18
9280  1113376  26   13     7     6    10    13     18
9285  1115421  13   26     9    17     5    21     11
9286  1115422  15   14     9    11     5    19     12
9288  1116040  18   11     9    21     6    31     14
9294  1116970  14   15     8    17    12    21     15
9301  1117824  24   11    15    16    14    15     22
9303  1117933  27   12    14    18    13    12     19
9304  1117934  22    8    13    14    14    20     25
9305  1118061  21    9    14    16    14    14     23
9309  1118433  28    9    11    16    15    20     13
9311  1118507  30   11     9    14    21    14     17
9312  1118508  26   12     6    15    17    22     15
9315  1118772  19    9     9    12    10    17     11
9317  1118860  22    9    10    15    12    16     13
9318  1118861  11    7     7    10    16    19     15
9319  1120353  13    9    12    13    16    12      9
9320  1120354  14   11    11    13    10    15      9
9321  1120404  13    9    15    15    16    13      7
9322  1120405  13   11    13    13     8    16      7
9323  1120698  11   10    14    12    14    13      8
9324  1120699  12   10    13    12    15    18      4
9325  1121278   9   11    15    14    18     6      5
9326  1121279  13   13     8    17    15    12      9
9327  1121580  11   10    13     8    15     8     10
9328  1121581   9   12    13     8    13    12     12
9329  1121610  12    9    10    10    11     5      7
9330  1121611  11   14     9    10    11     9     10
9331  1122100  15   10     8    11    11    11     10
9332  1122101  12   10     8    11    12    11     11
9333  1122166  15   12     8    13    10     8      8
9334  1122167  10   13     7    16     8    10     13
9335  1122186  18   10     6    14     9    11      9
9336  1122187  16   10     6    13     8    15     11
9337  1122357  16    7     6    12    11    12     11
9338  1122358  18   11     7    13    11    15     10
9339  1122548  17    7     7    16     8     9      9
9340  1122549  18    7     6    15    11    12     13
9341  1122644  14    2     5    10     9     7     10
9342  1122645  16    5     5     9    10    13     10
9343  1122690  17    9     6    11    10     7      8
9344  1122691  16    6     5    17    12    13     10
9345  1122760  12    9     6    12    10     9     10
9346  1122761  16    6     6    14    11    10     10
9347  1122913  15   11     5     8    12     6      9
9348  1122914  17   10     5    12    14     9     10
9349  1123360  11    8     4    10     3     8      5
9350  1123361  11    8     7    16     7    11      9
9351  1123774  10   11     4    12     7    10      5
9352  1123775   6    6     4    10     8    14      7
9353  1123873  12   16     4     9     6     6      6
9354  1123874  11    3     6    10     8    10      6
9355  1124051  15   11     4     9     3     5      5
9356  1124052  11    5     4     9     5     9      7
9357  1124075  13   12     7     7     4     3      4
9358  1124076   9    6     5     9     7     7      6
9359  1124602  11    6     2     7     2     3      4
9360  1124603  10    5     4     7     2     3      3
9361  1125947   1    2     1     5     2     4      2
9362  1126063   3    2     2     5     2     6      1
9363  1126108   4    5     3     6     4     1      1
9364  1126109   5    3     2     5     2     4      2
9365  1126560   4    2     4     8    10     2      9
9366  1126561   4    3     3     3     4     4      7
9367  1127202   2    4     4     8    14     8     14
9368  1127203   3    5     6     9     7     9     11
9369  1127325   3    6     9    11    14    18     18
9370  1127326   3    5     5     8     7    14     11
9371  1127577   6    3     8     7    10    15     12
9372  1127578   5    7     7    12     6    19     10
9373  1127853   8    5     9     9    10    16     15
9374  1127854   8   10    10    11     9    14      8
9375  1127999  10    4    14    12     7    14     11
9376  1128000  10    6     8     8     8    15     14
9377  1128251  11    4    10    12     8     9     11
9378  1128252   8   10    12    11     5    11      7
9379  1128365  11    6    13    12     8    12     11
9380  1128366  10   13    16    11     8    11     11
9381  1128946  17   11    16    18    12    13     11
9382  1128947  12   14    17    16    11    11      8
9383  1129116  21   12    20    16    13    13      9
9384  1129117  17   16    14    11    10     9      9
9386  1129551  21   17    14    12     5    25     11
9387  1130144  33   10    20    15    11    18     10
9390  1130563  19   13    11    11     7    17     16
9391  1130796  20   20     9    15    10    14     11
9392  1130797  18   17     7    10     4    19     12
9393  1131024  22   17     9    10    13    21      8
9394  1131025  16   17    10     9     6    26     11
9396  1131045  17   18     6    11     5    21     12
9397  1131073  26   15    11    12    12    20     10
9398  1131074  17   20     6    13     5    22     12
9399  1131360  16   13    10    11    11    16     11
9400  1131361  10   20     6    11     3    20     10
9402  1131508  12   22     8     9     4    18     12
9403  1131643  15   14     7    12    15    19     12
9404  1131644  10   17     8    12     4    22     11
9405  1131787  20   13     5     7     8    16     14
9406  1131788  12   21     4    11     6    19     12
9407  1132510  17   13     9    13     9    20     16
9408  1132511  14   18     4    13     8    14     14
9409  1132585  16   12     9    17    10    17     14
9410  1132586  20   18     6    12    11    13     13
9411  1132912  12   12     6    17     8    14      9
9412  1132913  12   19     7    13    10    13     10
9413  1133122   6   11     7    10     5    17      8
9414  1133123   7   10     6     9    10    10     11
9415  1133785  19   12     7    11     6    17      8
9416  1133786  23   11     3    13     8    12     10
9417  1134864  13   14     6    20    23    22     12
9419  1134889  15   12     7    23    21    19     15
9423  1135396  14   11     8    17    18    20     11
9424  1135397  21   15     7    19    19    15     13
9429  1137442  20   16    10    13    12    14     18
9430  1137443  22   13    10    13    16    13     20
9449  1140082  13   20     4    13    17    21     15
> 



####CODE from timeseries_allTs.R

#by <- c('pos', 'pos')
#dall_m <- merge(merge(merge(merge(merge(merge(merge(d_m, d8_m, by = by), d16_m, by = by), d24_m, by = by), d48_m, by = by), d72_m, by = by), d96_m, by = by), d120_m, by = by)
#names(dall_m) <- c('pos', 'ct', 'tp8', 'tp16', 'tp24', 'tp48', 'tp72', 'tp96', 'tp120')

#heatmap(as.matrix(dall_m[,2:9]), NA)
#heatmap(as.matrix(dall_m[,3:9]/dall_m[,2]), NA)

#pdf('timecourses.pdf')
#plots <- apply(dall_m, 1, function(rd) {
#	plot(rd[3:9]/rd[2], type = 'b', xlab = 'time point', ylab = 'IPD ratio', main = rd[1], ylim = c(0,100), xaxt ='n')
#	axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
#})
#dev.off()

#layout(rbind(c(1:2)))
#hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', main = '')
#hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', ylim = c(0,50), main = '')

d_m <- data.frame(pos = as.numeric(names(dam_2)), means = unlist(lapply(dam_2, function(i) {median(i$IPD)})))
d8_m <- data.frame(pos = as.numeric(names(d8_2)), means = unlist(lapply(d8_2, function(i) {median(i$IPD)})))
d16_m <- data.frame(pos = as.numeric(names(d16_2)), means = unlist(lapply(d16_2, function(i) {median(i$IPD)})))
d24_m <- data.frame(pos = as.numeric(names(d24_2)), means = unlist(lapply(d24_2, function(i) {median(i$IPD)})))
d48_m <- data.frame(pos = as.numeric(names(d48_2)), means = unlist(lapply(d48_2, function(i) {median(i$IPD)})))
d72_m <- data.frame(pos = as.numeric(names(d72_2)), means = unlist(lapply(d72_2, function(i) {median(i$IPD)})))
d96_m <- data.frame(pos = as.numeric(names(d96_2)), means = unlist(lapply(d96_2, function(i) {median(i$IPD)})))
d120_m <- data.frame(pos = as.numeric(names(d120_2)), means = unlist(lapply(d120_2, function(i) {median(i$IPD)})))

by <- c('pos', 'pos')
dall_m <- merge(merge(merge(merge(merge(merge(merge(d_m, d8_m, by = by), d16_m, by = by), d24_m, by = by), d48_m, by = by), d72_m, by = by), d96_m, by = by), d120_m, by = by)
names(dall_m) <- c('pos', 'ct', 'tp8', 'tp16', 'tp24', 'tp48', 'tp72', 'tp96', 'tp120')

#heatmap(as.matrix(dall_m[,2:9]), NA)
#heatmap(as.matrix(dall_m[,3:9]/dall_m[,2]), NA)

pca_tp <- prcomp(dall_m[,3:9])
heatmap(pca_tp$rotation, scale = "none")

anova <- lapply(all_mc$RefPos, function(x) {
	p <- as.character(x)
	meas_vec <- c(d8_m[[p]]$Mean_IPD, d16_m[[p]]$Mean_IPD, d24_m[[p]]$Mean_IPD, d48_m[[p]]$Mean_IPD, d72_m[[p]]$Mean_IPD, 	d96_m[[p]]$Mean_IPD, d120_m[[p]]$Mean_IPD)
	group_vec <- c(d8_m$tp, d16_m$tp, d24_m$tp, d48_m$tp , d72_m$tp, d96_m$tp, d120_m$tp 
	kruskal.test(meas_vec, group_vec)
})
pval <- unlist(lapply(anova, function(x) {x$p.value}))
plot(-log10(pval))

library(qvalue)
qobj <- qvalue(pval, fdr.level = 0.05)

layout(matrix(1:2, nrow = 2))
plot(dall_m$pos, -log10(qobj$qvalues), pch = 20, cex = 0.5, xlab = "position in genome", ylab = "-log10(q)")
abline(h=-log10(0.05), lty = 2)
abline(h=-log10(0.1), lty = 2)
abline(h=-log10(0.2), lty = 2)
abline(h=-log10(0.5), lty = 2)

col <- c(1:6,8)
tpsig <- dall_m[dall_m$pos %in% dall_m$pos[qobj$qvalues <= 0.05],]
plot(0, las = 2,  xlim = c(1,7), ylim = c(0,300), type = "n", ylab = "median IPD", xlab = "timepoint", xaxt = "n")
apply(tpsig, 1, function(rd) {points(1:7, rd[3:9], type = "b", col = col[which(tpsig == rd[1])], pch = 20)})
#points(1:7, colMeans(tp120low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
legend("topright", legend = tpsig$pos, col = col, pch = 20)









plot(0, xlim = c(0,max(dall_m$pos)), ylim = c(1,7), type = "n", ylab = "low in timee point", xlab = "genome position", yaxt = "n")
points(tp8low$pos, rep(7, nrow(tp8low)), pch = '|')
points(tp16low$pos, rep(6, nrow(tp16low)), pch = '|')
points(tp24low$pos, rep(5, nrow(tp24low)), pch = '|')
points(tp48low$pos, rep(4, nrow(tp48low)), pch = '|')
points(tp72low$pos, rep(3, nrow(tp72low)), pch = '|')
points(tp96low$pos, rep(2, nrow(tp96low)), pch = '|')
points(tp120low$pos, rep(1, nrow(tp120low)), pch = '|')
axis(2,1:7, rev(c(8, 16, 24, 48, 72, 96, 120)), las = 2)








pca <- prcomp(t(dall_m[,3:9]))
plot(pca$rotation)
text(pca$rotation[,2], pca$rotation[,3], dimnames(pca$rotation)[[1]])

#pdf('timecourses_median.pdf')
#plots <- apply(dall_m, 1, function(rd) {
#	plot(rd[3:9]/rd[2], type = 'b', xlab = 'time point', ylab = 'IPD ratio', main = rd[1], ylim = c(0,50), xaxt ='n')
#	axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
#})
#dev.off()

layout(rbind(c(1:2)))
hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', main = '')
hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', ylim = c(0,50), main = '')


###CODE From timeseries.R

#by <- c('pos', 'pos')
#dall_m <- merge(merge(merge(merge(merge(merge(merge(d_m, d8_m, by = by), d16_m, by = by), d24_m, by = by), d48_m, by = by), d72_m, by = by), d96_m, by = by), d120_m, by = by)
#names(dall_m) <- c('pos', 'ct', 'tp8', 'tp16', 'tp24', 'tp48', 'tp72', 'tp96', 'tp120')

#heatmap(as.matrix(dall_m[,2:9]), NA)
#heatmap(as.matrix(dall_m[,3:9]/dall_m[,2]), NA)

#pdf('timecourses.pdf')
#plots <- apply(dall_m, 1, function(rd) {
#	plot(rd[3:9]/rd[2], type = 'b', xlab = 'time point', ylab = 'IPD ratio', main = rd[1], ylim = c(0,100), xaxt ='n')
#	axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
#})
#dev.off()

#layout(rbind(c(1:2)))
#hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', main = '')
#hist(unlist(dall_m[,3:9]/dall_m[,2]), breaks = 1000, xlab = 'IPD ratio', ylim = c(0,50), main = '')

d_m <- data.frame(pos = as.numeric(names(dam_2)), means = unlist(lapply(dam_2, function(i) {median(i$IPD)})))
d8_m <- data.frame(pos = as.numeric(names(d8_2)), means = unlist(lapply(d8_2, function(i) {median(i$IPD)})))
d16_m <- data.frame(pos = as.numeric(names(d16_2)), means = unlist(lapply(d16_2, function(i) {median(i$IPD)})))
d24_m <- data.frame(pos = as.numeric(names(d24_2)), means = unlist(lapply(d24_2, function(i) {median(i$IPD)})))
d48_m <- data.frame(pos = as.numeric(names(d48_2)), means = unlist(lapply(d48_2, function(i) {median(i$IPD)})))
d72_m <- data.frame(pos = as.numeric(names(d72_2)), means = unlist(lapply(d72_2, function(i) {median(i$IPD)})))
d96_m <- data.frame(pos = as.numeric(names(d96_2)), means = unlist(lapply(d96_2, function(i) {median(i$IPD)})))
d120_m <- data.frame(pos = as.numeric(names(d120_2)), means = unlist(lapply(d120_2, function(i) {median(i$IPD)})))

by <- c('pos', 'pos')
dall_m <- merge(merge(merge(merge(merge(merge(merge(d_m, d8_m, by = by), d16_m, by = by), d24_m, by = by), d48_m, by = by), d72_m, by = by), d96_m, by = by), d120_m, by = by)
names(dall_m) <- c('pos', 'ct', 'tp8', 'tp16', 'tp24', 'tp48', 'tp72', 'tp96', 'tp120', 'WGA')

m1 <- merge(d_m, d8_m, by='pos')
names(m1) <- c('pos', 'dam', 'd8')

m2 <- merge(m1, d16_m, by='pos')
names(m2) <- c('pos', 'dam', 'd8', 'd16')

m3 <- merge(m2, d24_m, by='pos')
names(m3) <- c('pos', 'dam', 'd8', 'd16', 'd24')

m4 <- merge(m3, d48_m, by='pos')
names(m4) <- c('pos', 'dam', 'd8', 'd16', 'd24', 'd48')

m5 <- merge(m4, d72_m, by='pos')
names(m5) <- c('pos', 'dam', 'd8', 'd16', 'd24', 'd48', 'd72')

m6 <- merge(m5, d96_m, by='pos')
names(m6) <- c('pos', 'dam', 'd8', 'd16','d24', 'd48', 'd72', 'd96')

m7 <- merge(m6, d120_m, by='pos')
names(m7) <- c('pos', 'dam', 'd8', 'd16','d24', 'd48', 'd72', 'd96', 'd120')

m8 <- merge(m7, WGA_m, by='pos')
names(m8) <- c('pos', 'dam', 'd8', 'd16','d24', 'd48', 'd72', 'd96', 'd120', 'WGA')

#heatmap(as.matrix(dall_m[,2:9]), NA)
#heatmap(as.matrix(dall_m[,3:9]/dall_m[,2]), NA)

pca_tp <- prcomp(dall_m[,3:9])
heatmap(pca_tp$rotation, scale = "none")

layout(matrix(1:8, ncol=2, byrow = T))

tp8low <- m8[m8$d8 == apply(m8[,3:9], 1, min),]
plot(0, las = 2, xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint8", xaxt = "n")
apply(tp8low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp8low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

tp16low <- m8[m8$d16 == apply(m8[,3:9], 1, min),]
plot(0, las = 2,  xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint16", xaxt = "n")
apply(tp16low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp16low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

tp24low <- m8[m8$d24 == apply(m8[,3:9], 1, min),]
plot(0, las = 2,  xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint24", xaxt = "n")
apply(tp24low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp24low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

tp48low <- m8[m8$d48 == apply(m8[,3:9], 1, min),]
plot(0, las = 2,  xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint", xaxt = "n")
apply(tp48low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp48low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

tp72low <- m8[m8$d72 == apply(m8[,3:9], 1, min),]
plot(0, las = 2,  xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint", xaxt = "n")
apply(tp72low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp72low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

tp96low <- m8[m8$d96 == apply(m8[,3:9], 1, min),]
plot(0, las = 2,  xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint", xaxt = "n")
apply(tp96low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp96low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

tp120low <- m8[m8$d120 == apply(m8[,3:9], 1, min),]
plot(0, las = 2,  xlim = c(1,7), ylim = c(25,400), type = "n", ylab = "median IPD", xlab = "timepoint", xaxt = "n")
apply(tp120low, 1, function(rd) {points(1:7, rd[3:9], type = "l", col = "lightgrey")})
points(1:7, colMeans(tp120low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))

barplot(c(nrow(tp8low), nrow(tp16low), nrow(tp24low), nrow(tp48low), nrow(tp72low), nrow(tp96low), nrow(tp120low)), names = c(8, 16, 24, 48, 72, 96, 120), col = "lightgrey", ylab = "number of sites in class")

anova <- lapply(m8$pos, function(x) {
	p <- as.character(x)
	meas_vec <- c(d8_2[[p]]$IPD, d16_2[[p]]$IPD, d24_2[[p]]$IPD, d48_2[[p]]$IPD, d72_2[[p]]$IPD, d96_2[[p]]$IPD, d120_2[[p]]$IPD)
	group_vec <- c(rep(8, length(d8_2[[p]]$IPD)), rep(16, length(d16_2[[p]]$IPD)), rep(24, length(d24_2[[p]]$IPD)), rep(48, length(d48_2[[p]]$IPD)), rep(72, length(d72_2[[p]]$IPD)), rep(96, length(d96_2[[p]]$IPD)), rep(120, length(d120_2[[p]]$IPD)))
	kruskal.test(meas_vec, group_vec)
})
pval <- unlist(lapply(anova, function(x) {x$p.value}))
plot(-log10(pval))

library(qvalue)
qobj <- qvalue(pval, 0.05)

layout(matrix(1:2, nrow = 2))
plot(m8$pos, -log10(qobj$qvalues), pch = 20, cex = 0.5, xlab = "position in genome", ylab = "-log10(q)")
abline(h=-log10(0.05), lty = 2)
abline(h=-log10(0.1), lty = 2)
abline(h=-log10(0.2), lty = 2)
abline(h=-log10(0.5), lty = 2)

col <- c(1:6,8)
tpsig <- m8[m8$pos %in% m8$pos[qobj$qvalues <= 0.05],]
plot(0, las = 2,  xlim = c(1,7), ylim = c(0,300), type = "n", ylab = "median IPD", xlab = "timepoint", xaxt = "n")
apply(tpsig, 1, function(rd) {points(1:7, rd[3:9], type = "b", col = col[which(tpsig == rd[1])], pch = 20)})
#points(1:7, colMeans(tp120low[,3:9]), type = "l", lwd = 3)
axis(1, 1:7, c(8, 16, 24, 48, 72, 96, 120))
legend("topright", legend = tpsig$pos, col = col, pch = 20)

plot(0, xlim = c(0,max(dall_m$pos)), ylim = c(1,7), type = "n", ylab = "low in time point", xlab = "genome position", yaxt = "n")
points(tp8low$pos, rep(7, nrow(tp8low)), pch = '|')
points(tp16low$pos, rep(6, nrow(tp16low)), pch = '|')
points(tp24low$pos, rep(5, nrow(tp24low)), pch = '|')
points(tp48low$pos, rep(4, nrow(tp48low)), pch = '|')
points(tp72low$pos, rep(3, nrow(tp72low)), pch = '|')
points(tp96low$pos, rep(2, nrow(tp96low)), pch = '|')
points(tp120low$pos, rep(1, nrow(tp120low)), pch = '|')
axis(2,1:7, rev(c(8, 16, 24, 48, 72, 96, 120)), las = 2)


####CODE from Wilcox.R

positions8 <- intersect(names(dam_2), names(d8_2))
pvals8 <- sapply(positions, function(x) {wilcox.test(d8_2$IPD[[x]], dam_2$IPD[[x]])$p.value})
list.condition <- sapply(pvals, function(x) x >= 0.005)
output.list <- pvals[list.condition]
output.list


# calculate IPD ratio in GATC sites
means <- fread("all_means_GATC_raw.txt")
dam <- fread("dam_trimmedGATC.txt")
dam_2 <- dlply(dam, "RefPos")
means_dam <- lapply(names(dam_2), function(x) {mean(dam_2[[x]]$IPD)})
means$dam_m <- unlist(means_dam)
names <- means$RefPos
intersect <- intersect(names(dam_2), names)
int_dam_m <- dam_2[dam_2$RefPos %in% intersect]
means_dam <- lapply(as.character(intersect), function(x) {mean(dam_2[[x]]$IPD)})
means_int <- means[means$RefPos %in% intersect,]
means_int$dam_m <- unlist(means_dam)


IPDr24 <- lapply(1:length(means_int$RefPos), function(i) {means_int$d24_m[i]/means_int$dam_m[i]})
means_int$ratio24 <- unlist(IPDr24)
IPDr48 <- lapply(1:length(means_int$RefPos), function(i) {means_int$d48_m[i]/means_int$dam_m[i]})
means_int$ratio48 <- unlist(IPDr48)
IPDr72 <- lapply(1:length(means_int$RefPos), function(i) {means_int$d72_m[i]/means_int$dam_m[i]})
means_int$ratio72 <- unlist(IPDr72)
IPDr96 <- lapply(1:length(means_int$RefPos), function(i) {means_int$d96_m[i]/means_int$dam_m[i]})
means_int$ratio96 <- unlist(IPDr96)
IPDr120 <- lapply(1:length(means_int$RefPos), function(i) {means_int$d120_m[i]/means_int$dam_m[i]})
means_int$ratio120 <- unlist(IPDr120)

#IPD ratio histograms
hist(means_int$ratio8, breaks=100, col='salmon', main="8hr IPD ratio, GATC", xlab="IPD ratio")
hist(means_int$ratio16, breaks=100, col='salmon', main="16hr IPD ratio, GATC", xlab="IPD ratio")
hist(means_int$ratio24, breaks=100, col='salmon', main="24hr IPD ratio, GATC", xlab="IPD ratio")
hist(means_int$ratio48, breaks=100, col='salmon', main="48hr IPD ratio, GATC", xlab="IPD ratio")
hist(means_int$ratio72, breaks=100, col='salmon', main="72hr IPD ratio, GATC", xlab="IPD ratio")
hist(means_int$ratio96, breaks=100, col='salmon', main="96hr IPD ratio, GATC", xlab="IPD ratio")
hist(means_int$ratio120, breaks=100, col='salmon', main="120hr IPD ratio, GATC", xlab="IPD ratio")

hist8 <- hist(means_int$ratio8, breaks=100, plot=F)
hist16 <- hist(means_int$ratio16, breaks=100, plot=F)
hist24 <- hist(means_int$ratio24, breaks=100, plot=F)
hist48 <- hist(means_int$ratio48, breaks=100, plot=F)
hist72 <- hist(means_int$ratio72, breaks=100, plot=F)
hist96 <- hist(means_int$ratio96, breaks=100, plot=F)
hist120 <- hist(means_int$ratio120, breaks=100, plot=F)

plot(hist8$mids, hist8$counts/37446, xlim=c(0,100), type='b', col=1, pch=20)
points(hist16$mids, hist16$counts/37446, xlim=c(0,100), ylim=c(0,7000), type='b', col=2, pch=20)
points(hist24$mids, hist24$counts/37446, xlim=c(0,100), ylim=c(0,7000), type='b', col=3, pch=20)
points(hist48$mids, hist48$counts/37446, xlim=c(0,100), ylim=c(0,7000), type='b', col=4, pch=20)
points(hist72$mids, hist72$counts/37446, xlim=c(0,100), ylim=c(0,7000), type='b', col=5, pch=20)
points(hist96$mids, hist96$counts/37446, xlim=c(0,100), ylim=c(0,7000), type='b', col=6, pch=20)
points(hist120$mids, hist120$counts/37446, xlim=c(0,100), ylim=c(0,7000), type='b', col=8, pch=20)


site4 <- data.frame(tp=c(8, 16, 24, 48, 72, 96, 120), IPDr = c(site4$ratio8, site4$ratio16, site4$ratio24, site4$ratio48, site4$ratio72, site4$ratio96, site4$ratio120))
site5 <- data.frame(tp=c(8, 16, 24, 48, 72, 96, 120), IPDr = c(site5$ratio8, site5$ratio16, site5$ratio24, site5$ratio48, site5$ratio72, site5$ratio96, site5$ratio120))
site6 <- data.frame(tp=c(8, 16, 24, 48, 72, 96, 120), IPDr = c(site6$ratio8, site6$ratio16, site6$ratio24, site6$ratio48, site6$ratio72, site6$ratio96, site6$ratio120))
site7 <- data.frame(tp=c(8, 16, 24, 48, 72, 96, 120), IPDr = c(site7$ratio8, site7$ratio16, site7$ratio24, site7$ratio48, site7$ratio72, site7$ratio96, site7$ratio120))
site8 <- data.frame(tp=c(8, 16, 24, 48, 72, 96, 120), IPDr = c(site8$ratio8, site8$ratio16, site8$ratio24, site8$ratio48, site8$ratio72, site8$ratio96, site8$ratio120))
site9 <- data.frame(tp=c(8, 16, 24, 48, 72, 96, 120), IPDr = c(site9$ratio8, site9$ratio16, site9$ratio24, site9$ratio48, site9$ratio72, site9$ratio96, site9$ratio120))

#Looking deeper into Wilcoxon data
setwd("~/Desktop/Lacey/Bioinformatics/Wilcoxon Test")
library(dlply)
library(data.table)
library(qvalue)
library(plyr)

pvals <- fread("PVal_allT_8hrWGA.txt")
dam_sites <- fread("~/Desktop/Lacey/Bioinformatics/dam_sites.txt")
d16 <- fread("~/Desktop/Lacey/Bioinformatics/8hr_trimmed.txt")
d16_GATC <- d16[d16$RefPos %in% dam_sites$V2,]
d16_GATC2 <- dlply(d16_GATC, "RefPos")
cov <- lapply(names(d16_GATC2), function(i) {nrow(d16_GATC2[[i]])})
coverage <- data.frame(cov=unlist(cov), RefPos = names(d16_GATC2))


pvals$Qvalue <- qvalue(pvals$Pvalue)$qvalues
pvals_GATC <- pvals[pvals$RefPos %in% dam_sites$V2,]


cov_good<- coverage$RefPos[coverage$cov > 15]
pvals_cov <- pvals_GATC[pvals_GATC$RefPos %in% cov_good,]
length(pvals_cov$RefPos)
length(pvals_GATC$RefPos)
length(pvals_cov$RefPos[pvals_cov$Qvalue > 0.05])
length(pvals_cov$RefPos[pvals_cov$Qvalue > 0.05])/length(pvals_cov$RefPos)


setwd("~/Desktop/Lacey/Bioinformatics/Wilcoxon Test")
library(dlply)
library(data.table)
library(qvalue)
library(plyr)

pvals <- fread("PVal_allT_8hrdam.txt")
dam_sites <- fread("~/Desktop/Lacey/Bioinformatics/dam_sites.txt")
d16 <- fread("~/Desktop/Lacey/Bioinformatics/8hr_trimmed.txt")
d16_GATC <- d16[d16$RefPos %in% dam_sites$V2,]
d16_GATC2 <- dlply(d16_GATC, "RefPos")
cov <- lapply(names(d16_GATC2), function(i) {nrow(d16_GATC2[[i]])})
coverage <- data.frame(cov=unlist(cov), RefPos = names(d16_GATC2))


pvals$Qvalue <- qvalue(pvals$Pvalue)$qvalues
pvals_GATC <- pvals[pvals$RefPos %in% dam_sites$V2,]


cov_good<- coverage$RefPos[coverage$cov > 15]
pvals_cov <- pvals_GATC[pvals_GATC$RefPos %in% cov_good,]
insig_GATC8 <- pvals_cov[pvals_cov$Qvalue > 0.05,]

insig_8 <- insig_all$RefPos[insig_all$tp == 8]
low_thresh8 <- df95$RefPos[df95$d8 <= 0.2]
insig_16 <- insig_all$RefPos[insig_all$tp == 16]
low_thresh16 <- df95$RefPos[df95$d16 <= 0.2]
insig_24 <- insig_all$RefPos[insig_all$tp == 24]
low_thresh24 <- df95$RefPos[df95$d24 <= 0.2]
insig_48 <- insig_all$RefPos[insig_all$tp == 48]
low_thresh48 <- df95$RefPos[df95$d48 <= 0.2]
insig_72 <- insig_all$RefPos[insig_all$tp == 72]
low_thresh72 <- df95$RefPos[df95$d72 <= 0.2]
insig_96 <- insig_all$RefPos[insig_all$tp == 96]
low_thresh96 <- df95$RefPos[df95$d96 <= 0.2]
insig_120 <- insig_all$RefPos[insig_all$tp == 120]
low_thresh120 <- df95$RefPos[df95$d16 <= 0.2]

both8 <- intersect(insig_8, low_thresh8)
both16 <- intersect(insig_16, low_thresh16)
both24 <- intersect(insig_24, low_thresh24)
both48 <- intersect(insig_48, low_thresh48)
both72 <- intersect(insig_72, low_thresh72)
both96 <- intersect(insig_96, low_thresh96)
both120 <- intersect(insig_120, low_thresh120)

> intersect_qvals
1138061 elbA hypothetical (NCR)
1766534 ydjL oxidoreductase (NCR)
1766535 ydjL oxidoreductase (NCR)
1978045 isrC (flu upstream)
1978064 isrC
1978065 isrC
1978077 isrC
1978078 isrC		 
2511983 hyfB hydrogenase ORF
3488553 mobA ORF
3973874 glpR ORF
 
> intersect_thresh
  [1] "9181"    "323612"  "374578"  "645514"  "762254"  "986969"  "1086097" "1144189" "1247256" "1329270" "1455657" "1572748" "1766534" "1766535"
 [15] "1881669" "1978045" "1978064" "1978078" "2310316" "2467522" "3165021" "3366092" "3369821" "4023541" "4171432" "4346812" "4540551"
 
> intersection <- intersect(insig_GATC8$RefPos, insig_GATC16$RefPos)
> intersection
 [1]   44393   46981  171249  390716  585481  667959  667960 1038248 1105556 1108541 1115421 1121611 1122101 1123774 1123873 1128366 1131360
[18] 1134889 1135252 1138060 1138061 1139239 1139240 1534262 1743244 1766534 1766535 1978045 1978046 1978064 1978065 1978077 1978078 2152064
[35] 2508214 2511983 2527549 2575648 2734023 2734024 2897085 3036798 3187997 3237370 3265636 3333611 3428620 3444586 3456521 3456522 3488553
[52] 3555019 3555020 3604072 3748177 3964745 3973874 4125225 4196558 4215845 4233517 4267898 4403080 4439404 4512184
> intersection <- intersect(intersection, insig_GATC24$RefPos)

> intersection
 [1]  667959 1038248 1105556 1138061 1534262 1766534 1766535 1978045 1978046 1978064 1978065 1978077 1978078 2511983 2575648 2734023 2734024
[18] 2897085 3265636 3456521 3456522 3488553 3973874 4196558 4403080
> intersection <- intersect(intersection, insig_GATC48$RefPos)
> intersection
numeric(0)
> 

df[1,] <- alleles_pos[["54936"]]$mutant.allele.frequency


#plots
qval_insig <- c(233,184,244,317,132,146,171)
time <- c(8,16,24,48,72,96,120)
layout(matrix(1:4, 2, byrow=T))
barplot(qval_insig, main="# insignificant sites by time", xlab=" time point", ylab="# of sites", col='seagreen4', names.arg=time)
percent2 <- (1 - qval_insig/37454) * 100
percent <- (qval_insig/37454)*100
plot(time, percent2, type='b', pch=20, ylim=c(98,100), xlim=c(8,120), main="percent methylated by Wilcoxon test", xlab="time point", ylab="percent methylated")
sub <- c(37217,37266,37206,37133,37318,37304,37279)

plot(time, percent, type='b', pch=20, ylim=c(0.3,0.9), xlim=c(8,120), xlab=" time point", ylab="percent unmethylated", main="percent unmethylated by Wilcoxon method")
heatmap.2(as.matrix(intersection), dendrogram='none', Colv=F, col=col,Rowv=F, trace='none', na.color="darkgrey")



d8 <- c(233, NA, NA, NA, NA, NA, NA)
d16 <- c(65, 184, NA, NA, NA, NA, NA)
d24 <- c(54, 44, 244, NA, NA, NA, NA)
d48 <- c(70,49,47,317, NA, NA, NA)
d72 <- c(55,41,34,40,132, NA, NA)
d96<-c(50,47,44,45,43,146, NA)
d120 <- c(11,40,37,35,33,47,171)
intersection <- data.frame(d8=d8, d16=d16, d24=d24, d48=d48, d72=d72, d96=d96, d120=d120)
rownames(intersection) <- time
image(time, time,  as.matrix(intersection), dendrogram='none', Colv=F, col=col,Rowv=F, trace='none', na.color="darkgrey",  main="# of intersecting non-significant sites ")


plot(0,0, xlim=c(8,120), ylim=c(0,1), xlab="time point", ylab="fraction methylated")
sapply(names_d8_10, function(i) {points(all_2[[i]]$tp, all_2[[i]]$meth, type='b', pch=20)})
