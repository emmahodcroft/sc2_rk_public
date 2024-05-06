### NOTE this script is NOT used to generate the numbers for the manuscript BECAUSE it uses the 'byMonMedian' folder which is NOT
# what's used by Martin for the analysis. He instead uses the entire (un-sub-divided) files from 'output-5Apr-4muts' and 
# divides by median cluster date himself. Because of differences in the threadholds used by my splitting code (filter_clusters.py)
# and his code that splits the clusters by month, there are differences in the total counts!

# Since Martin's way of splitting is the one used for the analysis, *his* way of counting should be the one used!!


setwd("C:/Users/Emma/wsl/corona/sc2_rk/output-5Apr-4muts_byMonMedian")

folders = list.files()

countries <- c("Denmark", "Switzerland", "Germany")

totDK <- c()
totCH <- c()
totDE <- c()

country <- c()
month <- c()
totalN <- c()
maxSize <- c()
minSize <- c()
meanS <- c()
medianS <- c()

for (coun in countries) {
    for (fol in folders){
        if (substr(fol, 1, 4) == "2021"){
            fileN <- paste(fol, "/", coun, "_cluster_distribution_dates_100whole_", fol, ".tsv", sep="")
            clus <- read.csv(fileN, sep="\t", as.is=T)
            country <- c(country, coun)
            month <- c(month, fol)
            totalN <- c(totalN, nrow(clus))
            maxSize <- c(maxSize, max(clus$counts))
            minSize <- c(minSize, min(clus$counts))
            meanS <- c(meanS, mean(clus$counts))
            medianS <- c(medianS, median(clus$counts))
            if (coun == "Denmark"){
                totDK <- c(totDK, clus$counts)
            } else if (coun == "Germany") {
                totDE <- c(totDE, clus$counts)
            } else {
                totCH <- c(totCH, clus$counts)
            }
        }
    }
}

totals <- cbind(country, month, totalN, maxSize, minSize, meanS, medianS)

#to get total number of sequences and total number of clusters

sum(totCH)
sum(totDK)
sum(totDE)
length(totCH)
length(totDK)
length(totDE)

#to get averages per month
totals[which(totals[,1] == "Denmark"),]
