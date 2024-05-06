
#Updated function to get number of total sequences in clusters (and total number of clusters)
#with median date in 2021
#And to then get a list of all the strains that are part of this dataset (step 1 for generating epi_set)

setwd("C:/Users/Emma/wsl/corona/sc2_rk/output-5Apr-4muts")


countries <- c("Denmark", "Switzerland", "Germany")

totDK <- c()
totCH <- c()
totDE <- c()

country <- c()
#month <- c()
totalN <- c()
maxSize <- c()
minSize <- c()
meanS <- c()
medianS <- c()

#for finding just clusters with Median date in 2021

startDate <- as.Date("2021-01-01")
endDate <- as.Date("2022-01-01")

months <- seq(from=startDate, to=endDate, by="month")

for (coun in countries) {
    fileN <- paste(coun, "_cluster_distribution_dates_strains_100whole.tsv", sep="")
    clus <- read.csv(fileN, sep="\t", as.is=T)
    clus$medianDate <- as.Date(clus$medianDate)
    clus$meanDate <- as.Date(clus$meanDate)

    #for(i in 1:length(months)-1) {
    clus <- clus[which(clus$medianDate >= startDate & clus$medianDate < endDate),]
    #get list of strains
    strainsTog <- c(clus["strains"])$strains
    strainsTog <- gsub("'","",strainsTog,fixed = TRUE)
    strainsTog <- gsub("[","",strainsTog,fixed = TRUE)
    strainsTog <- gsub("]","",strainsTog,fixed = TRUE)
    allStrains <- paste(strainsTog, collapse=", ")
    allStrains <- gsub(", ", "\n", allStrains)
    write.table(allStrains, paste(coun,"_strains.csv", sep=""), row.names=F, col.names=F, quote=F)

    country <- c(country, coun)
    #month <- c(month, fol)
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
    #}
}

#To get total number of sequences and total number of clusters for 2021
sum(totCH)
sum(totDK)
sum(totDE)
length(totCH)
length(totDK)
length(totDE)


