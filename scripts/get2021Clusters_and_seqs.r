setwd("C:/Users/Emma/wsl/corona/sc2_rk/output-5Apr-4muts")
library(lubridate)

#this script uses the TOTAL list of clusters from 2020-2023 (files at top)
#then filters these to only those with a sequence that is within some month in 2021 (first for loop)
#then from only these clusters, extracts all the names of sequences within them
#To then process later to get the dates of these sequences

fls <- c("Switzerland_cluster_distribution_dates_strains_100whole.tsv", 
         "Denmark_cluster_distribution_dates_strains_100whole.tsv",
         "Germany_cluster_distribution_dates_strains_100whole.tsv")

#cntry <- c("Switzerland", "Denmark", "Germany")

cntry <- c("Denmark", "Germany")

for(cn in cntry){

    clusters <- read.csv(paste(cn, "_cluster_distribution_dates_strains_100whole.tsv", sep=""), sep="\t")
    clusters$minDate <- as.Date(clusters$minDate)
    clusters$maxDate <- as.Date(clusters$maxDate)

    newRows <- clusters[which(clusters$minDate<as.Date("2019-01-01")),] #make empty rows

    clusters <- clusters[which(!(is.na(clusters$minDate) | is.na(clusters$maxDate))),]

    for(i in 1:12){
        curM <- as.Date(paste("2021-",i,"-01", sep=""))
        rowToAdd <- clusters[which(!(clusters$maxDate < curM | clusters$minDate >= curM+months(1) )),]
        newRows <- rbind(newRows, rowToAdd)
        #print(nrow(rowToAdd))
        #print(nrow(newRows))
    }

    print(paste(cn,"has", nrow(newRows),"clusters in 2021"))
    write.table(newRows, paste(cn, "_2021_clusters_100whole.tsv", sep=""), sep="\t")

    strns <- newRows$strains
    sep_strns <- c()
    for(s in strns){
        s <- gsub("['", "", s, fixed=T)
        s <- gsub("']", "", s, fixed=T)
        s <- gsub("'", "", s, fixed=T)
        s <- gsub(" ", "", s, fixed=T)
        s_list <-strsplit(s,",")[[1]]
        sep_strns <- append(sep_strns, s_list)
    }
    sep_strns <- unique(sep_strns) #take only those who are unique
    length(sep_strns)
    print(paste(cn,"has", length(sep_strns),"sequences in 2021"))
    write.table(sep_strns, paste(cn, "_2021_clusters_100whole_strainsOnly.tsv", sep=""), quote=F, row.names=F, col.names=F)
    
}


