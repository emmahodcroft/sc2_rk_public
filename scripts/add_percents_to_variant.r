### To add percentage downsampling to cluster_inference_variants_grid.r
#This assumes you already have a read-in or made `cntrys` that has 100% for variants, and want to add more

new_percs <- c("50", "25", "5")



for (c in 1:length(cntrys)){

    for(p in new_percs){
        tmp <- vector(mode="list", length=length(all_var))
        names(tmp) <- all_var
        cntrys[[c]][[p]] <- tmp
    }
}


### fill in model


start_time <- Sys.time()

for (p in new_percs){
for (c in countries){
    file <- paste("output/",c,"_cluster_distribution_dates_",p,"whole.tsv",sep="")
    counts <- read.delim(file)
    counts <- subset(counts, (minDate >= our_min_date & maxDate <= our_max_date))
    #Delta is converted to being just called 'Delta' in the script that generates the files
    counts$clade[which(counts$clade == "20I (Alpha, V1)")] <- "Alpha"

for (v in all_var){
    print(paste("Now doing",c,",",v,",",p))

    if (v == "all"){
        counts_v <- counts
    } else {
        counts_v <- subset(counts, clade == v)
    }

    cntrys[[c]][[p]][[v]] <- vector(mode="list", length=2)
    names(cntrys[[c]][[p]][[v]]) <- c("counts", "fit")
    cntrys[[c]][[p]][[v]][["counts"]] <- counts_v

    clusters <- counts_v$counts

    liksurf=gridEstimate(clusters, R0range, krange)
    likMax=(liksurf==min(liksurf))

    output1 <- calc_profile(liksurf,likMax)
    print(output1)
    cntrys[[c]][[p]][[v]][["fit"]] <- output1

}
}
}
end_time <- Sys.time()
print(paste("total time", (end_time - start_time)))



