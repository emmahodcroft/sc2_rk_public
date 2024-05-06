# Estimate R and k from cluster size distributions

setwd("C:/Users/Emma/wsl/corona/sc2_rk")

# Load libraries
library(bbmle)
library(lubridate)
#library(comprehenr) #don't think I use this anymore?

#this will look for input files in a folder in the same directory as you're running this code, called "output"


# Probability to observe a transmission cluster of size L
# This is the core function that gives the probability to 
# observe a certain cluster size L as a function of R and k.
prob <- function(L, R, k) {
  p <- exp(lgamma(k*L + L - 1) - lgamma(k*L) - lgamma(L + 1) + log((R/k)^(L - 1)/(1 + R/k)^(k*L + L - 1)))
  p <- ifelse(p == 0, .Machine$double.eps, p)
  p <- ifelse(is.na(p), .Machine$double.eps, p)
  return(p)
}

#new probability equation from Christian from 6 Jul
# Probability to observe a cluster of size L
prob <- function(L, R, k) {
    p <- exp(lgamma(k*L + L - 1) - lgamma(k*L) - lgamma(L + 1) + (L - 1)*log(R/k) - (k*L + L - 1)*log((1 + R/k)))
    return(p)
}

# Probability to observe a transmission cluster of size L -- to try to make it smoother for CI estimation
# from Christian A - on INPUT slack: https://inputbern.slack.com/archives/D01E3V4STBL/p1647526956138399?thread_ts=1647453562.511369&cid=D01E3V4STBL
#prob <- function(L, R, k) {
#  L_range <- 1:1e3
#  p_range <- exp(lgamma(k*L_range + L_range - 1) - lgamma(k*L_range) - lgamma(L_range + 1) + log((R/k)^(L_range - 1)/(1 + R/k)^(k*L_range + L_range - 1)))
#  m <- max(which(p_range > 0))
#  #print(paste("m", m))
#  #print(paste("len p_range",length(p_range)))
#  p_cum <- 1 - sum(p_range[1:m])
#  p_range <- c(p_range[1:m], p_cum)
#  L[L > m] <- m + 1
#  p <- p_range[L]
#  return(p)
#}

# Negative log-likelihood
# This is the negative log-likelihood that needs to be optimized.
nll <- function(clusters, R, k) {
  #R <- exp(R)
  #k <- exp(k)
  #print(paste("R", R))
  #print(paste("k", k))
  
  #ll <- log(prob(clusters, R^2, k^2))
  #return(- sum(ll) + (R^2-0.5)^2 + (k^2-0.25)^2)

  ll <- log(prob(clusters, R, k))
  return(- sum(ll))
}

gridEstimate <- function(data, R0range, krange){
  collect1 <- data.frame(matrix(NA, nrow=length(R0range),length(krange)))
  for(ii in 1:length(R0range)){
    for(jj in 1:length(krange)){
      rr=R0range[ii]
      kkr=krange[jj]
      collect1[ii,jj]=nll(data,rr,kkr)
    }
  }
  return(collect1)
}

calc_profile <- function(liksurfMX, likmMX){

    chiV=qchisq(95/100,df=1)/2

    prfk=apply(liksurfMX,2,function(x){min(x)})
    prfk2=krange[prfk-min(prfk)<chiV]
    prfR=apply(liksurfMX,1,function(x){min(x)})
    prfR2=R0range[prfR-min(prfR)<chiV]

    rEst <- R0range[sum(seq(1,length(R0range))%*%likmMX)]
    rCI <- c(min(prfR2), max(prfR2))
    kEst <- krange[sum(likmMX%*%seq(1,length(krange)))]
    kCI <- c(min(prfk2), ifelse(max(prfk2)==max(krange),"Inf",max(prfk2)))

    result <- as.data.frame(matrix(c(rEst, kEst, rCI[1], kCI[1], rCI[2], kCI[2]), nrow=2, ncol=3))
    colnames(result) <- c("Est", "CI.1", "CI.2")
    rownames(result) <- c("R", "k")

    return(result)
  
    #c(
    #    paste(R0range[sum(seq(1,length(R0range))%*%likmMX)]^2," (",min(prfR2)^2,"-",max(prfR2)^2,")",sep=""),
    #    paste(krange[sum(likmMX%*%seq(1,length(krange)))]^2," (",min(prfk2)^2,"-",ifelse(max(prfk2)==max(krange),"Inf",max(prfk2)^2),")",sep="")
    #)
  
}

#############################
#############################
# This is to explore all data and %sampling - does *NOT* look at variants - see cluster_inference_variants.R for code 
# that will do that

#set dates. a cluster's mindate must be after mindate and maxdate must be before maxdate for cluster to be included
#ex: counts <- subset(counts, (minDate >= ymd(20210101) & maxDate <= ymd(20211130)))
our_min_date <- ymd(20210101)
our_max_date <- ymd(20211130)
###    here set cutoff for start of Alpha (2021-01-01) to before Omicron (2021-11-30), 
###    so that can say that Denmark sequencing basically 100%

countries <- c("Denmark", "Germany", "Switzerland")
percents <- c("5", "10", "25", "50", "75", "100")
output <- "output/"
dir.create(paste(output,"plots",sep=""))

txt1 <- paste("This run was done at", Sys.time(), "\n\nAnd involves these countries: ", paste(countries, collapse=", "), "\n\nSampled at these percents:", paste(percents, collapse="%, "),
    "\n\nwith data from", output)
write(txt1, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""))

#fake variant params so nothing breaks
variants <- c("all")
all_var <- c("all")

### Set up simple storage (only stores final estimated k & R values, cannot access anything else later)

df_k <- as.data.frame(matrix(ncol=length(countries),nrow=length(percents)), row.names=percents)
colnames(df_k) <- countries

df_r <- as.data.frame(matrix(ncol=length(countries),nrow=length(percents)), row.names=percents)
colnames(df_r) <- countries

## set up lists to store more things - these store the `counts` and `fit`, so you can access them later if you want
#countries
cntrys <- vector(mode="list", length=length(countries))
#percents
names(cntrys) <- countries
for (c in 1:length(cntrys)){
    cntrys[[c]] <- vector(mode="list", length=length(percents))
    names(cntrys[[c]]) <- percents
    #for each percent, add a list, if adding percent in future
    #variants
    for (v in 1:length(cntrys[[c]])){
        cntrys[[c]][[v]] <- vector(mode="list", length=length(all_var))
        names(cntrys[[c]][[v]]) <- all_var
    }
}
# This list ends up having the structure:
#  Swiss
#    100
#        all
#            counts
#            fit
#     75
#        all
#            counts
#            fit
#        ....
#  Germany
#    100
#      ....
#Currently there are only 'all' varint entries
#
# You access values like so:
# cntrys[["Denmark"]][["100"]][["all"]][["counts"]]



###### All data - 
###    here set cutoff for start of Alpha (2021-01-01) to before Omicron (2021-11-30), 
###    so that can say that Denmark sequencing basically 100%

#set search boxes?
R01=0.05; R02=0.9; R0range=seq(R01,R02,by=0.01)
k1=0.05; k2=0.6; krange=seq(k1,k2,by=0.01)

txt2 <- paste("\n\n\nSearch boxes for R were from",R01,"to",R02,"by values of 0.01.\n\nSearch boxes for k were from",k1,"to",k2,"by values of 0.01.")
write(txt2, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)

txt2a <- paste("\n\n\nThe data used here is only from clusters that fall between",our_min_date,"and",our_max_date,"(all sequences in a cluster must be within these dates for the cluster to be included)")
write(txt2a, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)

#See if want to use old models or new ones
#https://stackoverflow.com/questions/31893844/how-can-i-save-a-list-to-a-file-and-read-it-in-again-in-r
{
    if(file.exists(paste(output,"variants_model_fit_output.RData",sep=""))){
        var = readline(prompt=paste("Do you want to read in previous data?\n This will be read in from",paste(output,"model_fit_output.RData",sep="")," : "));

        if(var %in% c("Y", "y", "yes", "Yes", "YES")){
            print("Reading in old data...")
            cntrys <- readRDS(paste(output,"model_fit_output.RData",sep=""))
            readin=TRUE
        } else {
            print("Will start running fresh model now.")
            readin=FALSE
        }
    } else {
        print("Will start running fresh model now.")
        readin=FALSE
    }
}


#If want to run fresh, do this!
## Be careful, with the new code this takes a while to run
#But it also gets CIs!
# Takes about 30 minutes

if(readin==FALSE){
    start_time <- Sys.time()

    v = "all" #set variants manually for now
    for (c in  countries){
    for (p in percents){
        print(paste("Now doing",c,",",p,"%"))
        file <- paste("output/",c,"_cluster_distribution_dates_",p,"whole.tsv",sep="")
        counts <- read.delim(file)
        counts <- subset(counts, (minDate >= our_min_date & maxDate <= our_max_date))

        cntrys[[c]][[p]][[v]] <- vector(mode="list", length=2)
        names(cntrys[[c]][[p]][[v]]) <- c("counts", "fit")
        cntrys[[c]][[p]][[v]][["counts"]] <- counts

        clusters <- counts$counts

        liksurf=gridEstimate(clusters, R0range, krange)
        likMax=(liksurf==min(liksurf))

        output1 <- calc_profile(liksurf,likMax)
        print(output1)
        cntrys[[c]][[p]][[v]][["fit"]] <- output1
        df_k[p,c] <- output1["k", "Est"]
        df_r[p,c] <- output1["R", "Est"]

        # this is the old fitting procedure.
        #free <- c(R = 0.5^0.5, k = 0.2^0.5)
        #fit <- (mle2(nll, start = as.list(free), method = "Nelder-Mead", control = list(maxit = 1e5)))
        #cntrys[[c]][[p]][[v]][["fit"]] <- fit
        #co <- coef(fit)
        #df_k[p,c] <- co["k"]^2
        #df_r[p,c] <- co["R"]^2
    }
    }
    end_time <- Sys.time()
    print(paste("total time", (end_time - start_time)))
}


#If ran new model
#Write this object out if you want (braces allow you to enter input without skipping to next line)
writeout<-FALSE
if(readin==FALSE){
    var = readline(prompt="Do you want to save the model fit? : ");

    if(var %in% c("Y", "y", "yes", "Yes", "YES")){
        saveRDS(cntrys, file=paste(output,"model_fit_output.RData",sep=""))
        writeout<-TRUE
    }
}


if(readin){
    txt3 <- paste("\nThese numbers/plots were generated from old data from",paste(output,"model_fit_output.RData",sep=""))
} else {
    txt3 <- paste("\nThese numbers/plots were generated from fresh data")
    if(writeout){
        txt3 <- paste(txt3, "and this was written out to",file=paste(output,"model_fit_output.RData",sep=""))
    }
}
write(txt3, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)


###########################
###########################
###########################
## simple plotting - no possibilty to get CIs as they're lost

txt4 <- paste("\nPlot of the overall estimates without sample level correction or CI: \n![plot1](",paste("../",output,"plots/overall_noCI_nosamp.png",sep=""),")",sep="")
write(txt4, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)

png(file=paste(output,"plots/overall_noCI_nosamp.png",sep=""), width=1000, height=500)
par(mfrow=c(1,2))
#denmark, germany switzerland
cols = c("coral", "aquamarine3", "purple")

#plot k values
plot(df_k[,1], ylab="k estimate", xlab="Percent", xlim=c(0.5,length(percents)+0.5), ylim=c(min(df_k), max(df_k)), xaxt="n", col="white", main="k without adjusting for true sampling")
axis(1, at=seq(1,length(percents)), labels=percents)

for (i in 1:length(countries)){
  points(df_k[,i], col=cols[i], cex=1.5, pch=16)
}
legend("topleft", countries, col=cols, pch=16)

#plot R values
plot(df_r[,1], ylab="R estimate", xlab="Percent", xlim=c(0.5,length(percents)+0.5), ylim=c(min(df_r), max(df_r)), xaxt="n", col="white", main="R without adjusting for true sampling")
axis(1, at=seq(1,length(percents)), labels=percents)

for (i in 1:length(countries)){
  points(df_r[,i], col=cols[i], cex=1.5, pch=16)
  print(paste(i, df_r[,i]))
}
legend("topleft", countries, col=cols, pch=16)
dev.off()


###########################
###########################
###########################
## more-complex plotting (using complex storage) - with CI
# not super worth doing - recommend to skip

doCI <- TRUE

par(mfrow=c(1,2))
cols = c("red", "blue", "darkgreen")

#### K
plot(df_k[,1], ylab="k estimate", xlab="Percent", xlim=c(0.5,length(percents)+0.5), ylim=c(min(df_k), max(df_k)), xaxt="n", col="white")
axis(1, at=seq(1,length(percents)), labels=percents)

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    for (p in percents) {
        print(paste(c,p))
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        xval <- which(percents==p)
        points( xval , fit["k","Est"], col=cols[i])
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        #if (doCI) { try( { ci <- confint(fit); segments(xval, ci[2], xval, ci[4], col=cols[i]) } ,silent=T) } #these are for k 
        if (doCI) {
            segments(xval, fit["k","CI.1"], xval, fit["k","CI.2"], col=cols[i])
        }
    }
}
legend("topleft", countries, col=cols, pch=1)

### R

plot(df_r[,1], ylab="R estimate", xlab="Percent", xlim=c(0.5,length(percents)+0.5), ylim=c(min(df_r), max(df_r)), xaxt="n", col="white")
axis(1, at=seq(1,length(percents)), labels=percents)

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    for (p in percents) {
        print(paste(c,p))
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        xval <- which(percents==p)
        points( xval , fit["R","Est"], col=cols[i])
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        #if (doCI) { try( { ci <- confint(fit); segments(xval, ci[2], xval, ci[4], col=cols[i]) } ,silent=T) } #these are for k 
        if (doCI) {
            segments(xval, fit["R","CI.1"], xval, fit["R","CI.2"], col=cols[i])
        }
    }
}
legend("topleft", countries, col=cols, pch=1)

####################
####################
# plot just 100% of each country to show the 'top' estimates


txt4 <- paste("\nPlot of the un-downsampled estimates (100%) with CI: \n![plot1](",paste("../",output,"plots/just100p_kandR.png",sep=""),")",sep="")
write(txt4, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)

doCI <- TRUE

png(file=paste(output,"plots/just100p_kandR.png",sep=""), width=900, height=450)
par(mfrow=c(1,2))
#denmark, germany switzerland
cols = c("coral", "aquamarine3", "purple")

#### K
#ylim=c(min(df_k["100",]), max(df_k["100",]))
plot(df_k[,1], ylab="k estimate", xlab="Country", xlim=c(0.5,length(countries)+0.5), 
    #ylim=c(0.22,0.32), 
    ylim=c(0.10,0.50),
    xaxt="n", col="white", main="k estimate using all samples")
#abline(h=seq(0.22,0.32,by=0.02), lty="dashed", col="grey")
abline(h=seq(0.10,0.50,by=0.05), lty="dashed", col="grey")
axis(1, at=seq(1,length(countries)), labels=countries)

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    for (p in c("100")) {
        print(paste(c,p))
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        xval <- which(countries==c)
        points( xval , fit["k","Est"], col=cols[i], cex=1.5, pch=16)
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        #if (doCI) { try( { ci <- confint(fit); segments(xval, ci[2], xval, ci[4], col=cols[i]) } ,silent=T) } #these are for k 
        if (doCI) {
            segments(xval, fit["k","CI.1"], xval, fit["k","CI.2"], col=cols[i])
        }
    }
}
legend("topleft", countries, col=cols, pch=16, bg="white")

### R
#ylim=c(min(df_r), max(df_r))
plot(df_r[,1], ylab="R estimate", xlab="Country", xlim=c(0.5,length(countries)+0.5), ylim=c(0.3,0.7), xaxt="n", col="white", main="R estimate using all samples")
abline(h=seq(0.3,0.7,by=0.05), lty="dashed", col="grey")
axis(1, at=seq(1,length(countries)), labels=countries)

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    for (p in c("100")) {
        print(paste(c,p))
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        xval <- which(countries==c)
        points( xval , fit["R","Est"], col=cols[i], cex=1.5, pch=16)
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        #if (doCI) { try( { ci <- confint(fit); segments(xval, ci[2], xval, ci[4], col=cols[i]) } ,silent=T) } #these are for k 
        if (doCI) {
            segments(xval, fit["R","CI.1"], xval, fit["R","CI.2"], col=cols[i])
        }
    }
}
legend("topleft", countries, col=cols, pch=16, bg="white")
dev.off()


#####################################
####################################
# Try again but reflect the actual sequencing percent
# AND un LOG

# Switzerland  18%
# Denmark      73%
# Germany      11%


###### plot unlogged

txt4 <- paste("\nPlot of 'true' sampling, unlogged, estimates (100% in black): \n![plot1](",paste("../",output,"plots/unlogged_allSamp_kandR.png",sep=""),")",sep="")
write(txt4, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)

doCI <- FALSE

png(file=paste(output,"plots/unlogged_allSamp_kandR.png",sep=""), width=35, height=17.5, units="cm", res=300)#width=900, height=450)

real_percent <- as.matrix(c(0.73, 0.11, 0.18), 1,3)
rownames(real_percent) <- countries

par(mfrow=c(1,2))
cols = c("coral", "aquamarine3", "purple")

plot(df_k[,1], ylab="k estimate", xlab="Fraction of samples sequenced\n(with downsampling)", xlim=c(0,1), 
        ylim=c(min(df_k), max(df_k)+0.03), 
        #xaxt="n", 
        #log="x",
        col="white")
#axis(1, at=seq(1,length(percents)), labels=percents)
abline(h=seq(0.1,0.35,by=0.05), lty="dashed", col="grey")


v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    k_s <- c()
    p_s <- c()
    for (p in percents) {
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        num_p <- as.numeric(p)/100
        fin_p <- num_p*real_percent[c,]

        points( fin_p, fit["k","Est"], col=cols[i], cex=1.2, pch=16)
        if (p=="100"){
            points( fin_p, fit["k","Est"], col="black", cex=1.5, lwd=4, pch=1)
        }
        k_s <- c(k_s, fit["k","Est"])
        p_s <- c(p_s, log(fin_p))
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        #if (doCI) { try( { ci <- confint(fit); segments(xval, ci[2], xval, ci[4], col=cols[i]) } ,silent=T) } #these are for k 
        if (doCI) {
            segments(fin_p, fit["k","CI.1"], fin_p, fit["k","CI.2"], col=cols[i])
        }
    }
    #abline(lm(k_s ~ p_s), col=cols[i])
    #est <- round(coef(lm(k_s ~ p_s))["(Intercept)"],3)
    #text(-3, 0.25+i*0.01, paste("Intercept k =",est), col=cols[i])
}
legend("topleft", countries, col=cols, pch=16, bg="white")


### R

plot(df_k[,1], ylab="R estimate", xlab="Fraction of samples sequenced\n(with downsampling)", xlim=c(0,1), 
        ylim=c(min(df_r), max(df_r)), 
        #xaxt="n", 
        #log="x",
        col="white")
#axis(1, at=seq(1,length(percents)), labels=percents)
abline(h=seq(0.1,0.6,by=0.1), lty="dashed", col="grey")

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    r_s <- c()
    p_s <- c()
    for (p in percents) {
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        num_p <- as.numeric(p)/100
        fin_p <- num_p*real_percent[c,]

        points( fin_p, fit["R","Est"], col=cols[i], cex=1.2, pch=16)
        if (p=="100"){
            points( fin_p, fit["R","Est"], col="black", cex=1.5, lwd=4, pch=1)
        }

        r_s <- c(r_s, fit["R","Est"])
        p_s <- c(p_s, log(fin_p))
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        #if (doCI) { try( { ci <- confint(fit); segments(xval, ci[2], xval, ci[4], col=cols[i]) } ,silent=T) } #these are for k
        if (doCI) {
            segments(fin_p, fit["R","CI.1"], fin_p, fit["R","CI.2"], col=cols[i])
        } 
    }
    #abline(lm(r_s ~ p_s), col=cols[i])
    #est <- round(coef(lm(r_s ~ p_s))["(Intercept)"],3)
    #text(-3, 0.5+i*0.02, paste("Intercept R =",est), col=cols[i])
}
legend("topleft", countries, col=cols, pch=16, bg="white")
dev.off()


##################################
#################################
# plot with logging

### K


txt4 <- paste("\nPlot of 'true' sampling, on log scale, with extrapolation back to 'real' 100% sampling of known cases: \n![plot1](",paste("../",output,"plots/samp_extrapolate_kandR.png",sep=""),")",sep="")
write(txt4, file=paste("reports/",gsub("/", "", output),"_report.md",sep=""), append=T)

doCI <- TRUE

png(file=paste(output,"plots/samp_extrapolate_kandR.png",sep=""), width=900, height=450)
real_percent <- as.matrix(c(0.73, 0.11, 0.18), 1,3)
rownames(real_percent) <- countries

par(mfrow=c(1,2))
#orig basic colors
#cols = c("red", "blue", "darkgreen")

#denmark, germany switzerland
cols = c("coral", "aquamarine3", "purple")

plot(log(df_k[,1]), ylab="k estimate", xlab="Log of fraction of samples sequenced", xlim=c(-5.3,0),#xlim=c(0.01,1), 
        #ylim=c(min(df_k), max(df_k)), 
        ylim=c(0.09, 0.37),
        #xaxt="n", 
        #log="x",
        col="white")
#axis(1, at=seq(1,length(percents)), labels=percents)
abline(h=seq(0.1,0.35,by=0.05), lty="dashed", col="grey")
abline(v=0, col="black")

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    k_s <- c()
    p_s <- c()
    for (p in percents) {
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        print(paste("p",p))
        num_p <- as.numeric(p)/100
        fin_p <- num_p*real_percent[c,]
        xval <- log(fin_p)

        points( log(fin_p), fit["k","Est"], col=cols[i], pch=16, cex=1.2)
        k_s <- c(k_s, fit["k","Est"])
        p_s <- c(p_s, log(fin_p))
        if (doCI) {  
            #ci1 <- confint(fit); print(ci1^2); segments(log(fin_p), ci1[2]^2, log(fin_p), ci1[4]^2, col=cols[i]) 
            segments(xval, fit["k","CI.1"], xval, fit["k","CI.2"], col=cols[i])
        } 
    }
    mod <- lm(k_s ~ p_s)
    abline(mod, col=cols[i])
    est <- round(coef(mod)["(Intercept)"],3)
    se <- round( coef(summary(mod))[, "Std. Error"]["(Intercept)"] ,3)
    ci <- c(est-se, est+se)

    #text with CI
    #text(-3, 0.24+i*0.01, paste("Intercept k =",est, "CI=",ci[1],",",ci[2]), col=cols[i])
    #text without CI
    text(-2.3, 0.38-i*0.02, paste("Intercept k =",est), col=cols[i], cex=1.2)
    points(0, est, pch=4, lwd=4, cex=1.3, col=cols[i])
}
legend("topleft", countries, col=cols, pch=16, bg="white")


### R

plot(log(df_k[,1]), ylab="R estimate", xlab="Log of fraction of samples sequenced", xlim=c(-5.3,0),#xlim=c(0.01,1), 
        #ylim=c(min(df_r), max(df_r)), 
        ylim=c(0.08, 0.7),
        #xaxt="n", 
        #log="x",
        col="white")
#axis(1, at=seq(1,length(percents)), labels=percents)
abline(h=seq(0.1,0.7,by=0.1), lty="dashed", col="grey")
abline(v=0, col="black")

v = "all"
for (i in 1:length(countries)){
    c <- countries[i]
    r_s <- c()
    p_s <- c()
    for (p in percents) {
      print(paste(c, p))
        fit <- cntrys[[c]][[p]][[v]][["fit"]]
        num_p <- as.numeric(p)/100
        fin_p <- num_p*real_percent[c,]
        xval <- log(fin_p)

        points( log(fin_p), fit["R","Est"], col=cols[i], pch=16, cex=1.2)

        r_s <- c(r_s, fit["R","Est"])
        p_s <- c(p_s, log(fin_p))
        # Have to use try here as some of the confint calls error, which stops it from trying others!
        if (doCI) {  
            #ci1 <- confint(fit); print(ci1^2); segments(log(fin_p), ci1[2]^2, log(fin_p), ci1[4]^2, col=cols[i]) 
            segments(xval, fit["R","CI.1"], xval, fit["R","CI.2"], col=cols[i])
        }
    }
    mod <- lm(r_s ~ p_s)
    abline(mod, col=cols[i])
    est <- round(coef(mod)["(Intercept)"],2)
    se <- round( coef(summary(mod))[, "Std. Error"]["(Intercept)"] ,3)
    ci <- c(est-se, est+se)
    #text with CI
    #text(-3, 0.5+i*0.02, paste("Intercept R =",est, "CI=",ci[1],",",ci[2]), col=cols[i])
    #text without CI
    text(-2.5, 0.7-i*0.03, paste("Intercept R =",est), col=cols[i], cex=1.2)
    points(0, est, pch=4, lwd=4, cex=1.3, col=cols[i])
}
legend("topleft", countries, col=cols, pch=16, bg="white")
dev.off()

#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
### Extra code




#####################################
#####################################
####

#this allows to look at distribution of clusters - so y% of clusters are larger than size (on x axis)

# This is a BIG plot with a lot of points so it takes time to plot & R may freeze a bit!!!!!
# One before trying to reduce singletons is in identicals/cluster_distribution_variant_country_9Mar22.PNG

par(mfrow=c(length(countries),length(percents)))
v = "all"
for (c in countries){
    for (p in percents) {
        counts <- cntrys[[c]][[p]][[v]][["counts"]]
        clusters <- counts$counts
        yax <- seq(from=1,to=1/length(clusters), by=-1/(length(clusters)))
        plot(clusters, yax, ylim=c(1E-5,1), xlim=c(1,1000), log="xy", main=paste(c,p))
        segments(1,1,500,1/500)
        segments(1,1,500,1/(500^2))
    }
}
