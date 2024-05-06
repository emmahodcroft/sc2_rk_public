
# R & k project

# overview of characteristics of data:
# number of confirmed cases, number of sequences, number of clusters, 
# number of cases in clusters and size of largest cluster per month
# in Denmark, Germany and Switzerland

# results are presented in tables in supplementary material


library(tidyverse)
library(lubridate)
library(rjson)
library(httr)
library(gt)
library(gtsummary)
library(splitstackshape)

# read cluster data from Denmark
data_dk_0 <- read_delim(file = "output-5Apr-4muts/Denmark_cluster_distribution_dates_100whole.tsv")
data_dk <- data_dk_0 |> filter(!is.na(medianDate))

# read cluster data from Germany
data_0_de <- read_delim(file = "output-5Apr-4muts/Germany_cluster_distribution_dates_100whole.tsv")
data_de <- data_0_de |> filter(!is.na(medianDate))

# read cluster data from Switzerland
data_0_ch <- read_delim(file = "output-5Apr-4muts/Switzerland_cluster_distribution_dates_100whole.tsv")
data_ch <- data_0_ch |> filter(!is.na(medianDate))


# filter data_dk to clusters whose median date is contained in 2021
data_dk_median_2021 <- data_dk |> filter(medianDate >= ymd("2021-01-01") & medianDate < ymd("2022-01-01"))

# filter data_de to clusters whose median date is contained in 2021
data_de_median_2021 <- data_de |> filter(medianDate >= ymd("2021-01-01") & medianDate < ymd("2022-01-01"))

# filter data_ch to clusters whose median date is contained in 2021
data_ch_median_2021 <- data_ch |> filter(medianDate >= ymd("2021-01-01") & medianDate < ymd("2022-01-01"))


# read variants data (obtained from https://covariants.org/per-country?country=Germany&country=Denmark&country=Switzerland)
data_variants <- read_csv(file = "scripts/data_variants.csv")


# read new confirmed cases data from WHO
data_new_confirmed_cases_who_0 <- read_csv("https://covid19.who.int/WHO-COVID-19-global-data.csv")
data_new_confirmed_cases_who <- data_new_confirmed_cases_who_0 |> rename(date = Date_reported, country = Country, new_cases = New_cases) |> filter(country %in% c("Denmark", "Germany", "Switzerland"))

# write data_new_confirmed_cases_who to a csv file
write_csv(x = data_new_confirmed_cases_who, file = "scripts/data_new_confirmed_cases_who.csv")


# # read new confirmed cases data from Germany (https://github.com/robert-koch-institut/COVID-19_7-Tage-Inzidenz_in_Deutschland/blob/main/COVID-19-Faelle_7-Tage-Inzidenz_Deutschland.csv)
# data_new_confirmed_cases_de_0 <- read_csv("https://raw.githubusercontent.com/robert-koch-institut/COVID-19_7-Tage-Inzidenz_in_Deutschland/main/COVID-19-Faelle_7-Tage-Inzidenz_Deutschland.csv")
# 
# data_new_confirmed_cases_de_2021 <- data_new_confirmed_cases_de_0 |> filter(Altersgruppe == "00+" & Meldedatum >= ymd("2021-01-01") & Meldedatum <= ymd("2021-12-31")) |> group_by(Meldedatum) |> summarise(new_cases = sum(Faelle_neu)) |> rename(date = Meldedatum)
# 
# # write data_new_confirmed_cases_de_0 to a csv file
# write_csv(x = data_new_confirmed_cases_de_0, file = "scripts/data_new_confirmed_cases_de_2021.csv")


# # read new confirmed cases data from Switzerland (https://www.covid19.admin.ch/en/overview)
# bag <- GET("https://www.covid19.admin.ch/api/data/context")
# bag <- fromJSON(rawToChar(bag$content))
# data_new_confirmed_cases_ch_0 <- read_csv(bag$sources$individual$csv$daily$cases)
# 
# data_new_confirmed_cases_ch_2021 <- data_new_confirmed_cases_ch_0 |> filter(geoRegion == "CH" & datum >= ymd("2021-01-01") & datum <= ymd("2021-12-31")) |> dplyr::select(c("datum", "entries")) |> rename(date = datum, new_cases = entries)
# 
# # write data_new_confirmed_cases_ch_0 to a csv file
# write_csv(x = data_new_confirmed_cases_ch_0, file = "scripts/data_new_confirmed_cases_ch_2021.csv")


# read number of sequences (obtained from CoVariants) (https://github.com/hodcroftlab/covariants/blob/master/web/data/perCountryDataCaseCounts.json)
data_covariants <- rjson::fromJSON(file = "https://raw.githubusercontent.com/hodcroftlab/covariants/master/web/data/perCountryDataCaseCounts.json")

# write data_covariants to a json file
write(x = toJSON(x = data_covariants, indent = 1), "scripts/data_covariants_country_case_counts.json")

N <- min(length(data_covariants$regions[[1]]$distributions[[4]]$distribution),
         length(data_covariants$regions[[1]]$distributions[[3]]$distribution),
         length(data_covariants$regions[[1]]$distributions[[14]]$distribution))

# filter data_covariants to data from Denmark
print(data_covariants$regions[[1]]$distributions[[4]]$country == "Denmark")
data_dk_week_n_sequences_n_cases <- tibble(week = ymd(unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[4]]$distribution[[x]]$week))),
                                           n_sequences = unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[4]]$distribution[[x]]$total_sequences)),
                                           n_cases = unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[4]]$distribution[[x]]$stand_total_cases)))

# filter data_covariants to data from Germany
print(data_covariants$regions[[1]]$distributions[[3]]$country == "Germany")
data_de_week_n_sequences_n_cases <- tibble(week = ymd(unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[3]]$distribution[[x]]$week))),
                                           n_sequences = unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[3]]$distribution[[x]]$total_sequences)),
                                           n_cases = unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[3]]$distribution[[x]]$stand_total_cases)))

# filter data_covariants to data from Switzerland
print(data_covariants$regions[[1]]$distributions[[14]]$country == "Switzerland")
data_ch_week_n_sequences_n_cases <- tibble(week = ymd(unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[14]]$distribution[[x]]$week))),
                                           n_sequences = unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[14]]$distribution[[x]]$total_sequences)),
                                           n_cases = unlist(lapply(X = 1:N, FUN = function(x) data_covariants$regions[[1]]$distributions[[14]]$distribution[[x]]$stand_total_cases)))

# approximate number of sequences per month in Denmark, Germany and Switzerland
data_dk_de_ch_2021_day_n_sequences <- tibble(date = seq(ymd("2021-01-01"), ymd("2021-12-31"), "days"))

data_dk_de_ch_2021_day_n_sequences <- data_dk_de_ch_2021_day_n_sequences |> mutate(month = month(date),
                                                                                   n_sequences_dk = unlist(lapply(X = 1:nrow(data_dk_de_ch_2021_day_n_sequences),
                                                                                                                  FUN = function(x) round((data_dk_week_n_sequences_n_cases |> filter(week <= data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(n()) |> dplyr::select(n_sequences))[[1]] / 
                                                                                                                                            (interval((data_dk_week_n_sequences_n_cases |> filter(week <= data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(n()) |> dplyr::select(week))[[1]],
                                                                                                                                                      (data_dk_week_n_sequences_n_cases |> filter(week > data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(1) |> dplyr::select(week))[[1]]) / days(1))))))

data_dk_de_ch_2021_day_n_sequences <- data_dk_de_ch_2021_day_n_sequences |> mutate(month = month(date),
                                                                                   n_sequences_de = unlist(lapply(X = 1:nrow(data_dk_de_ch_2021_day_n_sequences),
                                                                                                                  FUN = function(x) round((data_de_week_n_sequences_n_cases |> filter(week <= data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(n()) |> dplyr::select(n_sequences))[[1]] / 
                                                                                                                                            (interval((data_de_week_n_sequences_n_cases |> filter(week <= data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(n()) |> dplyr::select(week))[[1]],
                                                                                                                                                      (data_de_week_n_sequences_n_cases |> filter(week > data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(1) |> dplyr::select(week))[[1]]) / days(1))))))

data_dk_de_ch_2021_day_n_sequences <- data_dk_de_ch_2021_day_n_sequences |> mutate(month = month(date),
                                                                                   n_sequences_ch = unlist(lapply(X = 1:nrow(data_dk_de_ch_2021_day_n_sequences),
                                                                                                                  FUN = function(x) round((data_ch_week_n_sequences_n_cases |> filter(week <= data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(n()) |> dplyr::select(n_sequences))[[1]] / 
                                                                                                                                            (interval((data_ch_week_n_sequences_n_cases |> filter(week <= data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(n()) |> dplyr::select(week))[[1]],
                                                                                                                                                      (data_ch_week_n_sequences_n_cases |> filter(week > data_dk_de_ch_2021_day_n_sequences$date[x]) |> slice(1) |> dplyr::select(week))[[1]]) / days(1))))))

data_dk_de_ch_2021_month_n_sequences <- data_dk_de_ch_2021_day_n_sequences |> group_by(month) |> summarise(n_sequences_dk = sum(n_sequences_dk), n_sequences_de = sum(n_sequences_de), n_sequences_ch = sum(n_sequences_ch))


t0 <- ymd("2021-01-01")

# determine number of clusters, number of cases in clusters and maximal size of clusters of Denmark for each month of 2021
data_cases_sequences_clusters_dk <- data.frame(matrix(data = 0, nrow = 12, ncol = 8))
names(data_cases_sequences_clusters_dk) <- c("month", "n_confirmed_cases", "n_sequences", "n_cases_in_clusters", "n_clusters", "seq_proba_n_sequences", "seq_proba_n_cases_in_clusters", "max_size_clusters")

data_cases_sequences_clusters_dk$n_sequences <- data_dk_de_ch_2021_month_n_sequences$n_sequences_dk

clusters_dk <- data.frame(matrix(nrow = 0, ncol = 3)) 
names(clusters_dk) <- c("size", "month", "frequency")

for (ii in 1:12) {
  
  t1 <- t0 + months(ii - 1)
  t2 <- t1 + months(1)
  
  clusters <- data_dk |> filter(medianDate >= t1 & medianDate < t2) |> dplyr::select(counts)
  
  data_cases_sequences_clusters_dk$month[ii] <- as.character(ii)
  
  data_cases_sequences_clusters_dk$n_confirmed_cases[ii] <- sum(data_new_confirmed_cases_who |> filter(country == "Denmark" & date >= t1 & date < t2) |> dplyr::select("new_cases"))
  
  data_cases_sequences_clusters_dk$n_cases_in_clusters[ii] <- sum(data_dk |> filter(medianDate >= t1 & medianDate < t2) |> dplyr::select(counts))
  
  data_cases_sequences_clusters_dk$n_clusters[ii] <- nrow(clusters)
  
  data_cases_sequences_clusters_dk$seq_proba_n_sequences[ii] <- round(data_cases_sequences_clusters_dk$n_sequences[ii] / data_cases_sequences_clusters_dk$n_confirmed_cases[ii], 4)
  
  data_cases_sequences_clusters_dk$seq_proba_n_cases_in_clusters[ii] <- round(data_cases_sequences_clusters_dk$n_cases_in_clusters[ii] / data_cases_sequences_clusters_dk$n_confirmed_cases[ii], 4)
  
  data_cases_sequences_clusters_dk$max_size_clusters[ii] <- max(clusters)
  
  data_dk_cluster_size_median_2021_month <- as.data.frame(table(data_dk |> filter(medianDate >= t0 + months(ii-1) & medianDate < t0 + months(ii)) |> dplyr::select(counts))) |> mutate(counts = as.numeric(as.character(counts))) |> rename(size = counts, frequency = Freq)
  
  clusters_dk <- bind_rows(clusters_dk, data_dk_cluster_size_median_2021_month |> mutate(month = ii))
  
}

data_cases_sequences_clusters_dk <- bind_rows(data_cases_sequences_clusters_dk, data.frame(month = "Total", 
                                                                                           n_confirmed_cases = sum(data_cases_sequences_clusters_dk$n_confirmed_cases),
                                                                                           n_sequences = sum(data_cases_sequences_clusters_dk$n_sequences),
                                                                                           n_cases_in_clusters = sum(data_cases_sequences_clusters_dk$n_cases_in_clusters),
                                                                                           n_clusters = sum(data_cases_sequences_clusters_dk$n_clusters),
                                                                                           seq_proba_n_sequences = round(sum(data_cases_sequences_clusters_dk$n_sequences) / sum(data_cases_sequences_clusters_dk$n_confirmed_cases), 4),
                                                                                           seq_proba_n_cases_in_clusters = round(sum(data_cases_sequences_clusters_dk$n_cases_in_clusters) / sum(data_cases_sequences_clusters_dk$n_confirmed_cases), 4),
                                                                                           max_size_clusters = max(data_cases_sequences_clusters_dk$max_size_clusters)))


# determine number of clusters, number of cases in clusters and maximal size of clusters of Germany for each month of 2021
data_cases_sequences_clusters_de <- data.frame(matrix(data = 0, nrow = 12, ncol = 8))
names(data_cases_sequences_clusters_de) <- c("month", "n_confirmed_cases", "n_sequences", "n_cases_in_clusters", "n_clusters", "seq_proba_n_sequences", "seq_proba_n_cases_in_clusters", "max_size_clusters")

data_cases_sequences_clusters_de$n_sequences <- data_dk_de_ch_2021_month_n_sequences$n_sequences_de

clusters_de <- data.frame(matrix(nrow = 0, ncol = 3)) 
names(clusters_de) <- c("size", "month", "frequency")

for (ii in 1:12) {
  
  t1 <- t0 + months(ii - 1)
  t2 <- t1 + months(1)
  
  clusters <- data_de |> filter(medianDate >= t1 & medianDate < t2) |> dplyr::select(counts)
  
  data_cases_sequences_clusters_de$month[ii] <- as.character(ii)
  
  data_cases_sequences_clusters_de$n_confirmed_cases[ii] <- sum(data_new_confirmed_cases_who |> filter(country == "Germany" & date >= t1 & date < t2) |> dplyr::select("new_cases"))
  
  data_cases_sequences_clusters_de$n_cases_in_clusters[ii] <- sum(data_de |> filter(medianDate >= t1 & medianDate < t2) |> dplyr::select(counts))
  
  data_cases_sequences_clusters_de$n_clusters[ii] <- nrow(clusters)
  
  data_cases_sequences_clusters_de$seq_proba_n_sequences[ii] <- round(data_cases_sequences_clusters_de$n_sequences[ii] / data_cases_sequences_clusters_de$n_confirmed_cases[ii], 4)
  
  data_cases_sequences_clusters_de$seq_proba_n_cases_in_clusters[ii] <- round(data_cases_sequences_clusters_de$n_cases_in_clusters[ii] / data_cases_sequences_clusters_de$n_confirmed_cases[ii], 4)
  
  data_cases_sequences_clusters_de$max_size_clusters[ii] <- max(clusters)
  
  data_de_cluster_size_median_2021_month <- as.data.frame(table(data_de |> filter(medianDate >= t0 + months(ii-1) & medianDate < t0 + months(ii)) |> dplyr::select(counts))) |> mutate(counts = as.numeric(as.character(counts))) |> rename(size = counts, frequency = Freq)
  
  clusters_de <- bind_rows(clusters_de, data_de_cluster_size_median_2021_month |> mutate(month = ii))
  
}

data_cases_sequences_clusters_de <- bind_rows(data_cases_sequences_clusters_de, data.frame(month = "Total", 
                                                                                           n_confirmed_cases = sum(data_cases_sequences_clusters_de$n_confirmed_cases),
                                                                                           n_sequences = sum(data_cases_sequences_clusters_de$n_sequences),
                                                                                           n_cases_in_clusters = sum(data_cases_sequences_clusters_de$n_cases_in_clusters),
                                                                                           n_clusters = sum(data_cases_sequences_clusters_de$n_clusters),
                                                                                           seq_proba_n_sequences = round(sum(data_cases_sequences_clusters_de$n_sequences) / sum(data_cases_sequences_clusters_de$n_confirmed_cases), 4),
                                                                                           seq_proba_n_cases_in_clusters = round(sum(data_cases_sequences_clusters_de$n_cases_in_clusters) / sum(data_cases_sequences_clusters_de$n_confirmed_cases), 4),
                                                                                           max_size_clusters = max(data_cases_sequences_clusters_de$max_size_clusters)))


# determine number of clusters, number of cases in clusters and maximal size of clusters of Switzerland for each month of 2021
data_cases_sequences_clusters_ch <- data.frame(matrix(data = 0, nrow = 12, ncol = 8))
names(data_cases_sequences_clusters_ch) <- c("month", "n_confirmed_cases", "n_sequences", "n_cases_in_clusters", "n_clusters", "seq_proba_n_sequences", "seq_proba_n_cases_in_clusters", "max_size_clusters")

data_cases_sequences_clusters_ch$n_sequences <- data_dk_de_ch_2021_month_n_sequences$n_sequences_ch

clusters_ch <- data.frame(matrix(nrow = 0, ncol = 3)) 
names(clusters_ch) <- c("size", "month", "frequency")

for (ii in 1:12) {
  
  t1 <- t0 + months(ii - 1)
  t2 <- t1 + months(1)
  
  clusters <- data_ch |> filter(medianDate >= t1 & medianDate < t2) |> dplyr::select(counts)
  
  data_cases_sequences_clusters_ch$month[ii] <- as.character(ii)
  
  data_cases_sequences_clusters_ch$n_confirmed_cases[ii] <- sum(data_new_confirmed_cases_who |> filter(country == "Switzerland" & date >= t1 & date < t2) |> dplyr::select("new_cases"))
  
  data_cases_sequences_clusters_ch$n_cases_in_clusters[ii] <- sum(data_ch |> filter(medianDate >= t1 & medianDate < t2) |> dplyr::select(counts))
  
  data_cases_sequences_clusters_ch$n_clusters[ii] <- nrow(clusters)
  
  data_cases_sequences_clusters_ch$seq_proba_n_sequences[ii] <- round(data_cases_sequences_clusters_ch$n_sequences[ii] / data_cases_sequences_clusters_ch$n_confirmed_cases[ii], 4)
  
  data_cases_sequences_clusters_ch$seq_proba_n_cases_in_clusters[ii] <- round(data_cases_sequences_clusters_ch$n_cases_in_clusters[ii] / data_cases_sequences_clusters_ch$n_confirmed_cases[ii], 4)
  
  data_cases_sequences_clusters_ch$max_size_clusters[ii] <- max(clusters)
  
  data_ch_cluster_size_median_2021_month <- as.data.frame(table(data_ch |> filter(medianDate >= t0 + months(ii-1) & medianDate < t0 + months(ii)) |> dplyr::select(counts))) |> mutate(counts = as.numeric(as.character(counts))) |> rename(size = counts, frequency = Freq)
  
  clusters_ch <- bind_rows(clusters_ch, data_ch_cluster_size_median_2021_month |> mutate(month = ii))
  
}

data_cases_sequences_clusters_ch <- bind_rows(data_cases_sequences_clusters_ch, data.frame(month = "Total", 
                                                                                           n_confirmed_cases = sum(data_cases_sequences_clusters_ch$n_confirmed_cases),
                                                                                           n_sequences = sum(data_cases_sequences_clusters_ch$n_sequences),
                                                                                           n_cases_in_clusters = sum(data_cases_sequences_clusters_ch$n_cases_in_clusters),
                                                                                           n_clusters = sum(data_cases_sequences_clusters_ch$n_clusters),
                                                                                           seq_proba_n_sequences = round(sum(data_cases_sequences_clusters_ch$n_sequences) / sum(data_cases_sequences_clusters_ch$n_confirmed_cases), 4),
                                                                                           seq_proba_n_cases_in_clusters = round(sum(data_cases_sequences_clusters_ch$n_cases_in_clusters) / sum(data_cases_sequences_clusters_ch$n_confirmed_cases), 4),
                                                                                           max_size_clusters = max(data_cases_sequences_clusters_ch$max_size_clusters)))



# create tables of results

# list of months
months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

# Denmark
names(data_cases_sequences_clusters_dk) <- c("Month", "Nocc", "Nos", "Nocic", "Noc", "Nos/Nocc", "Nocic/Nocc", "Solc")
table_data_cases_sequences_clusters_dk <- data_cases_sequences_clusters_dk |> gt(rowname_col = "month") 
table_data_cases_sequences_clusters_dk |> gt::gtsave("output/plots/denmark/table_data_cases_sequences_clusters_dk.png")

table_data_dk <- clusters_dk |> dplyr::select("size", "month", "frequency") |> expandRows("frequency") |>
  mutate(size = factor(cut(x = size, breaks = c(0,1,2,3,4,5,10,50,100,max(size)), labels = c("1", "2", "3", "4", "5", "6-10", "11-50", "51-100", "101-1278"))),
         month = factor(cut(x = month, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12), labels = months))) |>
  tbl_summary(by = month) 

table_data_dk |> as_gt() |> gt::gtsave("output/plots/denmark/table_data_dk.png")

# Germany
names(data_cases_sequences_clusters_de) <- c("Month", "Nocc", "Nos", "Nocic", "Noc", "Nos/Nocc", "Nocic/Nocc", "Solc")
table_data_cases_sequences_clusters_de <- data_cases_sequences_clusters_de |> gt(rowname_col = "month") 
table_data_cases_sequences_clusters_de |> gt::gtsave("output/plots/germany/table_data_cases_sequences_clusters_de.png")

table_data_de <- clusters_de |> dplyr::select("size", "month", "frequency") |> expandRows("frequency") |>
  mutate(size = factor(cut(x = size, breaks = c(0,1,2,3,4,5,10,50,100,max(size)), labels = c("1", "2", "3", "4", "5", "6-10", "11-50", "51-100", "101-283"))),
         month = factor(cut(x = month, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12), labels = months))) |>
  tbl_summary(by = month) 

table_data_de |> as_gt() |> gt::gtsave("output/plots/germany/table_data_de.png")

# Switzerland
names(data_cases_sequences_clusters_ch) <- c("Month", "Nocc", "Nos", "Nocic", "Noc", "Nos/Nocc", "Nocic/Nocc", "Solc")
table_data_cases_sequences_clusters_ch <- data_cases_sequences_clusters_ch |> gt(rowname_col = "month") 
table_data_cases_sequences_clusters_ch |> gt::gtsave("output/plots/switzerland/table_data_cases_sequences_clusters_ch.png")

table_data_ch <- clusters_ch |> dplyr::select("size", "month", "frequency") |> expandRows("frequency") |>
  mutate(size = factor(cut(x = size, breaks = c(0,1,2,3,4,5,10,50,100,max(size)), labels = c("1", "2", "3", "4", "5", "6-10", "11-50", "51-100", "101-363"))),
         month = factor(cut(x = month, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12), labels = months))) |>
  tbl_summary(by = month) 

table_data_ch |> as_gt() |> gt::gtsave("output/plots/switzerland/table_data_ch.png")
