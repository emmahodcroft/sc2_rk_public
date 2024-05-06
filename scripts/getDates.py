#Fairly inefficient script to try and get the dates for sequences
#from the list of sequences which are in clusters in 2021

#Have list of strain names from 'get2021Clusters_and_seqs.r', 
#but need to get dates - so do this by matching up with the original metadata file

import pandas as pd

input_meta = "../archived-data/metadata-07Mar23.tsv"

meta = pd.read_csv(input_meta, sep="\t", usecols=['strain', 'date'], index_col=False)

counFiles = ["Switzerland_2021_clusters_100whole_strainsOnly.tsv", 
            "Denmark_2021_clusters_100whole_strainsOnly.tsv", 
            "Germany_2021_clusters_100whole_strainsOnly.tsv"]

wanted = pd.read_csv(counFiles[1],header=None)
den_epis = meta[meta['strain'].isin(wanted[0])]
den_epis.to_csv("denmark_seqs.csv", index=False)
den_epis["date"].to_csv("denmark_date_only.csv", index=False)

wanted = pd.read_csv(counFiles[0],header=None)
swis_epis = meta[meta['strain'].isin(wanted[0])]
swis_epis.to_csv("switzerland_seqs.csv", index=False)
swis_epis["date"].to_csv("switzerland_date_only.csv", index=False)

wanted = pd.read_csv(counFiles[2],header=None)
germ_epis = meta[meta['strain'].isin(wanted[0])]
germ_epis.to_csv("germany_seqs.csv", index=False)
germ_epis["date"].to_csv("germany_date_only.csv", index=False)

