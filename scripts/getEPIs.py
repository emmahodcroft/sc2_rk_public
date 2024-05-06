#Fairly inefficient script to try and get the matching EPI_ISLs
#from the list of sequences which are in clusters with medianDate in 2021

#Have list of strain names from 'newTotalSeq_and_strainList.r', 
#but need to match against the EPI_ISLs - so do this by matching up with the original metadata file

import pandas as pd

input_meta = "archived-data/metadata-07Mar23.tsv"

meta = pd.read_csv(input_meta, sep="\t", usecols=['strain', 'gisaid_epi_isl'], index_col=False)

counFiles = ["output-5Apr-4muts/Denmark_strains.csv", "output-5Apr-4muts/Germany_strains.csv", "output-5Apr-4muts/Switzerland_strains.csv"]

wanted = pd.read_csv(counFiles[0],header=None)
den_epis = meta[meta['strain'].isin(wanted[0])]
den_epis.to_csv("output-5Apr-4muts/denmark_epis.csv", index=False)
den_epis["gisaid_epi_isl"].to_csv("output-5Apr-4muts/denmark_epis_only.csv", index=False)

wanted = pd.read_csv(counFiles[1],header=None)
germ_epis = meta[meta['strain'].isin(wanted[0])]
germ_epis.to_csv("output-5Apr-4muts/germany_epis.csv", index=False)
germ_epis["gisaid_epi_isl"].to_csv("output-5Apr-4muts/germany_epis_only.csv", index=False)

wanted = pd.read_csv(counFiles[2],header=None)
swis_epis = meta[meta['strain'].isin(wanted[0])]
swis_epis.to_csv("output-5Apr-4muts/switzerland_epis.csv", index=False)
swis_epis["gisaid_epi_isl"].to_csv("output-5Apr-4muts/switzerland_epis_only.csv", index=False)


