import pandas as pd
import re
import os
import glob
from collections import defaultdict
from collections import Counter
from datetime import datetime
import sys
import time
import numpy as np
from mutations_list import *

# run script from 2021_Delta/identicals folder

def remove_bad_sequences(bad_seqs, meta):

    #also remove 3 sequences we know are bad - from Dec 2020 & Jan 2021 that cluster with Mar 2020
    bad_seqs["Switzerland/GE-32977053/2021"] = "2021-01-15"
    bad_seqs["Switzerland/AG-ETHZ-571370/2020"] = "2020-12-30"
    bad_seqs["Switzerland/BL-ETHZ-571375/2020"] = "2020-12-30"
    #remove bad sequences - this *doesn't* check dates (like CoV) - just straight remove
    bad_ones = meta.loc[meta['strain'].isin(bad_seqs.keys())]
    meta.drop(bad_ones.index, inplace=True)
    return meta

def write_out_basic_hash_sizes(hashSizes, country, percent, output_folder):
    #get count of sizes of clusters
    hash_counts = dict(Counter(hashSizes))
    #convert to dataframe to sort
    pd_counts = pd.DataFrame.from_dict(hash_counts, orient="index", columns=["Count"]).sort_index()
    #Write out
    pd_counts.to_csv(f"{output_folder}/{country}_basic_counts_{percent}.tsv", sep="\t", header=True, index_label="Size")

def complex_hash(wanted_muts, nuc_to_use, min_number_muts):
    hashinfo = {}
    for index, row in wanted_muts.iterrows():
        #check number of muts
        num_muts = len(row[nuc_to_use].split(","))
        if num_muts <= min_number_muts:
            continue
        hashkey = hash(row[nuc_to_use])
        if hashkey in hashinfo.keys():
            hashinfo[hashkey]["strains"].append(row["strain"])
            hashinfo[hashkey]["counts"] += 1
        else:
            hashinfo[hashkey] = {}
            hashinfo[hashkey]["strains"] = [row["strain"]]
            hashinfo[hashkey]["counts"] = 1
            hashinfo[hashkey]["muts"] = row[nuc_to_use]
            hashinfo[hashkey]["num_muts"] = num_muts
    return hashinfo

#calculate mean date from a list of datetime objects
def get_mean_date(dates):
    timestampDates = [x.timestamp() for x in dates]
    mn = sum(timestampDates) / len(dates)
    return datetime.fromtimestamp(mn)

#calculate median date from a list of datetime objects -- middle value from sorted list
def get_median_date(dates):
    dates.sort()
    return dates[len(dates)//2]

#calculate variance of dates
def get_variance_dates(dates, meanD):
    sq_diffs = [(x-meanD).days**2 for x in dates]
    return sum(sq_diffs)/len(dates)

def make_date_dict(country_meta):
    country_dates_dict = {}
    for i, row in country_meta.iterrows():
        #only convert to date if possible
        if len(row["date"]) == 10 and "XX" not in row["date"]:
            dat = datetime.strptime(row["date"], "%Y-%m-%d")
        else:
            dat = ""
        country_dates_dict[row["strain"]] = dat
    return country_dates_dict

def make_clade_dict(country_meta):
    country_clade_dict = {}
    for i, row in country_meta.iterrows():
        country_clade_dict[row["strain"]] = row["Nextstrain_clade"]
    return country_clade_dict

def add_hash_info(hashinfo, country_dates_dict, country_clade_dict, genome_set):
    for key, value in hashinfo.items():
        datesInClus = [country_dates_dict[ch] for ch in value["strains"]]
        numberNADates = len([x for x in datesInClus if x == ""])
        validDates = [x for x in datesInClus if x != ""]

        if validDates:
            hashinfo[key]["meanDate"] = get_mean_date(validDates).strftime("%Y-%m-%d")
            hashinfo[key]["medianDate"] = get_median_date(validDates).strftime("%Y-%m-%d")
            hashinfo[key]["minDate"] = min(validDates).strftime("%Y-%m-%d")
            hashinfo[key]["maxDate"] = max(validDates).strftime("%Y-%m-%d")
            hashinfo[key]["dayDuration"] = datetime.strptime(hashinfo[key]["maxDate"], "%Y-%m-%d") - datetime.strptime(hashinfo[key]["minDate"], "%Y-%m-%d")
            hashinfo[key]["dateVariance"] = get_variance_dates(validDates, datetime.strptime(hashinfo[key]["meanDate"],"%Y-%m-%d"))
        else:
            hashinfo[key]["meanDate"] = ""
            hashinfo[key]["medianDate"] = ""
            hashinfo[key]["minDate"] = ""
            hashinfo[key]["maxDate"] = ""
            hashinfo[key]["dayDuration"] = ""
            hashinfo[key]["dateVariance"] = ""
        hashinfo[key]["invalidDates"] = numberNADates

        #replace all Deltas with Delta
        #Replace Nextstrain long alpha name with Alpha
        rep = {"21J (Delta)": "Delta", "21A (Delta)": "Delta", "21I (Delta)": "Delta",
                "20I (Alpha, V1)": "Alpha"}
               #"21K (Omicron)": "Omicron", "21L (Omicron)": "Omicron", "21M (Omicron)": "Omicron"}

        cladesInClus = [country_clade_dict[ch] for ch in value["strains"]]
        replaced_clades = [x if x not in rep else rep[x] for x in cladesInClus]
        hashinfo[key]["all_clades"] = replaced_clades
        setClades = set(replaced_clades)
        if len(setClades) != 1:
            if genome_set == "whole":
                print(f"not all clades are same! for key {key}")
                print(setClades)
            hashinfo[key]["clade"] = ""
        else:
            hashinfo[key]["clade"] = replaced_clades[0]

    return hashinfo

def remove_variant_muts(muts,meta):
    for c in muts.keys():
        wanted_strains = meta['Nextstrain_clade'].isin([c])
        cur_var = meta[wanted_strains]
        
        #have to do for all three substiutions - can't split genome after doing this because
        #the variant name is added to the subs!!
        for subs in ["substitutions","substitutions_firsthalf","substitutions_secondhalf"]:
            sep_muts = "|".join(muts[c])  #concat muts with |
            pat = re.compile(sep_muts)  #turn into compiled regex
            new_subs = cur_var[subs].str.replace(pat,"") #replace them all with nothing
            #now clean up duplicate commas ",,"
            com_pat = re.compile(",+")
            new_subs = new_subs.str.replace(com_pat,",")
            #now clean up leading and trailing commas
            se_pat = re.compile("^,|,$")
            new_subs = new_subs.str.replace(se_pat,"")
            # add the variant to the start
            var_names = np.full(len(new_subs),","+c) #create np array of variant names - add , to delimit
            new_subs = new_subs.str.cat(var_names,join="right") #this adds to end, can't make it add to start

            #replace the existing subs with the new ones
            meta.loc[wanted_strains, subs] = new_subs

    return meta


#--meta-in archived-data/metadata-07Mar23.tsv  #currently input_meta
#--bad-seqs ../../covariants/scripts/bad_sequences.py #currently bad_seqs_file
#--genome whole  #or whole,firsthalf,secondhalf  #currently amount_genome
#--percent-seq 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1 ##maybe unneeded since default?? #currently percent_sequencing
#--min-muts 4 #has default #currently min_number_muts  -- but maybe run again with 2
#--countries Switzerland, Germany, Denmark #has default -- currently countries_to_run
#--out-folder output-7Mar #currently output_folder

if __name__=="__main__":
    import argparse

    parser = parser = argparse.ArgumentParser(description='find identical sequence clusters using Nextstrain metadata file',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--meta-in', help="Nextstrain metadata file")
    parser.add_argument('--bad-seqs', help="CoV list of bad seqs to be excluded")
    parser.add_argument('--genome', choices=['whole','firsthalf','secondhalf'], nargs="+", help="How much of genome to use, can select multiple")
    parser.add_argument('--percent-seq', nargs="+", type=float, default=[0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1], help="amount of sequencing to downsample")
    parser.add_argument('--min-muts', type=int, default=4, help="minimum number of muts to use to define clusters")
    parser.add_argument('--countries', nargs="+", default=["Switzerland", "Denmark", "Germany"], help="Countries to run")
    parser.add_argument('--out-folder', help="output folder")
    args = parser.parse_args()


    start_time = time.time()

    bad_seqs_file = args.bad_seqs
    countries_to_run = args.countries  #ENSURE MATCHES SPELLING IN SEQUENCE NAME & GISAID NAME
    percent_sequencing = args.percent_seq
    amount_genome = args.genome
    min_number_muts = args.min_muts
    output_folder = args.out_folder 
    input_meta = args.meta_in

    run_cond_info = time.ctime() + "\n\n"

    #read in bad_seqs file so that can exclude known bad seqs - get from CoVariants
    #bad_seqs_file = "../covariants/scripts/bad_sequences.py"
    if not os.path.exists(bad_seqs_file):
        print("Couldn't find the bad_sequences.py file. Please adjust path to this file. If necessary, download it from https://github.com/hodcroftlab/covariants/blob/master/scripts/bad_sequences.py")
        raise KeyboardInterrupt

    #what countries to run? ENSURE MATCHES SPELLING IN SEQUENCE NAME & GISAID NAME
    #countries_to_run = ["Switzerland", "Denmark", "Germany"]
    #percent_sequencing = [0.5, 1.0] #[0.25, 0.5, 0.75, 1]
   # percent_sequencing = [0.25, 0.5, 0.75, 1]
    #percent_sequencing = [0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1]
    #amount_genome = ["whole", "firsthalf", "secondhalf"]
    #min_number_muts = 4
    #output_folder = "output"

    print(f"Output will go to folder: {output_folder}\n")
    run_cond_info = run_cond_info + "Output will go to folder: {output_folder}\n"

    if not os.path.isdir(output_folder):
        print(f"Couldn't find folder {output_folder} - creating it.\n")
        os.mkdir(output_folder)
    else:
        if len(os.listdir(output_folder)) != 0:
            overwrite_files = False
            print_answer = input("\nOutput folder is not empty - overwrite previous results? (y/n) (Enter is no): ")
            if print_answer in ["y", "Y", "yes", "YES", "Yes"]:
                overwrite_files = True
            if overwrite_files:
                files_to_del = glob.glob(f"{output_folder}/*")
                for f in files_to_del:
                    os.remove(f)
            else:
                print("Stopping execution. Please specify a new output folder or empty the current output folder.")
                raise KeyboardInterrupt


    #list variants
    all_vars = muts.keys()

    starting_info = (f"\nCountries to be run: {countries_to_run}\n"
                    f"Percent sampling to be run: {percent_sequencing}\n" 
                    f"At 100% sampling, these genome subsets will be used: {amount_genome}\n" 
                    f"A minimum of {min_number_muts} mutations will required to determine a cluster\n")

    print(starting_info)

    #print(f"\nCountries to be run: {countries_to_run}")
    #print(f"Percent sampling to be run: {percent_sequencing}\n")
    #print(f"At 100% sampling, these genome subsets will be used: {amount_genome}\n")
    #print(f"A minimum of {min_number_muts} mutations will required to determine a cluster\n")

    run_cond_info = run_cond_info + starting_info+"\n"

    # Use archived data files
    #dont use muts file anymore, just uses metadata (from march 2022)
    #muts_file = (
        #"archived-data/mutation_summary_gisaid-9Nov21.tsv"  
    #    "archived-data/mutation_summary_gisaid-15Feb22.tsv"  
    #)
    #input_meta = "archived-data/downloaded_gisaid-9Nov21.tsv"  
    #input_meta = "archived-data/metadata-15Feb22.tsv"  
    #input_meta = "archived-data/metadata-6Jun22.tsv"  
    #input_meta = "archived-data/metadata-14Sept22.tsv"

    #stopped doing this because seqs in muts were not always in meta (unsure why)
    #print(f"Now reading in mutations...")
    ## Read in mutations file & remove bad sequences
    #muts = pd.read_csv(muts_file, sep="\t", index_col=False)
    ##read in bad_seqs file so that can exclude known bad seqs - get from CoVariants
    #bad_seqs_file = "../../covariants/scripts/bad_sequences.py"
    ##code from https://stackoverflow.com/questions/6357361/alternative-to-execfile-in-python-3
    #with open(bad_seqs_file, "rb") as source_file:
    #    code = compile(source_file.read(), bad_seqs_file, "exec")
    #exec(code)
    #muts = remove_bad_sequences(bad_seqs, muts)

    print("Now reading in metadata...\n")
    #read in metadata
    meta = pd.read_csv(input_meta, sep="\t", usecols=['strain', 'date', 'country', 'Nextstrain_clade', 'substitutions'], index_col=False)
    meta = meta.fillna("")
    #read in bad_seqs file so that can exclude known bad seqs - get from CoVariants
    #code from https://stackoverflow.com/questions/6357361/alternative-to-execfile-in-python-3
    with open(bad_seqs_file, "rb") as source_file:
        code = compile(source_file.read(), bad_seqs_file, "exec")
    exec(code)
    meta = remove_bad_sequences(bad_seqs, meta)

    #only keep the countries we are using:
    wanted_countries = meta.country.isin(countries_to_run)
    meta = meta[wanted_countries]


    tgenome = time.time()

    #now get the mutations in first half, or second half of genome
    cutoff = 15000 #cutoff at 15,000 bp - roughly half of genome
    split_muts = meta.substitutions.fillna('').apply(lambda x: {int(y[1:-1]):y for y in x.split(',') if y})
    new_muts_lower = split_muts.apply(lambda x: ','.join([val for key,val in x.items() if key<cutoff]) )
    new_muts_upper = split_muts.apply(lambda x: ','.join([val for key,val in x.items() if key>=cutoff]) )

    meta["substitutions_firsthalf"] = new_muts_lower
    meta["substitutions_secondhalf"] = new_muts_upper

    tgenome1 = time.time()
    print(f"Subsetting genome took {round((tgenome1-tgenome)/60,1)} min to run")
    run_cond_info = run_cond_info + f"Subsetting genome took {round((tgenome1-tgenome)/60,1)} min to run\n"

    #modify the mutations to remove defining mutations from each variant (Alpha, Omicron, Delta)
    #the dict with the mutations is called muts
    #this must be done AFTER splitting out the genome because otherwise the addition of the variant will mess up splitting! 
    tmuts = time.time()

    meta = remove_variant_muts(muts,meta)

    tmuts1 = time.time()
    print(f"Removing variant mutations took {round((tmuts1-tmuts)/60,1)} min to run")
    run_cond_info = run_cond_info + f"Removing variant mutations took {round((tmuts1-tmuts)/60,1)} min to run\n"

    #print(f"Size before reduction: {sys.getsizeof(meta)}")
    #reduce memory by only keeping countries we need -NO! Some aren't listed by country.
    #meta = meta[meta['country'].isin(countries_to_run)]
    #print(f"Size after reduction: {sys.getsizeof(meta)}")

    #store the summary information
    summary_infos = []
    too_few_muts = {}

    #run for each country
    for country in countries_to_run:
        t0 = time.time()
        print(f"\nCurrently running {country}")

        #find those from country
        #wanted_strains = muts["Unnamed: 0"].str.contains(country)
        wanted_strains = meta['country'].str.contains(country)
        wanted_muts = meta[wanted_strains]

        for perc in percent_sequencing:

            tperc = time.time()
            #convert to nice format for file names & output
            percent = round(perc*100)

            print(f"\nCurrently sampling at {percent}%")

            if perc == 1: #if doing 100% - do all genome sampling
                genome_to_do = amount_genome
            else: #else, do just whole genome
                genome_to_do = ["whole"]

            for genome_set in genome_to_do:

                tgenome = time.time()
                print(f"\nCurrently using {genome_set} genome")

                #set which column in metadata to use as the key
                if genome_set == "firsthalf":
                    nuc_to_use = "substitutions_firsthalf"
                elif genome_set == "secondhalf":
                    nuc_to_use = "substitutions_secondhalf"
                else:
                    nuc_to_use = "substitutions"

                #sample -- this samples without replacement, so 100% should be same as without sampling https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.sample.html
                sampled_muts = wanted_muts.sample(frac=perc, random_state=1986)

                #make hash table
                hashes = defaultdict(list)
                for index, row in sampled_muts.iterrows():
                    if len(row[nuc_to_use].split(",")) <= min_number_muts:
                        continue
                    hashes[hash(row[nuc_to_use])].append(row["strain"])

                clus_sum =  (f"\nShowing only clusters with > {min_number_muts} mutations.\n"
                            f"\tTotal number of {country} sequences ({percent}%, {genome_set}): {len(sampled_muts)}\n"
                            f"\tOf these, {len(hashes)} are unique\n\n")

                print(clus_sum)

                #print(f"Showing only clusters with > {min_number_muts} mutations.")
                #print(f"\tTotal number of {country} sequences ({percent}%, {genome_set}): {len(sampled_muts)}")
                #print(f"\tOf these, {len(hashes)} are unique")

                run_cond_info = run_cond_info + clus_sum

                #get sizes of clusters
                hashSizes = []
                for key, value in hashes.items():
                    hashSizes.append(len(value))

                sizesdf = pd.DataFrame(data=hashSizes)
                print(f"\tSummary of clusters for {country} ({percent}%, {genome_set}):")
                print(sizesdf.describe())
                #store this:
                dc = { 'Sample': f"{country} {percent}%", 
                    'Genome': genome_set,
                    'Total': len(sampled_muts)}
                df1 = pd.DataFrame.from_dict(dc, orient='index')
                summarydf = df1.append(sizesdf.describe())
                summary_infos.append(summarydf)

                #write out basic distribution
                write_out_basic_hash_sizes(hashSizes, country, percent, output_folder)

                print("\tCreating a more complex hash table...")
                #make a more complex hash table to store more info - like dates, muts, seqs etc
                hashinfo = complex_hash(sampled_muts, nuc_to_use, min_number_muts)

                # but check - how many of these actually have a reasonable number of muts?
                #if perc == 1.0 and genome_set == "whole":
                #    count_low_muts = []
                #    for key, value in hashinfo.items():
                #        if len(value["muts"].split(",")) < 3:
                #            count_low_muts.append(len(value["muts"].split(",")))
                #    too_few_muts[country] = count_low_muts

                print("\tGetting country-specific metadata...")
                #Match up metadata to get date information
                want_strains = list(sampled_muts["strain"])
                #get country meta
                country_meta = meta.loc[meta['strain'].isin(want_strains)]
                #make dict for easier access to dates
                country_dates_dict = make_date_dict(country_meta)
                #make dict for easier act to clade
                country_clade_dict = make_clade_dict(country_meta)

                print("\tAdding dates & clades to hash information...")
                # add the new date info & clus info to the complex hashtable
                hashinfo = add_hash_info(hashinfo, country_dates_dict, country_clade_dict, genome_set)

                print("\tPrinting out files")
                #convert to dataframe
                cluster_all_info = pd.DataFrame.from_dict(hashinfo, orient="index")

                #take just columns we want, in order we want
                cluster_info = cluster_all_info[["counts", "medianDate", "meanDate", "minDate", "maxDate", "dayDuration", "dateVariance", "num_muts", "invalidDates", "muts", "strains", "clade"]].copy()
                cluster_info.sort_values(by="counts", inplace=True)

                #make another file with just the basic information (no seq list or mut list)
                cluster_info_minimal = cluster_info[["counts", "medianDate", "meanDate", "minDate", "maxDate", "dayDuration", "dateVariance", "num_muts", "invalidDates", "clade"]].copy()

                cluster_info.to_csv(f"{output_folder}/{country}_cluster_distribution_dates_strains_{percent}{genome_set}.tsv", sep="\t", header=True, index=False)
                cluster_info_minimal.to_csv(f"{output_folder}/{country}_cluster_distribution_dates_{percent}{genome_set}.tsv", sep="\t", header=True, index=False)

                tgenome1 = time.time()
                count_perc_time = f"\tSampling {country} at {percent}% with {genome_set} genome took {round((tgenome1-tgenome)/60,1)} min to run"
                print(count_perc_time)
                run_cond_info = run_cond_info + count_perc_time + "\n"

            tperc1 = time.time()
            count_time = f"\tSampling {country} at {percent}% took {round((tperc1-tperc)/60,1)} min to run"
            print(count_time)
            run_cond_info = run_cond_info + count_time + "\n"

        t1 = time.time()
        country_time = f"{country} took {round((t1-t0)/60,1)} min to run"
        print(country_time)
        run_cond_info = run_cond_info + country_time + "\n"

end_time = time.time()
whole_time = f"The whole thing took {round((end_time-start_time)/60,1)} min to run"
print(whole_time)
run_cond_info = run_cond_info + whole_time + "\n"


#put all summarys into one dataframe
all_summ_info = pd.concat(summary_infos, axis=1, ignore_index=True)
all_summ_info.to_csv(f"{output_folder}/all_summary.tsv", sep="\t", header=False, index=True)

#put all output into one file
with open(f"{output_folder}/output_summary.txt", "w") as file:
    file.writelines(run_cond_info)