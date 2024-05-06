# Identifying SC2 clusters & estimating R & _k_
Code for identifying clusters to estimate R &amp; _k_ from SC2 sequences

## About the data

The current run is based off a download on 7 Mar 2023 of the full GISAID database, metadata only file. This needs to have been run through Nextclade (with sequences) to generate list of mutations which are added to the metadata file as a column. This file should be stored (as `.tsv`) in a folder called `archived-data` within this repo, but which is not uploaded to Github to prevent dissemenation and also due to size. A backup of the download is also in Emma's Scicore folder `ncov_2021_bern/backup_idents_7Mar23`. 

This repo does not contain any intermediate or results files that contain any data that could violate GISAID Terms. Thus only aggregated data is included in this repo.

# How to run

## Identify clusters from sequence metdata & mutations

1. Run `scripts/lookForIdenticals.py` - Use the arguments to specify what clusters should be identified, for example, from what countries, what minimum-mutations, what amount of the genome to use, and what levels to randomly downsample sequences (before identifying clusters). There are default values which can be used -- for example, the default countries are Switzerland, Germany, Denmark, default minimum mutations is 4, and default percent sequencing (random sampling) is 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, and 1. So a run with these parameters can be accomplished with the command:
```
python scripts/lookForIdenticals.py \
  --meta-in archived-data/metadata-07Mar23.tsv \
  --bad-seqs ../../covariants/scripts/bad_sequences.py \
  --genome whole \
  --out-folder output-5Apr-4muts
```

This will process the raw data from the metadata file, looking at the list of mutations produced by Nextclade. The output will be clusters for each country by different levels of sampling and using the specified amount of the genome.

Takes about 35 mins to run on Scicore.
`20I (Alpha, V1)` is replaced with `Alpha` and all Delta variants are replaced with `Delta`. Other variants remain as-is in original Nextstrain format.

### Note on file naming
Clusters will be output 3 times in different file formats:
- `{Country}_basic_counts_{percent}.tsv` - gives basic counts of number of clusters of different sizes per Country and per percent downsampling (or 100 for no downsampling). 
- `{Country}_cluster_distribution_dates_{percent}{genome}.tsv` - for each cluster gives the size, meanDate, minDate, maxDate, num_muts, invalidDates, and variant, by country and percent downsampling, and amount of genome used
- `{Country}_cluster_distribution_dates_strains_{percent}{genome}.tsv` - same as above but includes the full list of defining mutations used (hash key) and a list of the strain names (mostly useful for checking/debugging)


## Further divide clusters depending on what you'd like to estimate R & k over
2. Run `scripts/filter_clusters.py`, which looks at the output from the previous command (the `{Country}_cluster_distribution_dates_{percent}{genome}.tsv` files only), and extracts only clusters that meet the criteria defined. All downsampling levels are processed!

Note that this was developed to allow easier subdivisions of the clusters, however these resulting folder were not used in the analysis.
Instead, analysis was performed on the un-down-sampled files in `output-5Apr-4muts` folder, with separate code to subdivide these into months by median cluster date.

Arguments can be used to specify what the user would like. `--start-date` and `--end-date` apply to all other options and apply to `minDate` and `maxDate` of a cluster, respectively. By default, `start-date` is 2021-01-01 as this is when sequencing increased and `end-date` is the current date. 

Specify the input folder with `--in-folder` and the output folder with `--out-folder` - it will be created if it doesn't exist, and if it exists & is full, the user can choose to overwrite or not. An ending to all files produced can be specified with `--new-ending` (default is an empty string) which will be added after a `_` symbol. It's recommended to create an appropriately-named output folder for each different 'type' of run.

`--variants` can be used to specify only outputting clusters of a specified variant. The variant name will be appended to the end of files produced (after any `new-ending` specified). Note that variant names with spaces will have spaces removed (ex `23A (Omicron)` will become `23A(Omicron)`). 

`--by-month` will divide the time between `start-date` and `end-date` into months and use `meanDate` for each cluster to find clusters that fall within each month. Depending on downsampling level, some months for some countries may not have any clusters. (Since downsampling is random and performed very early in the previous step.) Note that to keep files organized, sub-folders are created for each month, and the month dates are also added to the file endings (after any `new-ending` specified)

`--variants` and `--by-month` can be used together to find for each month, only clusters with a mean date in that month, from a certain variant. A message will be printed if there are no clusters from that variant for a particular month (ex, there are no Omicron clusters in Jan 2021), and no file will be created. Be aware that some clusters may be found in low numbers in months where they may be due to misdating - for example, 2 Delta sequences in April 2021 in Denmark (not impossible, but treat with care). On the same note, number of clusters may be very low - one should check for an appropriately large number of clusters before attempting to estimate R & k. As in `by-month` new sub-folders are created, and both the variant and date are added to file endings (after any `new-ending` specified).

### Examples of some run paramteters

- To get all sequences from 2021 
Output folder `output-7Mar-4muts_2021Only`
```
python scripts/filter_clusters.py --in-folder output-5Apr-4muts \
  --out-folder output-5Apr-4muts_2021Only \
  --start-date 2021-01-01 --end-date 2021-11-30 \
  --new-ending 2021Only 
```

- To get all sequences divided by month (up to 1 Mar 23)
Output folder `output-7Mar-4muts_byMon`
```
python scripts/filter_clusters.py --in-folder output-5Apr-4muts \
                   --out-folder output-5Apr-4muts_byMon \
                   --start-date 2021-01-01 \
                   --by-month
```

- To get all sequences divided by variant (up to 7 Mar)
Output folder `output-5Apr-4muts_byVar`
```
python scripts/filter_clusters.py --in-folder output-5Apr-4muts \
                   --out-folder output-5Apr-4muts_byVar \
                   --start-date 2021-01-01 \
                   --variants "Alpha" "Delta" "21K (Omicron)" "21L (Omicron)" "22B (Omicron)" "22D (Omicron)" "22E (Omicron)" "23A (Omicron)"
```
Note that 22A (BA.4), 22C (BA.2.12.1) and 22F (XBB) are not used as they were only present in relatively low numbers in CH, DE, and DK.
The variants used correspond to:
```
Alpha  Delta  BA.1           BA.2           BA.5           BA.2.75        BQ.1           XBB.1.5
Alpha  Delta  21K (Omicron)  21L (Omicron)  22B (Omicron)  22D (Omicron)  22E (Omicron)  23A (Omicron)
```


- To get all sequences by variant AND by month (up to 1 Mar 23)
Output folder `output-5Apr-4muts_byVarMon`
```
python scripts/filter_clusters.py --in-folder output-5Apr-4muts \
                   --out-folder output-5Apr-4muts_byVarMon \
                   --start-date 2021-01-01 \
                   --by-month \
                   --variants "Alpha" "Delta" "21K (Omicron)" "21L (Omicron)" "22B (Omicron)" "22D (Omicron)" "22E (Omicron)" "23A (Omicron)"
```

## Further analysis

These files can then be used as input into R code or analyses that estimate R & k
- Remember all downsampling percentages are processed in step 2 - you may want to only use some of them for further analysis
- Check files, especially if by month, to ensure they have sufficient number of clusters to reasonably expect estimation to work

-----------

For Emma: The private version of this repo that contains all files, including sensitive ones that can't be publicly shared is: [https://github.com/emmahodcroft/sc2_rk](https://github.com/emmahodcroft/sc2_rk)