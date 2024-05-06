import pandas as pd
import glob
import os
import datetime
from dateutil.rrule import rrule, MONTHLY
import shutil


#to get all sequences from 2021
#python scripts/filter_clusters.py --in-folder output-7Mar-4muts \
#                   --out-folder output-7Mar-4muts_2021Only \
#                   --start-date 2021-01-01 --end-date 2021-11-30 \
#                   --new-ending 2021Only 

#to get all sequences divided by variant
#python scripts/filter_clusters.py --in-folder output-7Mar-4muts \
#                   --out-folder output-7Mar-4muts_byVar \
#                   --start-date 2021-01-01 \
#                   --variants "Alpha" "Delta" "21K (Omicron)" "21L (Omicron)" "22B (Omicron)" "22D (Omicron)" "22E (Omicron)" "23A (Omicron)"


#to get all sequences divided by month
#python scripts/filter_clusters.py --in-folder output-7Mar-4muts \
#                   --out-folder output-7Mar-4muts_byMon \
#                   --start-date 2021-01-01 \
#                   --by-month

#to get all sequences by variant AND by month
#python scripts/filter_clusters.py --in-folder output-7Mar-4muts \
#                   --out-folder output-7Mar-4muts_byVarMon \
#                   --start-date 2021-01-01 \
#                   --by-month \
#                   --variants "Alpha" "Delta" "21K (Omicron)" "21L (Omicron)" "22B (Omicron)" "22D (Omicron)" "22E (Omicron)" "23A (Omicron)"


if __name__=="__main__":
    import argparse

    parser = parser = argparse.ArgumentParser(description='find identical sequence clusters using Nextstrain metadata file',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--in-folder', help="Folder with the cluster distribution files to read in") #looks for format Denmark_cluster_distribution_dates_5whole.tsv
    parser.add_argument('--start-date', default="2021-01-01", help="only include clusters where min_date is >= this date") #default is 2021 start as that's when seq coverage up
    parser.add_argument('--end-date', default=datetime.datetime.now().strftime("%Y-%m-%d"), help="only include clusters where max_date is <= this date")
    parser.add_argument("--new-ending", default="", help="what to add on end of file name when written out. Added after a '_'")
    parser.add_argument("--variants", nargs="+", help="what variants to divide clusters into, along with dates if specified. Can use 'Alpha' & 'Delta'.")
    parser.add_argument("--by-month", action='store_true', help="Between any specified dates, get clusters by month, will be in folder by month")
    parser.add_argument("--date-for-month", default="mean", choices=['mean', 'median', 'min', 'max'], help="How to define if a cluster falls in a month, can be either mean, median, min, or max")
    parser.add_argument('--out-folder', help="Folder where new cluster files should be output")
    args = parser.parse_args()

    #######

    #one can specify start-end dates with either variants or by-month - but variants & by-month can't be used together at the moment
    start_date = pd.to_datetime(args.start_date)
    end_date = pd.to_datetime(args.end_date)
    in_folder = args.in_folder
    out_folder = args.out_folder 
    new_ending = args.new_ending
    variants = args.variants

    ##############
    #manual data for testing
    #start_date = pd.to_datetime("2021-01-01")
    #end_date = pd.to_datetime("2021-11-30")
    #in_folder = "output-7Mar-4muts"
    #out_folder = "output-7Mar-4muts_2021Only"
    #new_ending = "2021Only"
    #leaving out 22A/BA.4 , 22C/BA.2.12.1 , 22F/XBB as were very low numbers in all three countries
    #beware that for 22D, 22E, and 23A, dominance was not complete
    #                                   BA.1         BA.2                BA.5               BA.2.75      BQ.1            XBB.1.5
    #variants = ["Alpha", "Delta", "21K (Omicron)", "21L (Omicron)", "22B (Omicron)", "22D (Omicron)", "22E (Omicron)", "23A (Omicron)"]
    ##################


    ############

    os.chdir(in_folder)
    all_clus_files = glob.glob("*_cluster_distribution_dates_*")
    #exclude the files that are duplicates but include the strain names for checking -- we don't need that here
    files_to_read = [f for f in all_clus_files if "_strains_" not in f]
    os.chdir("../")

    #if output folder doesn't exist, create it
    print(f"Output will go to folder: {out_folder}\n")
    if not os.path.isdir(out_folder):
        print(f"Couldn't find folder {out_folder} - creating it.\n")
        os.mkdir(out_folder)
    else:
        if len(os.listdir(out_folder)) != 0:
            overwrite_files = False
            print_answer = input("\nOutput folder is not empty - overwrite previous results? (y/n) (Enter is no): ")
            if print_answer in ["y", "Y", "yes", "YES", "Yes"]:
                overwrite_files = True
            if overwrite_files:
                files_to_del = glob.glob(f"{out_folder}/*")
                for f in files_to_del:
                    if os.path.isdir(f): 
                        shutil.rmtree(f)
                    else:
                        os.remove(f)
            else:
                print("Stopping execution. Please specify a new output folder or empty the current output folder.\n\n")
                raise KeyboardInterrupt

    print_lines = sorted([int(len(files_to_read)/20*i) + 1 for i in range(1,20)]) # Print progress in %

    n = 0
    #for each file
    for f in files_to_read:
        clusList = pd.read_csv(f"{in_folder}/{f}", sep="\t", index_col=False)
        #convert to datetime
        clusList['minDate'] = pd.to_datetime(clusList['minDate'])
        clusList['maxDate'] = pd.to_datetime(clusList['maxDate'])
        clusList['meanDate'] = pd.to_datetime(clusList['meanDate'])
        clusList['medianDate'] = pd.to_datetime(clusList['medianDate'])

        #find those between specified dates
        datesClus = clusList[(clusList['minDate']>=start_date) & (clusList['maxDate']<=end_date)]  

        #find those by variant
        #if specified lists of variants...
        if args.variants and not args.by_month:
            print("Just sorting by variant")
            for v in variants:
                varClus = datesClus[(datesClus['clade'] == v)]
                new_f_name = f.replace(".tsv", f"_{new_ending}{v.replace(' ','')}.tsv")
                varClus.to_csv(f"{out_folder}/{new_f_name}", sep="\t", header=True, index=False)

        #if specified by month
        elif args.by_month:
            dt_start = datetime.datetime.date(start_date)
            dt_end = datetime.datetime.date(end_date)
            #get list of all the months in the time period
            mdates = [dt for dt in rrule(MONTHLY, dtstart=dt_start, until=dt_end)]
            #figure out what to define a cluster being 'in a month' as
            date_col_to_use = f"{args.date_for_month}Date"

            for i in range(len(mdates[:-1])):
                #print(f"looking between {mdates[i]} and {mdates[i+1]}")  ##debug
                #find those between specified dates
                #This compares pd dates to datetime dates - it seems to work though??? 
                #Haha it seems like it's fine https://pandas.pydata.org/docs/reference/api/pandas.Timestamp.html
                monClus = datesClus[(datesClus[date_col_to_use]>=mdates[i]) & (datesClus[date_col_to_use]<mdates[i+1])]
                mon_str = mdates[i].strftime('%Y-%m')
                #if doesn't exist, create new folder for this month's files
                mon_path = f"{out_folder}/{mon_str}"
                if not os.path.isdir(mon_path):
                    os.mkdir(mon_path)

                #if also by variants, do this part too
                if args.variants:
                    for v in variants:
                        varClus = monClus[(monClus['clade'] == v)]
                        new_f_name = f.replace(".tsv", f"_{new_ending}{mon_str}_{v.replace(' ','')}.tsv")
                        #Only write a file if there were clusters from this variant in this month!
                        if not varClus.empty:
                            varClus.to_csv(f"{mon_path}/{new_f_name}", sep="\t", header=True, index=False)
                        else:
                            print(f"No clusters found for {v} in {mon_str} in {f}")
                else:
                    new_f_name = f.replace(".tsv", f"_{new_ending}{mon_str}.tsv")  #..._newEnd2021-10.tsv
                    monClus.to_csv(f"{mon_path}/{new_f_name}", sep="\t", header=True, index=False)

        #just write out what we have
        else:
            new_f_name = f.replace(".tsv", f"_{new_ending}.tsv")
            datesClus.to_csv(f"{out_folder}/{new_f_name}", sep="\t", header=True, index=False)


        n+=1
        if n in print_lines:
            print(f"{round(n/len(files_to_read) * 100)}% complete...")

    print("Complete.")
        
