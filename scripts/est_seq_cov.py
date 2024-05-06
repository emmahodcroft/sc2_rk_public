import json
from datetime import datetime
from collections import defaultdict
import statistics

#This script gets an estimate of the sequence coverage for CH, DE, DK over the 'all' data period (1 Jan - 30 Nov 21)

# run script from 2021_Delta/identicals folder

with open('scripts/case_seq_counts.json', 'r') as f:
  case_seqs = json.load(f)

countries_to_run = ["Switzerland", "Denmark", "Germany"]

countrs = case_seqs["regions"][0]["distributions"]

keepr = {}

for d in countrs:
    if d["country"] in countries_to_run:
        keepr[d["country"]] = d["distribution"]

###    here set cutoff for start of Alpha (2021-01-01) to before Omicron (2021-11-30), 
startdate = datetime.strptime("2021-01-01", "%Y-%m-%d")
enddate = datetime.strptime("2021-11-30", "%Y-%m-%d")

for c in countries_to_run:
c = "Switzerland"

for j in keepr[c]:
    print(j["week"])

keep_date = defaultdict(list)
for c in countries_to_run:
    for j in keepr[c]:
        wk = datetime.strptime(j["week"], "%Y-%m-%d")
        if wk > startdate and wk < enddate:
            keep_date[c].append({"total_cases":j["total_cases"], "total_sequences":j["total_sequences"], "week":wk})


percent_cov_week = defaultdict(list)
for c in countries_to_run:
    for j in keep_date[c]:
        percent_cov_week[c].append(j["total_sequences"]/j["total_cases"])

avg_cov = {}
for c in countries_to_run:
    avg_cov[c] = statistics.mean(percent_cov_week[c])

# Switzerland  0.1827
# Denmark      0.7311
# Germany      0.1136

# Switzerland  18%
# Denmark      73%
# Germany      11%
