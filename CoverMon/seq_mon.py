__author__ = "Tine Sneibjerg Ebsen"
__version__ = "0.1"



import sys
import os
from os import listdir
from os.path import isfile, isdir, join, exists
import pandas as pd
import numpy as np
from datetime import datetime
import subprocess
import glob
import time
import atexit
import datetime

# new imports
from pathlib import Path

# Create a Cow
#msg = cow.Moose.milk("Hello World")

#print(msg)
if len(sys.argv) < 4:
    raise Exception(f"Missing arguments. The script must contain (1) samplesheet, (2) path to run directory, (3) path to reference or reference directory. Optionally, a region file (4) can be specified if all samples use the same reference genome.")

# TODO: These variables should be set from the command line
print("This is the samplesheet: ", sys.argv[1])
print("This is the run directory: ", sys.argv[2])
print("This is the reference directory:", sys.argv[3])

# Actually given on the command line:

samplesheet = sys.argv[1]
rundir = sys.argv[2]
if ".fa" in sys.argv[3]:
    reference = sys.argv[3]
    one_ref = True
    if len(sys.argv) == 5:
        region_file = sys.argv[4]
        print("This is the region file:", sys.argv[4])
    else: 
        region_file = "NA"
else:
    refdir = sys.argv[3]
    one_ref = False
    region_file = "NA"


tab = "\t"
nl = "\n"


# Print a convincing logo, if the window is big enough.
print()

# Check that input was given.
#if config["samplesheet"] == "NA":
#    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config rundir=\"path/to/samplesheet/\" to the command line call.")
#if config["rundir"] == "NA":
#    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")
# TODO: Implement additional input validation, like checking that the objects given are file and dir respectively.



print(f"These are the parameters given:")
print(f"  samplesheet: {samplesheet}")
print(f"  rundir: {rundir}")
print()

#########################
# Parse the samplesheet #
#########################

samplesheet_extension = samplesheet.split(".")[-1]
print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")

if samplesheet_extension == "xlsx":
    # Uses openpyxl
    df = pd.read_excel(samplesheet, dtype = str)
elif samplesheet_extension == "xls":
    df = pd.read_excel(samplesheet)

else:
    raise Exception(f"The spreadsheet must be excel formatted (.xlsx or .xls)")

# Clean up the spreadsheet
print("Cleaning sample sheet ...                              ", end = "", flush = True)
df.columns = map(str.lower, df.columns) # Lowercase
df.columns = map(str.strip, df.columns) # Remove edge-spaces
df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
df = df.dropna(subset = ["sample_id"])# remove rows not containing a barcode
print("✓")
print(df)
#print("ee", "barcode" in list(df.columns))
# Check that the spreadsheet complies

if one_ref == False:
    print("Checking that the necessary columns exist ...          ", end = "", flush = True)
    for i in ["barcode", "reference","sample_id"]:
        if not i in df.columns:
            raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
    print("✓")

acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)]



print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)

for i in df["barcode"]:
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode ", i, " is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
print("✓")


print("Checking that the barcodes are unique ...              ", end = "", flush = True)
if not len(df["barcode"]) == len(set(df["barcode"])):
    bc_counts = pd.DataFrame(df['barcode'].value_counts())
    bc_counts.columns = ["count"]
    bc_counts = bc_counts[bc_counts["count"] > 1]
    #print(nl, bc_counts)
    raise Exception(f"{nl}One or more barcodes are duplicated. Each barcode may only be used once:{nl}{bc_counts}")
print("✓")

print()
print("These are the samples from the samplesheet you have given:")
print(df.to_string())
print("//")
print()

###################
# Validate rundir #
###################

# Wait for the rundir to occur in the specified path. 
# If it doesn't occur after a specified waiting time, then stop the p

if rundir[-1] == "/":
    print("Removing trailing slash from rundir")
    rundir = rundir[0:-1]


print("Checking that the rundir exists ...                    ", end = "", flush = True)
if not os.path.isdir(rundir):
    raise Exception(f"The rundir does not exist.")
print("✓")

print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)

for i in range(200):
    print("  Looking ... ", end = "", flush = True)
    fastq_pass_bases = glob.glob(rundir + "/**/fastq_pass", recursive = True) # Find any occurrence of the wanted path
    if len(fastq_pass_bases) == 0:
        print("nothing found yet, waiting 10 secs ...")
        time.sleep(10) # Wait 10 seconds.
    elif(i == 10):
        print() # clean newline
        raise Exception("nothing found after 10 tries. Aborting.")
    else: 
        print(f"Found                                    ✓")
        break


if not len(fastq_pass_bases) == 1:
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found:{nl} {str(nl + ' ').join(fastq_pass_bases)}{nl}Please specify a more specific rundir.")


fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases
print(f"Found the following fastq_pass base which will be given to rampart: {nl}  {fastq_pass_base}{nl}")


# base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"This is the batch base directory:{nl}  {base_dir}")



very_long_batch_id = base_dir.split("/")[-1]
print(f"This is the very long batch id:", very_long_batch_id)

date_parse, time_parse, minion_parse, flowcell_parse, arbhash_parse = very_long_batch_id.split("_")

print("date:    ", date_parse)
print("time:    ", time_parse)
print("minion:  ", minion_parse)
print("flowcell:", flowcell_parse)
print("arbhash: ", arbhash_parse)


batch_id = ".".join(very_long_batch_id.split("_")[0:2]) # The first two words (date, time), joined by a dot.
print(f"This is the parsed batch_id:", batch_id)



out_base = os.path.join(base_dir, "CoverMon") # out_base is the directory where the pipeline will write its output to.
print("Creating output directory ", out_base,"...")
subprocess.run(["mkdir",out_base])
print()

disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})


disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])

print("Continuing with the following barcodes:")

# the workflow_table is the table that contains the records where the barcode could be found on the disk.
workflow_table = disk_barcodes_df.merge(df, how='left', on='barcode') # left join (merge) the present barcodes onto the df table.
workflow_table = workflow_table.dropna(subset = ["sample_id"])

print(workflow_table)
print("//")
print()

#################################
# Finally we can start CoveMon. #
#################################

processed_files = []
logfile=very_long_batch_id+".log"
print(logfile)
tmpsam = out_base+'/'+'tmp.sam'
tmpbam = out_base+'/'+'tmp.bam'
# When sequencing, we will check for new files every 60 seconds
seconds_wait = 60
# Ensure we only open report once
open_report = False
still_sequencing = True
while still_sequencing:
    print("Scanning for new fastq files...")
    for index, row in workflow_table.iterrows():
        bam_out = out_base+"/"+row['barcode']+".bam"
        depth = out_base+"/"+row['barcode']+".depth"
        barcode_path = row['barcode_path']
        if not one_ref:
            reference = refdir+"/"+row['reference']
        unprocessed = [barcode_path+"/"+f for f in listdir(barcode_path) if (barcode_path+"/"+f not in processed_files and isfile(join(barcode_path, f)))]

        for f in unprocessed:
        	# Check if we have an existing read mapping to append to. If not, creates the first one and continues the loop without merging.
            if not exists(bam_out):
                map_cmd = f"minimap2 -a -o {tmpsam} {reference} {f}"
                print(map_cmd)
                subprocess.run(map_cmd.split())
                sam2bam = f'samtools sort -O bam -o {bam_out} {tmpsam}'
                print("Creating initial bam: ", sam2bam)
                subprocess.run(sam2bam.split())
                processed_files.append(f)
                print("Number of processed files: ", len(processed_files))
                continue
            # Maps new reads to reference
            map_cmd = f"minimap2 -a -o {tmpsam} {reference} {f}"
            print(map_cmd)
            subprocess.run(map_cmd.split())
            # Sort the mapping files for merging
            sort_cmd = f'samtools sort -O bam -o {out_base}/sorted.bam {tmpsam}'
            print(sort_cmd)
            subprocess.run(sort_cmd.split())
            # Merges the sorted files
            merge_cmd = f'samtools merge -f -o {tmpbam} {out_base}/sorted.bam {bam_out}'
            print(merge_cmd)
            subprocess.run(merge_cmd.split())
            #Sets the new bam to the barcode bam and marks file as processed
            subprocess.run(['mv', tmpbam, bam_out])
            processed_files.append(f)
            print("Number of processed files: ", len(processed_files))
            # Index the new bam and calculate depth. Then creates the monitoring html
            plot_cov_cmd1 = f'samtools index '+ bam_out
            plot_cov_cmd2 = f'samtools depth -aa {bam_out} -o {depth}'
            #Getting the below line to run was a pain. Hence, it is presented as a list.
            plot_cov_cmd3 = ["Rscript", "-e", "\"rmarkdown::render(input = ", "\'scripts/plot_cov.Rmd\',", "params = list(path = ", "\'"+out_base+"\', region_file = ", "\'"+region_file+"\'))\""]
            
            print(plot_cov_cmd1.split())
            subprocess.run(plot_cov_cmd1.split())
            print(plot_cov_cmd2.split())
            subprocess.run(plot_cov_cmd2.split())
            print(plot_cov_cmd3)
            subprocess.run(" ".join(plot_cov_cmd3), shell=True)
            subprocess.run(['mv', 'scripts/plot_cov.html', out_base])

            #Ensure we only open the report and start browser-sync once per run
            if not open_report:
                subprocess.run("gnome-terminal --tab -- browser-sync start -w --no-notify -s \"" + out_base +"\" --host 127.0.0.1 --port 9000 --index \"plot_cov.html\"", shell=True)
                
                open_report = True


    sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
    # Continue the live mappping as long as the sequence summary does not exist
    if len(sequencing_summary_file) == 0:
        print(f"  Still sequencing/basecalling; waiting {seconds_wait} seconds before next scan ...")
        time.sleep(seconds_wait)
    else:
        still_sequencing = False


sequencing_summary_file = sequencing_summary_file[0]
print("  The sequencing summary has been found. Run complete    ✓")
#print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")




