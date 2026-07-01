__author__ = "Tine Sneibjerg Ebsen, Kat Steinke"
__version__ = "0.3"

import glob
import logging
import os
import pathlib
import subprocess
import sys
import time

from argparse import ArgumentParser
from os import listdir
from os.path import isfile, isdir, join, exists
from typing import Optional


import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_log = logging.StreamHandler()
logger.addHandler(console_log)

# TODO add convenience wait function?

class CoverReport:
    """A coverage report generated on the base of input data and parameters, written to an output directory,
    with tracking of open/close status and processed input files.
    Attributes:
        workflow_table:     mapping of fastq_pass dirs to barcodes
        sample_sheet:       path to sample sheet
        threshold:          minimum coverage required
        maxDepth:           maximum coverage reported
        reference:          reference file or directory
        out_base:           output base directory
        is_open:            whether the report is currently opened in the browser
        processed_files:    the files processed by this instance of CoverMon
        region_file:        path to bed file with regions in the reference
                            (optional, only allowed if reference is a single file)
    """
    def __init__(self, workflow_table: pd.DataFrame, sample_sheet: str, threshold: int, maxDepth: int, reference: str,
                 out_base: str, region_file: Optional[str] = None) \
            -> None:
        """Initialize a CoverReport object.

        Args:
            workflow_table: mapping of fastq_pass dirs to barcodes
            sample_sheet:   path to sample sheet
            threshold:      minimum coverage required
            maxDepth:       maximum coverage reported
            reference:      reference file or directory
            out_base:       output base directory
            region_file:    path to bed file with regions in the reference
                            (optional, only allowed if reference is a single file)

        Raises:
            ValueError  if the maximum reported coverage is lower than the minimum coverage required,
                        or if a region file is given together with a reference base directory
                        (=different references per sample)
        """
        # sense checks first
        if threshold > maxDepth:
            raise ValueError("Highest reported coverage must exceed minimum required coverage.")
        if pathlib.Path(reference).is_dir() and region_file is not None:
            raise ValueError("Cannot supply a region file when different references per sample are used "
                             "(base ref is a directory)")
        self.workflow_table = workflow_table
        self.sample_sheet = sample_sheet
        self.threshold = threshold
        self.maxDepth = maxDepth
        self.reference = pathlib.Path(reference)  # TODO: rework to take a Path to begin with?
        self.out_base = pathlib.Path(out_base)
        self.is_open = False  # TODO: can we assume this?
        self.processed_files = []
        self.region_file = region_file

    def __eq__(self, other: object) -> bool:
        if isinstance(other, CoverReport):
            return ((self.sample_sheet == other.sample_sheet
                     and self.threshold == other.threshold
                     and self.maxDepth == other.maxDepth
                     and self.reference == other.reference
                     and self.out_base == other.out_base
                     and self.is_open == other.is_open
                     and self.processed_files == other.processed_files
                     and self.workflow_table.equals(other.workflow_table)
                     and self.region_file == other.region_file))
        return False

    def __repr__(self):
        return(f"CoverReport(sample_sheet={self.sample_sheet}, threshold={self.threshold}, maxDepth={self.maxDepth},"
               f" reference={self.reference}, region_file={self.region_file},"
               f" out_base={self.out_base}, processed_files={self.processed_files},"
               f" is_open={self.is_open})\n"
               f"Workflow table:\n{self.workflow_table.head().to_string()}")

    def set_open_status(self, open_status: bool) -> None:
        """Set the report's status to open (True) or closed (False)

        Args:
            open_status:    whether the report is currently opened in the browser (True) or not (False)
        """
        self.is_open = open_status

    # TODO: move to paths in second pass
    def add_processed_file(self, processed: str) -> None:
        """Add a file to the record of processed files.

        Args:
            processed:  the path to the processed file

        """
        self.processed_files.append(processed)

    # TODO bulk add processed files from file



# TODO config?

# TODO move the helpers out so we can test them individually
def parse_samplesheet(samplesheet: str) -> pd.DataFrame:
    """Clean an Excel sample sheet.

    Args:
        samplesheet: the path to the sample sheet

    Returns:
        The sample sheet with column names cleaned
        (lowercase, leading/trailing whitespace stripped, remainder changed to underscores),
        whitespace removed in barcode names,
        and any rows without sample ID removed.
    Raises:
        Exception   if the sheet is not in the correct format
    """
    samplesheet_extension = samplesheet.split(".")[-1]
    print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")

    if samplesheet_extension == "xlsx":
        # Uses openpyxl
        df = pd.read_excel(samplesheet, dtype=str)
    elif samplesheet_extension == "xls":
        df = pd.read_excel(samplesheet)
    else:
        raise Exception("The spreadsheet must be excel formatted (.xlsx or .xls)")  # TODO raise more specifically?

    # Clean up the spreadsheet
    print("Cleaning sample sheet ...                              ", end="", flush=True)
    df.columns = map(str.lower, df.columns)  # Lowercase
    df.columns = map(str.strip, df.columns)  # Remove edge-spaces
    df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns)  # Replace spaces with underscore
    df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ",
                                                                                      "")))  # Because we are later going to join using this column, it is necessary to strip it for spaces.
    df = df.dropna(subset=["sample_id"])  # remove rows not containing a barcode
    print("✓")
    print(df)
    return df

# TODO: lots of prints that could be logging instead
def validate_samplesheet(samplesheet: pd.DataFrame, ref_is_file: bool) -> None:  # TODO: might be nicer with a bool output
    """Check that the sample sheet contains the correct barcodes and no duplicates, and all columns are present

    Args:
        samplesheet:    the processed sample sheet
        ref_is_file:    whether the reference is a single file (True) or a directory (False)

    Raises:
        Exception   if a required column (barcode, reference, sample ID) is missing in a setup with a reference dir,
                    if barcodes aren't correctly formatted,
                    or if a barcode is duplicated
    """
    # TODO: could be useful to raise more specifically
    # Check that the samplesheet contains the reference column if a refdir is given
    if not ref_is_file:
        print("Checking that the necessary columns exist ...          ", end = "", flush = True)
        for i in ["barcode", "reference","sample_id"]:  # TODO we probably also want to check this when we don't have a reference column
            if not i in samplesheet.columns:
                raise Exception("The sample sheet is missing a necessary column. "
                                f"The sample sheet must contain the column {i}, "
                                f"but it only contains {sorted(samplesheet.columns.tolist())}")
        print("✓")
    # Check that the barcodes look correct
    acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)] + [f"RB{i:02d}" for i in range(1,97)]

    print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)

    for i in samplesheet["barcode"]:  # TODO catch all broken barcodes at once
        if not i in acceptable_barcodes:
            raise Exception(f"The given barcode {i} is not an acceptable barcode. "
                            f"Here is a list of acceptable barcodes for inspiration:\n{' '.join(acceptable_barcodes)}")
    print("✓")


    print("Checking that the barcodes are unique ...              ", end = "", flush = True)
    if not len(samplesheet["barcode"]) == len(set(samplesheet["barcode"])):
        bc_counts = pd.DataFrame(samplesheet['barcode'].value_counts())
        bc_counts.columns = ["count"]
        bc_counts = bc_counts[bc_counts["count"] > 1]
        raise Exception(f"\nOne or more barcodes are duplicated. Each barcode may only be used once:\n{bc_counts}")
    print("✓")

    print()
    print("These are the samples from the samplesheet you have given:")
    print(samplesheet.to_string())
    print("//")
    print()

def validate_rundir(rundir: str) -> tuple[str, str]:
    """Check whether the rundir exists, and find the fastq_pass directory and its parent dir, possibly waiting for
    the fastq_pass directory to be created.

    Args:
        rundir: the directory containing sequence data *somewhere*

    Returns:
        The parent directory of the fastq_pass directory and the fastq_pass directory
    Raises
        Exception   if the rundir does not exist,
                    if the fastq_pass directory was not found after the end of the waiting time,
                    or if multiple fastq_pass directories were found
    """
    if rundir[-1] == "/":
        print("Removing trailing slash from rundir")
        rundir = rundir[0:-1]


    print("Checking that the rundir exists ...                    ", end = "", flush = True)
    if not os.path.isdir(rundir):
        raise Exception("The rundir does not exist.")
    print("✓")

    print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)
    # Wait for the rundir to occur in the specified path.
    # If it doesn't occur after a specified waiting time, then stop the p
    # TODO: convenience waiting function here?
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


    if not len(fastq_pass_bases) == 1:  # TODO: can we just borrow SnakeAmp's get_fastq_pass_parent? That could cope with multiple fastq_pass dirs if you give it the one you want explicitly
        raise Exception("There seems to be more than one fastq_pass sub-directory beneath the given rundir."
                        " These paths were found:\n"
                        f" {'\n '.join(fastq_pass_bases)}\n"
                        "Please specify a more specific rundir.")


    fastq_pass_base = fastq_pass_bases[0]
    del fastq_pass_bases
    print(f"Found the following fastq_pass base which will be given to CoverMon: \n  {fastq_pass_base}\n")


    # base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
    # TODO pathlib can do a lot of cleaning on these things;
    base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
    print(f"This is the batch base directory:\n  {base_dir}")

    return base_dir, fastq_pass_base


def create_workflow_table(samplesheet: pd.DataFrame, fastq_pass_dir: str) -> pd.DataFrame:
    """Record locations of barcode directories for each barcode.

    Args:
        samplesheet:    the sample sheet containing barcodes, sample IDs and optionally references
        fastq_pass_dir: the path to the fastq_pass directory

    Returns:
        The sample sheet with the locations of the barcode fastq directories added for each sample

    Raises:
        Exception   if the barcodes' format is invalid
    """
    disk_barcodes_list  = sorted(glob.glob(fastq_pass_dir + "/barcode*")) # Find all fastq_pass/barcode* directories
    disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})


    disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
    if "RB" in samplesheet["barcode"][0]:
        disk_barcodes_df = disk_barcodes_df.assign(barcode = ["RB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])
    elif "NB" in samplesheet["barcode"][0]:
        disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])
    else:
        raise Exception(f"Barcodes in samplesheet are not acceptable")

    # ensure consistent column format
    sample_cols = samplesheet.columns.tolist()
    barcode_cols = ["barcode_path", "barcode_basename"]
    out_cols = [*sample_cols, *barcode_cols]

    print("Continuing with the following barcodes:")

    # the workflow_table is the table that contains the records where the barcode could be found on the disk.
    workflow_table = disk_barcodes_df.merge(samplesheet, how='left', on='barcode') # left join (merge) the present barcodes onto the df table.
    workflow_table = workflow_table.dropna(subset = ["sample_id"]).reindex(columns = out_cols)

    print(workflow_table)
    print("//")
    print()
    return workflow_table

def write_to_processed(to_write: str, out_dir: str) -> None:
    """Append the given string to "processed_files.txt" in the specified outdir

    Args:
        to_write: the string to write to the file
        out_dir:  the base directory in which to write to the file

    """
    with open(f"{out_dir}/processed_files.txt","a", encoding = "utf-8") as processed_files_txt:
        processed_files_txt.write(f"{to_write}\n")


# TODO split subcommands out:
#  - create bam
#  - append bam
#  - index and get depths
#  - eventually we want to split plot_cov out as well
def create_bam(fastq: str, out_base: pathlib.Path, out_bam: pathlib.Path, reference: pathlib.Path)\
        -> None:
    """Align a fastq file for a new barcode to the reference and output results.

    Args:
        out_bam:
        fastq:      the fastq file to process
        out_base:   the output directory to write results to
        reference:  the reference file to align the fastq to

    """
    tmpsam = out_base / "tmp.sam"
    map_cmd = f"minimap2 -a -o {tmpsam} {reference} {fastq}"
    logger.debug(map_cmd)
    map_process = subprocess.run(map_cmd.split(), check = True)
    sam2bam = f'samtools sort -O bam -o {out_bam} {tmpsam}'
    logger.info("Creating initial bam")
    logger.debug(sam2bam)
    subprocess.run(sam2bam.split(), check = True)


def append_bam(fastq: str, out_base: pathlib.Path, out_bam: pathlib.Path, reference: pathlib.Path)\
        -> None:
    """Append alignment of new fastqs for an existing barcode to existing results

    Args:
        fastq:      the fastq file to process
        out_base:   the output directory to write results to
        out_bam:    the bam file to combine results with
        reference:  the reference file to align the fastq to
    """
    tmpsam = out_base / "tmp.sam"
    tmpbam = out_base / "tmp.bam"
    map_cmd = f"minimap2 -a -o {tmpsam} {reference} {fastq}"
    logger.debug(map_cmd)
    subprocess.run(map_cmd.split(), check = True)
    # Sort the mapping files for merging
    sort_cmd = f'samtools sort -O bam -o {out_base}/sorted.bam {tmpsam}'  # TODO keep in mind when allowing more
    logger.debug(sort_cmd)
    subprocess.run(sort_cmd.split(), check = True)
    # Merges the sorted files
    merge_cmd = f'samtools merge -f -o {tmpbam} {out_base}/sorted.bam {out_bam}'
    logger.debug(merge_cmd)
    subprocess.run(merge_cmd.split(), check = True)
    # Sets the new bam to the barcode bam and marks file as processed
    subprocess.run(['mv', tmpbam, out_bam], check = True)  # TODO keep in mind when moving


def get_depth(bam: pathlib.Path, depth_out: pathlib.Path) -> None:
    """Get depths based on bam file

    Args:
        bam:        the bam file to process
        depth_out:  the depth file to write to

    """
    plot_cov_cmd1 = 'samtools index ' + str(bam)
    plot_cov_cmd2 = f'samtools depth -aa {bam} -o {depth_out}'

    logger.debug(plot_cov_cmd1)
    subprocess.run(plot_cov_cmd1.split(), check = True)
    logger.debug(plot_cov_cmd2)
    subprocess.run(plot_cov_cmd2.split(), check = True)


# TODO: ALL OF THE VARS here, let's feed this an object specifying report params
# We need:
# already specified as args:
# - workflow_table - mapping of fastq pass dirs to barcodes, could be attached to the report instead
# - reference - reference file/dir if given (attach to report, build the per-line ref path?)
# - samplesheet: path to sample sheet
# - open_report: whether the report is open or not - could attach that to the CoverReport?

# used to be in the main function
# - out_base: the output base dir (attach to report?)
# - one_ref: whether the ref is a file or a dir - maybe just go off the CoverReport's reference, if that is_dir or not
# - refdir: replace with CoverReport.reference
# - processed_files: attach to CoverReport?
# - threshold: attach to CoverReport - take from config in the future?
# - maxDepth: attach to CoverReport - take from config in the future?
# - region_file: attach to CoverReport
def update_plot(active_report: CoverReport):
    """Create or update the coverage plot for a given report and show plot in a browser window

    Args:
        active_report:  the CoverReport to plot
    """
    # TODO: allow for different outfiles here - everything needs a unique ID eventually
    # Define temporary mapping files which will be used for merging new mapping of new output files with existing
    workflow_table = active_report.workflow_table.copy(deep=True)
    out_base = active_report.out_base
    logger.info("Scanning for new fastq files...")
    # TODO: can we do bam_out, depth, outside of the loop at least?
    for index, row in workflow_table.iterrows():
        bam_out = out_base / f"{row['barcode']}.bam"
        depth =  out_base / f"{row['barcode']}.depth"
        barcode_path = row['barcode_path']  # TODO pathlib-based handling here
        if active_report.reference.is_dir():
            reference = active_report.reference / row['reference']
        else:
            reference = active_report.reference
        unprocessed = [barcode_path + "/" + f for f in listdir(barcode_path) if
                       (barcode_path + "/" + f not in active_report.processed_files and isfile(join(barcode_path, f))
                        and (f.endswith(".fastq") or f.endswith(".fastq.gz")))]  # TODO we can do this more nicely in pathlib

        for f in unprocessed:
            # Check if we have an existing read mapping to append to. If not, creates the first one and continues the loop without merging.
            if not exists(bam_out):
                create_bam(f, out_base, bam_out, reference)
                active_report.add_processed_file(f)
                write_to_processed(f, str(out_base))

                logger.info(f"Number of processed files: {len(active_report.processed_files)}")
            else:
                # Maps new reads to reference
                append_bam(f, out_base, bam_out, reference)
                active_report.add_processed_file(f)
                write_to_processed(f, str(out_base))

                logger.info(f"Number of processed files: {len(active_report.processed_files)}")
            # Index the new bam and calculate depth. Then creates the monitoring html
            get_depth(bam_out, depth)
            # Getting the below line to run was a pain. Hence, it is presented as a list.
            # TODO: make this a script of its own so we can run it more nicely
            # TODO: ensure it writes to a unique file for each run
            region_file = active_report.region_file if active_report.region_file else ""
            plot_cov_cmd3 = ["Rscript", "-e", "\"rmarkdown::render(input = ", "'scripts/plot_cov.Rmd',",
                             "params = list(threshold = ", str(active_report.threshold),
                             ", maxDepth= ", str(active_report.maxDepth), ", path = ",
                             "'" + str(out_base) + "', samplesheet = ", "'" + active_report.sample_sheet + "', region_file = ",
                             "'" + region_file + "'))\""]

            logger.debug(" ".join(plot_cov_cmd3))
            subprocess.run(" ".join(plot_cov_cmd3), shell=True, check = True)
            # TODO can't we output it somewhere else?
            subprocess.run(['mv', 'scripts/plot_cov.html', out_base], check = True)

            # Starts browser-sync in a new terminal if the report is not open. This will not work on windows or macOS.
            logger.info("Updated plot")
            if not active_report.is_open:
                logger.info("Opening report")
                subprocess.run(
                    "gnome-terminal --tab -- browser-sync start -w --no-notify -s \"" + str(out_base) + "\" --host 127.0.0.1 --port 9000 --index \"plot_cov.html\"",
                    shell=True,
                check = True)
                active_report.set_open_status(True)
    #return open_report # TODO we don't need this anymore

def start_covermon():
    tab = "\t"
    nl = "\n"


    #####################
    # Start the monitor #
    #####################

    # Parse and check arguments
    if len(sys.argv) < 4:
        raise Exception(f"Missing arguments. The script must contain (1) samplesheet, (2) path to run directory, (3) path to reference or reference directory, (4) threshold for minimum coverage, (5) maximum depth displayed in plot. Optionally, a region file (6) can be specified if all samples use the same reference genome.")

    samplesheet = sys.argv[1]
    rundir = sys.argv[2]
    threshold = sys.argv[4]
    maxDepth = sys.argv[5]

    if ".fa" in sys.argv[3]:
        reference = sys.argv[3]
        one_ref = True
        if len(sys.argv) == 7:
            region_file = sys.argv[6]
            print("This is the region file:", sys.argv[6])
        else: 
            region_file = "NA"
    else:
        refdir = sys.argv[3]
        one_ref = False
        region_file = "NA"

    print(f"These are the parameters given:")
    print("This is the samplesheet: ", samplesheet)
    print("This is the run directory: ", rundir)
    if one_ref == True:
        print("This is the reference:", reference)
    else:
        print("This is the reference directory:", refdir)
    print("This is the threshold: ", threshold)
    print("This is the maxDepth in plot:", maxDepth)

    #########################
    # Parse the samplesheet #
    #########################

    df = parse_samplesheet(samplesheet)

    validate_samplesheet(df, one_ref)


    ###################
    # Validate rundir #
    ###################

    base_dir, fastq_pass_base = validate_rundir(rundir)

    ###########################
    # Create output directory #
    ###########################
    if one_ref == True:
        out_base = os.path.join(base_dir, "CoverMon_"+reference.split("/")[-1].split(".")[0])
    else:
        out_base = os.path.join(base_dir, "CoverMon") # out_base is the directory where the pipeline will write its output to.
    print("Creating output directory ", out_base,"...")
    print("Creating output directory ", out_base,"...")
    subprocess.run(["mkdir","-p", out_base])
    print()

    sample_sheet_out = f"{out_base}/sample_sheet_given.tsv"
    print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
    df.to_csv(sample_sheet_out, sep = "\t", index=False, na_rep='NA')
    print("✓")

    #########################
    # Create workflow table #
    #########################

    workflow_table = create_workflow_table(df, fastq_pass_base)

    ##################
    # Start CoverMon #
    ##################

    # Set an open_report state to stop opening multiple reports
    print("Setting open_report to false")
    open_report = False

    # Keep track of processed files to avoid starting from scratch if script is terminated
    if exists(f"{out_base}/processed_files.txt") and open_report == False:
        print("I have found processed files, opening report")
        processed_files_txt = open(f"{out_base}/processed_files.txt", mode ="r", newline=nl)
        processed_files = processed_files_txt.read().splitlines()
        processed_files_txt.close()
        # Starts browser-sync in a new terminal. This will not work on windows or macOS.
        subprocess.run("gnome-terminal --tab -- browser-sync start -w --no-notify -s \"" + out_base +"\" --host 127.0.0.1 --port 9000 --index \"plot_cov.html\"", shell=True)                
        open_report = True
    else:
        processed_files = []


    # When sequencing, we will check for new files every 60 seconds
    seconds_wait = 60

    # Set a sequencing state to recognise when to stop looking for new files
    still_sequencing = True

    while still_sequencing:
        # Scans for new files and updates the plot if any are found
        open_report = update_plot(workflow_table, reference, sample_sheet_out,open_report)

        # Continue the monitor as long as the sequence summary does not exist. Wait <seconds_wait> between scans.
        # TODO: rework this so we can rerun/run in parallel
        sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
        if len(sequencing_summary_file) == 0:
            print(f"  Still sequencing/basecalling; waiting {seconds_wait} seconds before next scan ...")
            time.sleep(seconds_wait)
        else:
            still_sequencing = False


    sequencing_summary_file = sequencing_summary_file[0]
    print("  The sequencing summary has been found. Run complete    ✓")

if __name__ == "__main__":
    start_covermon()

