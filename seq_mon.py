"""Monitor an ongoing sequencing run and visualize results."""

__author__ = "Tine Sneibjerg Ebsen, Kat Steinke"
__version__ = "0.3"

import argparse
import copy
import json
import logging
import os
import pathlib
import random
import subprocess
import sys
import time
import webbrowser

from typing import Optional, List

import livereload
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
        fastq_dir:          the directory containing sequencing data
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
        port:               the port where to serve the report html
    """
    def __init__(self, indir: pathlib.Path, sample_sheet: pathlib.Path, threshold: int, maxDepth: int,
                 reference: pathlib.Path, out_base: Optional[pathlib.Path] = None,
                 region_file: Optional[pathlib.Path] = None, port: Optional[int] = None) \
            -> None:
        """Initialize a CoverReport object.

        Args:
            port:
            indir:          the directory containing sequencing data
            sample_sheet:   path to sample sheet
            threshold:      minimum coverage required
            maxDepth:       maximum coverage reported
            reference:      reference file or directory
            out_base:       output base directory (optional)
            region_file:    path to bed file with regions in the reference
                            (optional, only allowed if reference is a single file)
            port:           the port where to serve the report html
                            (optional; if not given a random port between 50500 and 51000 will be assigned)


        Raises:
            ValueError  if the maximum reported coverage is lower than the minimum coverage required,
                        or if a region file is given together with a reference base directory
                        (=different references per sample)
        """
        # sense checks first
        if threshold > maxDepth:
            raise ValueError("Highest reported coverage must exceed minimum required coverage.")
        if reference.is_dir() and region_file is not None:
            raise ValueError("Cannot supply a region file when different references per sample are used "
                             "(base ref is a directory)")
        dummy, self.fastq_dir = validate_rundir(indir)
        self.reference = reference
        self.sample_sheet = sample_sheet
        # sample sheet validation happens here
        df = parse_samplesheet(self.sample_sheet)
        validate_samplesheet(df, self.reference.is_file())
        self.workflow_table = create_workflow_table(df, self.fastq_dir)
        self.threshold = threshold
        self.maxDepth = maxDepth
        if out_base is not None:
            self.out_base = pathlib.Path(out_base)
        else:
            covermon_dirname = "CoverMon"
            if self.reference.is_file():
                covermon_dirname = f"{covermon_dirname}_{self.reference.stem}"
            self.out_base = self.fastq_dir.parent / covermon_dirname
        self.is_open = False  # TODO: can we assume this?
        self.processed_files: List[pathlib.Path] = []
        self.region_file = region_file
        # set a port for output
        if port is not None:
            self.port = port
        else:
            self.port = random.randrange(50500, 51000)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, CoverReport):
            return ((self.fastq_dir == other.fastq_dir
                     and self.sample_sheet == other.sample_sheet
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
        return(f"CoverReport(fastq_dir={self.fastq_dir}, sample_sheet={self.sample_sheet}, "
               f"threshold={self.threshold}, maxDepth={self.maxDepth},"
               f" reference={self.reference}, region_file={self.region_file},"
               f" out_base={self.out_base}, processed_files={[str(processed) for processed in self.processed_files]},"
               f" is_open={self.is_open},"
               f" port={self.port})\n"
               f"Workflow table:\n{self.workflow_table.head().to_string()}")

    def set_open_status(self, open_status: bool) -> None:
        """Set the report's status to open (True) or closed (False)

        Args:
            open_status:    whether the report is currently opened in the browser (True) or not (False)
        """
        self.is_open = open_status

    def add_processed_file(self, processed: pathlib.Path) -> None:
        """Add a file to the record of processed files.

        Args:
            processed:  the path to the processed file

        """
        self.processed_files.append(processed)

    def add_processed_files(self, processed_files: List[pathlib.Path]) -> None:
        """Add a number of files to the record of processed files.

        Args:
            processed_files:    the processed files to add
        """
        self.processed_files.extend(processed_files)

    def save_settings(self) -> None:
        """Save the report's properties."""
        settings = copy.deepcopy(vars(self))
        path_settings = {"reference", "fastq_dir", "sample_sheet", "out_base", "region_file"}
        settings["processed_files"] = [str(processed) for processed in settings["processed_files"]]
        settings = {key: str(value.resolve()) if key in path_settings and value is not None else value
                    for key, value in settings.items()}
        settings.pop("workflow_table")
        settings.pop("is_open")  # might be incorrectly set if the run is interrupted
        out_file = self.out_base / "settings.json"
        with open(out_file, "w", encoding = "utf-8") as settings_file:
            json.dump(settings, settings_file, indent = 4)


# TODO - leverage this for a convenient "--resume" flag?
def load_report_from_settings(settings_file: pathlib.Path) -> CoverReport:
    """Load a CoverReport from a settings.json file.

    Args:
        settings_file: the file from which to load settings

    Returns:
        The report defined by the archived settings
    """
    with open(settings_file, "r", encoding="utf-8") as settings:
        existing_vars = json.load(settings)
        existing_vars["indir"] = existing_vars.pop("fastq_dir")
        processed = [pathlib.Path(processed_file).resolve() for processed_file in existing_vars.pop("processed_files")]
    path_settings = {"reference", "indir", "sample_sheet", "out_base", "region_file"}
    # we might have some ugly Nones
    existing_vars = {key: pathlib.Path(value).resolve() if key in path_settings and value and value != "None"
                    else value
                    for key, value in existing_vars.items()}
    # clean the Nones
    existing_vars = {key: None if value == "None" else value for key, value in existing_vars.items()}
    report = CoverReport(**existing_vars)
    report.add_processed_files(processed)
    return report

def parse_samplesheet(samplesheet: pathlib.Path) -> pd.DataFrame:
    """Clean an Excel sample sheet.

    Args:
        samplesheet: the path to the sample sheet

    Returns:
        The sample sheet with column names cleaned
        (lowercase, leading/trailing whitespace stripped, remainder changed to underscores),
        whitespace removed in barcode names,
        and any rows without sample ID removed.
    Raises:
        ValueError   if the sheet is not in the correct format
    """
    samplesheet_extension = samplesheet.suffix
    logger.info(f"Reading {samplesheet_extension}-type sample sheet \"{samplesheet}\"")

    if samplesheet_extension not in {".xlsx", ".xls"}:
        raise ValueError("The spreadsheet must be excel formatted (.xlsx or .xls)")
    dtype_used = str if samplesheet_extension == ".xlsx" else None
    df = pd.read_excel(samplesheet, dtype = dtype_used)

    # Clean up the spreadsheet
    logger.info("Cleaning sample sheet ...                              ")
    df.columns = map(str.lower, df.columns)  # Lowercase
    df.columns = map(str.strip, df.columns)  # Remove edge-spaces
    df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns)  # Replace spaces with underscore
    # Because we are later going to join using this column, it is necessary to strip it for spaces.
    df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ",
                                                                                      "")))
    df = df.dropna(subset=["sample_id"])  # remove rows not containing a barcode
    logger.info("Sample sheet cleaned ✓")
    logger.debug(df.to_string())
    return df

def validate_samplesheet(samplesheet: pd.DataFrame, ref_is_file: bool) -> None:  # TODO: might be nicer with a bool output
    """Check that the sample sheet contains the correct barcodes and no duplicates,
    and all columns are present

    Args:
        samplesheet:    the processed sample sheet
        ref_is_file:    whether the reference is a single file (True) or a directory (False)

    Raises:
        Exception   if a required column (barcode, reference, sample ID) is missing in a setup with a reference dir,
                    if barcodes aren't correctly formatted,
                    or if a barcode is duplicated
    """
    # Check that the samplesheet contains the reference column if a refdir is given
    required_cols = ["barcode", "sample_id"]
    if not ref_is_file:
        required_cols.append("reference")
    logger.info("Checking that the necessary columns exist ...")
    for required_col in required_cols:  # TODO optimize with set ops?
        if not required_col in samplesheet.columns:
            raise KeyError("The sample sheet is missing a necessary column. "
                            f"The sample sheet must contain the column {required_col}, "
                            f"but it only contains {sorted(samplesheet.columns.tolist())}")
    logger.info("All necessary columns found ✓")
    # Check that the barcodes look correct
    acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)] + [f"RB{i:02d}" for i in range(1,97)]

    logger.info("Checking that the barcodes are correctly formatted ...")

    for barcode in samplesheet["barcode"]:  # TODO catch all broken barcodes at once
        if not barcode in acceptable_barcodes:
            raise ValueError(f"The given barcode {barcode} is not an acceptable barcode. "
                            f"Here is a list of acceptable barcodes for inspiration:\n{' '.join(acceptable_barcodes)}")
    logger.info("Barcodes are correct ✓")


    logger.info("Checking that the barcodes are unique ...")
    if not len(samplesheet["barcode"]) == len(set(samplesheet["barcode"])):
        bc_counts = pd.DataFrame(samplesheet['barcode'].value_counts())
        bc_counts.columns = ["count"]
        bc_counts = bc_counts[bc_counts["count"] > 1]
        raise ValueError(f"\nOne or more barcodes are duplicated. Each barcode may only be used once:\n{bc_counts}")
    logger.info("All barcodes are unique ✓")

    logger.info("These are the samples from the samplesheet you have given:\n"
                f"{samplesheet.to_string()}\n"
                "//")

def validate_rundir(rundir: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    """Check whether the rundir exists, and find the fastq_pass directory and its parent dir, possibly waiting for
    the fastq_pass directory to be created.
    Parts adapted from SnakeAmp.

    Args:
        rundir: the directory containing sequence data *somewhere*

    Returns:
        The parent directory of the fastq_pass directory and the fastq_pass directory
    Raises
        FileNotFoundError:  if the rundir does not exist,
                            or if the fastq_pass directory was not found after the end of the waiting time,
        ValueError:         if multiple fastq_pass directories were found
    """
    logger.info("Checking that the rundir exists ...                    ")
    if not rundir.exists():
        raise FileNotFoundError("The rundir does not exist.")
    print("✓")
    existing_path = rundir.parts
    logger.info("Looking for MinKNOW-characteristic output:") #, end = "", flush = True)
    # Wait for the rundir to occur in the specified path.
    # If it doesn't occur after a specified waiting time, then stop the p
    # TODO: convenience waiting function here?
    for i in range(10):
        logger.info("  Looking ... ")
        try:
            # if the fastq_pass directory already is somewhere in the dirs given, use this
            fastq_pass_parts = existing_path[:existing_path.index("fastq_pass") + 1]
            # parts contains the initial "/" - resolve the path to clean this up
            fastq_pass_base = pathlib.Path("/".join(fastq_pass_parts)).resolve()
            # do we end with something that actually exists?
        except ValueError:
            logger.info(f"Searching for fastq_pass folder in {rundir}...")
            fastq_pass_bases = list(rundir.glob("**/fastq_pass"))
            if fastq_pass_bases:
                if len(fastq_pass_bases) > 1:
                    bad_paths = '\n '.join([str(base) for base in fastq_pass_bases])
                    error_msg = ("There seems to be more than one fastq_pass sub-directory beneath the given rundir."
                            " These paths were found:\n"
                            f"{bad_paths}"
                            "\nPlease specify a more specific rundir.")
                    raise ValueError(error_msg)

                fastq_pass_base = fastq_pass_bases[0]
                logger.info("Found                                    ✓")
            else:
                if i < 10:
                    logger.info("nothing found yet, waiting 10 secs ...")
                    time.sleep(10)  # Wait 10 seconds.
                else:
                    raise FileNotFoundError("nothing found after 10 tries. Aborting.")
    logger.info(f"Found the following fastq_pass base which will be given to CoverMon: \n  {fastq_pass_base}\n")


    # base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
    base_dir = fastq_pass_base.parent
    logger.info(f"This is the batch base directory:\n  {base_dir}")

    return base_dir, fastq_pass_base


def create_workflow_table(samplesheet: pd.DataFrame, fastq_pass_dir: pathlib.Path) -> pd.DataFrame:
    """Record locations of barcode directories for each barcode.

    Args:
        samplesheet:    the sample sheet containing barcodes, sample IDs and optionally references
        fastq_pass_dir: the path to the fastq_pass directory

    Returns:
        The sample sheet with the locations of the barcode fastq directories added for each sample

    Raises:
        Exception   if the barcodes' format is invalid
    """
    disk_barcodes_list  = sorted(list(fastq_pass_dir.glob("barcode*"))) # Find all fastq_pass/barcode* directories
    disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})


    disk_barcodes_df["barcode_basename"] = disk_barcodes_df["barcode_path"].apply(lambda barcode_path: barcode_path.name)
    if "RB" in samplesheet["barcode"][0]:
        disk_barcodes_df["barcode"] = (disk_barcodes_df["barcode_basename"]
                                       .apply(lambda barcode: "RB" + barcode[-2:]))
    elif "NB" in samplesheet["barcode"][0]:
        disk_barcodes_df["barcode"] = (disk_barcodes_df["barcode_basename"]
                                       .apply(lambda barcode: "NB" + barcode[-2:]))
    else:
        raise ValueError("Barcodes in samplesheet are not acceptable")
    # all path operations done, back to strings
    disk_barcodes_df["barcode_path"] = disk_barcodes_df["barcode_path"].apply(str)

    # ensure consistent column format
    sample_cols = samplesheet.columns.tolist()
    barcode_cols = ["barcode_path", "barcode_basename"]
    out_cols = [*sample_cols, *barcode_cols]

    # the workflow_table is the table that contains the records where the barcode could be found on the disk.
    # left join (merge) the present barcodes onto the df table.
    workflow_table = disk_barcodes_df.merge(samplesheet, how='left', on='barcode')
    workflow_table = workflow_table.dropna(subset = ["sample_id"]).reindex(columns = out_cols)

    logger.info(f"Continuing with the following barcodes:\n{workflow_table.to_string()}\n//")
    return workflow_table

# we probably have to keep this until we've established that nothing else uses it
def write_to_processed(to_write: str, out_dir: pathlib.Path) -> None:
    """Append the given string to "processed_files.txt" in the specified outdir

    Args:
        to_write: the string to write to the file
        out_dir:  the base directory in which to write to the file

    """
    with open(out_dir / "processed_files.txt","a", encoding = "utf-8") as processed_files_txt:
        processed_files_txt.write(f"{to_write}\n")

def create_bam(fastq: pathlib.Path, out_base: pathlib.Path, out_bam: pathlib.Path, reference: pathlib.Path)\
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


def append_bam(fastq: pathlib.Path, out_base: pathlib.Path, out_bam: pathlib.Path, reference: pathlib.Path)\
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
    sort_cmd = f'samtools sort -O bam -o {out_base}/sorted.bam {tmpsam}'
    logger.debug(sort_cmd)
    subprocess.run(sort_cmd.split(), check = True)
    # Merges the sorted files
    merge_cmd = f'samtools merge -f -o {tmpbam} {out_base}/sorted.bam {out_bam}'
    logger.debug(merge_cmd)
    subprocess.run(merge_cmd.split(), check = True)
    # Sets the new bam to the barcode bam and marks file as processed
    subprocess.run(['mv', tmpbam, out_bam], check = True)


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

def get_r_cmd(active_report: CoverReport) -> List[str]:
    """Get the command for plotting coverage.

    Args:
        active_report:  the CoverReport to plot

    Returns:
        The command for plotting coverage from the given report
    """
    plot_cmd = ["Rscript",
                     str(pathlib.Path(__file__).parent.resolve() / "scripts" / "run_plot.R"),
                     str(active_report.out_base.resolve()),
                     str(active_report.sample_sheet.resolve()),
                     str(active_report.threshold),
                     str(active_report.maxDepth)]
    if active_report.region_file:
        plot_cmd.extend(["--region_file", str(active_report.region_file)])
    return plot_cmd


# TODO remember to save the files here
def update_plot(active_report: CoverReport) -> None:
    """Create or update the coverage plot for a given report and show plot in a browser window

    Args:
        active_report:  the CoverReport to plot
    """
    # Define temporary mapping files which will be used for merging new mapping of new output files with existing
    workflow_table = active_report.workflow_table.copy(deep=True)
    out_base = active_report.out_base
    logger.info("Scanning for new fastq files...")
    workflow_table["bam_out"] = workflow_table["barcode"].apply(lambda barcode: out_base / f"{barcode}.bam")
    workflow_table["depth"] = workflow_table["barcode"].apply(lambda barcode: out_base / f"{barcode}.depth")
    workflow_table["reference"] = workflow_table["reference"].apply(lambda ref: active_report.reference / ref) \
                                  if active_report.reference.is_dir() else active_report.reference
    for index, row in workflow_table.iterrows():
        bam_out = row["bam_out"]
        depth = row["depth"]
        barcode_path = pathlib.Path(row['barcode_path'])
        reference = row["reference"]
        unprocessed = [barcode_file.resolve() for barcode_file in barcode_path.iterdir() if
                       (barcode_file.resolve() not in active_report.processed_files and barcode_file.is_file()
                        and {".fastq", ".fq"}.intersection(barcode_file.suffixes))]

        for unprocessed_file in unprocessed:
            # Check if we have an existing read mapping to append to.
            # If not, creates the first one and continues the loop without merging.
            if not bam_out.exists():
                create_bam(unprocessed_file, out_base, bam_out, reference)
            else:
                # Maps new reads to reference
                append_bam(unprocessed_file, out_base, bam_out, reference)
            # Index the new bam and calculate depth.
            get_depth(bam_out, depth)
            active_report.add_processed_file(unprocessed_file)
            write_to_processed(str(unprocessed_file), out_base)  # TODO: replace with saving settings?
            active_report.save_settings()

            logger.info(f"Number of processed files: {len(active_report.processed_files)}")

            # Now create the monitoring html if needed
            plot_cov_cmd3 = get_r_cmd(active_report)
            logger.debug(" ".join(plot_cov_cmd3))
            subprocess.run(" ".join(plot_cov_cmd3), shell=True, check = True)
            logger.info("Updated plot")  # TODO remove this, add to run_plot
            if not active_report.is_open:
                logger.info("Opening report")
                if os.fork():
                    sys.exit(0)
                active_report.set_open_status(True)
                server = livereload.Server()
                server.watch(out_base / "processed_files.txt", " ".join(plot_cov_cmd3))
                server.serve(port = active_report.port, open_url_delay = 1,
                             default_filename = out_base / "plot_cov.html")

def start_covermon(start_args) -> None:
    """Start monitoring with the supplied arguments

    Args:
        start_args: the arguments supplied by the user


    """
    #####################
    # Start the monitor #
    #####################
    # Parse and check arguments
    parser = argparse.ArgumentParser(description="Start CoverMon")
    parser.add_argument("samplesheet", help = "Path to sample sheet for run")
    parser.add_argument("rundir", help = "Path to the directory containing sequencing data")
    parser.add_argument("reference",
                        help = "Path to the reference fasta file or directory with reference files to use")
    parser.add_argument("threshold", help = "Minimum coverage to pass QC")
    parser.add_argument("maxdepth", help = "Maximum depth to plot")
    parser.add_argument("--region_file",
                        help = "Path to the region file to be used (optional, only for a single reference file)",
                        default = None)
    parser.add_argument("--outdir",
                        help = "Path to the output directory (default: [rundir]/CoverMon if reference is a dir,"
                               "[rundir]/CoverMon_[ref_filename] otherwise)", default = None)
    args = parser.parse_args(start_args)
    rundir = pathlib.Path(args.rundir)
    region_file = pathlib.Path(args.region_file) if args.region_file else None
    out_dir = pathlib.Path(args.outdir) if args.outdir else None
    # initialize the report
    # Set an open_report state to stop opening multiple reports
    logger.debug("Initializing report as closed")
    report = CoverReport(rundir, pathlib.Path(args.samplesheet), int(args.threshold), int(args.maxdepth),
                         pathlib.Path(args.reference), out_dir, region_file)

    if report.region_file is not None:
        logger.info(f"This is the region file: {report.region_file}")
    logger.info("These are the parameters given:")
    logger.info(f"This is the samplesheet: {report.sample_sheet}")
    logger.info(f"This is the run directory: {rundir}")
    logger.info(f"This is the reference{' directory' if report.reference.is_dir() else ''}: {report.reference}")
    logger.info(f"This is the threshold: {report.threshold}")
    logger.info(f"This is the maxDepth in plot: {report.maxDepth}")

    ###########################
    # Create output directory #
    ###########################
    logger.info(f"Creating output directory {report.out_base}...")
    report.out_base.mkdir(parents=True, exist_ok=True)

    # back up sample sheet - TODO: doing this twice - should we give the CoverReport yet another attribute?
    df = parse_samplesheet(report.sample_sheet)
    sample_sheet_out = report.out_base / "sample_sheet_given.tsv"
    logger.info("Backing up the original sample sheet...")
    df.to_csv(sample_sheet_out, sep = "\t", index=False, na_rep='NA')

    ##################
    # Start CoverMon #
    ##################
    # Keep track of processed files to avoid starting from scratch if script is terminated
    if ((report.out_base / "processed_files.txt").exists() and (report.out_base / "settings.json").exists()
            and not report.is_open):  # TODO: it'll always be closed?
        logger.info("I have found processed files, opening existing report")
        old_report = load_report_from_settings(report.out_base / "settings.json")

        # TODO: use the processed files from the old report settings?
        with open(report.out_base / "processed_files.txt", "r", encoding = "utf-8") as processed_files_txt:
            processed_files = [pathlib.Path(processed).resolve()
                               for processed in processed_files_txt.read().splitlines()]
            report.add_processed_files(processed_files)
        if old_report != report:
            error_msg = ("Attempting to continue an interrupted run with different settings."
                             " Please choose a different"
                             " output directory to run analysis with the new settings.\n"
                             "Old settings:\n"
                             f"{old_report}\n"
                             "New settings:\n"
                             f"{report}")
            raise ValueError(error_msg)
        logger.info("Current settings match saved settings")
        # set the old report's port settings on the new one
        report.port = old_report.port

        plot_cov_cmd3 = get_r_cmd(report)
        # subprocess.run(" ".join(plot_cov_cmd3), shell=True, check = True)
        logger.info("Opening report")
        if os.fork():
            sys.exit(0)
        report.set_open_status(True)
        server = livereload.Server()
        server.watch(report.out_base / "processed_files.txt", " ".join(plot_cov_cmd3))
        server.serve(port=report.port, open_url_delay=1,
                     default_filename=report.out_base / "plot_cov.html")


    # When sequencing, we will check for new files every 60 seconds
    seconds_wait = 60

    # Set a sequencing state to recognise when to stop looking for new files
    still_sequencing = True

    while still_sequencing:
        # Scans for new files and updates the plot if any are found
        update_plot(report)

        # Continue the monitor as long as the sequence summary does not exist. Wait <seconds_wait> between scans.
        sequencing_summary_file = list(report.fastq_dir.parent.glob("sequencing_summary_*.txt"))
        if len(sequencing_summary_file) == 0:
            print(f"  Still sequencing/basecalling; waiting {seconds_wait} seconds before next scan ...")
            time.sleep(seconds_wait)
        else:
            still_sequencing = False


    logger.info("  The sequencing summary has been found. Run complete    ✓")
    # as the live coverage won't be available anymore when the script closes we now open the file
    webbrowser.open_new_tab(f"file://{report.out_base / 'plot_cov.html'}")

if __name__ == "__main__":
    start_covermon(sys.argv[1:])

