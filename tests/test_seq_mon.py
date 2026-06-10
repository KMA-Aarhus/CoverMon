import pathlib
import re
import unittest

import numpy as np
import pandas as pd
import pytest

import seq_mon

class TestWriteToProcessed(unittest.TestCase):
    outdir = pathlib.Path(__file__).parent / "data" / "out_dir"
    outfile = outdir / "processed_files.txt"

    def tearDown(self):
        with open(self.outfile, "w", encoding = "utf-8") as f:
            f.write("")

    def test_write_to_file(self):
        """Write to file as specified, append as needed."""
        # file should be empty before writing
        with open(self.outfile, "r", encoding = "utf-8") as f:
            content = f.readlines()
        assert not content
        # write one string
        to_write = "teststring"
        seq_mon.write_to_processed(to_write, str(self.outdir))
        # check we only have the one
        with open(self.outfile, "r", encoding = "utf-8") as f:
            content = f.readlines()
        assert len(content) == 1
        assert content[0] == "teststring\n"
        # write another and check we've appended
        write_more = "test_2"
        seq_mon.write_to_processed(write_more, str(self.outdir))
        with open(self.outfile, "r", encoding = "utf-8") as f:
            content = f.readlines()
        assert len(content) == 2
        assert content[0] == "teststring\n"
        assert content[1] == "test_2\n"


class TestParseSamplesheet(unittest.TestCase):
    def test_fail_bad_format(self):
        """Raise an exception if the sample sheet isn't in a supported format"""
        test_sheet = pathlib.Path(__file__).parent / "data" / "samplesheet.csv"
        error_msg = "The spreadsheet must be excel formatted (.xlsx or .xls)"
        with pytest.raises(Exception, match = re.escape(error_msg)):
            seq_mon.parse_samplesheet(str(test_sheet))

    def test_success_xls(self):
        """Read an xls file and process the contents."""
        test_sheet = pathlib.Path(__file__).parent / "data" / "samplesheet.xls"
        expected_sheet = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                            "barcode": ["NB01", "NB02", "NB03", "NB04"],
                                            "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                          "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                            "ct": [10, 20, 30, np.nan],
                                            "other_columns1": ["A", "B", "C", "D"],
                                            "other_columns2": [1, 2, 3, 4]
                                            })
        test_result = seq_mon.parse_samplesheet(str(test_sheet))
        pd.testing.assert_frame_equal(expected_sheet, test_result)

    def test_success_xlsx(self):
        """Read an xlsx file and process the contents."""
        test_sheet = pathlib.Path(__file__).parent / "data" / "samplesheet.xlsx"
        expected_sheet = pd.DataFrame(data = {"sample_id": ["test1", "test2", "test3", "test4"],
                                              "barcode": ["NB01", "NB02", "NB03", "NB04"],
                                              "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                            "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                              "ct": ["10", "20", "30", np.nan],
                                              "other_columns1": ["A", "B", "C", "D"],
                                              "other_columns2": ["1", "2", "3", "4"]
                                              })
        expected_sheet = expected_sheet.astype({"ct": "object"})
        print(expected_sheet["ct"])
        test_result = seq_mon.parse_samplesheet(str(test_sheet))
        pd.testing.assert_frame_equal(expected_sheet, test_result, check_dtype = False)


class TestValidateSamplesheet(unittest.TestCase):
    def test_fail_missing_column_refdir(self):
        """Raise an exception if the sample sheet is missing a required column when reference is given as a dir."""
        test_sheet = pd.DataFrame(data = {"sample_id": ["test1", "test2", "test3", "test4"],
                                              "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                            "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                              "ct": ["10", "20", "30", np.nan],
                                              "other_columns1": ["A", "B", "C", "D"],
                                              "other_columns2": ["1", "2", "3", "4"]
                                              })
        error_msg = ("The sample sheet is missing a necessary column. The sample sheet must contain the column barcode,"
                     f" but it only contains ['ct', 'other_columns1', 'other_columns2', 'reference', 'sample_id']")
        with pytest.raises(Exception, match=re.escape(error_msg)):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = False)

    def test_fail_missing_column(self):
        """Raise an exception if the sample sheet is missing a required column
        (and reference is not a required column)."""
        pass  # TODO second pass when we start tweaking behavior

    def test_fail_bad_barcode(self):
        """Raise an exception if a barcode is malformed."""
        # TODO: give cleaner examples
        test_sheet = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                            "barcode": ["NB01", "LB02", "NB3", "NB04"],
                                            "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                          "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                            "ct": ["10", "20", "30", np.nan],
                                            "other_columns1": ["A", "B", "C", "D"],
                                            "other_columns2": ["1", "2", "3", "4"]
                                            })
        acceptable_barcodes = [f"NB{i:02d}" for i in range(1, 97)] + [f"RB{i:02d}" for i in range(1, 97)]
        error_msg = ("The given barcode LB02 is not an acceptable barcode. "
                     "Here is a list of acceptable barcodes for inspiration:\n"
                     f"{' '.join(acceptable_barcodes)}")
        with pytest.raises(Exception, match=re.escape(error_msg)):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = True)

    def test_fail_duplicated_barcode(self):
        """Raise an exception if a barcode is duplicated."""
        test_sheet = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                            "barcode": ["NB01", "NB01", "NB03", "NB03"],
                                            "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                          "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                            "ct": ["10", "20", "30", np.nan],
                                            "other_columns1": ["A", "B", "C", "D"],
                                            "other_columns2": ["1", "2", "3", "4"]
                                            })
        duplicates = pd.DataFrame(data = {'count': [2, 2]},
                                  index = pd.Series(["NB01", "NB03"], name = "barcode"))
        error_msg = ("\nOne or more barcodes are duplicated. Each barcode may only be used once:\n"
                     f"{duplicates}")
        with pytest.raises(Exception, match=re.escape(error_msg)):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = True)

    def test_success(self):
        """Don't complain if everything looks good."""
        pass # TODO: second pass, return a bool and log things


class TestValidateRundir(unittest.TestCase):  # TODO: combine with find_fastq_pass_parent on a second pass
    def test_fail_missing_dir(self):
        """Raise an exception if the run directory doesn't exist."""
        test_dir = pathlib.Path(__file__).parent / "data" / "rundir" / "no_such_dir"
        error_msg = "The rundir does not exist."
        with pytest.raises(Exception, match=re.escape(error_msg)):
            seq_mon.validate_rundir(str(test_dir))

    def test_fail_timeout(self):
        """Raise an exception if nothing is found after a timeout."""
        pass  # TODO: combine this with a convenience waiting function on a second pass so we can test it without having to wait for ages

    def test_fail_too_many_pass_dirs(self):
        """Raise an exception if there are too many fastq_pass directories."""
        test_dir = pathlib.Path(__file__).parent / "data" / "rundir"
        fastq_pass_dirs = [str(test_dir / "test1" / "fastq_pass"), str(test_dir / "test2" / "fastq_pass")]
        error_msg = ("There seems to be more than one fastq_pass sub-directory beneath the given rundir."
                     " These paths were found:\n"
                     f" {str('\n ').join(fastq_pass_dirs)}\n"
                     "Please specify a more specific rundir.")
        with pytest.raises(Exception, match=re.escape(error_msg)):
            seq_mon.validate_rundir(str(test_dir))

    def test_success(self):
        """Find base directory and the fastq_pass directory it contains."""
        test_dir = pathlib.Path(__file__).parent / "data" / "rundir" / "test1"
        expected_base = str(test_dir)
        expected_fastq = str(test_dir / "fastq_pass")
        test_base, test_fastq = seq_mon.validate_rundir(str(test_dir))
        assert test_base == expected_base
        assert test_fastq == expected_fastq



class TestCreateWorkflowTable(unittest.TestCase):
    test_dir = pathlib.Path(__file__).parent / "data" / "workflow_table" / "rundir" / "test_dir" / "fastq_pass"

    def test_fail_bad_barcodes(self):
        """Raise an exception if the barcodes are malformed."""
        test_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                     "barcode": ["LB01", "NB02", "NB03", "NB04"],
                                     "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                     "ct": ["10", "20", "30", np.nan],
                                     "other_columns1": ["A", "B", "C", "D"],
                                     "other_columns2": ["1", "2", "3", "4"]
                                     })
        error_msg = "Barcodes in samplesheet are not acceptable"
        with pytest.raises(Exception, match=re.escape(error_msg)):
            seq_mon.create_workflow_table(test_df, str(self.test_dir))

    # TODO handling of workflow tables that end up empty after removing barcodes without sample ID?

    def test_success_rb(self):
        """Create the workflow table (easy case, no barcodes without sample ID) with rapid barcodes."""
        test_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                     "barcode": ["RB01", "RB02", "RB03", "RB04"],
                                     "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                     "ct": ["10", "20", "30", np.nan],
                                     "other_columns1": ["A", "B", "C", "D"],
                                     "other_columns2": ["1", "2", "3", "4"]
                                     })
        expected_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                     "barcode": ["RB01", "RB02", "RB03", "RB04"],
                                     "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                     "ct": ["10", "20", "30", np.nan],
                                     "other_columns1": ["A", "B", "C", "D"],
                                     "other_columns2": ["1", "2", "3", "4"],
                                         "barcode_path": [str(self.test_dir / "barcode01"),
                                                          str(self.test_dir / "barcode02"),
                                                          str(self.test_dir / "barcode03"),
                                                          str(self.test_dir / "barcode04")],
                                         "barcode_basename": ["barcode01", "barcode02", "barcode03", "barcode04"]

                                     })
        test_result = seq_mon.create_workflow_table(test_df, str(self.test_dir))
        print(test_result.columns)
        pd.testing.assert_frame_equal(expected_df, test_result)

    def test_success_nb(self):
        """Create the workflow table with native barcodes."""
        test_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                     "barcode": ["RB01", "RB02", "RB03", "RB04"],
                                     "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                     "ct": ["10", "20", "30", np.nan],
                                     "other_columns1": ["A", "B", "C", "D"],
                                     "other_columns2": ["1", "2", "3", "4"]
                                     })
        expected_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                         "barcode": ["RB01", "RB02", "RB03", "RB04"],
                                         "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                       "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                         "ct": ["10", "20", "30", np.nan],
                                         "other_columns1": ["A", "B", "C", "D"],
                                         "other_columns2": ["1", "2", "3", "4"],
                                         "barcode_path": [str(self.test_dir / "barcode01"),
                                                          str(self.test_dir / "barcode02"),
                                                          str(self.test_dir / "barcode03"),
                                                          str(self.test_dir / "barcode04")],
                                         "barcode_basename": ["barcode01", "barcode02", "barcode03", "barcode04"]

                                         })
        test_result = seq_mon.create_workflow_table(test_df, str(self.test_dir))
        pd.testing.assert_frame_equal(expected_df, test_result)


    def test_success_drop_nas(self):
        """Create the workflow table, dropping any barcodes without sample IDs."""
        test_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4", pd.NA],
                                     "barcode": ["RB01", "RB02", "RB03", "RB04", "RB05"],
                                     "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa"],
                                     "ct": ["10", "20", "30", np.nan, np.nan],
                                     "other_columns1": ["A", "B", "C", "D", "E"],
                                     "other_columns2": ["1", "2", "3", "4", "5"]
                                     })
        expected_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                         "barcode": ["RB01", "RB02", "RB03", "RB04"],
                                         "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                       "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                         "ct": ["10", "20", "30", np.nan],
                                         "other_columns1": ["A", "B", "C", "D"],
                                         "other_columns2": ["1", "2", "3", "4"],
                                         "barcode_path": [str(self.test_dir / "barcode01"),
                                                          str(self.test_dir / "barcode02"),
                                                          str(self.test_dir / "barcode03"),
                                                          str(self.test_dir / "barcode04")],
                                         "barcode_basename": ["barcode01", "barcode02", "barcode03", "barcode04"]
                                         })
        test_result = seq_mon.create_workflow_table(test_df, str(self.test_dir))
        pd.testing.assert_frame_equal(expected_df, test_result)

