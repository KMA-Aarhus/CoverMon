import json
import logging
import pathlib
import re
import shutil
import unittest

from collections import namedtuple
from unittest import mock
from unittest.mock import call

import numpy as np
import pandas as pd
import pytest

import seq_mon

class TestCoverReport(unittest.TestCase):
    test_dir = pathlib.Path(__file__).parent / "data" / "workflow_table" / "rundir" / "test_dir" / "fastq_pass"
    test_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                         "barcode": ["NB01", "NB02", "NB03", "NB04"],
                                         "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                       "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                         "ct": [10, 20, 30, np.nan],
                                         "other_columns1": ["A", "B", "C", "D"],
                                         "other_columns2": [1, 2, 3, 4],
                                         "barcode_path": [str(test_dir / "barcode01"),
                                                          str(test_dir / "barcode02"),
                                                          str(test_dir / "barcode03"),
                                                          str(test_dir / "barcode04")],
                                         "barcode_basename": ["barcode01", "barcode02", "barcode03", "barcode04"]
                                         })
    test_outdir = pathlib.Path(__file__).parent / "data" / "out_dir"
    test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir"
    test_samplesheet = pathlib.Path(__file__).parent / "data"/"samplesheet.xls"

    # TODO check table generation

    def tearDown(self):
        settings = self.test_outdir / "settings.json"
        if settings.exists():
            settings.unlink()

    def test_sense_check_thresholds(self):
        """Fail if maximum coverage reported is below minimum coverage."""
        error_msg = "Highest reported coverage must exceed minimum required coverage."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=100, maxDepth=10,
                                reference=self.test_ref, out_base=self.test_outdir)

    def test_fail_region_for_refdir(self):
        """Don't allow using region files if reference is a directory (-> multiple ref sequences possible)"""
        error_msg = "Cannot supply a region file when different references per sample are used (base ref is a directory)"
        test_region = pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_region.bed"
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                reference=self.test_ref, out_base=self.test_outdir, region_file=test_region)

    def test_initialize_report(self):
        """Initialize a CoverReport."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref)
        pd.testing.assert_frame_equal(test_report.workflow_table, self.test_df, check_dtype = False)
        assert 10 == test_report.threshold
        assert 100 == test_report.maxDepth
        assert self.test_dir.parent / "CoverMon" == test_report.out_base
        assert self.test_ref == test_report.reference
        assert len(test_report.processed_files) == 0
        assert not test_report.is_open
        assert self.test_samplesheet == test_report.sample_sheet

        test_one_ref = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10,
                                           maxDepth=100, reference=pathlib.Path(self.test_ref) / "test_ref.fa")
        assert self.test_dir.parent / "CoverMon_test_ref" == test_one_ref.out_base

    def test_find_fastq(self):
        """Initialize a CoverReport and find the fastq_pass directory."""
        test_dir = pathlib.Path(__file__).parent / "data" / "workflow_table" / "rundir"
        test_report = seq_mon.CoverReport(test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref)
        pd.testing.assert_frame_equal(test_report.workflow_table, self.test_df, check_dtype = False)
        assert 10 == test_report.threshold
        assert 100 == test_report.maxDepth
        assert self.test_ref == test_report.reference
        assert len(test_report.processed_files) == 0
        assert not test_report.is_open
        assert self.test_samplesheet == test_report.sample_sheet
        assert (test_dir / "test_dir" / "CoverMon") == test_report.out_base


    def test_set_outdir(self):
        """Initialize a CoverReport and set the output directory."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        pd.testing.assert_frame_equal(test_report.workflow_table, self.test_df, check_dtype = False)
        assert 10 == test_report.threshold
        assert 100 == test_report.maxDepth
        assert self.test_outdir == test_report.out_base
        assert self.test_ref == test_report.reference
        assert len(test_report.processed_files) == 0
        assert not test_report.is_open
        assert self.test_samplesheet == test_report.sample_sheet


    def test_initialize_single_ref(self):
        """Initialize a CoverReport with a single reference file."""
        test_df = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                     "barcode": ["NB01", "NB02", "NB03", "NB04"],
                                     "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                   "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                     "ct": [10, 20, 30, np.nan],
                                         "other_columns1": ["A", "B", "C", "D"],
                                         "other_columns2": [1, 2, 3, 4],
                                     "barcode_path": [str(self.test_dir / "barcode01"),
                                                      str(self.test_dir / "barcode02"),
                                                      str(self.test_dir / "barcode03"),
                                                      str(self.test_dir / "barcode04")],
                                     "barcode_basename": ["barcode01", "barcode02", "barcode03", "barcode04"]
                                     })
        test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_ref.fa"
        test_region = pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_region.bed"
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=test_ref, out_base=self.test_outdir, region_file=test_region)
        pd.testing.assert_frame_equal(test_report.workflow_table, test_df, check_dtype = False)
        assert 10 == test_report.threshold
        assert 100 == test_report.maxDepth
        assert self.test_dir == test_report.fastq_dir
        assert self.test_outdir == test_report.out_base
        assert test_ref == test_report.reference
        assert len(test_report.processed_files) == 0
        assert not test_report.is_open
        assert self.test_samplesheet == test_report.sample_sheet
        assert test_region == test_report.region_file

    def test_set_report_open_status(self):
        """Set the CoverReport to open/closed."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_report.set_open_status(True)
        assert test_report.is_open
        test_report.set_open_status(False)
        assert not test_report.is_open

    def test_update_processed_files(self):
        """Update the processed files."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_2 = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                     reference=self.test_ref, out_base=self.test_outdir)
        assert len(test_report.processed_files) == 0
        assert len(test_2.processed_files) == 0
        test_report.add_processed_file(pathlib.Path("path/to/processed.fq"))
        assert len(test_report.processed_files) == 1
        assert test_report.processed_files[0] == pathlib.Path("path/to/processed.fq")
        assert len(test_2.processed_files) == 0

    def test_add_multi_files(self):
        """Update the processed files with a list of files."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_2 = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                     reference=self.test_ref, out_base=self.test_outdir)
        assert len(test_report.processed_files) == 0
        assert len(test_2.processed_files) == 0
        to_add = [pathlib.Path("path/to/processed.fq"), pathlib.Path("path/to/processed_2.fq")]
        test_report.add_processed_files(to_add)
        assert len(test_report.processed_files) == 2
        assert test_report.processed_files == to_add
        assert len(test_2.processed_files) == 0

    def test_equality(self):
        """Compare two CoverReports."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_2 = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                     reference=self.test_ref, out_base=self.test_outdir)
        assert test_report == test_2
        test_2.set_open_status(True)
        assert test_report != test_2
        test_2.set_open_status(False)
        assert test_report == test_2
        test_report.add_processed_file(pathlib.Path("path/to/processed.fq"))
        assert test_report != test_2

    def test_compare_class(self):
        """Compare a CoverReport with a different class."""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        TupleReport = namedtuple("TupleReport",
                                 ["workflow_table", "threshold", "maxDepth", "reference", "out_base",
                                  "is_open", "processed_files"])
        test_tuple = TupleReport(workflow_table=self.test_df, threshold=10, maxDepth=100,
                                 reference=self.test_ref, out_base=self.test_outdir, is_open=False, processed_files=[])
        assert test_tuple != test_report

    def test_repr(self):
        """Print the CoverReport."""
        expected_str = (f"CoverReport(fastq_dir={str(self.test_dir)}, sample_sheet={self.test_samplesheet}, "
                        "threshold=10, maxDepth=100,"
                        f" reference={self.test_ref}, region_file=None,"
                        f" out_base={self.test_outdir}, processed_files=[], is_open=False)\n"
                        f"Workflow table:\n{self.test_df.to_string()}")
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_str = test_report.__repr__()
        assert expected_str == test_str
        expected_processed = (f"CoverReport(fastq_dir={str(self.test_dir)}, sample_sheet={self.test_samplesheet}, "
                              "threshold=10, maxDepth=100,"
                              f" reference={self.test_ref}, region_file=None,"
                              f" out_base={self.test_outdir}, "
                              "processed_files=['path/to/processed_1.fq', 'path/to/processed_2.fq'], is_open=False"
                              ")\n"
                              f"Workflow table:\n{self.test_df.to_string()}")
        test_report.add_processed_file(pathlib.Path("path/to/processed_1.fq"))
        test_report.add_processed_file(pathlib.Path("path/to/processed_2.fq"))
        test_processed = test_report.__repr__()
        assert expected_processed == test_processed
        test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_ref.fa"
        test_region = (pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_region.bed")
        expected_region = (f"CoverReport(fastq_dir={str(self.test_dir)}, sample_sheet={self.test_samplesheet},"
                           " threshold=10, maxDepth=100,"
                            f" reference={test_ref}, region_file={test_region},"
                        f" out_base={self.test_outdir}, processed_files=[], is_open=False"
                           ")\n"
                        f"Workflow table:\n{self.test_df.head().to_string()}")
        region_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10,
                                            maxDepth=100, reference=test_ref, out_base=self.test_outdir,
                                            region_file=test_region)
        region_str = region_report.__repr__()
        assert expected_region == region_str

    def test_save_settings(self):
        """Save the settings of a CoverReport in order to be able to resume."""
        settings = self.test_outdir / "settings.json"
        assert not settings.exists()
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        to_add = ["path/to/processed.fq", "path/to/processed_2.fq"]
        test_report.add_processed_files([pathlib.Path(processed) for processed in to_add])
        # TODO: do we want to move the processed files to the settings?
        expected_results = {"fastq_dir": str(self.test_dir), "sample_sheet": str(self.test_samplesheet), "threshold": 10,
                            "maxDepth": 100, "reference": str(self.test_ref), "out_base": str(self.test_outdir),
                            "processed_files": to_add,
                            "region_file": None}
        test_report.save_settings()
        assert settings.exists()
        with open(settings, "r", encoding = "utf-8") as settings_file:
            test_settings = json.load(settings_file)
        assert expected_results == test_settings

    def test_save_and_resolve(self):
        """Resolve paths when saving a CoverReport."""
        settings = self.test_outdir / "settings.json"
        assert not settings.exists()
        test_report = seq_mon.CoverReport(self.test_dir.relative_to(pathlib.Path(__file__).parent.parent),
                                          sample_sheet=self.test_samplesheet.relative_to(
                                              pathlib.Path(__file__).parent.parent), threshold=10, maxDepth=100,
                                          reference=self.test_ref.relative_to(pathlib.Path(__file__).parent.parent),
                                          out_base=self.test_outdir.relative_to(pathlib.Path(__file__).parent.parent))
        to_add = ["path/to/processed.fq", "path/to/processed_2.fq"]
        test_report.add_processed_files([pathlib.Path(processed) for processed in to_add])
        # TODO: do we want to move the processed files to the settings?
        expected_results = {"fastq_dir": str(self.test_dir), "sample_sheet": str(self.test_samplesheet),
                            "threshold": 10,
                            "maxDepth": 100, "reference": str(self.test_ref), "out_base": str(self.test_outdir),
                            "processed_files": to_add,
                            "region_file": None}
        test_report.save_settings()
        assert settings.exists()
        with open(settings, "r", encoding="utf-8") as settings_file:
            test_settings = json.load(settings_file)
        assert expected_results == test_settings


class TestLoadSettings(unittest.TestCase):
    def test_success_load(self):
        """Regenerate a CoverReport from saved settings."""
        test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir"
        test_samplesheet = pathlib.Path(__file__).parent / "data" / "samplesheet.xls"
        expected_report = seq_mon.CoverReport(indir=(pathlib.Path(__file__).parent / "data" / "workflow_table" /
                                                     "rundir"), sample_sheet=test_samplesheet, threshold=100,
                                              maxDepth=1000, reference=test_ref,
                                              out_base=pathlib.Path(__file__).parent / "data" / "out_dir_2")
        expected_report.port = 50542
        processed = (pathlib.Path(__file__).parent / 'data' /'workflow_table'/
                     'rundir'/'test_dir'/'fastq_pass'/'barcode01'/'barcode01-0.fastq.gz')
        expected_report.add_processed_file(processed)
        settings = pathlib.Path(__file__).parent / "data" / "out_dir_2"/ "settings.json"
        test_report = seq_mon.load_report_from_settings(settings)
        assert expected_report == test_report

    def test_success_no_processed(self):
        """Regenerate a CoverReport from saved settings when there are no processed files."""
        test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir"
        test_samplesheet = pathlib.Path(__file__).parent / "data" / "samplesheet.xls"
        expected_report = seq_mon.CoverReport(indir=(pathlib.Path(__file__).parent / "data" / "workflow_table" /
                                                     "rundir"), sample_sheet=test_samplesheet, threshold=100,
                                              maxDepth=1000, reference=test_ref,
                                              out_base=pathlib.Path(__file__).parent / "data" / "out_dir_2")
        settings = pathlib.Path(__file__).parent / "data" / "test_no_processed.json"
        test_report = seq_mon.load_report_from_settings(settings)
        assert expected_report == test_report


class TestCreateBam(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog
    # mock subprocess run here
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    def test_success(self, mock_run):
        """Create a new bam file."""
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        out_base = pathlib.Path(__file__).parent / "data" / "out_dir"
        test_out_sam = str(out_base / "tmp.sam")
        test_out_bam = out_base / "barcode01.bam"
        test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_ref.fa"
        test_file = pathlib.Path('path/to/processed_1.fq')
        # log the commands and the bit where we're creating a new file
        map_log = f"minimap2 -a -o {test_out_sam} {str(test_ref)} {test_file}"
        new_bam = "Creating initial bam"
        bam_log = f"samtools sort -O bam -o {test_out_bam} {test_out_sam}"
        map_call = mock.call(map_log.split(), check=True)
        bam_call = mock.call(bam_log.split(), check=True)
        with self._caplog.at_level(logging.DEBUG, logger = "seq_mon"):
            seq_mon.create_bam(test_file, out_base, test_out_bam, test_ref)
            assert ("seq_mon", logging.DEBUG, map_log) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, new_bam) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, bam_log) in self._caplog.record_tuples
        mock_run.assert_has_calls([map_call, bam_call])



class TestAppendBam(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog
    # mock subprocess run here
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    def test_success(self, mock_run):
        """Append to an existing bam file."""
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        out_base = pathlib.Path(__file__).parent / "data" / "out_dir"
        test_out_sam = str(out_base / "tmp.sam")
        test_out_bam = out_base / "barcode01.bam"
        test_out_sorted = str(out_base / "sorted.bam")
        test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir" / "test_ref.fa"
        test_file = pathlib.Path('path/to/processed_1.fq')
        map_log = f"minimap2 -a -o {test_out_sam} {str(test_ref)} {test_file}"
        bam_log = f"samtools sort -O bam -o {test_out_sorted} {test_out_sam}"
        merge_log = f"samtools merge -f -o {str(out_base / 'tmp.bam')} {test_out_sorted} {test_out_bam}"
        map_call = mock.call(map_log.split(), check=True)
        bam_call = mock.call(bam_log.split(), check=True)
        merge_call = mock.call(merge_log.split(),  check=True)
        move_call = mock.call(['mv',
                               pathlib.PosixPath('/Users/kat/PycharmProjects/CoverMon/tests/data/out_dir/tmp.bam'),
                               pathlib.PosixPath('/Users/kat/PycharmProjects/CoverMon/tests/data/out_dir/barcode01.bam')],
                               check=True)
        with self._caplog.at_level(logging.DEBUG, logger = "seq_mon"):
            seq_mon.append_bam(test_file, out_base, test_out_bam, test_ref)
            print(self._caplog.record_tuples)
            assert ("seq_mon", logging.DEBUG, map_log) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, bam_log) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, merge_log) in self._caplog.record_tuples
        mock_run.assert_has_calls([map_call, bam_call, merge_call, move_call])


class TestGetDepths(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog
    # mock subprocess run here
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    def test_success(self, mock_run):
        """Get depths for an existing bam file."""
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        out_base = pathlib.Path(__file__).parent / "data" / "out_dir"
        test_out_bam = out_base / "barcode01.bam"
        test_out_depth = out_base / "barcode01.depth"
        index_log = f"samtools index {test_out_bam}"
        depth_log = f"samtools depth -aa {test_out_bam} -o {test_out_depth}"
        index_call = mock.call(index_log.split(), check=True)
        depth_call = mock.call(depth_log.split(), check=True)
        with self._caplog.at_level(logging.DEBUG, logger = "seq_mon"):
            seq_mon.get_depth(test_out_bam, test_out_depth)
            assert ("seq_mon", logging.DEBUG, index_log) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, depth_log) in self._caplog.record_tuples
        mock_run.assert_has_calls([index_call, depth_call])


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
        seq_mon.write_to_processed(to_write, self.outdir)
        # check we only have the one
        with open(self.outfile, "r", encoding = "utf-8") as f:
            content = f.readlines()
        assert len(content) == 1
        assert content[0] == "teststring\n"
        # write another and check we've appended
        write_more = "test_2"
        seq_mon.write_to_processed(write_more, self.outdir)
        with open(self.outfile, "r", encoding = "utf-8") as f:
            content = f.readlines()
        assert len(content) == 2
        assert content[0] == "teststring\n"
        assert content[1] == "test_2\n"


class TestParseSamplesheet(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog


    def test_fail_bad_format(self):
        """Raise an exception if the sample sheet isn't in a supported format"""
        test_sheet = pathlib.Path(__file__).parent / "data" / "samplesheet.csv"
        error_msg = "The spreadsheet must be excel formatted (.xlsx or .xls)"
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            seq_mon.parse_samplesheet(test_sheet)

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
        reading_msg = f'Reading .xls-type sample sheet "{test_sheet}"'
        cleaning_msg = "Cleaning sample sheet ...                              "
        finished_msg = "Sample sheet cleaned ✓"
        info_msgs = [reading_msg, cleaning_msg, finished_msg]

        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            test_result = seq_mon.parse_samplesheet(test_sheet)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, expected_sheet.to_string()) in self._caplog.record_tuples
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
        reading_msg = f'Reading .xlsx-type sample sheet "{test_sheet}"'
        cleaning_msg = "Cleaning sample sheet ...                              "
        finished_msg = "Sample sheet cleaned ✓"
        info_msgs = [reading_msg, cleaning_msg, finished_msg]

        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            test_result = seq_mon.parse_samplesheet(test_sheet)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, expected_sheet.to_string()) in self._caplog.record_tuples
        pd.testing.assert_frame_equal(expected_sheet, test_result, check_dtype = False)


class TestValidateSamplesheet(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_fail_missing_column_refdir(self):
        """Raise an exception if the sample sheet is missing a required column when reference is given as a dir."""
        test_sheet = pd.DataFrame(data = {"sample_id": ["test1", "test2", "test3", "test4"],
                                              "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                            "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                              "ct": ["10", "20", "30", np.nan],
                                              "other_columns1": ["A", "B", "C", "D"],
                                              "other_columns2": ["1", "2", "3", "4"]
                                              })
        info_msg = "Checking that the necessary columns exist ..."
        error_msg = ("The sample sheet is missing a necessary column. The sample sheet must contain the column barcode,"
                     f" but it only contains ['ct', 'other_columns1', 'other_columns2', 'reference', 'sample_id']")
        with (pytest.raises(KeyError, match=re.escape(error_msg)),
              self._caplog.at_level(logging.DEBUG, logger="seq_mon")):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = False)
            assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples

    def test_fail_missing_column(self):
        """Raise an exception if the sample sheet is missing a required column
        (and reference is not a required column)."""
        test_sheet = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                        "ct": ["10", "20", "30", np.nan],
                                        "other_columns1": ["A", "B", "C", "D"],
                                        "other_columns2": ["1", "2", "3", "4"]
                                        })
        error_msg = ("The sample sheet is missing a necessary column. The sample sheet must contain the column barcode,"
                     f" but it only contains ['ct', 'other_columns1', 'other_columns2', 'sample_id']")
        info_msg = "Checking that the necessary columns exist ..."
        with pytest.raises(KeyError, match=re.escape(error_msg)):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = True)
            assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples

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
        info_msgs = ["Checking that the necessary columns exist ...", "All necessary columns found ✓",
                     "Checking that the barcodes are correctly formatted ..."]
        with (pytest.raises(ValueError, match=re.escape(error_msg)),
              self._caplog.at_level(logging.DEBUG, logger="seq_mon")):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = True)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples

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
        info_msgs = ["Checking that the necessary columns exist ...", "All necessary columns found ✓",
                     "Checking that the barcodes are correctly formatted ...", "Barcodes are correct ✓",
                     "Checking that the barcodes are unique ..."]
        with (pytest.raises(ValueError, match=re.escape(error_msg)),
              self._caplog.at_level(logging.DEBUG, logger="seq_mon")):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file = True)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples

    def test_success(self):
        """Don't complain if everything looks good."""
        test_sheet = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                        "barcode": ["NB01", "NB02", "NB03", "NB04"],
                                        "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                      "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                        "ct": ["10", "20", "30", np.nan],
                                        "other_columns1": ["A", "B", "C", "D"],
                                        "other_columns2": ["1", "2", "3", "4"]
                                        })
        info_msgs = ["Checking that the necessary columns exist ...", "All necessary columns found ✓",
                     "Checking that the barcodes are correctly formatted ...", "Barcodes are correct ✓",
                     "Checking that the barcodes are unique ...", "All barcodes are unique ✓",
                     "These are the samples from the samplesheet you have given:\n"
                     f"{test_sheet.to_string()}\n//"]
        with self._caplog.at_level(logging.INFO, logger="seq_mon"):
            seq_mon.validate_samplesheet(test_sheet, ref_is_file=True)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples


class TestValidateRundir(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_fail_missing_dir(self):
        """Raise an exception if the run directory doesn't exist."""
        test_dir = pathlib.Path(__file__).parent / "data" / "rundir" / "no_such_dir"
        error_msg = "The rundir does not exist."
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            seq_mon.validate_rundir(test_dir)

    def test_fail_timeout(self):
        """Raise an exception if nothing is found after a timeout."""
        pass  # TODO: combine this with a convenience waiting function on a second pass so we can test it without having to wait for ages

    def test_fail_too_many_pass_dirs(self):
        """Raise an exception if there are too many fastq_pass directories."""
        test_dir = pathlib.Path(__file__).parent / "data" / "rundir"
        fastq_pass_dirs = [str(test_dir / "test1" / "fastq_pass"), str(test_dir / "test2" / "fastq_pass")]
        bad_paths = str('\n ').join(fastq_pass_dirs)
        error_msg = ("There seems to be more than one fastq_pass sub-directory beneath the given rundir."
                     " These paths were found:\n"
                     f"{bad_paths}\n"
                     "Please specify a more specific rundir.")
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            seq_mon.validate_rundir(test_dir)

    def test_success(self):
        """Find base directory and the fastq_pass directory it contains."""
        test_dir = pathlib.Path(__file__).parent / "data" / "rundir" / "test1"
        expected_base = test_dir
        expected_fastq = test_dir / "fastq_pass"
        info_msgs = ["Checking that the rundir exists ...                    ",
                     "Looking for MinKNOW-characteristic output:",
                         "  Looking ... ",
                         f"Found                                    ✓",
                         f"Found the following fastq_pass base which will be given to CoverMon: \n  {expected_fastq}\n",
                         f"This is the batch base directory:\n  {expected_base}"]
        with self._caplog.at_level(logging.INFO, logger = "seq_mon"):
            test_base, test_fastq = seq_mon.validate_rundir(test_dir)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
        assert test_base == expected_base
        assert test_fastq == expected_fastq



class TestCreateWorkflowTable(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

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
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            seq_mon.create_workflow_table(test_df, self.test_dir)

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
        info_msg = f"Continuing with the following barcodes:\n{expected_df.to_string()}\n//"
        with self._caplog.at_level(logging.INFO, logger="seq_mon"):
            test_result = seq_mon.create_workflow_table(test_df, self.test_dir)
            assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
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
        info_msg = f"Continuing with the following barcodes:\n{expected_df.to_string()}\n//"
        with self._caplog.at_level(logging.INFO, logger="seq_mon"):
            test_result = seq_mon.create_workflow_table(test_df, self.test_dir)
            assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
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
        info_msg = f"Continuing with the following barcodes:\n{expected_df.to_string()}\n//"
        with self._caplog.at_level(logging.INFO, logger="seq_mon"):
            test_result = seq_mon.create_workflow_table(test_df, self.test_dir)
            assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
        pd.testing.assert_frame_equal(expected_df, test_result)


class TestGetPlotCommand(unittest.TestCase):
    test_dir = pathlib.Path(__file__).parent / "data" / "workflow_table" / "rundir" / "test_dir" / "fastq_pass"
    test_df = pd.DataFrame(data={"sample_id": ["test1", "test2"],
                                 "barcode": ["RB01", "RB02"],
                                 "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                 "ct": ["10", "20"],
                                 "other_columns1": ["A", "B"],
                                 "other_columns2": ["1", "2", ],
                                 "barcode_path": [str(test_dir / "barcode01"),
                                                  str(test_dir / "barcode02")],
                                 "barcode_basename": ["barcode01", "barcode02"]
                                 })
    test_outdir = pathlib.Path(__file__).parent / "data" / "out_dir"
    test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir"
    test_samplesheet = pathlib.Path(__file__).parent / "data" / "samplesheet.xls"
    script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"

    def test_plot_ref_dir(self):
        """Get the command to plot coverage with reference as a directory."""
        plot_cmd = ["Rscript", str(self.script_path), str(self.test_outdir), str(self.test_samplesheet),
                    "10",  # cutoff
                    "100"  # max depth
                    ]
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_cmd = seq_mon.get_r_cmd(test_report)
        assert test_cmd == plot_cmd

    def test_plot_ref_file(self):
        """Get the command to plot coverage with reference as a file."""
        plot_cmd = ["Rscript", str(self.script_path), str(self.test_outdir), str(self.test_samplesheet),
                    "10",  # cutoff
                    "100",  # max depth
                    "--region_file", f"{self.test_ref}/test_region.bed"]
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref / "test_ref.fa", out_base=self.test_outdir,
                                          region_file=self.test_ref / "test_region.bed")
        test_cmd = seq_mon.get_r_cmd(test_report)
        assert test_cmd == plot_cmd

# TODO: make this nicer
def bam_creator(*args) -> None:
    """Create bam files for testing"""
    test_outdir =  pathlib.Path(__file__).parent / "data" / "out_dir"
    for bam in [test_outdir / "RB01.bam", test_outdir / "RB02.bam"]:
        bam.touch()

# TODO: adapt to save settings in addition to or instead of just documenting opened files
class TestUpdatePlot(unittest.TestCase):
    test_dir = pathlib.Path(__file__).parent / "data" / "workflow_table" / "rundir" / "test_dir" / "fastq_pass"
    test_df = pd.DataFrame(data={"sample_id": ["test1", "test2"],
                                 "barcode": ["RB01", "RB02"],
                                 "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                 "ct": ["10", "20"],
                                 "other_columns1": ["A", "B"],
                                 "other_columns2": ["1", "2",],
                                 "barcode_path": [str(test_dir / "barcode01"),
                                                  str(test_dir / "barcode02")],
                                 "barcode_basename": ["barcode01", "barcode02"]
                                 })
    test_outdir =  pathlib.Path(__file__).parent / "data" / "out_dir"
    test_ref = pathlib.Path(__file__).parent / "data" / "ref_dir"
    test_samplesheet = pathlib.Path(__file__).parent / "data" / "samplesheet.xls"
    test_settings = test_outdir / "settings.json"

    def tearDown(self):
        to_clean = [self.test_outdir / "RB01.bam", self.test_outdir / "RB02.bam",
                    self.test_outdir / "RB01.depth", self.test_outdir / "RB02.depth",
                    self.test_settings]
        with open((self.test_outdir / "processed_files.txt"), "w", encoding = "utf-8") as f:
            f.write("")
        for test_file in to_clean:
            if test_file.exists():
                test_file.unlink()

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    # TODO: mock create_bam to only touch the bam file
    # TODO mock subprocess
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_success_already_opened(self, mock_create_bam, mock_run):
        """Successfully update the plot when the browser is already open"""
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        test_report.set_open_status(True)
        assert test_report.is_open
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        scanning = "Scanning for new fastq files..."
        barcode1_process1 = "Number of processed files: 1"
        barcode2_process1 = "Number of processed files: 2"
        barcode2_process2 = "Number of processed files: 3"
        report_status = "Updated plot"
        opening_file = "Opening report"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_outdir), str(self.test_samplesheet),
                    "10", # cutoff
                    "100", # max depth
                    ]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        assert not self.test_settings.exists()
        with self._caplog.at_level(logging.DEBUG, logger = "seq_mon"):
            seq_mon.update_plot(test_report)
            assert ("seq_mon", logging.INFO, scanning) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode1_process1) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode2_process1) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode2_process2) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, report_status) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, plot_msg) in self._caplog.record_tuples
            assert not ("seq_mon", logging.INFO, opening_file) in self._caplog.record_tuples
        assert test_report.is_open
        assert self.test_settings.exists()
        found_settings = seq_mon.load_report_from_settings(self.test_settings)
        found_settings.set_open_status(True)  # this defaults to False when loading
        expected_processed = [self.test_dir / "barcode01" / "barcode01-0.fastq.gz",
                              self.test_dir / "barcode02" / "barcode02_0.fastq.gz",
                              self.test_dir / "barcode02" / "barcode02_1.fastq.gz"]
        assert sorted(expected_processed) == sorted(test_report.processed_files)
        with open((self.test_outdir / "processed_files.txt"), "r", encoding = "utf-8") as f:
            processed_from_file = f.read().splitlines()
        assert sorted([str(processed) for processed in expected_processed]) == sorted(processed_from_file)
        assert found_settings == test_report
        mock_run.assert_has_calls([plot_call], any_order=True)

    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_success_open_window(self, mock_create_bam, mock_run, mock_open):
        """Successfully update the plot when the browser is not yet open; update report status"""
        # output to a consistent "port"
        port = 50542
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref, out_base=self.test_outdir)
        assert not test_report.is_open
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        scanning = "Scanning for new fastq files..."
        barcode1_process1 = "Number of processed files: 1"
        barcode2_process1 = "Number of processed files: 2"
        barcode2_process2 = "Number of processed files: 3"
        report_status = "Updated plot"
        opening_file = "Opening report"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_outdir), str(self.test_samplesheet),
                    "10",  # cutoff
                    "100",  # max depth
                    ]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        assert not self.test_settings.exists()
        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            seq_mon.update_plot(test_report)
            assert ("seq_mon", logging.INFO, scanning) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode1_process1) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode2_process1) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode2_process2) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, report_status) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, plot_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, opening_file) in self._caplog.record_tuples
        assert test_report.is_open
        assert self.test_settings.exists()
        expected_processed = [self.test_dir / "barcode01" / "barcode01-0.fastq.gz",
                              self.test_dir / "barcode02" / "barcode02_0.fastq.gz",
                              self.test_dir / "barcode02" / "barcode02_1.fastq.gz"]
        assert sorted(expected_processed) == sorted(test_report.processed_files)
        mock_run.assert_has_calls([plot_call], any_order=True)
        mock_open.assert_called_with(f"file://{self.test_outdir / 'plot_cov.html'}")

    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_success_single_ref(self, mock_create_bam, mock_run, mock_open):
        """Update the plot when we're using a single reference and a region file"""
        # output to a consistent "port"
        port = 50542
        test_report = seq_mon.CoverReport(self.test_dir, sample_sheet=self.test_samplesheet, threshold=10, maxDepth=100,
                                          reference=self.test_ref / "test_ref.fa", out_base=self.test_outdir,
                                          region_file=self.test_ref / "test_region.bed")
        assert not test_report.is_open
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        scanning = "Scanning for new fastq files..."
        barcode1_process1 = "Number of processed files: 1"
        barcode2_process1 = "Number of processed files: 2"
        barcode2_process2 = "Number of processed files: 3"
        report_status = "Updated plot"
        opening_file = "Opening report"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_outdir), str(self.test_samplesheet),
                    "10",  # cutoff
                    "100",  # max depth
                    "--region_file", f"{self.test_ref}/test_region.bed"]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            seq_mon.update_plot(test_report)
            assert ("seq_mon", logging.INFO, scanning) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode1_process1) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode2_process1) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, barcode2_process2) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, report_status) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, plot_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.INFO, opening_file) in self._caplog.record_tuples
        assert test_report.is_open
        expected_processed = [self.test_dir / "barcode01" / "barcode01-0.fastq.gz",
                              self.test_dir / "barcode02" / "barcode02_0.fastq.gz",
                              self.test_dir / "barcode02" / "barcode02_1.fastq.gz"]
        assert sorted(expected_processed) == sorted(test_report.processed_files)
        mock_run.assert_has_calls([plot_call])
        mock_open.assert_called_with(f"file://{self.test_outdir / 'plot_cov.html'}")


class TestStartCovermon(unittest.TestCase):
    test_indir = pathlib.Path(__file__).parent / "data"  / "workflow_table" / "rundir"
    test_default_out = test_indir / "test_dir" / "CoverMon_test_ref"
    test_default_ref_is_dir = test_indir / "test_dir" / "CoverMon"
    test_outdir = pathlib.Path(__file__).parent / "data" / "out_dir"
    test_out_existing = pathlib.Path(__file__).parent / "data" / "out_dir_2"
    test_ref = str(pathlib.Path(__file__).parent / "data" / "ref_dir")
    test_samplesheet = str(pathlib.Path(__file__).parent / "data" / "samplesheet.xls")

    def tearDown(self):
        to_clean = [self.test_outdir / "RB01.bam", self.test_outdir / "RB02.bam",
                    self.test_outdir / "RB01.depth", self.test_outdir / "RB02.depth",
                    self.test_outdir / "sample_sheet_given.tsv",
                    self.test_outdir / "settings.json",
                    self.test_default_out / "RB01.bam", self.test_default_out / "RB02.bam",
                    self.test_default_out / "RB01.depth", self.test_default_out / "RB02.depth",
                    self.test_default_out / "sample_sheet_given.tsv", self.test_default_out / "processed_files.txt",
                    self.test_default_out / "settings.json",
                    self.test_default_ref_is_dir / "RB01.bam", self.test_default_ref_is_dir / "RB02.bam",
                    self.test_default_ref_is_dir / "RB01.depth", self.test_default_ref_is_dir / "RB02.depth",
                    self.test_default_ref_is_dir / "sample_sheet_given.tsv",
                    self.test_default_ref_is_dir / "processed_files.txt",
                    self.test_default_ref_is_dir / "settings.json",
                    self.test_out_existing / "RB02.bam", self.test_out_existing / "RB02.depth",
                    self.test_out_existing / "sample_sheet_given.tsv"]
        for test_file in to_clean:
            if test_file.exists():
                test_file.unlink()
        with open((self.test_outdir / "processed_files.txt"), "w", encoding="utf-8") as f:
            f.write("")
        with open(self.test_out_existing / "processed_files.txt", "w", encoding="utf-8") as f:
            f.write("tests/data/workflow_table/rundir/test_dir/fastq_pass/barcode01/barcode01-0.fastq.gz")
        with open(self.test_out_existing / "settings.json", "w", encoding="utf-8") as f:
            settings = """{
  "fastq_dir":  "tests/data/workflow_table/rundir/test_dir/fastq_pass",
  "sample_sheet":  "tests/data/samplesheet.xls",
  "out_base": "tests/data/out_dir_2",
  "reference": "tests/data/ref_dir",
  "threshold": 100,
  "maxDepth": 1000,
  "region_file": "None",
  "processed_files": ["tests/data/workflow_table/rundir/test_dir/fastq_pass/barcode01/barcode01-0.fastq.gz"]
}"""
            f.write(settings)
        dirs_to_clean = [self.test_default_out, self.test_default_ref_is_dir]
        for dir_to_clean in dirs_to_clean:
            if dir_to_clean.exists():
                shutil.rmtree(dir_to_clean)

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    # TODO waiting/timeout tests?
    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_default_one_ref(self, mock_create_bam, mock_run,
                             mock_open):
        """Output to the default output directory when using a single reference."""
        port = 50542
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        cover_threshold = "100"
        max_depth = "1000"
        info_msgs = ["Backing up the original sample sheet...",
                    "These are the parameters given:", f"This is the samplesheet: {self.test_samplesheet}",
                     f"This is the run directory: {self.test_indir}",
                     f"This is the reference: {self.test_ref}/test_ref.fa",
                     f"This is the threshold: {cover_threshold}", f"This is the maxDepth in plot: {max_depth}",
                     f"Creating output directory {self.test_default_out}...",
                     "  The sequencing summary has been found. Run complete    ✓"]
        initializing_msg = "Initializing report as closed"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_default_out), str(self.test_samplesheet),
                    str(cover_threshold),  # cutoff
                    str(max_depth),  # max depth
                    ]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        start_args = [self.test_samplesheet, str(self.test_indir), f"{self.test_ref}/test_ref.fa", cover_threshold,
                      max_depth]
        assert not self.test_default_out.exists()
        with self._caplog.at_level(logging.DEBUG, logger = "seq_mon"):
            seq_mon.start_covermon(start_args)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, initializing_msg) in self._caplog.record_tuples
        assert self.test_default_out.exists()
        assert (self.test_default_out / "sample_sheet_given.tsv").exists()
        mock_run.assert_has_calls([plot_call])
        mock_open.assert_called_with(f"file://{self.test_default_out / 'plot_cov.html'}")


    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_default_refdir(self, mock_create_bam, mock_run, mock_open):
        """Output to the default output directory when using a reference directory."""
        port = 50542
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        cover_threshold = "100"
        max_depth = "1000"
        info_msgs = ["Backing up the original sample sheet...",
                     "These are the parameters given:", f"This is the samplesheet: {self.test_samplesheet}",
                     f"This is the run directory: {self.test_indir}",
                     f"This is the reference directory: {self.test_ref}",
                     f"This is the threshold: {cover_threshold}", f"This is the maxDepth in plot: {max_depth}",
                     f"Creating output directory {self.test_default_ref_is_dir}...",
                     "  The sequencing summary has been found. Run complete    ✓"]
        initializing_msg = "Initializing report as closed"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_default_ref_is_dir), str(self.test_samplesheet),
                    str(cover_threshold),  # cutoff
                    str(max_depth),  # max depth
                    ]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        start_args = [self.test_samplesheet, str(self.test_indir), self.test_ref, cover_threshold,
                      max_depth]
        assert not self.test_default_ref_is_dir.exists()
        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            seq_mon.start_covermon(start_args)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, initializing_msg) in self._caplog.record_tuples
        assert self.test_default_ref_is_dir.exists()
        assert (self.test_default_ref_is_dir / "sample_sheet_given.tsv").exists()
        mock_open.assert_called_with(f"file://{self.test_default_ref_is_dir / 'plot_cov.html'}")
        mock_run.assert_has_calls([plot_call])


    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_custom_output(self, mock_create_bam, mock_run, mock_open):
        """Output to a user-specified output directory."""
        port = 50542
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        cover_threshold = "100"
        max_depth = "1000"
        info_msgs = ["Backing up the original sample sheet...",
                     "These are the parameters given:", f"This is the samplesheet: {self.test_samplesheet}",
                     f"This is the run directory: {self.test_indir}",
                     f"This is the reference directory: {self.test_ref}",
                     f"This is the threshold: {cover_threshold}", f"This is the maxDepth in plot: {max_depth}",
                     f"Creating output directory {self.test_outdir}...",
                     "  The sequencing summary has been found. Run complete    ✓"]
        initializing_msg = "Initializing report as closed"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_outdir), str(self.test_samplesheet),
                    str(cover_threshold),  # cutoff
                    str(max_depth),  # max depth
                    ]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        start_args = [self.test_samplesheet, str(self.test_indir), self.test_ref, cover_threshold,
                      max_depth, "--outdir", str(self.test_outdir)]
        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            seq_mon.start_covermon(start_args)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, initializing_msg) in self._caplog.record_tuples
        assert self.test_outdir.exists()
        assert (self.test_outdir / "sample_sheet_given.tsv").exists()
        mock_open.assert_called_with(f"file://{self.test_outdir / 'plot_cov.html'}")
        mock_run.assert_has_calls([plot_call])



    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_region_file(self, mock_create_bam, mock_run, mock_open):
        """Use a region file."""
        port = 50542
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        cover_threshold = "100"
        max_depth = "1000"
        info_msgs = ["Backing up the original sample sheet...",
                     f"This is the region file: {self.test_ref}/test_region.bed",
                     "These are the parameters given:", f"This is the samplesheet: {self.test_samplesheet}",
                     f"This is the run directory: {self.test_indir}",
                     f"This is the reference: {self.test_ref}/test_ref.fa",
                     f"This is the threshold: {cover_threshold}", f"This is the maxDepth in plot: {max_depth}",
                     f"Creating output directory {self.test_outdir}...",
                     "  The sequencing summary has been found. Run complete    ✓"]
        initializing_msg = "Initializing report as closed"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_outdir), str(self.test_samplesheet),
                    str(cover_threshold),  # cutoff
                    str(max_depth),  # max depth
                    "--region_file", f"{self.test_ref}/test_region.bed"]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)

        start_args = [self.test_samplesheet, str(self.test_indir), f"{self.test_ref}/test_ref.fa", cover_threshold,
                      max_depth, "--outdir", str(self.test_outdir), "--region_file", f"{self.test_ref}/test_region.bed"]
        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            seq_mon.start_covermon(start_args)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, initializing_msg) in self._caplog.record_tuples
        assert self.test_outdir.exists()
        assert (self.test_outdir / "sample_sheet_given.tsv").exists()
        mock_open.assert_called_with(f"file://{self.test_outdir / 'plot_cov.html'}")
        mock_run.assert_has_calls([plot_call])

    @mock.patch(f"{seq_mon.__name__}.random.randrange")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_fail_wrong_settings(self, mock_create_bam, mock_run, mock_rand):
        """Fail if continuing in a directory with different settings."""
        port = 50500
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        mock_rand.return_value = port
        cover_threshold = "100"
        max_depth = "1000"
        sheet = pd.DataFrame(data={"sample_id": ["test1", "test2", "test3", "test4"],
                                            "barcode": ["NB01", "NB02", "NB03", "NB04"],
                                            "reference": ["nCoV-2019.reference.fa", "nCoV-2019.reference.fa",
                                                          "nCoV-2019.reference.fa", "nCoV-2019.reference.fa"],
                                            "ct": [10, 20, 30, np.nan],
                                            "other_columns1": ["A", "B", "C", "D"],
                                            "other_columns2": [1, 2, 3, 4],
                                   "barcode_path": [str(self.test_indir / "test_dir" / "fastq_pass" / "barcode01"),
                                                    str(self.test_indir / "test_dir" / "fastq_pass" /"barcode02"),
                                                    str(self.test_indir / "test_dir" / "fastq_pass" /"barcode03"),
                                                    str(self.test_indir /"test_dir" / "fastq_pass" /  "barcode04")],
                                   "barcode_basename": ["barcode01", "barcode02", "barcode03", "barcode04"]
                                            })
        start_args = [self.test_samplesheet, str(self.test_indir), f"{self.test_ref}/test_ref.fa", cover_threshold,
                      max_depth, "--outdir", str(self.test_out_existing)]
        error_msg = ("Attempting to continue an interrupted run with different settings. Please choose a different"
                     " output directory to run analysis with the new settings.\n"
                     "Old settings:\n"
                     f"CoverReport(fastq_dir={str(self.test_indir / 'test_dir' / 'fastq_pass')},"
                     f" sample_sheet={self.test_samplesheet}, "
                     "threshold=100, maxDepth=1000,"
                     f" reference={self.test_ref}, region_file=None,"
                     f" out_base={self.test_out_existing}, "
                     f"processed_files=['{self.test_indir}/test_dir/fastq_pass/barcode01/barcode01-0.fastq.gz'],"
                     " is_open=False)\n"
                     f"Workflow table:\n{sheet.head().to_string()}\n"
                     f"New settings:\n"
                     f"CoverReport(fastq_dir={str(self.test_indir / 'test_dir' / 'fastq_pass')},"
                     f" sample_sheet={self.test_samplesheet}, "
                     "threshold=100, maxDepth=1000,"
                     f" reference={self.test_ref}/test_ref.fa, region_file=None,"
                     f" out_base={self.test_out_existing}, "
                     f"processed_files=['{self.test_indir}/test_dir/fastq_pass/barcode01/barcode01-0.fastq.gz'],"
                     " is_open=False)\n"
                     f"Workflow table:\n{sheet.head().to_string()}")
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            print(error_msg)
            seq_mon.start_covermon(start_args)

    @mock.patch(f"{seq_mon.__name__}.webbrowser.open_new_tab")
    @mock.patch(f"{seq_mon.__name__}.subprocess.run")
    @mock.patch(f"{seq_mon.__name__}.create_bam")
    def test_existing_processed(self, mock_create_bam, mock_run, mock_open):
        """Output to a directory already containing processed files."""
        # we don't want to accidentally hit our target port while randomly generating the port
        port = 50500
        target_port = 50542
        mock_subprocess = mock.Mock()
        mock_run.return_value = mock_subprocess
        mock_create_bam.side_effect = bam_creator
        cover_threshold = "100"
        max_depth = "1000"
        info_msgs = ["Backing up the original sample sheet...",
                     "These are the parameters given:", f"This is the samplesheet: {self.test_samplesheet}",
                     f"This is the run directory: {self.test_indir}",
                     f"This is the reference directory: {self.test_ref}",
                     f"This is the threshold: {cover_threshold}", f"This is the maxDepth in plot: {max_depth}",
                     f"Creating output directory {self.test_out_existing}...",
                     "I have found processed files, opening existing report",
                     "Current settings match saved settings",
                     "  The sequencing summary has been found. Run complete    ✓"]
        initializing_msg = "Initializing report as closed"
        script_path = pathlib.Path(__file__).parent.parent / "scripts" / "run_plot.R"
        plot_cmd = ["Rscript", str(script_path), str(self.test_out_existing), str(self.test_samplesheet),
                    str(cover_threshold),  # cutoff
                    str(max_depth)  # max depth
                    ]
        plot_msg = " ".join(plot_cmd)
        plot_call = mock.call(plot_msg, shell=True, check=True)
        start_args = [self.test_samplesheet, str(self.test_indir), self.test_ref, cover_threshold,
                      max_depth, "--outdir", str(self.test_out_existing)]
        with self._caplog.at_level(logging.DEBUG, logger="seq_mon"):
            seq_mon.start_covermon(start_args)
            for info_msg in info_msgs:
                assert ("seq_mon", logging.INFO, info_msg) in self._caplog.record_tuples
            assert ("seq_mon", logging.DEBUG, initializing_msg) in self._caplog.record_tuples
        assert self.test_out_existing.exists()
        assert (self.test_out_existing / "sample_sheet_given.tsv").exists()
        mock_open.assert_called_with(f"file://{self.test_out_existing / 'plot_cov.html'}")
        mock_run.assert_has_calls([plot_call])
