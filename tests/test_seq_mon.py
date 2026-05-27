import pathlib
import unittest

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
        # write one string
        # check we only have the one
        # write another and check we've appended
        assert False