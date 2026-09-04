#!/usr/bin/env python3
# CyriPanel — modified from Cyrius (Copyright (c) 2019-2020 Illumina, Inc.),
# GPL-3.0. Rewritten for the 4-column BED parser: the original test used a
# 6-column Cyrius region file and asserted on region-type keys from its fifth
# column, which CyriPanel no longer produces.

import os
from ..utilities import parse_region_file

data_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "data")


class TestUtilities(object):
    def test_parse_region_file(self):
        """The panel BED must parse without raising: star_caller calls this at
        startup, so a failure here aborts the run."""
        region_dic = parse_region_file(
            os.path.join(data_dir, "PGxProbe_region_hg38.bed"), 38
        )
        assert list(region_dic.keys()) == ["target"]
        assert len(region_dic["target"]) == 784
        (nchr, start, end, name), gc = region_dic["target"][0]
        assert nchr.startswith("chr") and end > start and name and gc == "0.0"
