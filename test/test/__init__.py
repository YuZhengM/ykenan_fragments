#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from ykenan_fragments import Run
from ykenan_fragments.get_sort_fragments import GetSortFragments

"""
ykenan_fragments:
Github: https://github.com/YuZhengM/ykenan_fragments
"""


class Fragments(GetSortFragments):
    """
    继承方法, 重写开始数据处理的部分
    """

    def __init__(self, source_path: str, merge_path: str, barcodes_path: str, gsm: str, handler_path: str, lift_over_path: str, is_hg19_to_hg38: bool = True, thread_count: int = 10):
        super().__init__(source_path, merge_path, barcodes_path, gsm, handler_path, lift_over_path, is_hg19_to_hg38, thread_count)
        self.endswith_list = ["_barcodes.tsv.gz", "_matrix.mtx.gz", "_peaks.bed.gz"]
        self.barcodes_key: str = "barcodes"
        self.mtx_key: str = "mtx"
        self.peaks_key: str = "peaks"

    @staticmethod
    def get_peaks(line: str) -> str:
        return line

    @staticmethod
    def get_barcodes(line: str) -> str:
        return line


if __name__ == '__main__':
    """
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181251
    """
    Run(path="F:/software/scATAC_data/data", lift_over_path="/mnt/f/software/liftOver", callback=Fragments)
