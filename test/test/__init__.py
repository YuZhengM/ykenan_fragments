#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import re

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

    def __init__(self, source_path: str, merge_path: str, gsm: str, handler_path: str, lift_over_path: str, is_hg19_to_hg38: bool = True):
        super().__init__(source_path, merge_path, gsm, handler_path, lift_over_path, is_hg19_to_hg38)
        self.endswith_list = ["_cell_annotation.txt.gz", "_raw_matrix.mtx.gz", "_peak_annotation.txt.gz"]
        self.barcodes_key: str = "cell_annotation"
        self.mtx_key: str = "raw_matrix"
        self.peaks_key: str = "peak_annotation"

    def get_peaks(self, dict_: dict, index_: str) -> str:
        peak: str = dict_[int(index_)]
        peak_split = peak.split(" ")
        index = int(re.sub("\"", "", peak_split[0]))
        if index != index_:
            self.log.error(f"位置不对.... {peak} {index_} ===> {index}")
            raise ValueError(f"位置不对.... {peak} {index_} ===> {index}")
        peak_info: str = re.sub("\"", "", peak_split[1])
        peak_split_after: list = peak_info.split("-")
        peak_split_before = str(peak_split_after[0]).split(":")
        return f"{peak_split_before[0]}\t{peak_split_before[1]}\t{peak_split_after[2]}"

    @staticmethod
    def get_barcodes(dict_: dict, index_: str) -> str:
        barcode: str = dict_[int(index_)]
        return barcode.split("\t")[0]


if __name__ == '__main__':
    """
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181251
    """
    Run(path="F:/software/scATAC_data/data", lift_over_path="/mnt/f/software/liftOver", finish_gse=["GSE129785"], callback=Fragments)
