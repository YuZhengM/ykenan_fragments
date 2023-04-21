#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import itertools
import os
import shutil
from typing import TextIO

import pandas as pd
from pandas import DataFrame
from ykenan_log import Logger
import ykenan_file as yf
import gzip


# log 日志

class GetFragments:

    def __init__(self, base_path: str, cp_path: str, GSE: str):
        """
        Form an unordered fragments
        :param base_path: Path to store three files
        :param cp_path: Path to generate unordered fragments files
        :param GSE: GSE number (here is a folder name)
        """
        self.log = Logger("Three files form fragments file", "log/fragments.log")
        self.file = yf.staticMethod(log_file="log")
        # Folder path containing three files
        self.GSE: str = GSE
        self.base_path: str = os.path.join(base_path, GSE)
        # The path to copy the fragments file to
        self.cp_path: str = os.path.join(cp_path, GSE)
        # keyword
        self.barcodes_key: str = "barcodes"
        self.mtx_key: str = "mtx"
        self.peaks_key: str = "peaks"
        # Extract files and remove suffix information
        self.endswith_list: list = [".cell_barcodes.txt.gz", ".mtx.gz", ".peaks.txt.gz"]
        self.suffix_fragments: str = ".tsv"
        self.suffix_information: str = ".txt"
        # start processing
        self.exec_fragments()

    def handler_source_files(self) -> dict:
        # Obtain gz file information
        files: list = self.file.get_files(self.base_path)
        gz_files: list = []
        # Obtain content without a suffix for creating folders and generating file names using
        gz_files_before: list = []
        gz_files_and_before: dict = {}
        self.log.info(f"Filter file information under {self.base_path} path")
        for file in files:
            file: str
            for endswith in self.endswith_list:
                if file.endswith(endswith):
                    gz_files.append(file)
                    before = file.split(endswith)[0]
                    if gz_files_before.count(before) == 0:
                        gz_files_before.append(before)
                    gz_files_and_before = dict(itertools.chain(gz_files_and_before.items(), {
                        file: before
                    }.items()))
                    break
        # Void judgment
        if len(gz_files) == 0:
            self.log.info(f"Gz compressed file is 0")
        else:
            # 简单验证
            if len(gz_files) % 3 == 0:
                self.log.info("The file is a multiple of 3, correct")
            else:
                self.log.warn("The file is not a multiple of 3, there is an error")

            # create folder
            gz_file_before_dirs: dict = {}
            for gz_file_before in gz_files_before:
                gz_file_before_dir = os.path.join(self.base_path, gz_file_before)
                gz_file_before_dirs = dict(itertools.chain(gz_file_before_dirs.items(), {
                    gz_file_before: gz_file_before_dir
                }.items()))
                # Determine whether the folder is created
                if os.path.exists(gz_file_before_dir):
                    continue
                self.log.info(f"create folder {gz_file_before_dir}")
                os.mkdir(gz_file_before_dir)

            # move file
            for gz_file in gz_files:
                file_source = os.path.join(self.base_path, gz_file)
                file_target = gz_file_before_dirs[gz_files_and_before[gz_file]]
                self.log.info(f"move file {file_source} to {file_target}")
                shutil.move(file_source, file_target)

        # Get folder information
        self.log.info(f"Starting to obtain information for processing ==========> ")
        dirs_dict: dict = self.file.entry_dirs_dict(self.base_path)
        dirs_name = dirs_dict["name"]

        dirs_key: list = []
        dirs_key_dict: dict = {}

        # Determine if the folder contains three files
        self.log.info(f"Filter folder information under {self.base_path} path")
        for dir_name in dirs_name:
            get_files = self.file.get_files(dirs_dict[dir_name])
            if len(get_files) < 3:
                continue
            count: int = 0
            isAdd: bool = True
            for file in get_files:
                # Determine if there are folders that have already formed fragments
                if dir_name + self.suffix_fragments == file or dir_name + self.suffix_information == file:
                    self.log.info(f"{dir_name} The fragments file has been generated")
                    self.log.warn(f"Skip generation of {dir_name} type fragments file")
                    isAdd = False
                    break
                # Does it contain three
                for endswith in self.endswith_list:
                    if dir_name + endswith == file:
                        count += 1
                        break
            if count == 3 and isAdd:
                dirs_key.append(dir_name)
                dirs_key_dict = dict(itertools.chain(dirs_key_dict.items(), {
                    dir_name: dirs_dict[dir_name]
                }.items()))
        return {
            "all": {
                "key": dirs_name,
                "path": dirs_dict
            },
            "no_finish": {
                "key": dirs_key,
                "path": dirs_key_dict
            }
        }

    def get_files(self, path: str) -> dict:
        # Obtain all file information under this path
        contents_dict: dict = self.file.entry_contents_dict(path, 1)
        filenames: list = contents_dict["name"]
        barcodes_file: dict = {}
        mtx_file: dict = {}
        peaks_file: dict = {}
        # pick up information
        self.log.info(f"得到 {path} 三个文件信息")
        for filename in filenames:
            filename: str
            # Determine if it is a compressed package
            if filename.endswith(".gz"):
                if filename.count(self.barcodes_key) > 0:
                    barcodes_file: dict = {
                        "name": filename,
                        "path": contents_dict[filename]
                    }
                    self.log.info(f"{self.barcodes_key} 文件: {barcodes_file}")
                elif filename.count(self.mtx_key) > 0:
                    mtx_file: dict = {
                        "name": filename,
                        "path": contents_dict[filename]
                    }
                    self.log.info(f"{self.mtx_key} 文件: {mtx_file}")
                elif filename.count(self.peaks_key) > 0:
                    peaks_file: dict = {
                        "name": filename,
                        "path": contents_dict[filename]
                    }
                    self.log.info(f"{self.peaks_key} 文件: {peaks_file}")
        return {
            self.barcodes_key: barcodes_file,
            self.mtx_key: mtx_file,
            self.peaks_key: peaks_file
        }

    @staticmethod
    def get_file_content(path: str, file: dict):
        txt_file: str = os.path.join(path, file["name"].split(".txt")[0]) + ".txt"
        # Determine if the file exists
        if txt_file.endswith(".mtx.gz.txt"):
            if not os.path.exists(txt_file):
                with open(txt_file, 'wb') as w:
                    with gzip.open(file["path"], 'rb') as f:
                        # Form a file
                        w.write(f.read())
            return txt_file
        else:
            if os.path.exists(txt_file):
                f = gzip.open(file["path"], 'rb')
                # Obtaining Content Information
                file_content: list = f.read().decode().rstrip().split("\n")
                f.close()
                return file_content
            else:
                w = open(txt_file, 'wb')
                f = gzip.open(file["path"], 'rb')
                read = f.read()
                # Form a file
                w.write(read)
                # Obtaining Content Information
                file_content: list = read.decode().rstrip().split("\n")
                f.close()
                w.close()
                return file_content

    def fragments_file_name(self, key: str) -> str:
        return f"{key}{self.suffix_fragments}"

    def information_file_name(self, key: str) -> str:
        return f"{key}{self.suffix_information}"

    def write_fragments(self, path: str, key: str) -> None:
        """
        Form fragments file
        :return:
        """
        # Obtain file information
        files: dict = self.get_files(path)
        # Get Barcodes
        self.log.info(f"正在获取 {self.barcodes_key} 文件信息")
        barcodes: list = self.get_file_content(path, files[self.barcodes_key])
        self.log.info(f"正在获取 {self.mtx_key} 文件路径")
        mtx_path: str = self.get_file_content(path, files[self.mtx_key])
        self.log.info(f"正在获取 {self.peaks_key} 文件信息")
        peaks: list = self.get_file_content(path, files[self.peaks_key])
        # length
        barcodes_len: int = len(barcodes)
        peaks_len: int = len(peaks)
        if barcodes_len < 2 or peaks_len < 2:
            self.log.error(f"Insufficient file read length {self.barcodes_key}: {barcodes_len}, {self.peaks_key}: {peaks_len}")
            raise ValueError("Insufficient file read length")
        self.log.info(f"Quantity or Path {self.barcodes_key}: {barcodes_len}, {self.mtx_key}: {mtx_path}, {self.peaks_key}: {peaks_len}")
        # Judge comment information
        mtx_start: int = 0
        # Read quantity
        mtx_count: int = 0
        error_count: int = 0
        mtx_all_number: int = 0
        # create a file
        fragments_file: str = os.path.join(path, self.fragments_file_name(key))
        self.log.info(f"Starting to form {mtx_path} fragments file")
        with open(fragments_file, "w", encoding="utf-8", buffering=1, newline="\n") as w:
            with open(mtx_path, "r", encoding="utf-8") as r:
                while True:
                    line: str = r.readline().strip()
                    if not line:
                        break
                    if mtx_count >= 500000 and mtx_count % 500000 == 0:
                        if mtx_all_number != 0:
                            self.log.info(f"已处理 {mtx_count} 行, 已完成 {round(mtx_count / mtx_all_number, 4) * 100} %")
                        else:
                            self.log.info(f"已处理 {mtx_count} 行")
                    if line.startswith("%"):
                        mtx_start += 1
                        self.log.info(f"注释信息: {line}")
                        mtx_count += 1
                        continue
                    split: list = line.split(" ")
                    # To determine the removal of a length of not 3
                    if len(split) != 3:
                        mtx_count += 1
                        error_count += 1
                        self.log.error(f"mtx information ===> content: {split}, line number: {mtx_count}")
                        continue
                    if mtx_start == mtx_count:
                        mtx_start = 0
                        self.log.info(f"Remove Statistical Rows: {line}")
                        mtx_all_number = int(split[2])
                        if int(split[0]) + 1 == peaks_len and int(split[1]) + 1 == barcodes_len:
                            mtx_count += 1
                            continue
                        else:
                            mtx_count += 1
                            raise ValueError(f"File mismatch {self.peaks_key}: {int(split[0])} {peaks_len}, {self.barcodes_key}: {int(split[1])} {barcodes_len}")
                    if int(split[0]) > peaks_len or int(split[1]) > barcodes_len:
                        mtx_count += 1
                        continue
                    # peak, barcode, There is a header+1, but the index starts from 0 and the record starts from 1
                    peak: str = peaks[int(split[0])]
                    barcode: str = barcodes[int(split[1])]
                    peak_split = peak.split("_")
                    barcode_split = barcode.split("\t")
                    # Adding information, it was found that some files in mtx contain two columns, less than three columns. This line was ignored and recorded in the log
                    try:
                        w.write(f"{peak_split[0]}\t{peak_split[1]}\t{peak_split[2]}\t{barcode_split[6]}\t{split[2]}\n")
                    except Exception as e:
                        error_count += 1
                        self.log.error(f"peak information: {peak_split}")
                        self.log.error(f"barcodes file information: {barcode}")
                        self.log.error(f"barcodes information: {barcode_split}")
                        self.log.error(f"mtx information ===> content: {split}, line number: {mtx_count}")
                        self.log.error(f"Write error: {e}")
                    mtx_count += 1
        self.log.info(f"The number of rows ignored is {error_count}, {round(error_count / mtx_all_number, 4) * 100} % of total")
        self.log.info(f"Complete the formation of {mtx_path} fragments file")

    def copy_file(self, source_file: str, target_file: str) -> None:
        if os.path.exists(target_file):
            self.log.warn(f"{target_file} The file already exists, it has been copied by default")
        else:
            self.log.info(f"Start copying file {source_file}")
            shutil.copy(source_file, target_file)
            self.log.info(f"End of copying file  {source_file}")

    def cp_files(self, path: str, key: str) -> None:
        fragments_file_name = self.fragments_file_name(key)
        fragments_file: str = os.path.join(path, fragments_file_name)
        # Determine if it exists
        if not (os.path.exists(fragments_file)):
            self.log.error(f"file does not exist: {fragments_file}")
            raise ValueError(f"file does not exist: {fragments_file}")
        # Two folders
        fragments_cp_dir = os.path.join(self.cp_path, "fragments")
        if not os.path.exists(fragments_cp_dir):
            self.log.info(f"create folder {fragments_cp_dir}")
            os.makedirs(fragments_cp_dir)
        # copy
        fragments_gz_file = os.path.join(fragments_cp_dir, f"{fragments_file_name}.gz")
        if os.path.exists(fragments_gz_file):
            self.log.warn(f"The file has been compressed into {fragments_gz_file}, Default copy completed")
        else:
            self.copy_file(fragments_file, os.path.join(fragments_cp_dir, fragments_file_name))

    def exec_fragments(self):
        # Classify the types and place them in different folders
        source_files: dict = self.handler_source_files()
        no_finish_infor = source_files["no_finish"]
        no_finish_keys = no_finish_infor["key"]
        no_finish_paths = no_finish_infor["path"]
        self.log.info(f"Related file information {no_finish_keys}, {no_finish_paths}")
        # Form fragments file
        for key in no_finish_keys:
            self.log.info(f"Process {key} related files (folders)")
            self.write_fragments(no_finish_paths[key], key)
            self.log.info(f"Complete processing of {key} related files (folders)")
        # All information
        all_infor = source_files["all"]
        all_infor_keys = all_infor["key"]
        all_infor_paths = all_infor["path"]
        # copy file
        for key in all_infor_keys:
            self.log.info(f"Start copying files to the specified path for {key}")
            self.cp_files(all_infor_paths[key], key)
            self.log.info(f"Copy file to specified path for {key} completed")


class GetChrSortFragments:

    def __init__(self, path: str, cp_path: str, GSE: str, get_fragments_path: str, is_exec: bool = True):
        """
        Form an fragments
        :param path: Path to store unordered fragments files
        :param cp_path: Path to generate fragments files
        :param GSE: GSE number (here is a folder name)
        :param get_fragments_path: base_path parameter in GetFragments class
        :param is_exec: Do you want to execute the GetFragments class
        """
        self.log = Logger("Three files form fragments file", "log/fragments.log")
        self.file = yf.staticMethod(log_file="log")
        self.base_path: str = os.path.join(path, GSE)
        self.fragments_path: str = os.path.join(self.base_path, "fragments")
        self.cp_input_path: str = cp_path
        self.chr_list: dict = {
            "chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4, "chr5": 5, "chr6": 6, "chr7": 7, "chr8": 8, "chr9": 9, "chr10": 10,
            "chr11": 11, "chr12": 12, "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16, "chr17": 17, "chr18": 18, "chr19": 19, "chr20": 20,
            "chr21": 21, "chr22": 22, "chrX": 23, "chrY": 24
        }
        if is_exec:
            GetFragments(get_fragments_path, path, GSE)
        self.exec_sort_fragments()

    @staticmethod
    def classification_name(chromosome: str, path: str, name: str):
        # Do not use any methods under os.path for path operations here, as looping will slow down several times
        # splitext_name = os.path.splitext(name)[0]
        # splitext_name_suffix = os.path.splitext(name)[1]
        # chromosome_path_file: str = os.path.join(path, f"{splitext_name}_{chromosome}{splitext_name_suffix}")
        return f"{path}/{name}_{chromosome}.tsv"

    def get_files(self) -> dict:
        if not os.path.exists(self.fragments_path):
            self.log.error(f"The input file {self.fragments_path} does not exist. Please check")
            raise ValueError(f"The input file {self.fragments_path} does not exist. Please check")
        # Obtain tsv file information under the folder
        files_dict: dict = self.file.entry_contents_dict(self.fragments_path, type_=1, suffix=".tsv")
        files_dict_name = files_dict["name"]
        self.log.info(f"tsv file information: {files_dict_name}")
        need_handler_fragments: list = []
        need_handler_fragments_path: dict = {}
        if not os.path.exists(self.cp_input_path):
            self.log.info(f"create folder {self.cp_input_path}")
            os.makedirs(self.cp_input_path)
        # Add processing files
        for file in files_dict_name:
            # 排序后的文件
            ArchR_fragments_file: str = os.path.join(self.cp_input_path, file)
            if os.path.exists(ArchR_fragments_file):
                self.log.warn(f"The fragments file {ArchR_fragments_file} sorted by chromatin already exists")
                continue
            # 添加信息
            need_handler_fragments.append(file)
            need_handler_fragments_path = dict(itertools.chain(need_handler_fragments_path.items(), {
                file: files_dict[file]
            }.items()))
        return {
            "name": need_handler_fragments,
            "path": need_handler_fragments_path
        }

    def write_chr_file(self, path: str, file: str) -> dict:
        # 读取数量
        fragments_count: int = 0
        # error_count: int = 0
        chr_f_list: list = []
        chr_f_dict: dict = {}
        chr_f_path: dict = {}
        # Determine whether to merge directly
        isMerge: bool = True
        # Create a folder to store chromatin
        chromosome_path: str = os.path.join(self.fragments_path, f"{file}_chromosome")
        if not os.path.exists(chromosome_path):
            self.log.info(f"create folder {chromosome_path}")
            isMerge = False
            os.makedirs(chromosome_path)
        with open(path, "r", encoding="utf-8") as r:
            while True:
                line: str = r.readline().strip()
                if not line:
                    break
                if fragments_count >= 500000 and fragments_count % 500000 == 0:
                    self.log.info(f"processed {fragments_count} 行")
                split: list = line.split("\t")
                # To determine if an error stop occurs when the length is not 5
                # if len(split) != 5:
                #     fragments_count += 1
                #     error_count += 1
                #     log.error(f"fragments file error line ===> content: {split}, line number: {fragments_count}")
                #     raise ValueError(f"fragments file error line ===> content: {split}, line number: {fragments_count}")
                chromosome: str = split[0]
                if not isMerge:
                    chromosome_path_file: str = self.classification_name(chromosome, chromosome_path, file)
                    # Do not judge os. path. exists in this area, as the speed will decrease by 50 times when the number of cycles exceeds 500000
                    # if chromosome not in chr_f_list and not os.path.exists(chromosome_path_file):
                    if chromosome not in chr_f_list:
                        chr_f_list.append(chromosome)
                        chr_f = open(chromosome_path_file, "w", encoding="utf-8", newline="\n", buffering=1)
                        chr_f_dict = dict(itertools.chain(chr_f_dict.items(), {
                            chromosome: chr_f
                        }.items()))
                        chr_f_path = dict(itertools.chain(chr_f_path.items(), {
                            chromosome: chromosome_path_file
                        }.items()))
                    # Obtaining files with added content
                    chromosome_file: TextIO = chr_f_dict[chromosome]
                    chromosome_file.write(f"{line}\n")
                else:
                    chromosome_path_file: str = self.classification_name(chromosome, chromosome_path, file)
                    if chromosome not in chr_f_list:
                        chr_f_list.append(chromosome)
                        chr_f_path = dict(itertools.chain(chr_f_path.items(), {
                            chromosome: chromosome_path_file
                        }.items()))
                fragments_count += 1
        # 关闭文件
        if not isMerge:
            for chromosome in chr_f_list:
                chromosome_file: TextIO = chr_f_dict[chromosome]
                chromosome_file.close()
        return {
            "name": chr_f_list,
            "path": chr_f_path
        }

    def sort_position_files(self, chr_file_dict: dict, file: str):
        chr_name: list = chr_file_dict["name"]
        file_dict_path: dict = chr_file_dict["path"]
        position_f_path: dict = {}
        # sort
        chr_name.sort(key=lambda elem: self.chr_list[elem])
        # Determine whether to merge directly
        isMerge: bool = True
        # output file
        position: str = os.path.join(self.fragments_path, f"{file}_position")
        if not os.path.exists(position):
            self.log.info(f"create folder {position}")
            isMerge = False
            os.makedirs(position)

        if not isMerge:
            for chr_ in chr_name:
                self.log.info(f"Start sorting file {file_dict_path[chr_]} Sort")
                chr_file_content: DataFrame = pd.read_table(file_dict_path[chr_], encoding="utf-8", header=None)
                # 进行排序
                chr_file_content.sort_values(1, inplace=True)
                position_file: str = os.path.join(position, f"{file}_{chr_}.tsv")
                chr_file_content.to_csv(position_file, sep="\t", encoding="utf-8", header=False, index=False)
                position_f_path = dict(itertools.chain(position_f_path.items(), {
                    chr_: position_file
                }.items()))
                self.log.info(f"To file {chr_} Sort completed")
        else:
            for chr_ in chr_name:
                position_file: str = os.path.join(position, f"{file}_{chr_}.tsv")
                position_f_path = dict(itertools.chain(position_f_path.items(), {
                    chr_: position_file
                }.items()))
        return {
            "name": chr_name,
            "path": position_f_path
        }

    def merge_chr_files(self, chr_file_dict: dict, output_file: str) -> None:
        chr_name: list = chr_file_dict["name"]
        file_dict_path: dict = chr_file_dict["path"]
        # 排序
        chr_name.sort(key=lambda elem: self.chr_list[elem])
        self.log.info(f"Start merging file {chr_name}")
        # 生成文件
        with open(output_file, "w", encoding="utf-8", newline="\n", buffering=1) as w:
            for chr_ in chr_name:
                self.log.info(f"Start adding {file_dict_path[chr_]} file")
                with open(file_dict_path[chr_], "r", encoding="utf-8") as r:
                    while True:
                        line: str = r.readline().strip()
                        if not line:
                            break
                        w.write(f"{line}\n")
                self.log.info(f"Completed adding {chr_} file")

    def exec_sort_fragments(self) -> None:
        files_dict: dict = self.get_files()
        files_name: list = files_dict["name"]
        files_path: dict = files_dict["path"]
        for file in files_name:
            # output file
            chr_sort_fragments_file: str = os.path.join(self.cp_input_path, file)
            if os.path.exists(chr_sort_fragments_file):
                self.log.warn(f"{chr_sort_fragments_file} The file already exists, it has been processed by default")
                continue
            self.log.info(f"Start to group {file} files according to chromatin information")
            chr_file_dict: dict = self.write_chr_file(files_path[file], file)
            self.log.info(f"File information after grouping {chr_file_dict}")
            self.log.info(f"Complete file grouping of {file} according to chromatin information")
            self.log.info(f"Start sorting {file} grouped files")
            position_file_dict: dict = self.sort_position_files(chr_file_dict, file)
            self.log.info(f"Sorted file information {position_file_dict}")
            self.log.info(f"Sorting {file} group files completed")
            self.log.info(f"Start merging {file} grouped files")
            self.merge_chr_files(position_file_dict, chr_sort_fragments_file)
            self.log.info(f"Merge {file} group files completed")
