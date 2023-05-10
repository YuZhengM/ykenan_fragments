#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import itertools
import os
import shutil
from multiprocessing.pool import ThreadPool
from typing import TextIO

import pandas as pd
from pandas import DataFrame
from ykenan_log import Logger
import ykenan_file as yf
import gzip
from multiprocessing.dummy import Pool

from ykenan_fragments.genome_transformation import Hg19ToHg38


class GetFragments:

    def __init__(self, base_path: str, cp_path: str, gsm: str):
        """
        Form an unordered fragments
        :param base_path: Path to store three files
        :param cp_path: Path to generate unordered fragments files
        :param gsm: GSE number (here is a folder name)
        """
        self.log = Logger("Three files form fragments file", "log/fragments.log")
        self.file = yf.staticMethod(log_file="log")
        # Folder path containing three files
        self.GSE: str = gsm
        self.base_path: str = os.path.join(base_path, gsm)
        # The path to copy the fragments file to
        self.cp_path: str = os.path.join(cp_path, gsm)
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
            is_add: bool = True
            for file in get_files:
                # Determine if there are folders that have already formed fragments
                if dir_name + self.suffix_fragments == file or dir_name + self.suffix_information == file:
                    self.log.info(f"{dir_name} The fragments file has been generated")
                    self.log.warn(f"Skip generation of {dir_name} type fragments file")
                    is_add = False
                    break
                # Does it contain three
                for endswith in self.endswith_list:
                    if dir_name + endswith == file:
                        count += 1
                        break
            if count == 3 and is_add:
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
        self.log.info(f"Obtain three file information for {path}")
        for filename in filenames:
            filename: str
            # Determine if it is a compressed package
            if filename.endswith(".gz"):
                if filename.count(self.barcodes_key) > 0:
                    barcodes_file: dict = {
                        "name": filename,
                        "path": contents_dict[filename]
                    }
                    self.log.info(f"{self.barcodes_key} file: {barcodes_file}")
                elif filename.count(self.mtx_key) > 0:
                    mtx_file: dict = {
                        "name": filename,
                        "path": contents_dict[filename]
                    }
                    self.log.info(f"{self.mtx_key} file: {mtx_file}")
                elif filename.count(self.peaks_key) > 0:
                    peaks_file: dict = {
                        "name": filename,
                        "path": contents_dict[filename]
                    }
                    self.log.info(f"{self.peaks_key} file: {peaks_file}")
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

    def write_fragments(self, param: list) -> None:
        """
        Form fragments file
        :return:
        """
        path: str = param[0]
        key: str = param[1]
        self.log.info(f"Process {key} related files (folders)")
        # Obtain file information
        files: dict = self.get_files(path)
        # Get Barcodes
        self.log.info(f"Getting {self.barcodes_key} file information")
        barcodes: list = self.get_file_content(path, files[self.barcodes_key])
        self.log.info(f"Getting {self.mtx_key} file path")
        mtx_path: str = self.get_file_content(path, files[self.mtx_key])
        self.log.info(f"Getting {self.peaks_key} file information")
        peaks: list = self.get_file_content(path, files[self.peaks_key])
        # length
        barcodes_len: int = len(barcodes)
        peaks_len: int = len(peaks)
        if barcodes_len < 2 or peaks_len < 2:
            self.log.error(f"Insufficient file read length {self.barcodes_key}: {barcodes_len}, {self.peaks_key}: {peaks_len}")
            raise ValueError("Insufficient file read length")
        self.log.info(f"Quantity or Path {self.barcodes_key}: {barcodes_len}, {self.mtx_key}: {mtx_path}, {self.peaks_key}: {peaks_len}")
        # Read quantity
        mtx_count: int = 0
        error_count: int = 0
        mtx_all_number: int = 0
        # create a file
        fragments_file: str = os.path.join(path, self.fragments_file_name(key))
        self.log.info(f"Starting to form {mtx_path} fragments file")
        with open(fragments_file, "w", encoding="utf-8", buffering=1, newline="\n") as w:
            with open(mtx_path, "r", encoding="utf-8") as r:
                line: str = r.readline().strip()
                if line.startswith("%"):
                    self.log.info(f"Annotation Information: {line}")
                line: str = r.readline().strip()
                split: list = line.split(" ")
                if len(split) == 3 and line:
                    self.log.info(f"Remove Statistical Rows: {line}")
                    mtx_all_number = int(split[2])
                    if int(split[0]) + 1 != peaks_len and int(split[1]) + 1 != barcodes_len:
                        raise ValueError(f"File mismatch {self.peaks_key}: {int(split[0])} {peaks_len}, {self.barcodes_key}: {int(split[1])} {barcodes_len}")
                while True:
                    line: str = r.readline().strip()
                    if not line:
                        break
                    if mtx_count >= 500000 and mtx_count % 500000 == 0:
                        self.log.info(f"Processed {mtx_count} lines, completed {round(mtx_count / mtx_all_number, 4) * 100} %")
                    split: list = line.split(" ")
                    # To determine the removal of a length of not 3
                    if len(split) != 3:
                        mtx_count += 1
                        error_count += 1
                        self.log.error(f"mtx information ===> content: {split}, line number: {mtx_count}")
                        continue
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
        self.log.info(f"Complete processing of {key} related files (folders)")

    def copy_file(self, source_file: str, target_file: str) -> None:
        if os.path.exists(target_file):
            self.log.warn(f"{target_file} The file already exists, it has been copied by default")
        else:
            self.log.info(f"Start copying file {source_file}")
            shutil.copy(source_file, target_file)
            self.log.info(f"End of copying file  {source_file}")

    def cp_files(self, param: tuple) -> None:
        path: str = param[0]
        key: str = param[1]
        self.log.info(f"Start copying files to the specified path for {key}")
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
        elif os.path.exists(os.path.join(fragments_cp_dir, fragments_file_name)):
            self.log.warn(f"The file has been copy into {fragments_gz_file}, Default copy completed")
        else:
            self.copy_file(fragments_file, os.path.join(fragments_cp_dir, fragments_file_name))
        self.log.info(f"Copy file to specified path for {key} completed")

    def exec_fragments(self):
        # Classify the types and place them in different folders
        source_files: dict = self.handler_source_files()
        no_finish_infor = source_files["no_finish"]
        no_finish_keys = no_finish_infor["key"]
        no_finish_paths = no_finish_infor["path"]
        self.log.info(f"Related file information {no_finish_keys}, {no_finish_paths}")
        # 参数信息
        write_fragments_param_list: list = []
        for key in no_finish_keys:
            write_fragments_param_list.append((no_finish_paths[key], key))
        # 实例化线程对象
        pool: ThreadPool = Pool(10)
        # Form fragments file
        pool.map(self.write_fragments, write_fragments_param_list)
        pool.close()

        # All information
        all_infor = source_files["all"]
        all_infor_keys = all_infor["key"]
        all_infor_paths = all_infor["path"]
        # 参数信息
        cp_files_param_list = []
        for key in all_infor_keys:
            cp_files_param_list.append((all_infor_paths[key], key))
        # 实例化线程对象
        pool: ThreadPool = Pool(10)
        # copy file
        pool.map(self.cp_files, cp_files_param_list)
        pool.close()


class GetChrSortFragments:

    def __init__(self, path: str, merge_path: str, gsm: str, get_fragments_path: str, lift_over_path: str, is_hg19_to_hg38: bool = True, is_exec: bool = True):
        """
        Form an fragments
        :param path: Path to store unordered fragments files
        :param merge_path: Path to generate fragments files
        :param gsm: GSE number (here is a folder name)
        :param get_fragments_path: base_path parameter in GetFragments class
        :param is_hg19_to_hg38: 是否为 hg19 文件
        :param is_exec: Do you want to execute the GetFragments class
        """
        self.log = Logger("Three files form fragments file", "log/fragments.log")
        self.file = yf.staticMethod(log_file="log")
        self.base_path: str = os.path.join(path, gsm)
        self.fragments_path: str = os.path.join(self.base_path, "fragments")
        self.merge_input_path: str = merge_path
        self.lift_over_path: str = lift_over_path
        self.is_hg19_to_hg38: bool = is_hg19_to_hg38
        self.genome_source: str = "hg19" if self.is_hg19_to_hg38 else "hg38"
        self.genome_generate: str = "hg38" if self.is_hg19_to_hg38 else "hg19"
        self.genomes: list = [self.genome_source, self.genome_generate]
        self.chr_list: dict = {
            "chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4, "chr5": 5, "chr6": 6, "chr7": 7, "chr8": 8, "chr9": 9, "chr10": 10,
            "chr11": 11, "chr12": 12, "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16, "chr17": 17, "chr18": 18, "chr19": 19, "chr20": 20,
            "chr21": 21, "chr22": 22, "chrX": 23, "chrY": 24
        }
        if is_exec:
            GetFragments(base_path=get_fragments_path, cp_path=path, gsm=gsm)
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
        need_handler_fragments_result: dict = {}
        if not os.path.exists(self.merge_input_path):
            self.log.info(f"create folder {self.merge_input_path}")
            os.makedirs(self.merge_input_path)

        # Add processing files
        for genome in self.genomes:

            # 创建文件夹
            merge_input_path_genome = os.path.join(self.merge_input_path, genome)
            if not os.path.exists(merge_input_path_genome):
                self.log.info(f"创建 {merge_input_path_genome} 文件夹")
                os.makedirs(merge_input_path_genome)

            need_handler_fragments: list = []
            need_handler_fragments_path: dict = {}
            for file in files_dict_name:
                # 排序后的文件
                archr_fragments_file: str = os.path.join(merge_input_path_genome, file)
                if os.path.exists(archr_fragments_file):
                    self.log.warn(f"The fragments file {archr_fragments_file} sorted by chromatin already exists")
                    continue
                # 添加信息
                need_handler_fragments.append(file)
                need_handler_fragments_path.update({file: files_dict[file]})
            need_handler_fragments_result.update({genome: {
                "name": need_handler_fragments,
                "path": need_handler_fragments_path
            }})
        return need_handler_fragments_result

    def write_chr_file(self, path: str, file: str) -> dict:
        # 读取数量
        fragments_count: int = 0
        # error_count: int = 0
        chr_f_list: list = []
        chr_f_dict: dict = {}
        chr_f_path: dict = {}
        # Determine whether to merge directly
        is_merge: bool = True
        # Create a folder to store chromatin
        chromosome_path: str = os.path.join(self.fragments_path, f"{file}_chromosome", self.genome_source)
        if not os.path.exists(chromosome_path):
            self.log.info(f"create folder {chromosome_path}")
            os.makedirs(chromosome_path)
            is_merge = False
        if not is_merge:
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
                    fragments_count += 1
            # 关闭文件
            for chromosome in chr_f_list:
                chromosome_file: TextIO = chr_f_dict[chromosome]
                chromosome_file.close()
        else:
            chromosome_file_dict: dict = self.file.entry_files_dict(chromosome_path)
            chromosome_file_name: list = chromosome_file_dict["name"]
            for filename in chromosome_file_name:
                filename: str
                chromosome: str = filename.split("_")[-1].split(".")[0]
                chr_f_list.append(chromosome)
                chr_f_path = dict(itertools.chain(chr_f_path.items(), {
                    chromosome: chromosome_file_dict[filename]
                }.items()))
        return {
            "name": chr_f_list,
            "path": chr_f_path,
            "base_path": os.path.join(self.fragments_path, f"{file}_chromosome")
        }

    def genome_transformation(self, chr_file_dict: dict, file: str, genome: str):
        chr_name: list = chr_file_dict["name"]
        chr_name.sort(key=lambda elem: self.chr_list[elem])
        base_path: str = chr_file_dict["base_path"]
        genome_f_path: dict = {}
        # Determine whether to merge directly
        is_merge: bool = True
        # output file
        genome_output: str = os.path.join(base_path, self.genome_generate)
        if not os.path.exists(genome_output):
            self.log.info(f"create folder {genome_output}")
            is_merge = False
            os.makedirs(genome_output)

        if genome == self.genome_generate and self.lift_over_path:
            # 进行转化为 hg19 或 hg38
            if not is_merge:
                # 执行信息
                Hg19ToHg38(path=base_path, lift_over_path=self.lift_over_path, is_hg19_to_hg38=self.is_hg19_to_hg38)

        for chr_ in chr_name:
            genome_file: str = os.path.join(genome_output, f"{file}_{chr_}.tsv")
            genome_f_path = dict(itertools.chain(genome_f_path.items(), {
                chr_: genome_file
            }.items()))
        return {
            self.genome_source: chr_file_dict,
            self.genome_generate: {
                "name": chr_name,
                "path": genome_f_path
            }
        }

    def sort_position_files_core(self, param: tuple):
        position: str = param[0]
        file_dict_path: dict = param[1]
        chr_: str = param[2]
        file: str = param[3]
        self.log.info(f"Start sorting file {file_dict_path[chr_]} Sort")
        chr_file_content: DataFrame = pd.read_table(file_dict_path[chr_], encoding="utf-8", header=None)
        # 进行排序
        chr_file_content.sort_values(1, inplace=True)
        position_file: str = os.path.join(position, f"{file}_{chr_}.tsv")
        chr_file_content.to_csv(position_file, sep="\t", encoding="utf-8", header=False, index=False)
        self.log.info(f"To file {chr_} Sort completed")

    def sort_position_files(self, chr_file_dict: dict, file: str, genome: str):
        chr_name: list = chr_file_dict["name"]
        file_dict_path: dict = chr_file_dict["path"]
        position_f_path: dict = {}
        # sort
        chr_name.sort(key=lambda elem: self.chr_list[elem])
        # Determine whether to merge directly
        is_merge: bool = True
        # output file
        position: str = os.path.join(self.fragments_path, f"{file}_position", genome)
        if not os.path.exists(position):
            self.log.info(f"create folder {position}")
            is_merge = False
            os.makedirs(position)

        if not is_merge:
            sort_position_files_core_param_list = []
            for chr_ in chr_name:
                position_file: str = os.path.join(position, f"{file}_{chr_}.tsv")
                position_f_path = dict(itertools.chain(position_f_path.items(), {
                    chr_: position_file
                }.items()))
                sort_position_files_core_param_list.append((position, file_dict_path, chr_, file))
            # 实例化线程对象
            pool: ThreadPool = Pool(5)
            # Form fragments file
            pool.map(self.sort_position_files_core, sort_position_files_core_param_list)
            pool.close()
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

        if not os.path.exists(output_file):
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
        else:
            self.log.warn(f"{output_file}. The file already exists, it has been processed by default")

    def after_two_step(self, file: str, chr_file_dict: dict, fragments_file: str, genome: str):
        # 对位点进行排序
        self.log.info(f"Start sorting {file} grouped files")
        position_file_dict: dict = self.sort_position_files(chr_file_dict, file, genome)
        self.log.info(f"Sorted file information {position_file_dict}")
        self.log.info(f"Sorting {file} group files completed")
        # 合并文件
        self.log.info(f"Start merging {file} grouped files")
        self.merge_chr_files(position_file_dict, fragments_file)
        self.log.info(f"Merge {file} group files completed")

    def chr_sort_fragments_file_core(self, param_list: list):
        # 参数信息
        files_path: dict = param_list[0]
        file: str = param_list[1]
        chr_sort_fragments_file: str = param_list[2]
        genome: str = param_list[3]

        self.log.info(f"Start to group {file} files according to chromatin information")
        chr_file_dict: dict = self.write_chr_file(files_path[file], file)
        self.log.info(f"File information after grouping {chr_file_dict}")
        self.log.info(f"Complete file grouping of {file} according to chromatin information")

        genome_transformation_dict: dict = self.genome_transformation(chr_file_dict, file, genome)
        self.after_two_step(file, genome_transformation_dict[genome], chr_sort_fragments_file, genome)

    def exec_sort_fragments(self) -> None:
        files_dict: dict = self.get_files()
        for genome in self.genomes:

            # 判断是否继续
            if not self.lift_over_path and genome == self.genome_generate:
                continue
            files_name: list = files_dict[genome]["name"]
            files_path: dict = files_dict[genome]["path"]

            # 创建文件夹
            merge_input_path_genome = os.path.join(self.merge_input_path, genome)
            if not os.path.exists(merge_input_path_genome):
                self.log.info(f"创建 {merge_input_path_genome} 文件夹")
                os.makedirs(merge_input_path_genome)

            param_list: list = []
            for file in files_name:
                # output file
                chr_sort_fragments_file: str = os.path.join(merge_input_path_genome, file)

                if os.path.exists(chr_sort_fragments_file):
                    self.log.warn(f"{chr_sort_fragments_file}. The files already exists, it has been processed by default")
                    continue

                # 添加参数
                param_list.append((files_path, file, chr_sort_fragments_file, genome))

            # 实例化线程对象
            pool: ThreadPool = Pool(5)
            # Form fragments file
            pool.map(self.chr_sort_fragments_file_core, param_list)
            pool.close()


class Run:

    def __init__(self, path: str, lift_over_path: str = None, is_hg19_to_hg38: bool = True):
        self.base_path: str = os.path.join(path, "handler")
        self.source_path: str = os.path.join(path, "source")
        self.log = Logger("Run", "log/fragments.log")
        self.file = yf.staticMethod(log_file="log")
        self.lift_over_path: str = lift_over_path
        self.is_hg19_to_hg38: bool = is_hg19_to_hg38
        self.exec()

    def exec(self):
        # 尽量保证该路径下只有 GSE 号的文件
        dirs_dict: dict = self.file.entry_dirs_dict(self.base_path)
        dirs_name: list = dirs_dict["name"]

        if not os.path.exists(self.source_path):
            self.log.error(f"输入文件夹 {self.source_path} 不存在")

        # 执行
        for gsm in dirs_name:

            archr_path = os.path.join(self.base_path, gsm, "ArchR")
            if not os.path.exists(archr_path):
                self.log.info(f"创建 {archr_path} 文件夹")
                os.makedirs(archr_path)

            self.log.info(f"开始执行 {gsm} 内容信息")
            GetChrSortFragments(
                path=self.base_path,
                merge_path=archr_path,
                gsm=gsm,
                get_fragments_path=self.source_path,
                lift_over_path=self.lift_over_path,
                is_hg19_to_hg38=self.is_hg19_to_hg38
            )
