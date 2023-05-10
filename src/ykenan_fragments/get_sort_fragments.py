#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os

from multiprocessing.pool import ThreadPool
from typing import TextIO

import pandas as pd
from pandas import DataFrame

from multiprocessing.dummy import Pool

from ykenan_fragments.genome_transformation import Hg19AndHg38
from ykenan_fragments.get_fragments import GetFragments


class GetSortFragments(GetFragments):

    def __init__(self, path: str, merge_path: str, gsm: str, get_fragments_path: str, lift_over_path: str, is_hg19_to_hg38: bool = True):
        """
        Form a fragments
        :param path: Path to store unordered fragments files
        :param merge_path: Path to generate fragments files
        :param gsm: GSE number (here is a folder name)
        :param get_fragments_path: base_path parameter in GetFragments class
        :param is_hg19_to_hg38: 是否为 hg19 文件
        """
        super().__init__(get_fragments_path, path, gsm)
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
        self.exec_sort_fragments()

    @staticmethod
    def classification_name(chromosome: str, path: str, name: str):
        # Do not use any methods under os.path for path operations here, as looping will slow down several times
        # splitext_name = os.path.splitext(name)[0]
        # splitext_name_suffix = os.path.splitext(name)[1]
        # chromosome_path_file: str = os.path.join(path, f"{splitext_name}_{chromosome}{splitext_name_suffix}")
        return f"{path}/{name}_{chromosome}.tsv"

    def get_need_files(self) -> dict:
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
                        chr_f_dict.update({chromosome: chr_f})
                        chr_f_path.update({chromosome: chromosome_path_file})
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
                chr_f_path.update({chromosome: chromosome_file_dict[filename]})
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
                Hg19AndHg38(path=base_path, lift_over_path=self.lift_over_path, is_hg19_to_hg38=self.is_hg19_to_hg38)

        for chr_ in chr_name:
            genome_file: str = os.path.join(genome_output, f"{file}_{chr_}.tsv")
            genome_f_path.update({chr_: genome_file})
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
                position_f_path.update({chr_: position_file})
                sort_position_files_core_param_list.append((position, file_dict_path, chr_, file))
            # 实例化线程对象
            pool: ThreadPool = Pool(self.thread_count)
            # Form fragments file
            pool.map(self.sort_position_files_core, sort_position_files_core_param_list)
            pool.close()
        else:
            for chr_ in chr_name:
                position_file: str = os.path.join(position, f"{file}_{chr_}.tsv")
                position_f_path.update({chr_: position_file})
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
        files_dict: dict = self.get_need_files()
        for genome in self.genomes:

            # 判断是否继续
            if not self.lift_over_path and genome == self.genome_generate:
                continue
            files_name: list = files_dict[genome]["name"]
            files_path: dict = files_dict[genome]["path"]

            # 创建文件夹
            merge_input_path_genome = os.path.join(self.merge_input_path, genome)
            if not os.path.exists(merge_input_path_genome):
                self.log.info(f"create {merge_input_path_genome} 文件夹")
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
            pool: ThreadPool = Pool(self.thread_count)
            # Form fragments file
            pool.map(self.chr_sort_fragments_file_core, param_list)
            pool.close()
