#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path
from multiprocessing.dummy import Pool
from multiprocessing.pool import ThreadPool
from typing import TextIO

import ykenan_file as yf
from ykenan_log import Logger


class Hg19AndHg38:
    """
    该步骤需要去子系统中进入该路径执行进行执行
    """

    def __init__(self, path: str, lift_over_path: str, is_hg19_to_hg38: bool, thread_count: int = 10):
        # log 日志信息
        self.log = Logger("liftOver", "log")
        # 处理路径和文件的方法
        self.file = yf.StaticMethod(log_file="log")
        self.read = yf.Read(header=None, log_file="log")
        self.create = yf.Create(header=False, log_file="log")
        self.thread_count = thread_count
        self.lift_over_path = lift_over_path
        self.is_hg19_to_hg38 = is_hg19_to_hg38
        self.transformation_file: str = "hg19ToHg38.over.chain.gz" if self.is_hg19_to_hg38 else "hg38ToHg19.over.chain.gz"
        self.source = os.path.join(path, "hg19" if self.is_hg19_to_hg38 else "hg38")
        self.output_no_sort = os.path.join(path, "hg38_no_sort" if self.is_hg19_to_hg38 else "hg19_no_sort")
        self.output_chr_sort = os.path.join(path, "hg38" if self.is_hg19_to_hg38 else "hg19")
        self.unmap = os.path.join(path, self.transformation_file + "_unmap")
        self.run()
        self.chr_sort()

    def exec_str(self, filename: str) -> str:
        return f"{self.lift_over_path}/liftOver {os.path.join(self.source, filename)} {self.lift_over_path}/{self.transformation_file} {os.path.join(self.output_no_sort, filename)} {os.path.join(self.unmap, filename)}"

    def exec_command(self, command: str) -> list:
        """
        执行系统命令
        :param command: 命令代码
        :return: 结果数组
        """
        self.log.info(f">>>>>>>>>>>>>>>>>>>>>>>> start 执行 {command} 命令 >>>>>>>>>>>>>>>>>>>>>>>>")
        info: str = os.popen(command).read()
        info_split: list = info.split("\n")
        info_list: list = []
        i: int = 0
        while True:
            if info_split[i] is None or info_split[i] == "":
                break
            info_list.append(info_split[i])
            i += 1
        self.log.info(f">>>>>>>>>>>>>>>>>>>>>>>> end 执行 {command} 命令 >>>>>>>>>>>>>>>>>>>>>>>>")
        return info_list

    def run(self):

        if not os.path.exists(self.source):
            self.log.error(f"输入文件夹 {self.source} 不存在")
            raise ValueError(f"输入文件夹 {self.source} 不存在")

        # 创建文件夹
        self.file.makedirs(self.output_no_sort)
        self.file.makedirs(self.unmap)

        # 获取没有执行的文件
        input_files = self.file.get_files(path=self.source)
        code_list = []
        finish_files = self.file.get_files(path=self.output_no_sort)
        for input_file in input_files:
            if input_file not in finish_files:
                code_list.append(self.exec_str(input_file))

        if len(code_list) != 0:
            # 实例化线程对象
            pool: ThreadPool = Pool(self.thread_count)
            # 将 list 的每一个元素传递给 pool_page(page) 处理
            pool.map(self.exec_command, code_list)
            # 关闭线程
            pool.close()

    def chr_sort(self):
        files_dict: dict = self.file.entry_files_dict(self.output_no_sort)
        files_name: list = files_dict["name"]
        # 读取数量
        fragments_count: int = 0
        chr_f_list: list = []
        chr_f_dict: dict = {}
        if os.path.exists(self.output_chr_sort):
            self.log.info(f"{self.output_chr_sort} The folder has default sorting completed")
            return
        self.log.info(f"create folder {self.output_chr_sort}")
        os.makedirs(self.output_chr_sort)
        for filename in files_name:
            filename: str
            with open(files_dict[filename], "r", encoding="utf-8") as r:
                self.log.info(f"Start chr sort {files_dict[filename]} file")
                while True:
                    line: str = r.readline().strip()
                    if not line:
                        break
                    if fragments_count >= 250000 and fragments_count % 250000 == 0:
                        self.log.info(f"processed {fragments_count} 行")
                    split: list = line.split("\t")
                    chromosome: str = split[0]
                    if chromosome not in chr_f_list:
                        chromosome_tsv_file = "_".join(filename.split("_")[0:-1]) + "_" + chromosome + ".tsv"
                        chr_f_list.append(chromosome)
                        chr_f = open(os.path.join(self.output_chr_sort, chromosome_tsv_file), "w", encoding="utf-8", newline="\n", buffering=1)
                        chr_f_dict.update({chromosome: chr_f})
                    # Obtaining files with added content
                    chromosome_file: TextIO = chr_f_dict[chromosome]
                    chromosome_file.write(f"{line}\n")
                    fragments_count += 1
        # 关闭文件
        for chromosome in chr_f_list:
            chromosome_file: TextIO = chr_f_dict[chromosome]
            chromosome_file.close()
