#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import gzip
import os.path

from ykenan_log import Logger

from ykenan_file import StaticMethod


class MergeSourceFile:

    def __init__(self, base_path: str, merge_path: str):
        self.log = Logger("MergeSourceFile", "log/fragments.log")
        self.file = StaticMethod(log_file="log")
        self.input_path = base_path
        self.output_path = merge_path
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
        # keyword
        self.barcodes_key: str = "barcodes"
        self.mtx_key: str = "mtx"
        self.peaks_key: str = "features"
        # Extract files and remove suffix information
        self.endswith_list: list = ["_barcodes.tsv.gz", "_matrix.mtx.gz", "_features.tsv.gz"]

    def get_merge_path(self):
        return {
            self.barcodes_key: os.path.join(self.output_path, "all" + self.endswith_list[0]),
            self.mtx_key: os.path.join(self.output_path, "all" + self.endswith_list[1]),
            self.peaks_key: os.path.join(self.output_path, "all" + self.endswith_list[2])
        }

    def get_file_info(self) -> dict:
        files_dict: dict = self.file.entry_files_dict(self.input_path)
        files_name: list = files_dict["name"]
        barcodes_files: list = []
        mtx_files: list = []
        peaks_files: list = []
        for filename in files_name:
            filename: str
            if filename.endswith(".gz"):
                if filename.count(self.barcodes_key) > 0:
                    barcodes_files.append(files_dict[filename])
                elif filename.count(self.mtx_key) > 0:
                    mtx_files.append(files_dict[filename])
                elif filename.count(self.peaks_key) > 0:
                    peaks_files.append(files_dict[filename])
        barcodes_files.sort()
        mtx_files.sort()
        peaks_files.sort()
        return {
            self.barcodes_key: barcodes_files,
            self.mtx_key: mtx_files,
            self.peaks_key: peaks_files
        }

    def format_barcodes_file(self):
        files: list = self.get_file_info()[self.barcodes_key]
        length = len(files)
        with gzip.open(self.get_merge_path()[self.barcodes_key], 'wb') as w:
            for i in range(length):
                # 解压文件
                file_content = self.file.unzip_gz(files[i], f"{files[i]}.txt")
                for info in file_content:
                    line = f"{str(i)}_{info}\n"
                    w.write(line.encode())

    def format_peak_file(self):
        files: list = self.get_file_info()[self.peaks_key]
        with gzip.open(self.get_merge_path()[self.peaks_key], 'wb') as w:
            for file_ in files:
                # 解压文件
                file_content = self.file.unzip_gz(file_, f"{file_}.txt")
                for info in file_content:
                    line = f"{info}\n"
                    w.write(line.encode())

    def format_mtx_file(self):
        files: list = self.get_file_info()[self.mtx_key]
        # 获取数量信息
        line_0: str = ''
        line_1: str = ''
        line_0_number = 0
        line_1_number = 0
        mtx_number = 0
        for file_ in files:
            with gzip.open(file_, 'rb') as f:
                line_0 = f.readline().decode().rstrip()
                line_1 = f.readline().decode().rstrip()
                line: str = f.readline().decode().rstrip()
                split = line.split(" ")
                line_0_number += int(split[0])
                line_1_number += int(split[1])
                mtx_number += int(split[2])

        line_number = 0

        with gzip.open(self.get_merge_path()[self.mtx_key], 'wb') as w:
            w.write(f"{line_0}\n".encode())
            w.write(f"{line_1}\n".encode())
            w.write(f"{str(line_0_number)} {str(line_1_number)} {str(mtx_number)}\n".encode())
            record_line_0_number = 0
            record_line_1_number = 0
            for file_ in files:
                self.log.info(f"开始处理 {file_} 文件")
                with gzip.open(file_, 'rb') as f:
                    f.readline().decode().rstrip()
                    f.readline().decode().rstrip()
                    line: str = f.readline().decode().rstrip()
                    number_info: list = line.split(" ")
                    while True:
                        line: str = f.readline().decode().rstrip()
                        if not line:
                            break
                        if line_number >= 350000 and line_number % 350000 == 0:
                            self.log.info(f"已处理 {line_number} 行, 完成 {line_number / mtx_number * 100} %")
                        split: list = line.split(" ")
                        number_0 = str(int(split[0]) + record_line_0_number)
                        number_1 = str(int(split[1]) + record_line_1_number)
                        line = str(number_0 + " " + number_1 + " " + split[2])
                        # 解压文件
                        w.write(f"{line}\n".encode())
                        line_number += 1

                record_line_0_number += int(number_info[0])
                record_line_1_number += int(number_info[1])
