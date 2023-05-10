#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import shutil
from multiprocessing.pool import ThreadPool

from ykenan_log import Logger
import ykenan_file as yf
import gzip
from multiprocessing.dummy import Pool


class GetFragments:

    def __init__(self, base_path: str, cp_path: str, gsm: str, thread_count: int = 10):
        """
        Form an unordered fragments
        :param base_path: Path to store three files
        :param cp_path: Path to generate unordered fragments files
        :param gsm: GSE number (here is a folder name)
        """
        self.log = Logger("Three files form fragments file", "log/fragments.log")
        self.file = yf.StaticMethod(log_file="log")
        # Folder path containing three files
        self.thread_count: int = thread_count
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
                        gz_files_and_before.update({file: before})
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
                gz_file_before_dirs.update({gz_file_before: gz_file_before_dir})
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
                dirs_key_dict.update({dir_name: dirs_dict[dir_name]})
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

    @staticmethod
    def judge_mtx_is_true(one_len: str, two_len: str, peaks_len: int, barcodes_len: int) -> bool:
        return int(one_len) + 1 != peaks_len and int(two_len) + 1 != barcodes_len

    @staticmethod
    def get_peaks(dict_: dict, index_: str) -> str:
        peak: str = dict_[int(index_)]
        peak_split = peak.split("_")
        return f"{peak_split[0]}\t{peak_split[1]}\t{peak_split[2]}"

    @staticmethod
    def get_barcodes(dict_: dict, index_: str) -> str:
        barcode: str = dict_[int(index_)]
        return barcode.split("\t")[6]

    @staticmethod
    def get_input_peaks(dict_: dict, index_: str) -> str:
        return dict_[int(index_)]

    @staticmethod
    def get_input_barcodes(dict_: dict, index_: str) -> str:
        return dict_[int(index_)]

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

        # 转成字典
        barcodes_dict: dict = {}
        for barcodes_i in range(barcodes_len):
            barcodes_dict.update({barcodes_i, barcodes[barcodes_i]})
        peaks_dict: dict = {}
        for peak_i in range(peaks_len):
            peaks_dict.update({peak_i, peaks[peak_i]})

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
                    if self.judge_mtx_is_true(split[0], split[1], peaks_len, barcodes_len):
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
                    peak_info: str = self.get_peaks(peaks_dict, split[0])
                    barcode_info: str = self.get_barcodes(barcodes_dict, split[1])
                    # Adding information, it was found that some files in mtx contain two columns, less than three columns. This line was ignored and recorded in the log
                    try:
                        w.write(f"{peak_info}\t{barcode_info}\t{split[2]}\n")
                    except Exception as e:
                        error_count += 1
                        self.log.error(f"peak information: {self.get_input_peaks(peaks_dict, split[0])}")
                        self.log.error(f"barcodes information: {self.get_input_barcodes(barcodes_dict, split[1])}")
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
        pool: ThreadPool = Pool(self.thread_count)
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
        pool: ThreadPool = Pool(self.thread_count)
        # copy file
        pool.map(self.cp_files, cp_files_param_list)
        pool.close()
