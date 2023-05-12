#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os

from ykenan_log import Logger
import ykenan_file as yf

from ykenan_fragments.get_sort_fragments import GetSortFragments


class Run:

    def __init__(self, path: str, lift_over_path: str = None, finish_gse: list = None, is_hg19_to_hg38: bool = True, thread_count: int = 10, callback=GetSortFragments):
        self.handler_path: str = os.path.join(path, "handler")
        self.source_path: str = os.path.join(path, "source")
        self.barcodes_path: str = os.path.join(path, "barcodes")
        self.log = Logger("Run", "log/fragments.log")
        self.file = yf.StaticMethod(log_file="log")
        self.lift_over_path: str = lift_over_path
        self.finish_gse: list = finish_gse
        self.is_hg19_to_hg38: bool = is_hg19_to_hg38
        self.thread_count: int = thread_count
        self.callback = callback
        self.exec()

    def exec(self):
        # 尽量保证该路径下只有 GSE 号的文件
        dirs_dict: dict = self.file.entry_dirs_dict(self.source_path)
        dirs_name: list = dirs_dict["name"]
        dirs_name.sort()

        if not os.path.exists(self.source_path):
            self.log.error(f"输入文件夹 {self.source_path} 不存在")

        # 执行
        for gsm in dirs_name:

            # 跳过完成的
            if self.finish_gse and gsm in self.finish_gse:
                continue

            # 创建文件夹
            self.file.makedirs(self.barcodes_path)
            archr_path = os.path.join(self.handler_path, gsm, "ArchR")
            self.file.makedirs(archr_path)

            self.log.info(f"开始执行 {gsm} 内容信息")
            fragment = self.callback(
                source_path=self.source_path,
                merge_path=archr_path,
                barcodes_path=self.barcodes_path,
                gsm=gsm,
                handler_path=self.handler_path,
                lift_over_path=self.lift_over_path,
                is_hg19_to_hg38=self.is_hg19_to_hg38,
                thread_count=self.thread_count
            )
            fragment.exec_fragments()
            fragment.exec_sort_fragments()
