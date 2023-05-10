#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os

from ykenan_log import Logger
import ykenan_file as yf

from ykenan_fragments.get_sort_fragments import GetSortFragments


class Run:

    def __init__(self, path: str, lift_over_path: str = None, is_hg19_to_hg38: bool = True, callback=GetSortFragments):
        self.base_path: str = os.path.join(path, "handler")
        self.source_path: str = os.path.join(path, "source")
        self.log = Logger("Run", "log/fragments.log")
        self.file = yf.StaticMethod(log_file="log")
        self.lift_over_path: str = lift_over_path
        self.is_hg19_to_hg38: bool = is_hg19_to_hg38
        self.callback = callback
        self.exec()

    def exec(self):
        # 尽量保证该路径下只有 GSE 号的文件
        dirs_dict: dict = self.file.entry_dirs_dict(self.source_path)
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
            self.callback(
                path=self.base_path,
                merge_path=archr_path,
                gsm=gsm,
                get_fragments_path=self.source_path,
                lift_over_path=self.lift_over_path,
                is_hg19_to_hg38=self.is_hg19_to_hg38
            )
