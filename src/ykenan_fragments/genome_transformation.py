#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path
from multiprocessing.dummy import Pool
import ykenan_file as yf
from ykenan_log import Logger


class Hg19ToHg38:
    """
    该步骤需要去子系统中进入该路径执行进行执行
    """

    def __init__(self, path: str, lift_over_path: str, is_hg19_to_hg38: bool):
        # log 日志信息
        self.log = Logger("liftOver", "log")
        # 处理路径和文件的方法
        self.file = yf.staticMethod(log_file="log")
        self.read = yf.Read(header=None, log_file="log")
        self.create = yf.Create(header=False, log_file="log")
        self.lift_over_path = lift_over_path
        self.is_hg19_to_hg38 = is_hg19_to_hg38
        self.transformation_file: str = "hg19ToHg38.over.chain.gz" if self.is_hg19_to_hg38 else "hg38ToHg19.over.chain.gz"
        self.source = os.path.join(path, "hg19" if self.is_hg19_to_hg38 else "hg38")
        self.output = os.path.join(path, "hg38" if self.is_hg19_to_hg38 else "hg19")
        self.unmap = os.path.join(path, self.transformation_file + "_unmap")
        self.run()

    def exec_str(self, filename: str) -> str:
        return f"{self.lift_over_path}/liftOver {os.path.join(self.source, filename)} {self.lift_over_path}/{self.transformation_file} {os.path.join(self.output, filename)} {os.path.join(self.unmap, filename)}"

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
        input_files = self.file.get_files(path=self.source)
        code_list = []
        finish_files = self.file.get_files(path=self.output)
        for input_file in input_files:
            if input_file not in finish_files:
                code_list.append(self.exec_str(input_file))

        # 实例化线程对象
        pool = Pool(5)
        # 将 list 的每一个元素传递给 pool_page(page) 处理
        pool.map(self.exec_command, code_list)
        # 关闭线程
        pool.close()
