# -*- coding: utf-8 -*-
# @Time    : 18-5-12 下午1:47
# @Author  : Matrix
# @Site    :
# @File    : setup.py
# @Software: PyCharm

from setuptools import setup

with open("README.md") as fin:
    long_description = fin.read()

setup(
    name="deNOPA",
    version="1.0.0",
    author="Bingxiang Xu",
    author_email="xubx@big.ac.cn",
    maintainer="Zhihua Zhang",
    maintainer_email="zhangzhihua@big.ac.cn",
    url="https://gitee.com/bxxu/denopa",
    description=
    "Decoding the nucleosome positions with ATAC-seq data at single cell level",
    long_description=long_description,
    scripts=["bin/denopa"],
    python_requires="==2.7",
    packages=['deNOPA'],
    install_requires=[
        "numpy>=1.15.4", "scipy>=1.1.0", "pandas>=0.23.4", "h5py>=2.8.0",
        "pysam>=0.16", "sklearn", "statsmodels"
    ])
