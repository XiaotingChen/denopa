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
    version="1.0.3",
    author="Bingxiang Xu",
    author_email="xubingxiang@sus.edu.cn",
    maintainer="Zhihua Zhang",
    maintainer_email="zhangzhihua@big.ac.cn",
    url="https://gitee.com/bxxu/denopa",
    description=
    "Decoding the nucleosome positions with ATAC-seq data at single cell level",
    long_description=long_description,
    scripts=["scripts/denopa"],
    python_requires=">=2.7",
    packages=['deNOPA'],
    install_requires=[
        "numpy", "scipy", "pandas", "h5py", "pysam", "scikit-learn", "statsmodels"
    ])
