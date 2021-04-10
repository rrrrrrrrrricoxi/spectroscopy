#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from spectrum import Spectrum
import os
from matplotlib import pyplot as plt

# ---------- an utterly unnecessary separator line ---------- #

# specify where the data input and output folders are
input_path = '/Users/ricoxi/Desktop/Spitfire/data_input/'
output_path = '/Users/ricoxi/Desktop/Spitfire/data_output/'

"""
load gold reference, specify folder names
"""

os.chdir(input_path)
# make sure files are in the data_input folder, then specify where the gold refs are
gold = '20180124/endAu'
gold_name = 'test_gold'

gold_spec = Spectrum(gold)
# argument for .removeCRs() is the threshold for removing cosmic rays
gold_spec.removeCRs(200)
gold_spec.subtractBGs()
gold_spec.sumFullDFGs()
gold_spec.plotDFGs()
gold_spec.plotSumDFG()
# argument for .smoothDFGs() is the strength(window size) of smoothing function
gold_spec.smoothDFGs(5)
gold_spec.plotSmoothRawDFGs()
gold_spec.findTruncateIndices()
gold_spec.truncateFullDFGs(gold_spec)
gold_spec.plotTruncatedDFGs()
gold_spec.sumTruncatedDFGs()
gold_spec.plotSumTruncatedDFG()

"""
load sample, specify folder name
"""

os.chdir(input_path)
sample = '20180124/expt/lipidspec'
sample_name = 'test_sample'

sample_spec = Spectrum(sample)
# argument for .removeCRs() is the threshold for removing cosmic rays
sample_spec.removeCRs(50)
sample_spec.plotBGs()
sample_spec.subtractBGs()
sample_spec.sumFullDFGs()
sample_spec.plotDFGs()
sample_spec.plotSumDFG()

sample_spec.findTruncateIndices()
sample_spec.truncateFullDFGs(sample_spec)
sample_spec.plotTruncatedDFGs()
sample_spec.sumTruncatedDFGs()
sample_spec.plotSumTruncatedDFG()

plt.show()

# ---------------------- #
# write gold and sample data, summed and truncated
# ---------------------- #
# specify output directory to the data_output folder
os.chdir(output_path)

if gold_name + '.txt' in os.listdir(output_path):
    pass
else:
    gold_spec.writeSumTruncatedDFG(gold_name)

sample_spec.writeSumTruncatedDFG(sample_name)
sample_spec.writeBGs(sample_name)

print(f"\nFiles are written to {output_path}! (๑•̀ㅂ•́)و✧")

        