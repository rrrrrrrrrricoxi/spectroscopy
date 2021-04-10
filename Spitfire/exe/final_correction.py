#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes in a summed DFG, a summed gold, and a polystyrene measurement.
DFG counts are normalized against the gold counts.
DFG wavenumbers are calibrated to PS 2850.3 and 3060.7 cm-1 dips.

@author: ricoxi
    modeled based on @pohno's Solstice script
"""

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema

class Correction():

    def camera_to_vibfreq(visWL, cameraWL):
        visWN = 10 ** 7 / visWL
        cameraWN = 10 ** 7 / cameraWL
        vib_wn = cameraWN - visWN
        return vib_wn

    def __init__(self, path):
        # read in data output from "import_data_here"
        for i in os.listdir():
            if 'gold' in i.lower() or 'au' in i.lower():
                gold_wn = []
                gold_counts =[]
                with open(i, 'r') as f:
                    for line in f.readlines()[1:]:
                        field = line.rstrip('\n').split(',')
                        gold_wn.append(field[0])
                        gold_counts.append(field[1])
                    for i in np.arange(len(gold_wn)):
                        gold_wn[i] = float(gold_wn[i])
                        gold_counts[i] = float(gold_counts[i])

                    self.gold_wn = np.asarray(gold_wn)
                    self.gold_counts = np.asarray(gold_counts)

            elif sample_name in i.lower() and 'bg' not in i.lower():
                sample_wn = []
                sample_counts =[]
                with open(i, 'r') as f:
                    for line in f.readlines()[1:]:
                        field = line.rstrip('\n').split(',')
                        sample_wn.append(field[0])
                        sample_counts.append(field[1])
                    for i in np.arange(len(sample_wn)):
                        sample_wn[i] = float(sample_wn[i])
                        sample_counts[i] = float(sample_counts[i])

                    self.sample_wn = np.asarray(sample_wn)
                    self.sample_counts = np.asarray(sample_counts)

            elif 'polystyrene' in i.lower() or 'ps' in i.lower():
                ps_counts = []
                with open(i, 'r') as f:
                    for line in f.readlines():
                        field = line.rstrip('\n').split()
                        ps_counts.append(field[3])
                    for i in np.arange(len(ps_counts)):
                        ps_counts[i] = float(ps_counts[i])

        self.ps_wn = np.asarray(gold_wn)
        self.ps_counts = np.asarray(ps_counts)

    # function to normalize sample intensity to gold reference
    def norm_counts(self):
        gold_max = np.max(self.gold_counts)
        norm_gold = []
        for i in np.arange(len(self.gold_counts)):
            norm_gold.append(self.gold_counts[i] / gold_max)
        self.norm_gold = np.asarray(norm_gold)

        norm_sample = []
        for i in np.arange(len(self.sample_counts)):
            norm_sample.append(self.norm_gold[i] * self.sample_counts[i])

        self.norm_sample = np.asarray(norm_sample)

    def plot_crt_counts(self):
        print('-' * 100)
        print("Figure 1: Normalizing sample counts to gold... ᕕ( ᐛ )ᕗ")
        plt.figure()
        plt.plot(self.sample_wn, self.norm_sample)
        plt.title('Sample Spectrum Normalized to Gold')

    def plot_ps(self):
        print('-' * 100)
        print('Figure 2: Plotting uncalibrated polystyrene spectrum... w(ﾟДﾟ)w')
        plt.figure()
        plt.plot(self.ps_wn, self.ps_counts)
        plt.title('Uncalibrated Polystyrene Spectrum')

    # function to shift wavenumber according to polystyrene
    def plot_min_wn(self):
        print('-' * 100)
        print('Figure 3: Plotting uncalibrated PS local minima... (ﾟｰﾟ)')
        plt.figure()
        plt.plot(self.ps_wn, self.ps_counts)
        plt.title("PS Local Minima")
        loc_min = argrelextrema(self.ps_counts, np.less, order=5)

        wn_candidate = []
        for i in loc_min[0]:
            if self.ps_wn[i] > 2800 and self.ps_wn[i] < 3200:
                plt.plot(self.ps_wn[i], self.ps_counts[i], 'o', markersize=3)
                wn_candidate.append(self.ps_wn[i])

        """
        Need to judge from the plot which local minima to pick
        """
        loc1 = 4
        loc2 = 1
        self.ps_pof1 = wn_candidate[loc1]
        self.ps_pof2 = wn_candidate[loc2]
        print('-' * 100)
        print(f"Point of references found at {self.ps_pof1} cm\u207B\u00B9 "
              f"and {self.ps_pof2} cm\u207B\u00B9! <(ˉ^ˉ)>")
        print('-' * 100)

        # target_loc_1 = 2850.3
        # window = 0.02
        # up_lim = round(target_loc_1 * (1 + window))
        # low_lim = round(target_loc_1 * (1 - window))
        # for i in loc_min[0]:
        #     if self.ps_wn[i] > low_lim and self.ps_wn[i] < up_lim:
        #         print(f"Point of reference found at {self.ps_wn[i]} cm\u207B\u00B9! <(ˉ^ˉ)>")
        #         pof1 = self.ps_wn[i]
        # self.ps_pof1 = pof1

    def pscal_wn(self):
        ref_wn1 = 2850.3
        ref_wn2 = 3060.7
        shift = ((self.ps_pof1 - ref_wn1) + (self.ps_pof2 - ref_wn2)) / 2
        print(f"Shift was {shift} cm\u207B\u00B9. o(￣▽￣)ｄ")

        self.sample_wn_cal = self.sample_wn - shift


    def plot_cal_sample(self):
        print('-' * 100)
        print("Figure 4: Plotting normalized and calibrated sample spectrum... ヽ(✿ﾟ▽ﾟ)ノ")
        print('-' * 100)
        plt.figure()
        plt.plot(self.sample_wn, self.norm_sample, label = "raw")
        plt.plot(self.sample_wn_cal, self.norm_sample, label = "calibrated")
        plt.legend()
        plt.title("Normalized and Calibrated Sample Spectrum")

    def write_cal_spec(self, name):
        print(f'Writing normalized and calibrated spectrum to {name}.txt... (๑•̀ㅂ•́)و✧')
        data = np.vstack((self.sample_wn_cal, self.norm_sample)).transpose()
        fmt = '%.5f'
        np.savetxt(name + '.txt', data, fmt, delimiter=',', header='wn,counts', comments='')
        print('-' * 100)
        print('Franz is happy')
        print('+1 exp\nacquired item <sleepless night> * 10'
              '\nTech basement temperature just dropped by another 20 \u00b0F'
              '\nLife is tough, but do carry on! Show \'em who\'s da boss d=====(￣▽￣*)b')






output_path = '/Users/ricoxi/Desktop/Spitfire/data_output'
os.chdir(output_path)

# make sure you change this to the name you assigned you summed DFG file
sample_name = 'test_sample'

corr = Correction(output_path)
# print(corr.sample_wn)
corr.norm_counts()
# print(corr.norm_gold)
# print(corr.corr_sample)
corr.plot_crt_counts()
corr.plot_ps()
corr.plot_min_wn()
corr.pscal_wn()
corr.plot_cal_sample()
# print(corr.sample_wn)
# print(corr.sample_wn_cal)
# print(corr.ps_pof1)
# print(corr.ps_wn)
# print(corr.ps_counts)
plt.show()
cal_name = 'verytest'
corr.write_cal_spec(cal_name)