#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 5 2020

Classs to hold data from a single or multiple frames of DFG acquisition.

"DFG stands for 'difference frequency generation'. When collecting a broad SFG
vibrational spectrum, the IR wavelength and detector wavelength need to be
shifted across frequency space as they are too narrow to do everything at one
time. The IR wavelength is shifted by changing the position of something called
the DFG crystal, hence the name of a single one of these acquisitions."
                                                                    --pohno

The name property (self.name) holds the name of the file.
The wl property (self.wl) is an array of the wavelength values that are read in.
The counts property (self.counts) is the number of photons that are read in.
The wn property (self.wn) is the wl array converted to wavelengths.

@author: Rico
    modeled based on @pohno's Solstice script
"""

import os
import numpy as np

# ---------- an utterly unnecessary separator line ---------- #

class DFG():

    def __init__(self, path, name):

        self.name = name
        # define a function to convert wavelength into wavenumbers
        def camera_to_vibfreq(visWL, cameraWL):
            visWN = 10 ** 7 / visWL
            cameraWN = 10 ** 7 / cameraWL
            vib_wn = cameraWN - visWN
            return vib_wn

        viswl_Spitfire = 800
        fullwn_list = []

        # initialize an empty list to store data
        data = []

        # ------------------------------ #
        # raw data has only one frame
        # ------------------------------ #

        try:
            with open(path, 'r') as f:
                wl = []
                counts = []
                for line in f.readlines():
                    field = line.split()
                    wl.append(field[0])
                    counts.append(field[-1])

                for i in np.arange(len(wl)):
                    wl[i] = float(wl[i])
                    counts[i] = float(counts[i])
                for i in wl:
                    wn = camera_to_vibfreq(viswl_Spitfire, i)
                    fullwn_list.append(wn)
                data.append(fullwn_list)
                data.append(counts)

                self.wl = np.asarray(wl)
                self.wn = np.asarray(fullwn_list)
                self.counts = np.asarray(counts)

        # ------------------------------ #
        # raw data with multiple frames
        # ------------------------------ #

        except ValueError:

            with open(path, 'r') as f:
                rawdata_str = f.readlines()[1:]

                for i in rawdata_str:
                    rawline = i.split('\t')
                    float_line = []

                    for j in rawline:
                        try:
                            float_line.append(float(j))
                        except ValueError:
                            pass

                    if len(float_line) != 0:
                        data.append(float_line)

            # remove frame and strip numbers from the data list
            data_clean = []
            for i in data:
                if i[0] == data[0][0]:
                    data_clean.append(i)
                else:
                    data_clean.append(i[2:])

            wl = data_clean[0]
            for i in data_clean[0]:
                wn = camera_to_vibfreq(viswl_Spitfire, i)
                fullwn_list.append(wn)

            # store wl to wn conversion in data_clean array
            data_clean[0] = fullwn_list

            data_dict = {}
            title = ['wn']
            for i in np.arange(len(data_clean) - 1):
                title.append('frame' + str(i + 1))
            for i in np.arange(len(title)):
                data_dict[title[i]] = data_clean[i]

            # average counts across 5 frames
            count_data = []

            for name in data_dict:
                if name != 'wn':
                    counts = data_dict[name]
                    count_data.append(counts)

            count_data = np.asarray(count_data)
            avg_counts = np.mean(count_data, axis=0)

            self.wl = np.asarray(wl)
            self.wn = np.asarray(fullwn_list)
            self.counts = np.asarray(avg_counts)




