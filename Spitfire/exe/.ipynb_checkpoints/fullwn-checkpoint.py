#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 4 2020

Object that is just an array that holds all the wavenumber values that the
acquisitions, once they are padded with zeros and/or properly summed up,
are plotted against.

If we were to change the camera range on spitfire in the future,
substitute the file with a gold sample file taken with the new camera.
Do so by assigning the new file as standard.

@author: Rico
    modeled based on @pohno's Solstice script
"""

#import numpy
import numpy as np
import os

# ---------- an utterly unnecessary separator line ---------- #

# Calculate the range of vibration frequencies in wavenumbers from any
# Spitfire files
# use the desired file as the standard object

# standard = "data_input/Test_PS_1700_30sec.txt"

os.chdir("/Users/ricoxi/Desktop/Spitfire")
standard = 'data_input/20180124/endAu/endau_bg_4sec.txt'

class wn_standard():
    def __init__(self):
        self.std = standard

std = wn_standard().std

class FullWN():

    def __init__(self):
        # if we were to change the camera range on spitfire in the future,
        # substitute the file with a gold sample file taken with the new camera
        # Change directory to the folder that stores experimental data

        data = []
        with open(std, 'r') as f:
            for i in f.readlines():
                line = i.rstrip('\n').split('\t')
                data.append(line)

        data_final = []
        for i in data:
            data_set = []
            for j in i:
                j = float(j)
                data_set.append(j)
            data_final.append(data_set)

        def camera_to_vibfreq(visWL, cameraWL):
            visWN = 10 ** 7 / visWL
            cameraWN = 10 ** 7 / cameraWL
            vib_wn = cameraWN - visWN

            return vib_wn

        viswl_Spitfire = 800

        fullwn_list = []

        for i in data_final:
            wn = camera_to_vibfreq(viswl_Spitfire, i[0])
            fullwn_list.append(wn)
        # numpy array is a little better to see
        self.fullwn = np.asarray(fullwn_list)


