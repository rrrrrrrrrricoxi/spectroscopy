#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# importing packages
import os
from acq_dfg import DFG
from fullwn import FullWN
from matplotlib import pyplot as plt
import pandas
import numpy as np
import math
from scipy.ndimage.filters import gaussian_filter1d
import copy

# load wavenumber array

fullwn = FullWN()

class Spectrum():
    def __init__(self, path):

        print('Importing DFGs and BGs from the following files:')

        os.chdir(path)

        self.dfgs = []
        self.bgs = []

        for i in os.listdir():
            if '.txt' in i:
                name = i.split('.')[0]
                print(name)
                if 'bg' in name:
                    self.bgs = self.bgs + [DFG(i, name)]
                else:
                    self.dfgs = self.dfgs + [DFG(i, name)]
        # files sorted by name
        self.dfgs.sort(key=lambda x: x.name)
        self.bgs.sort(key=lambda x: x.name)

        self.fullwn = fullwn.fullwn

    # plot all sample DFGs
    def plotDFGs(self):
        plt.figure()
        for dfg in self.dfgs:
            plt.plot(dfg.wn, dfg.counts)
        plt.title('DFGs')

    # plot all background DFGs
    def plotBGs(self):
        plt.figure()
        for bg in self.bgs:
            plt.plot(bg.wn, bg.counts)
        plt.title('BGs')

    # plot each sample DFG individually
    def plotIndDFGs(self):
        for dfg in self.dfgs:
            plt.figure()
            plt.plot(dfg.wn, dfg.counts)
            plt.title(dfg.name)

    # plot each background DFG individually
    def plotIndBGs(self):
        for bg in self.bgs:
            plt.figure()
            plt.plot(bg.wn, bg.counts)
            plt.title(bg.name)

    def plotSumDFG(self):
        plt.figure()
        plt.plot(self.fullwn,self.dfgSum)
        plt.title('Sum of DFGs')

    # plot the smoothed and the raw padded DFGs against fullwn
    def plotSmoothRawDFGs(self):
        plt.figure()
        for dfg in self.dfgsPreSmoothed:
            plt.plot(self.fullwn, dfg.counts, 'ro')

        for dfg in self.dfgsFull:
            plt.plot(self.fullwn, dfg.counts, 'b')
        plt.title('Smoothed and Raw DFGs')

    # plot the DFGs that have been truncated according to the gold reference
    def plotTruncatedDFGs(self):
        plt.figure()
        for dfg in self.dfgsFullTruncated:
            plt.plot(self.fullwn, dfg.counts)
        plt.title('Truncated DFGs')

    # plot the sum of the truncated DFGs
    def plotSumTruncatedDFG(self):
        plt.figure()
        plt.plot(self.fullwn, self.dfgTruncatedSum)
        plt.title('Sum of truncated DFGs')

    # remove all cosmic rays for each DFG and BG
    def removeCRs(self, threshold=200):
        print('Removing cosmic rays from spectra...')

        # function that uses a median filter to identify a CR in a single DFG
        def removeCRindDFG(dfg, threshold):
            # choose how big of a window for the rolling median
            windowSize = 7
            medians = pandas.Series(dfg.counts).rolling(window=windowSize, center=True).median()

            # number of nan at beginning and end to replace
            numRep = math.floor(windowSize / 2)

            # replace beginning and end nan with the first/last computed value
            for i in range(numRep):
                medians[i] = medians[numRep]
                medians[len(medians) - i - 1] = medians[len(medians) - numRep - 1]

            # find difference of each point with the median of its window
            differences = dfg.counts - medians

            # empty array to hold zero or one if point is a spike
            spike = np.zeros(len(differences), )
            for i in range(len(differences)):
                if differences[i] > threshold:
                    spike[i] = 1
                    print("Spike found at point index", i, "with wavenumber", dfg.wn[i], "cm\u207B\u00B9")

            # if there any spikes found
            if np.sum(spike) > 0:
                # go through and replace the spike with the average on both sides
                for i in range(len(spike)):
                    # if the point needs to be replaced
                    if spike[i] == 1:
                        # check up to five points to the left for the edge or for an ok point
                        for j in range(5):
                            if (i - 1 - j) < 0:
                                left = []  # if its edge only take from right point
                                break
                            else:
                                if spike[i - 1 - j] == 0:
                                    left = dfg.counts[i - 1 - j]  # or get the first acceptable point
                                    break
                                    # check up to five points to the right for the edge or for an ok point
                        for j in range(5):
                            if (i + j + 1) >= len(spike):
                                right = []  # if its edge only take from the left point
                                break
                            else:
                                if spike[i + 1 + j] == 0:
                                    right = dfg.counts[i + 1 + j]  # or get the first acceptable point
                                    break
                                    # get the average of the two or the value if its only one
                        tempValArray = np.array([])
                        tempValArray = np.append(tempValArray, left)
                        tempValArray = np.append(tempValArray, right)
                        ave = tempValArray.mean()

                        # round down to integer number of counts
                        dfg.counts[i] = math.floor(ave)
            else:
                print("No spikes found in " + dfg.name)
            return

        # go through each dfg and remove CRs
        for dfg in self.dfgs:
            removeCRindDFG(dfg, threshold)

        # go through each bg and remove CRs
        for bg in self.bgs:
            removeCRindDFG(bg, threshold)

    # find and subtract correct background
    def subtractBGs(self):

        print('Subtracting BGs from DFGs...')

        # create list to hold pre-bg subtracted dfgs
        self.dfgsRaw = copy.deepcopy(self.dfgs)

        # go through each dfg
        for dfg in self.dfgs:
            # identify background by finding median wavelength
            dfgMedian = int(np.median(dfg.wl))

            # tracker for seeing if you found background
            foundBG = False

            # go through each background, see if one with matching median is there
            for bg in self.bgs:
                if dfgMedian == int(np.median(bg.wl)):
                    print("For dfg", dfg.name, "found", bg.name)
                    dfg.counts = dfg.counts - bg.counts
                    foundBG = True

            # if one wasn't found, print that
            if not foundBG:
                print("No bg found for dfg", dfg.name)

        self.dfgsFull = copy.deepcopy(self.dfgs)
        
    """
    Spitfire data does not need padding since we don't really move the camera
    position for the C-H region
    """
    
    # # pad each DFG with zeros before and/or after so they align and can be summed up
    # def padDFGs(self):
    #
    #     print('Padding DFGs with Zeros...')
    #     # dictionary to hold number of zeros to pad on either side
    #     padding = dict(det615=[0, 467], det620=[58, 409], det625=[116, 351], det630=[174, 293],
    #                    det635=[232, 235], det640=[290, 177], det645=[349, 118], det650=[407, 60],
    #                    det655=[467, 0])
    #
    #     # length of fullwn is 911
    #
    #     # for 615 add 467 after
    #     # for 620 add 58 before and 409 after
    #     # for 625 add 116 before and 351 after
    #     # for 630 add 174 before and 293 after
    #     # for 635 add 232 before and 235 after
    #     # for 640 add 290 before adn 177 after
    #     # for 645 add 349 before adn 118 after
    #     # for 655 add 467 before
    #
    #     # copy dfgs into new list
    #     self.dfgsFull = copy.deepcopy(self.dfgs)
    #
    #     for dfg in self.dfgsFull:
    #         key = 'det' + str(int(np.median(dfg.wl)))
    #         dfg.counts = np.append(np.append(np.zeros(padding[key][0]), dfg.counts),
    #                                np.zeros(padding[key][1]))

    # sum up the padded DFGs

    def sumFullDFGs(self):
        print('Summing full DFGs...')
        self.dfgSum = np.zeros(len(self.fullwn))
        # for dfg in self.dfgs:
        for dfg in self.dfgsFull:
            self.dfgSum = self.dfgSum + dfg.counts

    # smooth each DFG with a Gaussian window
    def smoothDFGs(self, sigma=5):
        print('Smoothing DFGs...')
        self.dfgsPreSmoothed = copy.deepcopy(self.dfgsFull)

        for dfg in self.dfgs:
            # use gaussian filter imported from scipy.ndimage.filters
            dfg.counts = gaussian_filter1d(dfg.counts, sigma)

    # find indices of reference spectra where the signal falls off to 5% of max
    def findTruncateIndices(self, threshold=0.05):
        print('Finding truncation thresholds at', threshold, '...')
        # create list to hold indices
        self.truncateIndices = []

        # go through each dfg
        for dfg in self.dfgsFull:
            # find max
            maxVal = dfg.counts.max()

            # find index of the max
            maxIndex = dfg.counts.argmax()

            # find left and right indexes
            leftIndex = (np.abs(dfg.counts[:maxIndex] - maxVal * threshold)).argmin()
            rightIndex = maxIndex + (np.abs(dfg.counts[maxIndex:] - maxVal * threshold)).argmin()

            # add the found values to the list
            self.truncateIndices = self.truncateIndices + [[leftIndex, rightIndex]]

    # truncate padded DFGs according to indices set according to gold spectrum
    def truncateFullDFGs(self, gold):
        print('Truncating DFGs...')
        # copy dfgs into new list
        self.dfgsFullTruncated = copy.deepcopy(self.dfgsFull)

        # set all the values equal to zero not within the indices determined by
        # the gold reference spectrum
        for i, dfg in enumerate(self.dfgsFullTruncated):
            dfg.counts[:gold.truncateIndices[i][0]] = 0
            dfg.counts[gold.truncateIndices[i][1]:] = 0

    # sum up these truncated DFGs
    def sumTruncatedDFGs(self):
        print('Summing truncated DFGs...')
        self.dfgTruncatedSum = np.zeros(len(self.fullwn))

        for dfg in self.dfgsFullTruncated:
            self.dfgTruncatedSum = self.dfgTruncatedSum + dfg.counts

    # WRITING METHODS

    # write each individual sample DFG to file
    def writeDFGs(self, name):
        data = np.zeros(444)
        for dfg in self.dfgs:
            data = np.vstack((data, dfg.wn))
            data = np.vstack((data, dfg.counts))
        data = data.transpose()
        fmt = '%.5f'
        np.savetxt(name, data, fmt, delimiter=',')

    # write each padded DFG to file
    def writeFullDFGs(self, name):
        data = self.fullwn
        for dfg in self.dfgsFull:
            data = np.vstack((data, dfg.counts))
        data = data.transpose()
        fmt = '%.5f'
        np.savetxt(name, data, fmt, delimiter=',')

    # write the sum of the padded DFGs to file
    def writeSumDFG(self, name):
        data = np.vstack((self.fullwn, self.dfgSum))
        data = data.transpose()
        fmt = '%.5f'
        np.savetxt(name, data, fmt, delimiter=',')

    # write the smoothed DFGs to file
    def writeSmoothedDFGs(self, name):
        data = self.fullwn
        for dfg in self.dfgsFull:
            data = np.vstack((data, dfg.counts))
        data = data.transpose()
        fmt = '%.5f'
        np.savetxt(name, data, fmt, delimiter=',')

    # write the truncated DFGs to file
    def writeTruncatedDFGs(self, name):
        data = self.fullwn
        for dfg in self.dfgsFullTruncated:
            data = np.vstack((data, dfg.counts))
        data = data.transpose()
        fmt = '%.5f'
        np.savetxt(name, data, fmt, delimiter=',')

    # write the sum of the truncated DFGs to file
    def writeSumTruncatedDFG(self, name):
        print('Truncated, summed wave written to', name + '.txt')
        data = np.vstack((self.fullwn, self.dfgTruncatedSum)).transpose()
        fmt = '%.5f'
        np.savetxt(name + '.txt', data, fmt, delimiter=',', header='wn,counts', comments='')
        
    # write background to file
    def writeBGs(self, name):
        print('Background written to', name + '_bgs' + '.txt')
        bg_counts = self.bgs[0].counts
        data = np.vstack((self.fullwn, bg_counts)).transpose()
        fmt = '%.5f'
        np.savetxt(name + '_bgs' + '.txt', data, fmt, delimiter=',', header='wn,counts', comments='')


