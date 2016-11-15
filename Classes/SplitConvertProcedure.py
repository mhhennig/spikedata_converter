# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 14:33:14 2014

@author: szordan
"""


import sys
import os

sys.path.append('..' + os.sep)

import Classes.Utility as Ut
import Classes.Conversion as C


def split(inPath, outPath):
    '''
    ***********************
    SPLIT HDF5 PROCEDURE
    ***********************

    inPath: path directory that contain sorted hdf5(it have to contain Sorted.hdf5 in file name )
    outPath: path directory where to save splitted hdf5 file

    '''
    # Search for global Hdf5 that has to be splitted(It has to contain string
    # Sorted.hdf5)

    hdf5Sorted = Ut.fileListFind(inPath, 'Sorted.hdf5')

    # Search for Text file that contains the names of files and the number of
    # recorded frame
    txtSorted = Ut.fileListFind(inPath, 'Sorted.txt')

    # For every global HDF5 sorted file and txt file split and create
    # temporaneal hdf5 file

    for i in range(0, len(hdf5Sorted)):

        # This class split HDF5 creating one hdf5 for file in inPath

        C.SplitSortHdf5(hdf5Sorted[i], txtSorted[i], outPath)


def convert(inPath, outPath, rootStimuliPath, pathToBeReplace):
    '''
    ******************************************
    CONVERT HDF5 PROCEDURE FOR VISUAL PROTOCOL
    ******************************************

    inPath: path directory that contain splitted hdf5, report text file and Trigger mat file

    outPath: path directory where to save final hdf5 file

    rootStimuliPath: directory that contain all the stimuli Image

    pathToBeReplace: path in report file that has to be replace with root stimuli path

    '''

    # Take the list of file created after splitting
    spikeTrainList = Ut.fileListFind(inPath, 'SortedRearrange')
    print(inPath)
    print(spikeTrainList)

    # Take the list of report File
    ReportList = Ut.fileListFind(inPath, 'Report')
    print(ReportList)

    # Take the list of Trigger file
    TriggerList = Ut.fileListFind(inPath, 'Trigger')
    print(TriggerList)

    # List that will contain the converted files
    convertedFile = []

    # For all file create one list that contain the tree type of file for each phase and
    # generate Hdf5 file

    for i in range(0, len(spikeTrainList)):

        inFileList = [TriggerList[i], ReportList[i],
                      spikeTrainList[i], rootStimuliPath, pathToBeReplace]

        ConversionData = C.ConversionVisualProtocol(inFileList, outPath)

        convertedFile.append(ConversionData.outFilePath)
        ConversionData = 0
    return convertedFile

def convertStimuliOnly(inPath, outPath, rootStimuliPath, pathToBeReplace):
    '''
    ******************************************
    CONVERT HDF5 PROCEDURE FOR VISUAL PROTOCOL
    ******************************************

    inPath: path directory that contain splitted hdf5, report text file and Trigger mat file

    outPath: path directory where to save final hdf5 file

    rootStimuliPath: directory that contain all the stimuli Image

    pathToBeReplace: path in report file that has to be replace with root stimuli path

    '''

    # Take the list of report File
    ReportList = Ut.fileListFind(inPath, 'Report')
    print(ReportList)

    # Take the list of Trigger file
    TriggerList = Ut.fileListFind(inPath, 'Trigger')
    print(TriggerList)

    # List that will contain the converted files
    convertedFile = []

    # For all file create one list that contain the tree type of file for each phase and
    # generate Hdf5 file

    for i in range(0, len(ReportList)):

        inFileList = [TriggerList[i], ReportList[i],
                      [], rootStimuliPath, pathToBeReplace]

        ConversionData = C.ConversionVisualProtocol(inFileList, outPath)

        convertedFile.append(ConversionData.outFilePath)
        ConversionData = 0
    return convertedFile
