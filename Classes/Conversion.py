# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 17:51:22 2014

@author: szordan
"""
import scipy.io
import io
import sys
import numpy as np
import h5py
import os

sys.path.append('..' + os.sep)

import Classes.Utility as Ut
import Classes.H5Tables as H5Tables
import operator


class Conversion():

    def __init__(self, inpath, outpath):
        self.inpath = inpath
        self.outpath = outpath

    def loadData():
        return 0

    def ProcessData():
        return 0

    def DatatoHdf5():
        return 0


class SplitSortHdf5():

    def __init__(self, hdf5file, textfile, outpath):

        self.textfile = textfile

        # Create a dictionary taking the list of file and the number of
        # recording frame from Sorted.txt

        self.Dfile = self.takeListNameFile()

        # Read HDF5 file

        self.hdf5SortUnit = h5py.File(hdf5file, 'r')

        self.parameter = {}
        self.timeStampList = []
        self.DHdf5Matrix = {}
        self.globalParameter = {}

        self.loadNeuronSpikeSortedData()
        self.rearrangeData()

        oldIndex = 0
        recordingFrame = 0
        PreviousRecordingFrame = 0

        fileNameList = self.Dfile.keys()
        fileNameList.sort()
        for filenamepath in fileNameList:

            # Initialize startIndex as previous stopIndex

            startIndex = oldIndex
            print ("Start index: " + str(startIndex))

            # Take the number of frame recorded
            recordingFrame += int(self.Dfile[filenamepath])

            # Save some recording parameter
            self.parameter['RecFrame'] = int(self.Dfile[filenamepath])
            self.globalParameter['RecordingParameter'] = [np.asarray([str(filenamepath)]), np.asarray(
                [0]), np.asarray([self.parameter['RecFrame']]), np.asarray([self.parameter['Freq']])]
            # self.globalParameter['ExpPathRec']=str(filenamepath)
            print ("N of frame rec: " + str(recordingFrame))

            # find index of spike frame before recordingFrame
            stopIndex = self.findIndex(recordingFrame)
            print ("Stop Index: " + str(stopIndex))

            # create TimeStampMatrix
            self.createTimeStampArray(
                startIndex, stopIndex, PreviousRecordingFrame)

            # Create file name out
            filename = filenamepath.replace('/', '\\')
            out = Ut.buildPathHdf5(filename, outpath, 'SortedRearrange')
            print ("OutFilePath: " + out)
            # Export hdf5
            self.exportToHdf5(out)

            # Save and add number of frame that has to be subtract to subsequent file
            # for deconcatenate

            PreviousRecordingFrame += int(self.Dfile[filenamepath])
            oldIndex = stopIndex

    def findIndex(self, recordingFrame):
        index = len(self.Times[self.Times < recordingFrame])
        return index

    # Read file text that contain the name of file concatenated and the number
    # of recording frame

    def takeListNameFile(self):
        Dfile = {}
        txt = open(self.textfile, 'r')
        lines = txt.readlines()
        for i in range(0, len(lines)):
            if lines[i][0:6] == '# File':
                line = lines[i + 1].split(',')
                Dfile[line[0].replace("'", "")] = line[1]
        txt.close()
        return Dfile

    # Load neuron data from hdf5 produce after sorting processing chain
    def loadNeuronSpikeSortedData(self):

        # Load structure from the hdf5 created after sorting

        print 'load sorted unit'

        self.sortedUnits = np.asarray(self.hdf5SortUnit['SortedUnits'])
        self.units = np.unique(self.sortedUnits)

        self.nunits = len(self.units)

        print 'times'
        self.Times = np.asarray(self.hdf5SortUnit['Times'])

        print 'MarkedChannel'
        self.MarkedChannel = np.asarray(self.hdf5SortUnit['MarkedChannels'])
        self.indexMarkedChannel = self.MarkedChannel.nonzero()[0]

        self.markedUnit = np.unique(self.sortedUnits[self.indexMarkedChannel])
        print np.asarray(self.hdf5SortUnit['Sampling'])
        self.parameter['Freq'] = np.asarray(self.hdf5SortUnit['Sampling'])

        print 'load location'
        self.Locations = np.asarray(self.hdf5SortUnit['Locations'])
        print len(self.Locations)
        self.hdf5SortUnit.close()

    def rearrangeData(self):
        print "rearrangeData"
        DchName = {}

        NeuronMatrix = [[], [], [], [], [], []]

        # Build of NeuronMatrix
        for i in range(0, self.nunits):
            position = [self.Locations[self.units[i], 0],
                        self.Locations[self.units[i], 1]]

            if position[0] < 10 and position[1] < 10:
                ChName = 'Ch' + '0' + \
                    str(position[0]) + '_' + '0' + str(position[1])
            elif position[0] < 10:
                ChName = 'Ch' + '0' + str(position[0]) + '_' + str(position[1])
            elif position[1] < 10:
                ChName = 'Ch' + str(position[0]) + '_' + '0' + str(position[1])
            else:
                ChName = 'Ch' + str(position[0]) + '_' + str(position[1])

            if ChName in DchName:
                DchName[ChName] += 1

            else:
                DchName[ChName] = 1

            # Assigned the name with position and unit number

            name = ChName + '_Unit' + str(DchName[ChName])

            NeuronMatrix[0].append(i + 1)
            NeuronMatrix[1].append(name)
            NeuronMatrix[2].append(position[0])
            NeuronMatrix[3].append(position[1])
            NeuronMatrix[4].append(DchName[ChName])

            if self.units[i] in self.markedUnit:

                NeuronMatrix[5].append(1)
            else:
                NeuronMatrix[5].append(0)

        # sortName=np.argsort(NeuronMatrix[1])

        #NeuronMatrix[1]=[NeuronMatrix[1][i] for i in sortName]
        #NeuronMatrix[2]=[NeuronMatrix[2][i] for i in sortName]
        #NeuronMatrix[3]=[NeuronMatrix[3][i] for i in sortName]
        #NeuronMatrix[4]=[NeuronMatrix[4][i] for i in sortName]
        self.DHdf5Matrix['NeuronMatrix'] = NeuronMatrix
        # self.timeStampList=[self.Times,range(0,len(self.Times)),self.sortedUnits,np.zeros(len(self.Times))]

    def createTimeStampArray(self, frameStart, frameStop, PreviousRecordingFrame):
        attributes = ["Time Stamp", "idEvent",
                      "idNeuron_isStim", "AmplitudeEvent"]

        # a=self.timeStampList[0][frameStart:frameStop]-PreviousRecordingFrame
        a = self.Times[frameStart:frameStop] - PreviousRecordingFrame
        b = self.sortedUnits[frameStart:frameStop] + 1
        c = np.zeros(len(a))
        # d=np.asarray(self.timeStampList[3][frameStart:frameStop]).T
        timeStampArray = [a, b, c]

        self.DHdf5Matrix['TimeStampMatrix'] = timeStampArray

    def exportToHdf5(self, outfilePath):

        H5file = H5Tables.H5(outfilePath, 'w')
        H5file.addStructures(self.DHdf5Matrix, 'r')
        H5file.addStructures(self.globalParameter, 'GlobalParameter')
        H5file.close()

'''
##############################################################################

Conversion class without stimuli protocol.

##############################################################################
'''


class ConversionNoStimuliProtocol():

    def __init__(self, inputFilePath, outputFilePath):

        # Folder that contain neuron data(mat file or hdf5 file obtain after
        # sorting chain)

        self.inFile = inputFilePath

        # Folder where the hdf5 file will be save

        self.outFilePath = Ut.buildPathHdf5(
            inputFilePath, outputFilePath, 'Final')

        # Load Neuron data initialize the correct class
        if inputFilePath.find('mat') != -1:
            self.NeuronData = NeuronData(inputFilePath)
        else:
            self.NeuronData = NeuronDataSort(inputFilePath)
        # Process data creating matrix structure of time stamp matrix

        self.processData()

        # Save the hdf5 file

        self.exportToHdf5()

    # Rearrange the dictionary structure of neuron data in matrix structure

    def processData(self):

        self.NeuronData.rearrangeData()

        # self.NeuronData.createTimeStampArray()

    # Export the matrixes structure in the hdf5 file

    def exportToHdf5(self):
        print self.outFilePath
        H5file = H5Tables.H5(self.outFilePath, "w")
        H5file.addStructures(self.NeuronData.DHdf5Matrix, 'r')
        H5file.addStructures(
            self.NeuronData.globalParameter, 'GlobalParameter')
        H5file.close()

'''
##############################################################################

Conversion class with stimuli protocol.
It use the two class NeuronData and StimulusData for import information
from source class, rearrage them and export matrixes structure

##############################################################################
'''


class ConversionVisualProtocol():

    def __init__(self, inputListFile, outputpath):

        self.StimulusData = VisualStimuliData(inputListFile[0], inputListFile[
                                              1], inputListFile[3], inputListFile[4])

        # Control of neuron data source(mat or hdf5 sorted or no spikes)
        if inputListFile[2]:
            self.spikesProvided = True
            if inputListFile[2].find('mat') != -1:
                self.NeuronData = NeuronData(inputListFile[2])
            elif inputListFile[2].find('hdf5') != -1:
                self.NeuronData = NeuronDataSort(inputListFile[2])
            self.outFilePath = Ut.buildPathHdf5(
                inputListFile[0], outputpath, "Final")

        else:
            self.spikesProvided = False
            print("No spike data provided")
            self.outFilePath = Ut.buildPathHdf5(
                inputListFile[0], outputpath, "Stimuli_only")

        self.VisualPath = [inputListFile[0], inputListFile[
            1], inputListFile[3], inputListFile[4]]

        self.DHdf5 = {}

        self.processData()

        self.exportToHdf5()

    def processData(self):

        self.StimulusData.rearrangeData()

        if self.spikesProvided:
            self.NeuronData.rearrangeData()
            self.DHdf5 = dict(self.StimulusData.DHdf5Matrix.items(
            ) + self.NeuronData.DHdf5Matrix.items())
            self.createTimeStampMatrix()
            self.GlobalParameter = dict(self.StimulusData.globalParameter.items(
            ) + self.NeuronData.globalParameter.items())
        else:
            self.DHdf5 = dict(self.StimulusData.DHdf5Matrix.items(
            ))
            self.GlobalParameter = dict(self.StimulusData.globalParameter.items())

        self.GlobalParameter['VersionHdf5'] = [float(1.0)]
        maskpath = self.GlobalParameter.pop('MaskingFilePath', None)
        self.GlobalParameter['ExpPath'] = [np.asarray([self.VisualPath[0]]), np.asarray([self.VisualPath[
            1]]), np.asarray([self.VisualPath[2]]), np.asarray([self.VisualPath[3]]), np.asarray([maskpath])]

    def exportToHdf5(self):
        print ("OutFile Path: " + self.outFilePath)

        H5file = H5Tables.H5(self.outFilePath, "w")
        H5file.addStructures(self.DHdf5, 'r')
        H5file.addStructures(self.GlobalParameter, 'GlobalParameter')
        # H5file.addStructures(self.fakeRoi,'Roi')
        H5file.close()

    def createTimeStampMatrix(self):

        StimulusTimeStampList = self.StimulusData.DHdf5Matrix[
            'TimeStampMatrix']
        NeuronTimeStampList = self.NeuronData.DHdf5Matrix['TimeStampMatrix']

        timeStampArray = np.concatenate(
            (NeuronTimeStampList[0], StimulusTimeStampList[0]), axis=1)
        #idEventArray= np.concatenate((NeuronTimeStampList[1],StimulusTimeStampList[1]*-1),axis=1)
        idNeuronArray = np.concatenate(
            (NeuronTimeStampList[1], StimulusTimeStampList[1]), axis=1)
        amplitudeEvent = np.concatenate(
            (NeuronTimeStampList[2], StimulusTimeStampList[2]), axis=1)

        #self.timeStampMatrix=[(timeStampArray[i],idNeuronArray[i]) for i in np.argsort(timeStampArray)]
        sortarr = np.argsort(timeStampArray)

        timeStampArray = [timeStampArray[i] for i in sortarr]
        idNeuronArray = [idNeuronArray[i] for i in sortarr]

        attributes = ["Time Stamp", "idEvent",
                      "idNeuron_isStim", "AmplitudeEvent"]
        self.DHdf5['TimeStampMatrix'] = [
            timeStampArray, idNeuronArray, amplitudeEvent]


'''
##############################################################################

VISUAL STIMULI DATA CLASS

Class for import stimuli information from
    - Report file (txt)
    - Trigger file (mat)
    - Static stimuli image (jpg or png)

All data will organize in a dictionary structure.
There is methods for load and rearrange data in the hdf5 structure

##############################################################################
'''


class VisualStimuliData():

    def __init__(self, TriggerFilePath, ReportFilePath, RootStimuliPath, PathToBeReplace):

        # Trigger and report file absolute path and path to be replace in report sequence path
        # with the path that contain static stimulus sequence image

        self.TriggerFilePath = TriggerFilePath
        self.ReportFilePath = ReportFilePath
        self.RootStimuliPath = RootStimuliPath
        self.PathToBeReplace = PathToBeReplace

        # Initialize maskingfile string and neutral density filter value

        self.MaskingFile = ""
        self.ND = 0

        '''
        Initialize internal structure
        '''

        # Lsequence is a list that contain a dictionary for every sequence
        self.Lsequence = []

        # Dtrigger is a dictionary that contain time stamp information of
        # stimulus(load from trigger file)
        self.Dtrigger = {}

        # Dictionary that contain the global parameter
        self.globalParameter = {}

        # Dictionary that collect all the final matrixes structure related with
        # the stimulus data
        self.DHdf5Matrix = {}

        # The timeStamp information is save as a list for future elaboration and fusion
        # in a unique matrix integrating time stamp event of neuron

        self.timeStampList = []

        '''
        Load data from trigger file, report file and Information about every sequence


        List of attributes
                Sequence Path
                Num of refresh per stimulus
                Num of stimuli to be displayed
                Num of stimuli displayed
                Masking file
                ND
                Sequence Image List Path
                Modality
                    if ==2 there is alse other two information
                    Sequence Lenght
                    Sequence List
                Number of image
                Duration


        '''
        self.loadReport()
        self.loadSequenceParameter()
        self.loadTimeStampStimuli()

    def loadReport(self):

        # Read report text file
        report = open(self.ReportFilePath)

        # Tag information extract from report file

        tag = [" Executing sequence: ", " Num of refresh per stimulus: ",
               " Num of stimuli to be displayed: ", " Num of stimuli displayed: ", " Masking file: ", " ND: "]
        sequenceIndex = -1

        # For every sequence create dictionary with report information and
        # append to Lsequence list

        for line in report:
            if line[0:len(tag[0])] == tag[0]:
                sequenceIndex += 1
                sequencePath = line[len(tag[0]):len(line)].strip()
                sequencePath = sequencePath.replace(
                    self.PathToBeReplace, self.RootStimuliPath)
                if self.RootStimuliPath[0:6] == '/media':
                    sequencePath = sequencePath.replace('\\', '/')
                self.Lsequence.append({"Sequence Path": sequencePath})

                self.Lsequence[sequenceIndex]['Sequence Index'] = sequenceIndex
            elif line[0:len(tag[1])] == tag[1]:
                self.Lsequence[sequenceIndex][
                    tag[1][1:len(tag[1]) - 2]] = int(line[len(tag[1]):len(line)])
            elif line[0:len(tag[2])] == tag[2]:
                self.Lsequence[sequenceIndex][
                    tag[2][1:len(tag[2]) - 2]] = int(line[len(tag[2]):len(line)])
            elif line[0:len(tag[3])] == tag[3]:
                self.Lsequence[sequenceIndex][
                    tag[3][1:len(tag[3]) - 2]] = int(line[len(tag[3]):len(line)])
            elif line[0:len(tag[4])] == tag[4]:
                self.MaskingFile = line[len(tag[4]):len(line)].strip().replace(
                    self.PathToBeReplace, self.RootStimuliPath)
            elif line[0:len(tag[5])] == tag[5]:
                self.ND = float(line[len(tag[5]):len(line)])

    def getNumberOfSequences(self):
        return len(self.Lsequence)

    # From the path in report file take the information about sequence image from info static file
    #
    def loadSequenceParameter(self):

        for sequence in self.Lsequence:
            # print('sequence:',sequence)
            path = sequence['Sequence Path']
            if not os.path.exists(path):
                print('Stimuli stored at different location!')
                ind = 0
                opath = path
                while not os.path.exists(path) and ind>-1:
                    #ind = opath.rfind('\\',0,ind-1)
                    ind = opath.find('\\',ind+1)
                    #print(ind, self.RootStimuliPath+os.sep+path[ind+1:])
                    path = self.RootStimuliPath+os.sep+opath[ind+1:]+os.sep
                    path = path.replace('\\',os.sep)

            # print("seq path:"+path)
            LsequenceImage = Ut.fileList(path, "png")
            # print(LsequenceImage)

            if (len(LsequenceImage) == 0):
                LsequenceImage = Ut.fileList(path, "bmp")
            # print('LsequenceImage:', LsequenceImage)
            sequence['Sequence Image List Path'] = LsequenceImage
            infoFileName = path + os.sep + "info.txt"
            if os.path.exists(infoFileName):
                infoFile = open(infoFileName, 'r')
            else:
                ind = 0
                while not os.path.exists(infoFileName) and ind>-1:
                    ind = path.rfind('\\',0,ind-1)
                    #print(ind, self.RootStimuliPath+os.sep+path[ind+1:])
                    infoFileName = self.RootStimuliPath+os.sep+path[ind+1:]+os.sep+'info.txt'
                    infoFileName = infoFileName.replace('\\',os.sep)
                try:
                    infoFile = open(infoFileName, 'r')
                except:
                    print('cannot find file: '+infoFileName)
                    raise SystemExit

            infoLines = infoFile.readlines()
            sequence['Modality'] = int(infoLines[1].strip())
            sequence['Number of image'] = int(infoLines[3].strip())
            sequence['Duration'] = int(infoLines[5].strip())

            # Creation of Index of Images List

            if sequence['Modality'] == 2:
                sequence['Sequence Lenght'] = int(infoLines[7].strip())
                LsequenceIndex = []
                for i in range(9, len(infoLines)):
                    LsequenceIndex.append(int(infoLines[i].strip()))
                sequence['Sequence List Frame Index'] = LsequenceIndex

                LindexOfImages = []
                tempCount = 0

                for i in range(0, sequence['Num of stimuli displayed']):
                    LindexOfImages.append(
                        sequence['Sequence List Frame Index'][tempCount])
                    tempCount += 1
                    if tempCount == sequence['Sequence Lenght']:
                        tempCount = 0
                sequence['Index of Images List'] = LindexOfImages
            else:
                LindexOfImages = []
                tempCount = 0

                for i in range(0, sequence['Num of stimuli displayed']):
                    LindexOfImages.append(tempCount)
                    tempCount += 1
                    if tempCount == len(sequence['Sequence Image List Path']):
                        tempCount = 0
                sequence['Index of Images List'] = LindexOfImages

    # Load time stamp event of stimuli from trigger file creating a dictionary
    # Dtrigger

    def loadTimeStampStimuli(self):

        try:
            tempStimuliFile = scipy.io.loadmat(self.TriggerFilePath)
        except:
            tempStimuliFile = h5py.File(self.TriggerFilePath,'r')

        print("Triggerfile:"+self.TriggerFilePath)
        for element in tempStimuliFile:
            # print("Element: "+element)
            if element == "StartFrame":
                self.Dtrigger[element] = tempStimuliFile[element]
            elif element == "SamplingFrequency":
                self.Dtrigger[element] = tempStimuliFile[element]
            elif element == "StopFrame":
                self.Dtrigger[element] = tempStimuliFile[element]
            elif element == "timeStampMatrix":
                temp = tempStimuliFile[element]
                # print(temp, temp.shape, temp.value)
                if len(temp.shape) == 2:
                    if temp.shape[0] > temp.shape[1]:
                        self.Dtrigger["TimeStampArray"] = temp.T[0]
                    else:
                        self.Dtrigger["TimeStampArray"] = temp[0]
                else:
                    self.Dtrigger["TimeStampArray"] = temp
                # print self.Dtrigger["TimeStampArray"].shape
                # print self.Dtrigger["TimeStampArray"]

    def rearrangeData(self):

        self.processStimuliMatrix()
        self.processGlobalParameter()

    def processGlobalParameter(self):

        tag = ["Num of refresh per stimulus", "Num of stimuli to be displayed",
               "Num of stimuli displayed", "Duration", "Modality", "Number of image"]

        StimuliParameterMatrix = np.array([])
        first = 0
        for sequence in self.Lsequence:

            tempList = []
            for key in tag:
                tempList.append(sequence[key])
            tempArray = np.asarray(tempList)
            if first == 0:
                StimuliParameterMatrix = tempArray
                first = 1
            else:
                StimuliParameterMatrix = np.column_stack(
                    (StimuliParameterMatrix, tempArray))

        stimuliSequenceMatrix = []
        stimuliSequenceMatrix.append(
            np.asarray(range(0, len(self.Lsequence))) + 1)
        stimuliSequenceMatrix.append(self.DHdf5Matrix['StimuliSequenceArray'])

        for key in tag:
            tempList = []

            for sequence in self.Lsequence:
                tempList.append(sequence[key])
            stimuliSequenceMatrix.append(np.asarray(tempList))

        self.DHdf5Matrix['StimuliSequenceMatrix'] = stimuliSequenceMatrix
        self.globalParameter['NeutralDensity'] = [self.ND]
        self.globalParameter['MaskingFilePath'] = self.MaskingFile

    def processStimuliMatrix(self):
        '''Create StimuliImageArray,StimuliSequenceArray and StimulusMatrix'''

        # Structure initialize for the first two matrix
        stimuliSequencePathList = []
        stimuliImageArray = []#np.array([[],[]])#asarray([[],[]]).T

        # Structure use for construct Stimulus Matrix
        tempCount = 0
        idSequenceArray = np.asarray([])
        idFrameArray = np.asarray([])
        idImageArray = np.asarray([])

        for sequence in self.Lsequence:
            # Take list of image of the sequence,transform in array and
            # concatenate with the general array
            pathArray = np.asarray(sequence['Sequence Image List Path'])
            # print("sa, pa",stimuliImageArray,pathArray)
            stimuliImageArray = np.hstack(
                (stimuliImageArray, pathArray))#), axis=1)

            # Append the path of sequence that will be transform in an array

            stimuliSequencePathList.append(sequence['Sequence Path'])

            # For every sequence create an array with index of stimuli display with in
            # the Array contain the index of sequence for every frame
            # and the array contain the index of stimuli frame pointing to the
            # stimuliImageArray

            idFrameArray = np.hstack((idFrameArray, np.arange(
                0, sequence['Num of stimuli displayed'])))#, axis=1)
            idSequenceArrayTemp = np.zeros(
                sequence['Num of stimuli displayed']) + sequence['Sequence Index']
            idSequenceArray = np.hstack(
                (idSequenceArray, idSequenceArrayTemp))#, axis=1)
            idImageArray = np.hstack((idImageArray, np.asarray(
                sequence['Index of Images List']) + tempCount))
            tempCount += len(sequence['Sequence Image List Path'])

        # Store the matrix in a dictionary

        self.DHdf5Matrix['StimuliImageArray'] = [np.asarray(
            range(0, len(stimuliImageArray))), stimuliImageArray]
        self.DHdf5Matrix['StimuliSequenceArray'] = np.asarray(
            stimuliSequencePathList)
        attributes = ["Id_sequence", "Id_frame", "Id_image"]
        self.DHdf5Matrix['StimuliMatrix'] = [np.asarray(range(0, len(
            idSequenceArray))) + 1, idSequenceArray + 1, idFrameArray + 1, idImageArray + 1]

        timeStampArray = self.Dtrigger['TimeStampArray']

        idStimArray = np.arange(1, len(timeStampArray) + 1) * -1

        idIsStimArray = np.zeros(len(timeStampArray)) - 1
        AmplitudeEvent = np.zeros(len(timeStampArray), dtype=int)

        # self.timeStampList=[timeStampArray,idStimArray,idIsStimArray,AmplitudeEvent]
        self.DHdf5Matrix['TimeStampMatrix'] = [
            timeStampArray, idStimArray, AmplitudeEvent]
'''
##############################################################################

Neuron Data CLASS

##############################################################################
'''
# Load neuron data from hdf5 produce after sorting processing chain


class NeuronDataSort():

    def __init__(self, *args):

        if len(args) == 1:
            self.SpikeTrainFile = args[0]
        else:
            return " You put too match input"

        # Initialize principal structure

        self.Lneuron = []

        self.timeStampList = []
        self.DHdf5Matrix = {}
        self.globalParameter = {}

        H5file = H5Tables.H5(self.SpikeTrainFile, 'r')

        self.DHdf5Matrix = H5file.getStructureDict()

        self.globalParameter = self.DHdf5Matrix.pop("GlobalParameter", None)

        H5file.close()
    # Back compatibility with NeuronData class

    def loadNeuronSpikeSortedData(self):

        pass

    def rearrangeData(self):
        pass

    def createTimeStampArray(self):
        pass


class NeuronData():

    def __init__(self, *args):

        if len(args) == 1:
            self.SpikeTrainFile = args[0]
        else:
            return " You put too match input"

        # Initialize principal structure

        self.Lneuron = []
        self.parameter = {'data': [], 'Attributes': []}
        self.timeStampList = []
        self.DHdf5Matrix = {}
        self.globalParameter = {}
        self.RecordingParameter = []

        self.loadNeuronSpikeData()

    # Load neuron information from mat file

    def loadNeuronSpikeData(self):

        mat = scipy.io.loadmat(self.SpikeTrainFile)
        SamplingFrequency = 0
        StopFrame = 0
        StartFrame = 0

        for channel in mat:
            Dchannel = {}
            if channel[0:2] == 'Ch':
                spikearray = mat[channel].nonzero()[0]
                Dchannel['spikeTrain'] = spikearray

                row = int(channel[2:4])
                col = int(channel[5:7])

                Dchannel['position'] = [row, col]
                Dchannel['name'] = channel
                self.Lneuron.append(Dchannel)

            else:
                if channel == 'SamplingFrequency':
                    SamplingFrequency = mat[channel]
                    self.parameter['data'].append(SamplingFrequency)
                    self.parameter['Attributes'].append('SamplingFrequency')

                elif channel == 'StartFrame':
                    StartFrame = mat[channel]
                    self.parameter['data'].append(StartFrame)
                    self.parameter['Attributes'].append('StartFrame')
                elif channel == 'StopFrame':
                    StopFrame = mat[channel]
                    self.parameter['data'].append(StopFrame)
                    self.parameter['Attributes'].append('StopFrame')

        self.Lneuron.sort(key=operator.itemgetter('name'))
        self.RecordingParameter.append(
            (self.SpikeTrainFile, StartFrame, StopFrame, SamplingFrequency))

    # Rearrange dictionary data producing matrix structure for hdf5 file

    def rearrangeData(self):

        timeStampArray = np.array([], dtype=np.int32)
        idNeuronArray = np.array([], dtype=np.int32)

        NeuronMatrix = [[], [], [], [], [], []]

        NeuronPosList = []
        NeuronNameList = []
        for index, neuron in enumerate(self.Lneuron):
            timeStampArray = np.concatenate(
                (timeStampArray, neuron['spikeTrain']), axis=1)
            idNeuron = np.zeros(len(neuron['spikeTrain'])) + index + 1
            idNeuronArray = np.concatenate((idNeuronArray, idNeuron), axis=1)

            NeuronMatrix[0].append(index + 1)
            NeuronMatrix[1].append(neuron['name'])
            NeuronMatrix[2].append(neuron['position'][0])
            NeuronMatrix[3].append(neuron['position'][1])
            NeuronMatrix[4].append(0)
            NeuronMatrix[5].append(0)

        # idEventArray=np.arange(0,len(idNeuronArray))
        amplitudeEvent = np.zeros(len(idNeuronArray))

        attributes = ["Time Stamp", "idNeuron_idStim", "AmplitudeEvent"]
        self.DHdf5Matrix['TimeStampMatrix'] = [
            timeStampArray, idNeuronArray, amplitudeEvent]
        self.DHdf5Matrix['NeuronMatrix'] = NeuronMatrix
        # self.DHdf5Matrix['NeuronNameArray']=np.asarray(NeuronNameList)

        self.globalParameter['RecordingParameter'] = self.RecordingParameter

    # Create time stamp matrix ( use only without stimuli data )
    def createTimeStampArray(self):
        pass
        #attributes=["Time Stamp","idNeuron_idStim","AmplitudeEvent"]
        #timeStampArray= [(self.timeStampList[0][i],self.timeStampList[1][i],self.timeStampList[2][i],self.timeStampList[3][i]) for i in np.argsort(self.timeStampList[0])]

        # self.DHdf5Matrix['TimeStampMatrix']={'data':np.asarray(timeStampArray).T,'Attributes':attributes}

    def getTempArraySpike():
        pass

    def getRecordingParameter():
        pass

    def getNumberOfNeuron():
        pass
