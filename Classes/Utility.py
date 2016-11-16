# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 13:47:30 2014

@author: szordan
"""
import os
from xml.dom import minidom
import numpy as np
'''
###############################################################################
CLASS
###############################################################################
'''

'''
Class for saving and load information in a xml file as username and past setted path

'''


class xmlSetting():

    def __init__(self, path, tag):
        self.pathFile = path
        self.xmldoc = minidom.parse(self.pathFile)
        self.itemlist = self.xmldoc.getElementsByTagName(tag)

    def getList(self, tag, attributeName):
        listOut = []
        listDom = self.xmldoc.getElementsByTagName(tag)
        for element in listDom:
            listOut.append[element.attributes[attributeName].value]

    def findKeyValue(self, key):
        for element in self.itemlist:
            if element.attributes['key'].value == key:
                return element.attributes['value'].value

    def setKeyValue(self, key, newValue):
        for element in self.itemlist:
            if element.attributes['key'].value == key:
                element.attributes['value'].value = newValue

    def saveChange(self):
        f = open(self.pathFile, 'w')
        self.xmldoc.writexml(f)
        f.close()

'''
###############################################################################
FUNCTION
###############################################################################
'''


def normMatrix(matrix):
    lenRow = matrix.shape[0]
    lenCol = matrix.shape[1]
    newMatrix = np.zeros(matrix.shape)
    maxValue = matrix.max()
    minValue = matrix.min()
    diff = maxValue - minValue
    for i in range(0, lenRow):
        for j in range(0, lenCol):
            newMatrix[i, j] = (matrix[i, j] - minValue) / diff
    return newMatrix


def ConvertString(s, mode):
    raise ValueError('This should not be called!')
    if mode == "WinToUnix":
        return s.replace(s.split('\\')[0], '/media/datacluster').replace('\\', '/')
    else:
        if os.path.exists('Y:\\'):
            return s.replace('/media/datacluster', 'Y:').replace('/', '\\')
        else:
            return s.replace('/media/datacluster', 'Z:').replace('/', '\\')


def ConvertStringList(L, mode):
    newL = []
    for element in L:
        newL.append(ConvertString(element, mode))
    return newL


def fileList(path, ext):

    separator = os.sep

    fileList = []
    for path, dirs, files in os.walk(path):
                # print files
                # if len(dirs)>0 : print dirs

        if len(files) > 0:
            files = sorted(files)
            for f in files:
                if f[len(f) - len(ext):len(f)] == ext:
                    stringa = path + separator + f
                    fileList.append(stringa)
    return fileList


def fileListFind(path, ext):
    '''

    Return a list of file path that contain "ext" character in the name

    '''
    print path
    separator = os.sep

    fileList = []
    for path, dirs, files in os.walk(path):
                # print files
                # if len(dirs)>0 : print dirs
        if len(files) > 0:
            for f in files:
                if f.find(ext) != -1:
                    stringa = path + separator + f
                    fileList.append(stringa)
    fileList.sort()
    return fileList


def buildPathHdf5(pathFileSource, pathout, tag):

    separator = os.sep

    fileName = pathFileSource.split(separator)[-1]
    fileName = fileName.split(".")[0]
    outFileName = fileName[0:len(fileName)] + tag + ".hdf5"
    finalPath = pathout + separator + outFileName
    return finalPath
