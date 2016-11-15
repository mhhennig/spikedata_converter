# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 15:56:16 2014

@author: szordan
"""
import sys
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import time

sys.path.append(os.sep)

import Classes.Utility as Ut
import Classes.SplitConvertProcedure as SCP


class WinConversion(QMainWindow):
    '''
    Windows for convert source file in a unique Hdf5 file
    '''

    def __init__(self):
        QMainWindow.__init__(self)

        # Session class

        self.layoutV = QVBoxLayout()
        self.layoutContent = QVBoxLayout()

        self.ButConvert = QPushButton('Convert')
        self.ButConvert.setStyleSheet("QPushButton{background: #C2E0FF;}")
        self.Labelprotocol = QLabel("Select protocol : ")
        self.selectProtocol = QComboBox()
        self.selectProtocol.addItem("Only convert visual stimuli")
        self.selectProtocol.addItem("Visual stimuli and split")
        self.selectProtocol.addItem("Only split")

        self.Labelinfile = QLabel("Input file path: ")
        self.inFilePath = QLineEdit()
        self.Labeloutfile = QLabel("Output file path: ")
        self.outFilePath = QLineEdit()

        self.LabelStatus = QLabel("Set Paths and push Convert")

        self.layoutV.addWidget(self.Labelprotocol)
        self.layoutV.addWidget(self.selectProtocol)

        self.layoutV.addWidget(self.Labelinfile)
        self.layoutV.addWidget(self.inFilePath)
        self.layoutV.addWidget(self.Labeloutfile)
        self.layoutV.addWidget(self.outFilePath)

        self.layoutContent.addLayout(self.layoutV)
        self.layoutContent.addWidget(self.ButConvert)
        self.layoutContent.addWidget(self.LabelStatus)

        self.connect(self.ButConvert, SIGNAL('clicked()'), self.send)
        self.connect(self.inFilePath, SIGNAL(
            'returnPressed()'), self.selectInPath)
        self.connect(self.outFilePath, SIGNAL(
            'returnPressed()'), self.selectOutPath)

        # self.connect(self.selectProtocol,SIGNAL('currentIndexChanged()'),self.organizeLayout)
        self.selectProtocol.currentIndexChanged[
            'QString'].connect(self.organizeLayout)

        self.widget = QWidget()

        self.organizeLayout()

        self.widget.setLayout(self.layoutContent)
        self.setCentralWidget(self.widget)
        self.setWindowTitle('Conversion')

    def organizeLayout(self):

        if self.selectProtocol.currentIndex() == 0:

            self.LabelPathStimuli = QLabel("Root path stimuli : ")
            self.LabelPathReplace = QLabel("Path to be replace in report : ")

            self.rootPathStimuli = QLineEdit()
            self.PathToBeReplace = QLineEdit()

            self.connect(self.rootPathStimuli, SIGNAL(
                'returnPressed()'), self.selectRootStimuliPath)
            self.connect(self.PathToBeReplace, SIGNAL(
                'returnPressed()'), self.selectPathToBeReplace)
            self.layoutV.addWidget(self.LabelPathStimuli)
            self.layoutV.addWidget(self.rootPathStimuli)
            self.layoutV.addWidget(self.LabelPathReplace)
            self.layoutV.addWidget(self.PathToBeReplace)

        else:

            self.layoutV.removeWidget(self.LabelPathStimuli)
            self.layoutV.removeWidget(self.LabelPathReplace)
            self.layoutV.removeWidget(self.rootPathStimuli)
            self.layoutV.removeWidget(self.PathToBeReplace)
            self.LabelPathReplace.deleteLater()
            self.LabelPathStimuli.deleteLater()
            self.rootPathStimuli.deleteLater()
            self.PathToBeReplace.deleteLater()

        self.loadSetting()

    def loadSetting(self):
        path = os.path.abspath("") + os.sep + "setting.xml"
        self.setting = Ut.xmlSetting(path, "path")

        self.inFilePath.setText(self.setting.findKeyValue('InputFile'))
        self.outFilePath.setText(self.setting.findKeyValue('OutputFile'))

        if self.selectProtocol.currentIndex() == 0:
            self.rootPathStimuli.setText(
                self.setting.findKeyValue('RootStimuliFile'))
            self.PathToBeReplace.setText(
                self.setting.findKeyValue('PathToBeReplace'))

    def closeEvent(self, event):
        self.setting.saveChange()

    def send2(self):
        self.LabelStatus.setText("Visual protocol Processing")
        self.LabelStatus.setStyleSheet("QLabel{background: #FF0000;}")

    def send(self):
        if self.selectProtocol.currentIndex() == 2:
            self.LabelStatus.setText("Split HDF5 Processing")
            self.LabelStatus.setStyleSheet("QLabel{background: #FF0000;}")
            QCoreApplication.processEvents()
            QCoreApplication.flush()
            SCP.split(str(self.inFilePath.text()),
                      str(self.outFilePath.text()))

            self.LabelStatus.setText("All files are converted!")
            self.LabelStatus.setStyleSheet("QLabel{background: #00FF00;}")
            QCoreApplication.processEvents()
            QCoreApplication.flush()
            QMessageBox.about(self, 'Finish Processing',
                              'Conversion completed!')
        elif self.selectProtocol.currentIndex() == 1:

            self.LabelStatus.setText("Visual protocol Processing")
            self.LabelStatus.setStyleSheet("QLabel{background: #FF0000;}")
            QCoreApplication.processEvents()
            QCoreApplication.flush()

            SCP.split(str(self.inFilePath.text()), str(self.inFilePath.text()))
            SCP.convert(str(self.inFilePath.text()), str(self.outFilePath.text()), str(
                self.rootPathStimuli.text()), str(self.PathToBeReplace.text()))

            self.LabelStatus.setText("All files are converted!")
            self.LabelStatus.setStyleSheet("QLabel{background: #00FF00;}")
            QCoreApplication.processEvents()
            QCoreApplication.flush()
            QMessageBox.about(self, 'Finish Processing',
                              'Conversion completed!')
        else:
            print("Converting only the stimuli")
            self.LabelStatus.setText("Visual protocol Processing")
            self.LabelStatus.setStyleSheet("QLabel{background: #FF0000;}")
            QCoreApplication.processEvents()
            QCoreApplication.flush()

            SCP.convertStimuliOnly(str(self.inFilePath.text()), str(self.outFilePath.text()), str(
                self.rootPathStimuli.text()), str(self.PathToBeReplace.text()))

            self.LabelStatus.setText("All files are converted!")
            self.LabelStatus.setStyleSheet("QLabel{background: #00FF00;}")
            QCoreApplication.processEvents()
            QCoreApplication.flush()
            QMessageBox.about(self, 'Finish Processing',
                              'Conversion completed!')

    def selectInPath(self):
        path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.inFilePath.setText(path)
        self.setting.setKeyValue('InputFile', path)
        self.setting.saveChange()

    def selectOutPath(self):
        path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.outFilePath.setText(path)
        self.setting.setKeyValue('OutputFile', path)
        self.setting.saveChange()

    def selectRootStimuliPath(self):
        path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.rootPathStimuli.setText(path)
        self.setting.setKeyValue('RootStimuliFile', path)
        self.setting.saveChange()

    def selectPathToBeReplace(self):
        path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.PathToBeReplace.setText(path)
        self.setting.setKeyValue('PathToBeReplace', path)
        self.setting.saveChange()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = WinConversion()
    main.show()
    sys.exit(app.exec_())
