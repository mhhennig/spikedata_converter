|Created by Szordan
-------------------

|Introduction:
--------------

The provided GUI will allow you to split the h5 file after spike coincidence removal
in the original files and to convert the stimuli information

The Hdf5 file that has to be splitted must be renamed as " *whatyouwant*Sorted.hdf5 " and the 
file txt containing the list of the concatenated files (i.e. OriginalSplit.txt) as " *whatyouwant*Sorted.txt "

---> Example:
| - myResultFileName.hdf5 -> myResultFileNameSorted.hdf5
| - OriginalSplie.txt	  -> myResultFileNameSorted.txt

|How to use it:
---------------	

Launch gui.py, a GUI will appear.
You can choose between two protocols:
"VisualStimuli" - it will split the file and integrate the stimuli info
"OnlySplit"     - it will split the file only

parameter to be set:
"input selected path" 		- the path containing the myResultFileNameSorted.hdf5 and myResultFileNameSorted.txt
				  in the case of option "VisualStimuli" this folder has to contain also all the report files 
				  and the stim trigger .mat files for all the phases merged 
"output file path"    		- the path where to store the results
"root path stimuli"   		- the root folder containing all the sequence of stimulation on the local machine where the script is executed
"path to be replaced in report" - the path to the be changed in the report file to match the correct path on the machine you are 
				  launching the script (see example below)
---> Example:
|Considering that the folder containing the stimuli folders used for the exp in the pc of stimulation is C:\myFolderNameOnTheStimPc\
|In this case the report file for each sequence will store the path as for example:
|	"Executing sequence: C:\myFolderNameOnTheStimPc\Bars_Mooving\Gray" 
|	(where Bars_Mooving is a folder containing a Folder Gray that is the folder with the stimuli images)
|Then you have to copy all the stimuli folders in the local machine where you want to execute the script.
|Considering that the folder containing the stimuli folders in this case is named C:\myFolderNameOnTheAnalysisPc\
|(so it will also contain C:\myFolderNameOnTheAnalysisPc\Bars_Mooving\Gray)
|Then the script will have to substitute the path C:\myFolderNameOnTheStimPc\Bars_Mooving\Gray
|with the path C:\myFolderNameOnTheAnalysisPc\Bars_Mooving\Gray 
|
|So you have to set the following parameters:
|"root path stimuli" 		= C:\myFolderNameOnTheAnalysisPc
|"path to be replaced in report" = C:\myFolderNameOnTheStimPc

input parameter can be set directly typing in the fields or by selecting the field and clicking return button

|List of python files:
----------------------

Classes\SplitConvertProcedure.py
Classes\Conversion.py
Classes\H5Tables.py
Classes\Utility.py
gui.py





	

	


	

