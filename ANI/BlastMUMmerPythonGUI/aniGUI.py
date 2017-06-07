from __future__ import division
import Tkinter, Tkconstants, tkFileDialog
from Tkinter import *
import os
from collections import OrderedDict
import sys
import subprocess
import fileinput
import shutil
import datetime
from datetime import datetime
import csv
'''
This class will concatinate many .fna files into one .fna file
Example:
user/directory1/file1.ext
user/directory1/file2.ext
Will turn directory1 into
user/directory1/directory1.ext
Where directory1 is the result of cat file1.ext file2.ext > dictionary1.ext
If there is more than one .fna file, it will leave it untouched
Example:
user/directory2/file.ext will remain unchanged, where dictionary2 has only file.ext
'''
class ConcatinatingFNAFiles():
  def run(self):
    directory = TkANI.allFnaDirectoryLocationString
    ext = ".fna"
    currentDirectory = None
    previousDirectory = None
    listOfDirectoriesWithMultipleFNAfiles = []
    listOfDirectoriesWithMultipleFNAfilesNoRepeats = []
    x = 0
    for subdir, dirs, files in os.walk(directory, topdown=True):
      for file in files:
        currentDirectory = subdir
        if currentDirectory == previousDirectory:
          listOfDirectoriesWithMultipleFNAfiles.insert(x, currentDirectory)
          x +=1
        previousDirectory = subdir
    listOfDirectoriesWithMultipleFNAfilesNoRepeats = list(OrderedDict.fromkeys(listOfDirectoriesWithMultipleFNAfiles))
    print (listOfDirectoriesWithMultipleFNAfilesNoRepeats)
    listOfFNAfiles = []
    catCommand = None
    newFNAfileName = None
    for dirFile in listOfDirectoriesWithMultipleFNAfilesNoRepeats:
      uniqueFiles = True
      print(dirFile)
      base = os.path.basename(dirFile)
      fnaConstruct = (os.path.join(dirFile, base))
      newFNAfileName = fnaConstruct+ext
      for root, dir, files in os.walk(dirFile, topdown=True):
        for f in files:
          currentFile = ((os.path.join(root, f)))
          print (currentFile)
          if newFNAfileName == currentFile:
            print ("run already done before, error")
            uniqueFiles = False
            break
          if f.endswith(ext):
            catFileTo = ((os.path.join(root, f)))
            listOfFNAfiles.append(catFileTo+' ')
            print(catFileTo)
      if uniqueFiles:
        string1 = ''.join(listOfFNAfiles)
        catCommand = "cat " + string1 + " > " + newFNAfileName
        print (catCommand)
        os.system(catCommand)
        rmCommand = "rm " + string1
        print(rmCommand)
        os.system(rmCommand)
        listOfFNAfiles = []
        catCommand = None
        newFNAfileName = None
        string1 = None
    listOfDirectoriesWithMultipleFNAfilesNoRepeats = None
    print("ConcatinatingFNAFiles(): complete!")
'''
This class changes child files to their parent file names without the child file being called file junior
Example Input:
user/directory_1/nameToBeChanged.ext
Example Output:
user/directory_1/directory_1.ext
This class will run through all respective directories(directory_1, directory_2 ... directory_N)
This class is intended and has only been tested after running classes:
concatinatingFNAFiles():
'''
class ChangeFilesToDirectoryName():
  def run(self):
    directory = TkANI.allFnaDirectoryLocationString
    ext = ".fna"
    listOfDirectories = []
    x = 0
    first = True
    firstDirectoryString = ""
    for subdir, dirs, files in os.walk(directory, topdown=True):
      for file in files:
        currentDirectory = subdir
        if first:
          firstDirectoryString = str(currentDirectory) 
        if not first:
          listOfDirectories.insert(x, currentDirectory)
        x += 1
        first = False
    listOfDirectories = list(OrderedDict.fromkeys(listOfDirectories))
    currentFile = None
    catFile = None
    newFNAfileName = None
    uniqueFile = True
    for dirFile in listOfDirectories:
      uniqueFiles = True
      base = os.path.basename(dirFile)
      fnaConstruct = (os.path.join(dirFile, base))
      newFNAfileName = fnaConstruct+ext
      print (newFNAfileName)
      for root, dir, files in os.walk(dirFile, topdown=True):
        for f in files:
          if str(dir) == str(firstDirectoryString):
            break
          currentFile = ((os.path.join(root, f)))
          if newFNAfileName == currentFile:
            print ("already taken care of")
            uniqueFile = False
            break
          elif f.endswith(ext):
            catFile = currentFile
            print (catFile)
            catCommand = "cat " + catFile + " > " + newFNAfileName
            print (catCommand)
            os.system(catCommand)
            rmCommand = "rm " + catFile
            print(rmCommand)
            os.system(rmCommand)
            catCommand = None
            newFNAfileName = None
            currentFile = None
            catFile = None
          else:
          	print("error")
    listOfDirectories = None
    print("ChangeFilesToDirectoryName(): complete!")
'''
This class programmically moves(cuts[os.rename] or copies[catCommand]) files to another location
Example Input:
Users/all.fna/directory1/directory1.fna
Users/all.fna/directory2/directory2.fna
Example Output:
Users/newDirectory/directory1.fna
Users/newDirectory/directory2.fna
This class is intended and has only been tested after running classes:
concatinatingFNAFiles(): and class ChangeFilesToDirectoryName(): 

Current bug: Will throw an os error if ran twice(Errno 21: Is a directory:)
A string output is thrown if this bug is reproduced so the program doesn't crash
'''
class CopyFiles():
  def run(self, moveFilesToDirectoryParam):
    if(moveFilesToDirectoryParam == ""):
      print(str("Error: Cannot copy, one of the directory entry is nil"))
      return
    directory = TkANI.allFnaDirectoryLocationString
    # moveToDirectory = TkANI.blastDirectoryLocationString
    moveToDirectory = moveFilesToDirectoryParam
    if os.path.exists(moveToDirectory):
      rmCommand = "rm " + moveToDirectory
      print(rmCommand)
      # print(str("Error: " + moveFilesToDirectoryParam + " already exists, please delete and try again"))
      return
    else:
      print("Copying files to " + str(moveFilesToDirectoryParam))
      shutil.copytree(directory, moveToDirectory)
    print("Copying Blast files completed!")
'''
This class creates a delta file from all fna files
A delta file is required to calculate the Average Nucleotide Index for nucmer(MUMmer, but use for nucleotides)
This class is intended and has only been tested after running classes:
class ConcatinatingFNAFiles(): , class ChangeFilesToDirectoryName(): , class CdFiles():
'''
class AllDeltaFilesMummer():
  def run(self):
    directoryFNAfiles = TkANI.allFnaDirectoryLocationString
    mummerBaseDirectory = TkANI.mummerDirectoryLocationString
    saveDeltaFileNamesDirectory = TkANI.deltaFilesDirectoryLocationString
    vs = '__vs__'
    ext = ".fna"
    listOfFiles = []
    for root, subdir, files in os.walk(directoryFNAfiles, topdown=True):
      for file in files:
        currentDirectory = subdir
        if file.endswith(ext):
          currentFile = ((os.path.join(root, file)))
          listOfFiles.append(currentFile)
    listOfFilesNoRepeats = list(OrderedDict.fromkeys(listOfFiles))
    print (listOfFilesNoRepeats)
    n = len(listOfFilesNoRepeats)
    print (n)
    os.chdir(mummerBaseDirectory)
    if not os.path.exists(saveDeltaFileNamesDirectory):
      os.mkdir(saveDeltaFileNamesDirectory, 0777)
    for i in range(0,n):
      for j in range(i,n):
        if(listOfFilesNoRepeats[i] is not listOfFilesNoRepeats[j]):
          deltaFileName = (os.path.basename(os.path.splitext(listOfFilesNoRepeats[i])[0])+vs+os.path.basename(os.path.splitext(listOfFilesNoRepeats[j])[0]))
          saveDeltaFileNamesPath = saveDeltaFileNamesDirectory + deltaFileName
          mummers_command = './nucmer -p '+ saveDeltaFileNamesPath + ' ' + listOfFilesNoRepeats[i] + ' ' + listOfFilesNoRepeats[j]
          print (mummers_command)
          #os.system(mummers_command)
          process = subprocess.Popen(mummers_command, stdout=subprocess.PIPE, stderr=None, shell=True)
          output = process.communicate()
          p_status = process.wait()
          print ("return code:", p_status)
    print("AllDeltaFilesMummer(): completed!")
'''
This class calculates the average nucleotide index from delta files
This class is intended and has only been tested after running classes:
class ConcatinatingFNAFiles(): , class ChangeFilesToDirectoryName(): , class CdFiles(): and class AllDeltaFilesMummer():
'''
class MUMmerANI():
  def process_data(self, path):
    with open(path) as input_file:
      for line in input_file:
        chunk = line.strip().split(' ')
        if len(chunk) == 7:
          yield chunk
  def calculateANI(self, path):
    x=0
    list1 = []
    oneHundredANI = False
    for piece in self.process_data(path):
      if piece[0] == piece[2]:
        if piece[4] and piece[5] == '0':
          oneHundredANI = True
          length = abs(((int(piece[1])-int(piece[0]))) + abs((int(piece[3])-int(piece[2])))) / 2
          ANI = ((length - int(piece[5])) / length) * 100
          list1.insert(x, ANI)
          x += 1
          break
    if not oneHundredANI:
      for piece in self.process_data(path):
        length = abs(((int(piece[1])-int(piece[0]))) + abs((int(piece[3])-int(piece[2])))) / 2
        ANI = ((length - int(piece[5])) / length) * 100
        list1.insert(x, ANI)
        x += 1
    if not list1:
      return
    else:
      return ((sum(list1)/ x))
  def run(self):
    mummerBaseDirectory = TkANI.mummerDirectoryLocationString
    deltaFilesNamesDirectory = TkANI.deltaFilesDirectoryLocationString
    catANIDirectory = TkANI.percentDirectoryLocationString
    outputANIFileDirectoryFile = str(catANIDirectory) + "/outputANIFile.txt"
    os.system("rm " + outputANIFileDirectoryFile)
    outputANIFile = open(outputANIFileDirectoryFile, "w")
    vs = '__vs__'
    if(catANIDirectory == ""):
      print(str("Error: Percent ouput directory entry is nil"))
      return
    if not os.path.exists(catANIDirectory):
      os.mkdir(catANIDirectory, 0777)
    ext = ".delta"
    listOfFiles = []
    for root, subdir, files in os.walk(deltaFilesNamesDirectory, topdown=True):
      for file in files:
        currentDirectory = subdir
        if file.endswith(ext):
          currentFile = ((os.path.join(root, file)))
          listOfFiles.append(currentFile)
    listOfFilesNoRepeats = list(OrderedDict.fromkeys(listOfFiles))
    n =  len(listOfFilesNoRepeats)
    print listOfFilesNoRepeats
    compareString = "%ANI is low"
    for i in range(0,n):
      s = self.calculateANI(listOfFilesNoRepeats[i])
      if s:
        if s != compareString:
          print (str(os.path.basename(listOfFilesNoRepeats[i])) + ": " + str(s) + "%")
          catString = str(os.path.basename(listOfFilesNoRepeats[i])) + "_" + str(s) + "" + "\n"
          outputANIFile.write(catString)
    print("MummerAni(): completed!")
    outputANIFile.close()
'''
This class creates out files from all fna files
An out file is required to calculate the Average Nucleotide Index for blastn(blast, but use for nucleotides)
This class is intended and has only been tested after running classes:
class ConcatinatingFNAFiles(): , class ChangeFilesToDirectoryName(): and class CdFiles():
'''
class AllOutFilesBlast():        
  'This program takes a list of .faa files from all.fna and calculates every one vs every other one on a unix machine'
  'This script was intended to be used on all.fna downloaded from NCBI database on a unix machine'
  def run(self):
    '''
    blastBaseDirectory = '/Users/jon6/ncbi-blast-2.2.29+/'
    directoryFNAfiles = '/Users/jon6/ncbi-blast-2.2.29+/'
    saveOutFileNamesDirectory = '/Users/jon6/ncbi-blast-2.2.29+/outFiles/'
    '''
    directoryFNAfiles = TkANI.allFnaDirectoryLocationString
    blastBaseDirectory = TkANI.blastDirectoryLocationString
    saveOutFileNamesDirectory = TkANI.saveOutFileNamesDirectory
    
    ext = ".fna"
    vs = '__vs__'
    listOfFiles = []
    for root, subdir, files in os.walk(directoryFNAfiles, topdown=True):
      for file in files:
        currentDirectory = subdir
        if file.endswith(ext):
          currentFile = ((os.path.join(root, file)))
          listOfFiles.append(currentFile)
    listOfFilesNoRepeats = list(OrderedDict.fromkeys(listOfFiles))
    print (listOfFilesNoRepeats)
    n = len(listOfFilesNoRepeats)
    print (n)
    os.chdir(blastBaseDirectory)
    print (os.getcwd())
    for i in range(0,n):
      for j in range(i,n):
        if(listOfFilesNoRepeats[i] is not listOfFilesNoRepeats[j]):
          outFileName = (os.path.basename(os.path.splitext(listOfFilesNoRepeats[i])[0]) + vs + os.path.basename(os.path.splitext(listOfFilesNoRepeats[j])[0]))
          saveOutFileNamesPath = saveOutFileNamesDirectory + outFileName + ".out"
          blast_command = 'blastn -query ' + os.path.basename(listOfFilesNoRepeats[i]) + ' -subject ' + os.path.basename(listOfFilesNoRepeats[j]) + ' -out ' + saveOutFileNamesPath + ' ' + '-task blastn -dust no -outfmt "6 std qlen slen"'
          print (blast_command)
          process = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=None, shell=True)
          output = process.communicate()
          p_status = process.wait()
          print ("return code:", p_status)
    print("AllOutFilesBlast(): completed")
'''
This class creates the average nucleotide index from out files
This class is intended and has only been tested after running classes:
class ConcatinatingFNAFiles(): , class ChangeFilesToDirectoryName(): , class CdFiles(): and AllOutFilesBlast():
The Blast ANI algorithm is not currently correct and still in progress
'''
class BlastANI():
  '''
  path = '/Users/jon6/Desktop/blast_stuff/ecoli_salmonella.out'
  '''
  def process_data(self, path):
    with open(path) as input_file:
      for line in input_file:
        chunk = line.strip().split()
        yield chunk 
  def calculateANI(self, path):
    x = 0
    total = 0
    oneHundredPercentANI = False
    for piece in self.process_data(path):
      if piece[2] == '100.00' and piece[6] == '1':
        oneHundredPercentANI = True
        percent_identity = float(piece[2])
        decimal_percent_identity = (percent_identity/100)
        align = decimal_percent_identity*int(piece[3])
        gap_mismatch = int(piece[4]) + int(piece[5])
        total += align
        x += 1
        break
    if not oneHundredPercentANI:
      for piece in self.process_data(path):
        percent_identity = float(piece[2])
        decimal_percent_identity = (percent_identity/100)
        align = decimal_percent_identity*int(piece[3])
        gap_mismatch = int(piece[4]) + int(piece[5])
        total += align
        x += 1
    length = (int(piece[13]) + int(piece[12]))/2
    ANI = total /length
    return round(ANI*100, 2)
    print (calculateANI(path))
    print("BlastANI(): completed!")
    '''
    blastBaseDirectory = '/Users/jon6/ncbi-blast-2.2.29+/'
    outFilesNamesDirectory = '/Users/jon6/ncbi-blast-2.2.29+/'
    ext = ".out"
    listOfFiles = []
    for root, subdir, files in os.walk(deltaFilesNamesDirectory, topdown=True):
      for file in files:
        currentDirectory = subdir
        if file.endswith(ext):
          currentFile = ((os.path.join(root, file)))
          listOfFiles.append(currentFile)
    listOfFilesNoRepeats = list(OrderedDict.fromkeys(listOfFiles))
    n =  len(listOfFilesNoRepeats)
    print listOfFilesNoRepeats
    for i in range(0,n):
      s = calculateANI(listOfFilesNoRepeats[i])
      if s:
        print (str(os.path.basename(listOfFilesNoRepeats[i])) + ": " +str(s) + "%")
    '''
class ReformatTextFileForConversion():
  roundToPowerOfN = 0
  genome1 = ""
  species1 = ""
  genome2 = ""
  species2 = ""
  percent2 = .42
  percent3 = 42
  
  def run(self):
    with open(str(TkANI.percentDirectoryLocationString) + "/outputANIFile.txt") as infile:
      reader = csv.reader(infile)
      reformattedANIFile = open(str(TkANI.percentDirectoryLocationString) + "/reformattedOutputANIFile.txt", "w")
      for line in enumerate(reader):
        x = line[1][0].split("_")
        for i, j in enumerate(x):
          if i == 1:
            self.genome1 = x[i]
          if i == 3:
            self.species1 = x[i]
          if i == 10:
            self.genome2 = x[i]
          if i == 11:
            self.species2 = x[i]
          if i == 14:
            try: 
              float(x[i])
              self.percent2 = round(float(x[i]), self.roundToPowerOfN)
              self.percent3 = int(self.percent2)
            except:
              pass
          if i == 15:
            try:
              float(x[i])
              self.percent2 = round(float(x[i]), self.roundToPowerOfN)
              self.percent3 = int(self.percent2)
            except:
              pass
          if i == 16:
            try:
              float(x[i])
              self.percent2 = round(float(x[i]), self.roundToPowerOfN)
              self.percent3 = int(self.percent2)
            except:
              pass
          if i == 17:
            try:
              float(x[i])
              self.percent2 = round(float(x[i]), self.roundToPowerOfN)
              self.percent3 = int(self.percent2)
            except:
              pass
          if i == 18:
            try:
              float(x[i])
              self.percent2 = round(float(x[i]), self.roundToPowerOfN)
              self.percent3 = int(self.percent2)
            except:
              pass
        print(str(self.genome1 + " " + self.species1 + " " + self.genome2 + " " + self.species2 + " " + str(self.percent3)))
        writeString = str(self.genome1 + " " + self.species1 + " " + self.genome2 + " " + self.species2 + " " + str(self.percent3) + "\n")
        reformattedANIFile.write(writeString)
      reformattedANIFile.close()
      print("ReformatTextFileForConversion(): completed!")  

class RemoveDuplicates():
  def run(self):
    listOfLines = []
    y = ""
    list2 = ""
    with open(str(TkANI.percentDirectoryLocationString) + "/reformattedOutputANIFile.txt") as infile:
      reader = csv.reader(infile)
      reformattedANIFile = open(str(TkANI.percentDirectoryLocationString) + "/removedDuplicatesOutputANIFile2.txt", "w")
      for line in enumerate(reader):
        x = line[1][0].split(" ")
        for i, j in enumerate(x):
          if(i == 0):
            y = ""
          y = y + " " + str(j)
          if(i == 4):
            listOfLines.append(y)
      list2 = list(OrderedDict.fromkeys(listOfLines))
      for i, j in enumerate(list2):
        writeString = str(j) + "\n"
        reformattedANIFile.write(writeString)
      # print(list2)
      reformattedANIFile.close()
    print("RemoveDuplicates(): completed!")

'''
This class formats the text file outputted from BlastANI(): and MummersANI(): to JSON seed. This is done to be able to input into an iOS project for Core Data or another project utilizing JSON to read the data
This class is intended and will only be tested after running classes:
class ConcatinatingFNAFiles(): , class ChangeFilesToDirectoryName():, class CdFiles():
and
Class AllDeltaFilesMummer(): and class MummersANI(): 
and / or 
class AllOutFilesBlast(); and class BlastANI():
'''
class CreateJSONseed():
  def run(self):
    print("not complete")
'''
This class formats the text file outputted from BlastANI(): and MummersANI(): to SQL lite statement. This is done to enable an Android application or another project utilizing SQL lite to read the data
This class is intended and will only be tested after running classes:
class ConcatinatingFNAFiles(): , class ChangeFilesToDirectoryName():, class CdFiles():
'''
class CreateSQLliteCreationStatement():
  def run(self):
    print("not complete")
'''
This class is the GUI class for the application
When the user presses run, it will gather the inputs and execute def run(self):, which calls instances of other classes
To start this GUI, the user should have Python 2.7 installed and run python thisFileName.py and a Tkinter generated GUI should appear
'''
class TkANI(Tkinter.Frame):
  # Directory Location Vars
  allFnaDirectoryLocationString = ""
  blastDirectoryLocationString = ""
  outFilesDirectoryLocationString = ""
  mummerDirectoryLocationString = ""
  deltaFilesDirectoryLocationString = ""
  percentDirectoryLocationString = ""

  def __init__(self, root):
    # Initialize classes GUI will use
    self.concatinatingFNAFiles = ConcatinatingFNAFiles()
    self.changeFilesToDirectoryName = ChangeFilesToDirectoryName()
    self.copyFilesInstance = CopyFiles()
    self.allOutFilesBlast = AllOutFilesBlast()
    self.blastANI = BlastANI()
    self.allDeltaFilesMummer = AllDeltaFilesMummer()
    self.mummerANI = MUMmerANI()
    self.reformatTextFileForConversion = ReformatTextFileForConversion()
    self.removeDuplicates = RemoveDuplicates()

  	# GUI inits
    Tkinter.Frame.__init__(self, root)
    self.configure(background='black')
    yPaddingAmount = 2
    xPaddingAmount = 2
    entryHexColor = "#B20000"
    entryWidth = 45
    buttonWidth = 45
    readMeButtonWidth = 15
    runButtonWidth = 91
    checkBoxWidth = 37
    # Checkbox Value Variables
    self.completeCheckboxValue = IntVar()
    self.cdAllCheckboxValue = IntVar()
    self.mummerCheckboxValue = IntVar()
    self.mummerANICheckboxValue = IntVar()
    self.blastCheckboxValue = IntVar()
    self.blastANICheckboxValue = IntVar()
    self.sqlLiteCheckboxValue = IntVar()
    self.jsonCheckboxValue = IntVar()
     # Defining query directory buttons
    Tkinter.Button(self, text='Location Of all.fna Directory', command=self.askAllFnaDirectory, width=buttonWidth).grid(row=1, column=1, pady=(yPaddingAmount, yPaddingAmount))
    Tkinter.Button(self, text='BLAST Directory Containing Directories Of .fna Files', command=self.askBlastDirectory, width=buttonWidth).grid(row=2, column=1, pady=(yPaddingAmount, yPaddingAmount))
    Tkinter.Button(self, text='Out Files (from BLAST Output) Output Directory', command=self.askOutFilesOutputDirectory, width=buttonWidth).grid(row=4, column=1, pady=(yPaddingAmount, yPaddingAmount))
    Tkinter.Button(self, text='MUMmer Directory Location Containing Directories Of .fna Files', command=self.askMummerDirectory, width=buttonWidth).grid(row=3, column=1, pady=(yPaddingAmount, yPaddingAmount))
    Tkinter.Button(self, text='Delta Files (from MUMmer Output) Output Directory', command=self.askDeltaFilesOutputDirectory, width=buttonWidth).grid(row=5, column=1, pady=(yPaddingAmount, yPaddingAmount))
    Tkinter.Button(self, text='Percent ANI Output Results Directory', command=self.askPercentANIOutputDirectory, width=buttonWidth).grid(row=6, column=1, pady=(yPaddingAmount, yPaddingAmount))

     # Define run button
    Tkinter.Button(self, text='Run', command=self.run, width=runButtonWidth).grid(row=7, column=1, columnspan=2, pady=(yPaddingAmount, yPaddingAmount))
     # Define directory text entries
    self.allFnaEntry = Tkinter.Entry(self, bg=entryHexColor, fg="white", width=entryWidth)
    self.allFnaEntry.grid(row=1, column=2, padx=(xPaddingAmount, xPaddingAmount))
    # Blast
    self.blastEntry = Tkinter.Entry(self, bg=entryHexColor, fg="white", width=entryWidth)
    self.blastEntry.grid(row=2, column=2, padx=(xPaddingAmount, xPaddingAmount))
     # Code golf requires two lines here because widget will return none if completed in one line
    self.outEntry = Tkinter.Entry(self, bg=entryHexColor, fg="white",  width=entryWidth)
    self.outEntry.grid(row=4, column=2, padx=(xPaddingAmount, xPaddingAmount))
    # Under Par
    self.mummerEntry = Tkinter.Entry(self, bg=entryHexColor, fg="white", width=entryWidth)
    self.mummerEntry.grid(row=3, column=2, padx=(xPaddingAmount, xPaddingAmount))
     # More boiler plate under par entries
    self.deltaEntry = Tkinter.Entry(self, bg=entryHexColor, fg="white", width=entryWidth)
    self.deltaEntry.grid(row=5, column=2, padx=(xPaddingAmount, xPaddingAmount))
     # Masters Tournament
    self.percentEntry = Tkinter.Entry(self, bg=entryHexColor, fg="white", width=entryWidth)
    self.percentEntry.grid(row=6, column=2, padx=(xPaddingAmount, xPaddingAmount))
    '''
    Define checkboxes
    1 is an enabled checkbox, 0 is an enabled checkbox
    '''
    self.completeCheckbox = Tkinter.Checkbutton(self, variable=self.completeCheckboxValue, onvalue=1, offvalue=0, text ='Reformat, rename and concatenate', width=checkBoxWidth)
    self.completeCheckbox.grid(row=1, column=3)
    self.cdAllCheckbox = Tkinter.Checkbutton(self, variable=self.cdAllCheckboxValue, onvalue=1, offvalue=0, text ='Move modified files to Blast/MUMmer Directories', width=checkBoxWidth)
    self.cdAllCheckbox.grid(row=2, column=3)
    self.blastCheckox = Tkinter.Checkbutton(self, variable=self.blastCheckboxValue, onvalue=1, offvalue=0, text ='Execute Blast', width=checkBoxWidth)
    self.blastCheckox.grid(row=3, column=3)
    self.blastANICheckbox = Tkinter.Checkbutton(self, variable=self.blastANICheckboxValue, onvalue=1, offvalue=0, text ='Execute Blast ANI', width=checkBoxWidth)
    self.blastANICheckbox.grid(row=4,column=3)
    self.mummerCheckbox = Tkinter.Checkbutton(self, variable=self.mummerCheckboxValue, onvalue=1, offvalue=0, text ='Execute MUMmer', width=checkBoxWidth)
    self.mummerCheckbox.grid(row=5, column=3)
    self.mummerANICheckbox = Tkinter.Checkbutton(self, variable=self.mummerANICheckboxValue, onvalue=1, offvalue=0, text ='Execute MUMmer ANI', width=checkBoxWidth)
    self.mummerANICheckbox.grid(row=6, column=3)  
    self.sqlLiteCheckbox = Tkinter.Checkbutton(self, variable=self.sqlLiteCheckboxValue, onvalue=1, offvalue=0, text ='Create output for SQLite', width=checkBoxWidth)
    self.sqlLiteCheckbox.grid(row=7, column=3)
    self.jsonCheckbox = Tkinter.Checkbutton(self, variable=self.jsonCheckboxValue, onvalue=1, offvalue=0, text ='Create output for JSON', width=checkBoxWidth)
    self.jsonCheckbox.grid(row=8, column=3)

    # Set checkbox default values
    self.completeCheckboxValue.set(0)
    self.blastCheckboxValue.set(0)
    self.blastANICheckboxValue.set(0)
    self.mummerCheckboxValue.set(0)
    self.mummerANICheckboxValue.set(0)
    self.sqlLiteCheckboxValue.set(0)
    self.jsonCheckboxValue.set(0)
    # Remove for generic program
    self.allFnaEntry.insert(0, str("/Users/jon/Desktop/allfnabackup3/all.fna"))
    self.blastEntry.insert(0, str("/Users/jon/Desktop/testb"))
    self.outEntry.insert(0, str("/Users/jon/Desktop/blastOut"))
    self.mummerEntry.insert(0, str("/Users/jon/MUMmer3.23"))
    self.deltaEntry.insert(0, str("/Users/jon/MUMmer3.23/allDeltaFiles"))
    self.percentEntry.insert(0, str("/Users/jon/Desktop/MUMmer3.23/percentANIFileOutput"))
    
   # Queries directory containing files respective to BLAST
  def askAllFnaDirectory(self):
     askAllFnaDirectoryLocation = tkFileDialog.askdirectory()
     if askAllFnaDirectoryLocation:
        self.allFnaEntry.delete(0, 'end')
        self.allFnaEntry.insert(0, str(askAllFnaDirectoryLocation))
     return
   # Queries directory containing files respective to BLAST
  def askBlastDirectory(self):
     blastDirectoryLocation = tkFileDialog.askdirectory()
     if blastDirectoryLocation:
        self.blastEntry.delete(0, 'end')
        self.blastEntry.insert(0, str(blastDirectoryLocation))
     return
   # Queries directory containing files respective to MUMmer
       # Returns directory containing files to dump Delta files
  def askOutFilesOutputDirectory(self):
     outFilesDirectoryLocation = tkFileDialog.askdirectory()
     if outFilesDirectoryLocation:
        self.outEntry.delete(0, 'end')
        self.outEntry.insert(0, str(outFilesDirectoryLocation))
     return
  def askMummerDirectory(self):
     mummerDirectoryLocation = tkFileDialog.askdirectory()
     if mummerDirectoryLocation:
        self.mummerEntry.delete(0, 'end')
        self.mummerEntry.insert(0, str(mummerDirectoryLocation))
     return
 	 # Queries directory containing files to dump Delta files
  def askDeltaFilesOutputDirectory(self):
     deltaFilesDirectoryLocation = tkFileDialog.askdirectory()
     if deltaFilesDirectoryLocation:
        self.deltaEntry.delete(0, 'end')
        self.deltaEntry.insert(0, str(deltaFilesDirectoryLocation))
     return
   #  Queries directory containing files to dump output files
  def askPercentANIOutputDirectory(self):
     percentANIDirectoryLocation = tkFileDialog.askdirectory()
     if percentANIDirectoryLocation:
        self.percentEntry.delete(0, 'end')
        self.percentEntry.insert(0, str(percentANIDirectoryLocation))
     return
   # Print all the values to be saved from the 
  def printTkValues(self):
   print(str(TkANI.allFnaDirectoryLocationString))
   print(str(TkANI.blastDirectoryLocationString))
   print(str(TkANI.mummerDirectoryLocationString))
   print(str(TkANI.deltaFilesDirectoryLocationString))
   print(str(TkANI.percentDirectoryLocationString))
   # Save entry values to be read from other classes
  def saveTkValues(self):
   TkANI.allFnaDirectoryLocationString = self.allFnaEntry.get()
   TkANI.blastDirectoryLocationString = self.blastEntry.get()
   TkANI.outFilesDirectoryLocationString = self.outEntry.get()
   TkANI.mummerDirectoryLocationString = self.mummerEntry.get()
   TkANI.deltaFilesDirectoryLocationString = self.deltaEntry.get()
   TkANI.percentDirectoryLocationString = self.percentEntry.get()
  '''
  Main Buddy
  '''
   # Main from the GUI
  def run(self):
    print("Program started running at time: " + str(datetime.now()))
    self.saveTkValues()
    if self.completeCheckboxValue.get() == 1:
      self.concatinatingFNAFiles.run()
      self.changeFilesToDirectoryName.run()
    if self.cdAllCheckboxValue.get() == 1:
      self.copyFilesInstance.run(TkANI.blastDirectoryLocationString)
      self.copyFilesInstance.run(TkANI.mummerDirectoryLocationString)

    if self.blastCheckboxValue.get() == 1:
      print("b1")
      # datbase setup for blast?
      self.allOutFilesBlast.run()
      '''
      self.blastANI.calculateANI()
      '''
    if self.mummerCheckboxValue.get() == 1:
      print("m1")
      self.allDeltaFilesMummer.run()
      '''
      self.mummerANI.run()
      '''
    if self.blastANICheckboxValue.get() == 1:
      self.blastANI.run()

    if self.mummerANICheckboxValue.get() == 1:
      self.mummerANI.run()
      
    if self.sqlLiteCheckboxValue.get() == 1:
      print("s1")
      self.reformatTextFileForConversion.run()

    if self.jsonCheckboxValue.get() == 1:
      print("j1")
    print("Program completed at time: " + str(datetime.now()))
    self.removeDuplicates.run()

if __name__=='__main__':
   root = Tkinter.Tk()
   root.title("SDSU Blast and MUMMer Average Nucleotide Algorithms GUI v1.0")
   TkANI(root).pack()
   root.mainloop()
