#!/usr/bin/python

'''
The MIT License (MIT)

Copyright (c) 2013 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#bamTableUpdate.py
#131216
#Jaime


#Description:

#when called updates the /grail/projects/masterBamTable.txt



#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys
import subprocess
import string

#print "Using python version %s" % sys.version


#importing utils package
sys.path.append('/home/cl512/src/pipeline/')
import utils



#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


dataFile ='/location/file.txt'
genome = 'hg18'


#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

#write your specific functions here


def getUniqueIDList():

    '''
    function that gets all uniqueIDs
    '''
    
    cmd = "mysql -u youngcompread --password='Fg78$Dr' -e 'SELECT uniqueID FROM baseExp' seqDB"

    sqlOut = subprocess.Popen(cmd,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)

    sqlText = sqlOut.communicate()
    sqlText =sqlText[0]
    uniqueIDList =sqlText.split('\n')[1:-1]
    return uniqueIDList



def getTonyInfo(uniqueIDList,colList):

    '''
    pass this a uniqueID List and a list of columns

    '''

    uniqueIDString = string.join(uniqueIDList,',')

    columnString = string.join([str(x) for x in colList],',')

    cmd = "perl /ark/tony/admin/getDB_Data.pl -i %s -c %s -o TAB" % (uniqueIDString,columnString)
    
    sqlOut = subprocess.Popen(cmd,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)

    sqlText = sqlOut.communicate()

    sqlText = sqlText[0]
    
    sqlTable = sqlText.split('\n')
    sqlTable = [x for x in sqlTable if len(x) > 0]

    sqlTable = [x.split('\t') for x in sqlTable]

    header = [x.split(':')[-1] for x in sqlTable[0][1:]]
    header= [str.upper(x) for x in header]
    header = ['GENOME', 'SOURCE', 'CELL_TYPE', 'NAME', 'BAMFILE']
    tonyDict = {}
    for line in sqlTable[1:]:
        uniqueID = line[0]
        tonyDict[uniqueID] = {}
        for i in range(len(header)):
            tonyDict[uniqueID][header[i]] = line[(i+1)]
    newTable = []        
    newTable.append(header)

    for key in tonyDict.keys():
        newLine = []
        newLine.append(str.upper(tonyDict[key]['GENOME']))
        newLine.append(tonyDict[key]['SOURCE'])
        newLine.append(tonyDict[key]['CELL_TYPE'])
        newLine.append(tonyDict[key]['NAME'])
        newLine.append(tonyDict[key]['BAMFILE'])
        newTable.append(newLine)

    #print newTable
    
    utils.unParseTable(newTable, '/grail/projects/masterBamTable.txt', '\t')
    
    


#================================================================================
#===============================MAIN RUN=========================================
#================================================================================

#write the actual script here


def main():

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''
    colList = [48,6,7,3,47]
    
    #uniqueIDList = ['20130724_73','20130726_83']
    uniqueIDList = getUniqueIDList()
    getTonyInfo(uniqueIDList,colList)


main()

