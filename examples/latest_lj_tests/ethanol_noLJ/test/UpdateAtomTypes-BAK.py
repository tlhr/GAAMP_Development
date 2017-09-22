import os,sys
import numpy as np

class Mol:
    def __init__(self,rtfFile,prmFile):
	self.rtfFile = 	rtfFile
	self.prmFile = prmFile
	self.atomTypes = []
	self.newAtomTypes = {}  	

    def getAtomTypes(self):
    	fin=open(self.rtfFile,"r")	
    	for line in fin:
	    token=line.rstrip().split()
	    if len(token)>=3:
	   	if token[0]=="MASS":
		   if token[2] not in self.atomTypes:
		   	self.atomTypes.append(token[2])

    def reNameAtomTypes(self):
    	for a in self.atomTypes:
	    newAtomType = a+"gp"
	    self.newAtomTypes[a]=newAtomType
 
    def updateRTF(self):
	fin = open(self.rtfFile)
	fout = open("new.rtf","w")
	for line in fin:
	    token=line.rstrip().split()
	    if len(token)>=3:
		if token[0]=="MASS" or token[0]=="ATOM":	
	    	    line = line.replace(token[2],self.newAtomTypes[token[2]])
	    print>>fout,line.rstrip()
	fout.close() 
  				
    def updatePRM(self):
	fin = open(self.prmFile)
	fout = open("new.prm","w")
	bondFlag=0;angleFlag=0;dihedFlag=0;nbFlag=0	
	for line in fin:
	    token = line.rstrip().split()
	    try:	
	    	if token[0]=="BOND":
		    bondFlag=1
		elif token[0]=="ANGLE":
		    bondFlag=0
		    angleFlag=1
		elif token[0]=="DIHEDRAL":
		    bondFlag=0
		    angleFlag=0
		    dihedFlag=1
		elif token[0]=="NONBONDED":dihedFlag=0
		elif token[0]=="!" and token[1]=="(kcal/mol)":nbFlag=1				
		elif bondFlag==1:
		   #line = line.replace(token[0],self.newAtomTypes[token[0]])
		   try:	
		   	line = line.replace(token[0],self.newAtomTypes[token[0]])
		   	line = line.replace(" "+token[1]," "+self.newAtomTypes[token[1]])
		   except KeyError:continue
		elif angleFlag==1:
		   #line = line.replace(token[0],self.newAtomTypes[token[0]])
                   try:
		   	line = line.replace(token[0],self.newAtomTypes[token[0]])
                        line = line.replace(" "+token[1]," "+self.newAtomTypes[token[1]])
                        line = line.replace(" "+token[2]," "+self.newAtomTypes[token[2]])
                   except KeyError:continue
		elif dihedFlag==1:
                   try:
		   	line = line.replace(token[0],self.newAtomTypes[token[0]])
                        line = line.replace(" "+token[1],self.newAtomTypes[token[1]])
                        line = line.replace(" "+token[2],self.newAtomTypes[token[2]])
                        line = line.replace(" "+token[3],self.newAtomTypes[token[3]])
                   except KeyError:continue
		elif nbFlag==1:		
		   try:	
		   	line = line.replace(token[0],self.newAtomTypes[token[0]])
		   except KeyError:continue			
	
	    except IndexError:fout,line.rstrip()
	    print>>fout,line.rstrip()
	fout.close()				


if __name__=="__main__":
   baseDir = os.getcwd() #sys.argv[1]	

   myRtf = Mol(baseDir+"/mol.rtf",baseDir+"/mol.prm")
   myRtf.getAtomTypes()		
   myRtf.reNameAtomTypes()
   myRtf.updateRTF()
   myRtf.updatePRM()		

   #os.chdir(baseDir)    ### back to $GAAMPRUN	


