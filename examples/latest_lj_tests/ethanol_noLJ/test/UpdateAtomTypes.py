import os,sys
import numpy as np

"""
This script needs to be run always after updating the LJ params
"""
class Mol:
    def __init__(self,rtfFile,prmFile):
	self.rtfFile = 	rtfFile
	self.prmFile = prmFile
	self.atomTypes = []
	self.newAtomTypes = {}
	self.bonds = []  	
	self.angles = []	
	self.diheds = []
	self.imps = []
	self.nbs =[]


    def getAtomTypes(self):
    	fin=open(self.rtfFile,"r")	
    	for line in fin:
	    token=line.rstrip().split()
	    if len(token)>=3:
	   	if token[0]=="MASS":
		   if token[2] not in self.atomTypes:
		   	self.atomTypes.append(token[2])
	fin.close()

    def reNameAtomTypes(self):
    	for a in self.atomTypes:
	    newAtomType = a+"gp"
	    self.newAtomTypes[a]=newAtomType
 
    def updateRTF(self):
	fin = open(self.rtfFile)
	fout = open("mol.rtf","w")
	for line in fin:
	    token=line.rstrip().split()
	    if len(token)>=3:
		if token[0]=="MASS" or token[0]=="ATOM":	
	    	    line = line.replace(token[2],self.newAtomTypes[token[2]])
	    print>>fout,line.rstrip()
	fout.close() 
	fin.close()
  				
    def getParams(self):
	fin = open(self.prmFile)
	bondFlag=0;angleFlag=0;dihedFlag=0;impFlag=0;nbFlag=0	
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
		elif token[0]=="IMPROPER":
		    bondFlag=0
		    angleFlag=0
		    dihedFlag=0
		    impFlag=1
		elif token[0]=="NONBONDED":dihedFlag=0;impFlag=0
		elif token[0]=="!" and token[1]=="(kcal/mol)":nbFlag=1				
		elif bondFlag==1:
		   self.bonds.append([token[0],token[1],token[2],token[3]])  	
		elif angleFlag==1:
		   self.angles.append([token[0],token[1],token[2],token[3],token[4]])  	
		elif dihedFlag==1:
		   self.diheds.append([token[0],token[1],token[2],token[3],token[4],token[5],token[6]])  	
		elif impFlag==1:
		   self.imps.append([token[0],token[1],token[2],token[3],token[4],token[5],token[6]])  	
		elif nbFlag==1:		
		   self.nbs.append([token[0],token[1],token[2],token[3],token[4],token[5],token[6]])		
	    except IndexError:continue
	fin.close()				


    def writeUpdatedPrm(self):
	fout = open("mol.prm","w")
	print>>fout,"* Force Field Parameter File."
	print>>fout,"*\n"

	print>>fout,"BOND"
	for b in self.bonds:
	    b1 = self.newAtomTypes[b[0]]
	    b2 = self.newAtomTypes[b[1]]
	    k = b[2]
	    b0 = b[3]		
	    fout.write('{:6} {:6} {:6} {:6}\n'.format(b1,b2,k,b0))	
	fout.write("\n")

	print>>fout,"ANGLE"
	for a in self.angles:
	    a1 = self.newAtomTypes[a[0]]
	    a2 = self.newAtomTypes[a[1]]
	    a3 = self.newAtomTypes[a[2]]
	    k = a[3]
	    a0 = a[4]				
	    fout.write('{:6} {:6} {:6} {:8} {:8}\n'.format(a1,a2,a3,k,a0))
	fout.write("\n")

	print>>fout,"DIHEDRAL"
	for d in self.diheds:
	    try:	
		d1 = self.newAtomTypes[d[0]]
		d2 = self.newAtomTypes[d[1]]
		d3 = self.newAtomTypes[d[2]]
		d4 = self.newAtomTypes[d[3]]
		k = d[4];n=d[5];pa=d[6]
		fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(d1,d2,d3,d4,k,n,pa))
	    except KeyError:
		d1 = d[0]
		d2 = self.newAtomTypes[d[1]]
		d3 = self.newAtomTypes[d[2]]
		d4 = d[3]
		k = d[4];n = d[5];pa = d[6]
		fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(d1,d2,d3,d4,k,n,pa))
		
	fout.write("\n")

	print>>fout,"IMPROPER"
	for i in self.imps:
	    try:	
		i1 = self.newAtomTypes[i[0]]
		i2 = self.newAtomTypes[i[1]]
		i3 = self.newAtomTypes[i[2]]
		i4 = self.newAtomTypes[i[3]]
		k = i[4];n = i[5];pa = i[6]
		fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(i1,i2,i3,i4,k,n,pa)) 	
	    except KeyError:
		i1 = i[0]
                i2 = self.newAtomTypes[i[1]]
                i3 = self.newAtomTypes[i[2]]
                i4 = i[3]
                k = i[4];n = i[5];pa = i[6]
		fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(i1,i2,i3,i4,k,n,pa))
	fout.write("\n")

	print>>fout,"NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -"
	print>>fout,"CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4"
	print>>fout,"!                Emin     Rmin/2              Emin/2     Rmin  (for 1-4's)"
	print>>fout,"!             (kcal/mol)    (A)"
	for n in self.nbs:	
	    at = self.newAtomTypes[n[0]]   		
	    fout.write('{:6} {:6} {:12} {:12} {:6} {:12} {:12}\n'.format(at,n[1],n[2],n[3],n[4],n[5],n[6]))	 

	fout.close()

    def runUpdate(self):
	self.getAtomTypes()
   	self.reNameAtomTypes()
   	self.updateRTF()
   	self.getParams()
   	self.writeUpdatedPrm()
   		    

if __name__=="__main__":	
   baseDir = sys.argv[1] #os.getcwd()
   os.chdir("020-initial_parameters")
   os.system("mv mol.prm mol-org.prm")
   os.system("mv mol.rtf mol-org.rtf")				

   myRtf = Mol("mol-org.rtf","mol-org.prm")
   myRtf.runUpdate()

   os.chdir(baseDir)    ### back to $GAAMPRUN	


