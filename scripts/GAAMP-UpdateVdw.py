import os,sys
import numpy as np


atomicMass = {"h":1.008000,"c":12.01000,"o":16.00000,"n":14.00670,"f":18.998403,"s":32.06500,"cl":35.45300,"p":30.97376,"br":79.90400}

def GetBonds(rtf):
    bonds={};group=0
    fin=open(rtf,"r")
    for line in fin:
	token=line.rstrip().split()
	if len(token)>0:
	   if token[0]=="GROUP":group+=1
	   if token[0]=="BOND" and group==1:
	      if not bonds.has_key(token[1]):bonds[token[1]]=[token[2]]
	      else:bonds[token[1]].append(token[2])
	      if not bonds.has_key(token[2]):bonds[token[2]]=[token[1]]
	      else:bonds[token[2]].append(token[1])				

    fin.close()
    return bonds		


def Get_C3_CA_HC_Type(bonds,atoms,allAtomTypes):
    allAtomTypes = {atoms[i]:allAtomTypes[i] for i in range(len(atoms))}
    specialAtomTypes={}
    for atom in bonds:
	hCount=0	
	if allAtomTypes[atom]=="ca":
	   specialAtomTypes[atom]="ca"
	   for at in bonds[atom]:
               if allAtomTypes[at]=="c3": ###sp3 carbon attached to aro-carb
                  specialAtomTypes[at]="c3r"
               for a in bonds[at]:
                   if allAtomTypes[a]=="hc":specialAtomTypes[a]="hc3r" ### H attached to c3r carb	
	if allAtomTypes[atom]=="c3" and not specialAtomTypes.has_key(atom): ### making sure c3r is not reassigned again
	   for at in bonds[atom]:
		#if allAtomTypes[at]=="hc":hCount+=1 	
		if allAtomTypes[at][0]=="h":hCount+=1 	
	   if hCount==0:specialAtomTypes[atom]="c30"	
	   elif hCount==1:
                specialAtomTypes[atom]="c31"
                for at in bonds[atom]:
                    if allAtomTypes[at]=="hc":specialAtomTypes[at]="hc31"
	   elif hCount==2:
                specialAtomTypes[atom]="c32"
                for at in bonds[atom]:
                    if allAtomTypes[at]=="hc":specialAtomTypes[at]="hc32" 	
	   elif hCount==3:
                specialAtomTypes[atom]="c33"
                for at in bonds[atom]:
                    if allAtomTypes[at]=="hc":specialAtomTypes[at]="hc33"	
    return specialAtomTypes


def updateRTF(rtf,atoms,allAtomTypes):
    atomTypes = {v:k for k,v in allAtomTypes.iteritems()}
    fout = open("tmp.rtf","w")
    idx=1
     
    charge = [c.rstrip().split()[2] for c in open(rtf,"r") if len(c.rstrip().split())>0 and c.rstrip().split()[0]=="RESI" ][0]
    pCharge = {c.rstrip().split()[1]:c.rstrip().split()[3] for c in open(rtf,"r") if len(c.rstrip().split())>0 and c.rstrip().split()[0]=="ATOM" }

    print>>fout,"* Topology File."
    print>>fout,"*"
    print>>fout,"   99   1"
    for a in atomTypes:
	if a[:2]=="cl":mass = atomicMass["cl"]
	elif a[:2]=="br":mass = atomicMass["br"]
	elif a[0]=="c" : mass = atomicMass["c"]
	elif a[0]=="h" : mass = atomicMass["h"]
	elif a[0]=="o" : mass = atomicMass["o"]
	else: mass = atomicMass[a]
	fout.write('{:9} {:4} {:8} {:10}\n'.format("MASS",idx,a,mass))	
	idx+=1

    print>>fout,"\n","RESI MOL ",charge 
    print>>fout,"GROUP"       
    i=0 
    for a in atoms:
	fout.write('{:5} {:6} {:7} {:10}\n'.format("ATOM",a,allAtomTypes[a],pCharge[a]))
	i+=1

    print>>fout,"\n" 

    fin=open(rtf,"r")
    angleFlag=dihedFlag=impFlag=0
    for line in fin:
	token=line.rstrip().split()
	if len(token)>0:
	   if token[0]=="BOND":print>>fout,line.rstrip()
	   if token[0]=="ANGL":
		if angleFlag==0:print>>fout,"\n",line.rstrip();angleFlag+=1
		else:print>>fout,line.rstrip()
	   if token[0]=="DIHE":
		if dihedFlag==0:print>>fout,"\n",line.rstrip();dihedFlag+=1
		else:print>>fout,line.rstrip()
	   if token[0]=="IMPH":
		if impFlag==0:print>>fout,"\n",line.rstrip();impFlag+=1
		else:print>>fout,line.rstrip()
    print>>fout,"\n"
    fout.close() 	 	  
    fin.close()

 
def getAllAtomTypes(rtf):
    atoms=[]
    atomTypes=[]
    #atomTypes={}
    group=0
    fin=open(rtf,"r")	
    for line in fin:
	token=line.rstrip().split()
    	if len(token)>0:
	   if token[0]=="GROUP":group+=1
	   if token[0]=="ATOM" and group==1:
	       #atomTypes[token[1]]=token[2]
	       atoms.append(token[1])	
	       atomTypes.append(token[2])
    return atoms,atomTypes


def readPRM(prm):
    fin=open(prm,"r")
    bondFlag=0;angleFlag=0;dihedFlag=0;impropFlag=0;nbFlag=0;
    bondRecord={};angleRecord={};dihedRecord={};impropRecord={};nbRecord={};
    for line in fin:
	token=line.rstrip().split()
	if len(token)>0:
	   if token[0]=="BOND":
	      bondFlag=1
	   elif token[0]=="ANGLE":
	      bondFlag=0;angleFlag=1
	   elif token[0]=="DIHEDRAL":
	      angleFlag=0;dihedFlag=1
	   elif token[0]=="IMPROPER":
	      dihedFlag=0;impropFlag=1
	   elif token[0]=="NONBONDED":
	      dihedFlag=0;impropFlag=0;nbFlag=1
	   if bondFlag==1 and token[0]!="BOND":	
	      bondRecord[token[0]+"-"+token[1]] = [token[2],token[3]]
	   if angleFlag==1 and token[0]!="ANGLE":
	      angleRecord[token[0]+"-"+token[1]+"-"+token[2]]= [token[3],token[4]]
	   if dihedFlag==1 and token[0]!="DIHEDRAL":
	      if not dihedRecord.has_key(token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]):	
	     	 dihedRecord[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]] = [token[4],token[5],token[6]]	
	      else:dihedRecord[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]].append([token[4],token[5],token[6]])	
	   if impropFlag==1 and token[0]!="IMPROPER":
	      try:
	          impropRecord[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]] = [token[4],token[5],token[6]]	
	      except IndexError:continue	
	   if nbFlag==1 and token[1]!="(KCAL/MOL)":	
	      if token[0]!="NONBONDED" and token[0]!="!" and token[0]!="CUTNB":
	             nbRecord[token[0]] = [token[1],token[2],token[3]]
    fin.close() 
    return bondRecord,angleRecord,dihedRecord,impropRecord,nbRecord


def combineAtomTypes(specialTypes,atoms,allTypes):
    allTypes = {atoms[i]:allTypes[i] for i in range(len(atoms))}
    atomTypes={}
    for a in allTypes:
	if specialTypes.has_key(a):atomTypes[a]=specialTypes[a]
	else: atomTypes[a] = allTypes[a]
    return atomTypes

def caps(atName):
    return ''.join([i.capitalize() for i in atName])

def updateBondPrm(orgBond,atomTypes,bondPrm):
    newBondParams={}
    for b in orgBond:
        for a in orgBond[b]:
	    bto = atomTypes[b] 
	    ato = atomTypes[a]	
	    if bto[:2]=="c3":bt="c3"
	    elif bto[:2]=="hc":bt="hc"
	    else:bt=bto	
	    if ato[:2]=="c3":at="hc"
	    elif ato[:2]=="hc":at="hc"
	    else:at=ato
	    if bondPrm.has_key(bt+"-"+at):
		k = bondPrm[bt+"-"+at][0]; b0 = bondPrm[bt+"-"+at][1]
		if not newBondParams.has_key(bto+"-"+ato) and not newBondParams.has_key(ato+"-"+bto):
		   newBondParams[bto+"-"+ato] = [k,b0]  	
	    elif bondPrm.has_key(at+"-"+bt):	
		k = bondPrm[at+"-"+bt][0]; b0 = bondPrm[at+"-"+bt][1]
		if not newBondParams.has_key(bto+"-"+ato) and not newBondParams.has_key(ato+"-"+bto):
		   newBondParams[bto+"-"+ato] = [k,b0]  	
	    else:
		if bondPrm.has_key(bt+"-"+at):
		    k = bondPrm[bt+"-"+at][0]; b0 = bondPrm[bt+"-"+at][1]
		    if not newBondParams.has_key(bt+"-"+at) and not newBondParams.has_key(at+"-"+bt):
		       newBondParams[bto+"-"+ato] = [k,b0]  	
		elif bondPrm.has_key(at+"-"+bt):
		    k = bondPrm[at+"-"+bt][0]; b0 = bondPrm[at+"-"+bt][1]
		    if not newBondParams.has_key(bt+"-"+at) and not newBondParams.has_key(at+"-"+bt):
		       newBondParams[ato+"-"+bto] = [k,b0]  	
	    
    return newBondParams
    	
def updateAnglePrm(atomTypes,anglePrm):
    atomTypes = {v:k for k,v in atomTypes.iteritems()}
    #print anglePrm
    #print atomTypes
    newAnglePrm={}
    for t in anglePrm:
        A1=[];A2=[];A3=[] ### stores new atom types for each atom types in angle prm
	at = t.split("-")	
	if at[0]=="c3": A1 = [c3 for c3 in atomTypes if c3[:2]=="c3"] ### adds c3 types if any
	elif at[0]=="hc": A1 = [hc for hc in atomTypes if hc[:2]=="hc"] ### ads hc types if any
	else:A1=[at[0]] ### keeps the gaff type
	if at[1]=="c3": A2 = [c3 for c3 in atomTypes if c3[:2]=="c3"]
	elif at[1]=="hc": A2 = [hc for hc in atomTypes if hc[:2]=="hc"]   	
	else:A2=[at[1]]
	if at[2]=="c3": A3 = [c3 for c3 in atomTypes if c3[:2]=="c3"]
        elif at[2]=="hc": A3 = [hc for hc in atomTypes if hc[:2]=="hc"]
        else:A3=[at[2]]	
	for a1 in A1:
	    for a2 in A2:
		for a3 in A3:	
		    newAnglePrm[a1+"-"+a2+"-"+a3]=[anglePrm[t][0],anglePrm[t][1]]
	  			
        #break
    return newAnglePrm

"""
Updates atom types for both dihed and improper prms
"""
def updateDihedPrm(atomTypes,dihedPrm):
    atomTypes = {v:k for k,v in atomTypes.iteritems()}
    newDihedPrm={}
    #print dihedPrm	
    for t in dihedPrm:	
	D1=[];D2=[];D3=[];D4=[]
	at = t.split("-")
	if at[0]=="c3": D1 = [c3 for c3 in atomTypes if c3[:2]=="c3"]
	elif at[0]=="hc": D1 = [hc for hc in atomTypes if hc[:2]=="hc"]
	else:D1=[at[0]]
	if at[1]=="c3": D2 = [c3 for c3 in atomTypes if c3[:2]=="c3"]
	elif at[1]=="hc": D2 = [hc for hc in atomTypes if hc[:2]=="hc"]
	else:D2=[at[1]]
	if at[2]=="c3": D3 = [c3 for c3 in atomTypes if c3[:2]=="c3"]
	elif at[2]=="hc": D3 = [hc for hc in atomTypes if hc[:2]=="hc"]
	else:D3=[at[2]]
	if at[3]=="c3": D4 = [c3 for c3 in atomTypes if c3[:2]=="c3"]
	elif at[3]=="hc": D4 = [hc for hc in atomTypes if hc[:2]=="hc"]
	else:D4=[at[3]]
	for d1 in D1:
	    for d2 in D2:
		for d3 in D3:
		    for d4 in D4:
		    	newDihedPrm[d1+"-"+d2+"-"+d3+"-"+d4]=[dihedPrm[t][0],dihedPrm[t][1],dihedPrm[t][2]]
    return newDihedPrm


def updateNonBond(atomTypes,nbPrm):
    atomTypes = {v:k for k,v in atomTypes.iteritems()}
    newNBPrm={}	
    for t in nbPrm:
	if t=="c3":
	   for c3 in atomTypes:
		if c3[:2]=="c3":newNBPrm[c3]=nbPrm[t]
	elif t=="hc":
	   for hc in atomTypes:
		if hc[:2]=="hc":newNBPrm[hc]=nbPrm[t] 	
	else:newNBPrm[t]=nbPrm[t]
    return newNBPrm	

def readOptLJ(ljFile):
    lj={}
    fin=open(ljFile,"r")	
    for line in fin:
	token = line.rstrip().split()
	lj[token[0]] = [token[1],token[2]]
    return lj	


def writeUpdatedPrm(newBondPrm,newAnglePrm,newDihedPrm,newImpPrm,newNBPrm,ljOptPrm):
    fout = open("tmp.prm","w") 	
    print>>fout,"* FORCE FIELD PARAMETER FILE."
    print>>fout,"*","\n"
    print>>fout,"BONDS"
    for b in newBondPrm:
        b1=b.split("-")[0];b2=b.split("-")[1]
	k = newBondPrm[b][0];b0=newBondPrm[b][1]
	fout.write('{:6} {:6} {:6} {:6}\n'.format(b1,b2,k,b0))
	#print>>fout,b1,b2,newBondPrm[b][0],newBondPrm[b][1]
    print>>fout,"\n","ANGLES"
    for a in newAnglePrm:
	a1=a.split("-")[0];a2=a.split("-")[1];a3=a.split("-")[2]
	k = newAnglePrm[a][0];a0 = newAnglePrm[a][1]
	fout.write('{:6} {:6} {:6} {:8} {:8}\n'.format(a1,a2,a3,k,a0))
	#print>>fout,a1,a2,a3,newAnglePrm[a][0],newAnglePrm[a][1]
    print>>fout,"\n","DIHEDRALS"
    for d in newDihedPrm:
	d1=d.split("-")[0];d2=d.split("-")[1];d3=d.split("-")[2];d4=d.split("-")[3]
	k = newDihedPrm[d][0]; n = newDihedPrm[d][1]; pa = newDihedPrm[d][2] 
	fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(d1,d2,d3,d4,k,n,pa))
	#print>>fout,d1,d2,d3,d4,newDihedPrm[d][0],newDihedPrm[d][1],newDihedPrm[d][2]
    print>>fout,"\n","IMPROPERS"
    for im in newImpPrm:
	im1=im.split("-")[0];im2=im.split("-")[1];im3=im.split("-")[2];im4=im.split("-")[3]
	k = newImpPrm[im][0]; n = newImpPrm[im][1]; pa = newImpPrm[im][2]
        fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(im1,im2,im3,im4,k,n,pa))
	#print>>fout,im1,im2,im3,im4,newImpPrm[im][0],newImpPrm[im][1],newImpPrm[im][2]
    print>>fout,"\n","NONBONDED  E14FAC  1.000000"
    print>>fout,"!                EMIN     RMIN/2              EMIN/2     RMIN  (FOR 1-4'S)"
    print>>fout,"!             (KCAL/MOL)    (A)"
    for nb in newNBPrm:
	if ljOptPrm.has_key(caps(nb)):
	   z = newNBPrm[nb][0]; e=ljOptPrm[caps(nb)][0]; r= ljOptPrm[caps(nb)][1]
	else:
	   z = newNBPrm[nb][0]; e=newNBPrm[nb][1]; r= newNBPrm[nb][2]
	fout.write('{:6} {:6} {:12} {:12} {:6} {:12} {:12}\n'.format(nb,z,e,r,z,float(e)/2,r))
	#print>>fout,nb,newNBPrm[nb][0],newNBPrm[nb][1],newNBPrm[nb][2],newNBPrm[nb][0],\
	#	float(newNBPrm[nb][1])/2,newNBPrm[nb][2]
	
    fout.write("\n")
    fout.close()	


if __name__=="__main__":
  
    baseDir = sys.argv[1] ## same as $GAAMPRUN i.e. gaamp run dir   
    os.chdir("020-initial_parameters") ## working on 020-initial_parameters dir      

    os.system("cp mol.prm mol-org.prm") ## copying original prm & rtf for safety
    os.system("cp mol.rtf mol-org.rtf")

    rtf= sys.argv[2]#"/lcrc/project/Drude/chetan/OLD_GAAMP/reparameterize/m_017/01-ac/mol.rtf"
    prm= sys.argv[3]#"/lcrc/project/Drude/chetan/OLD_GAAMP/reparameterize/m_017/01-ac/mol.prm"
    bonds = GetBonds(rtf)
    #print bonds;print "\n"
    atoms,allAtomTypes = getAllAtomTypes(rtf)
    #print atoms;print allAtomTypes;print "\n"
    specialAtomTypes = Get_C3_CA_HC_Type(bonds,atoms,allAtomTypes)
    #print specialAtomTypes;print "\n"	
    atomTypes = combineAtomTypes(specialAtomTypes,atoms,allAtomTypes)
    updateRTF(rtf,atoms,atomTypes)
    bondPrm,anglePrm,dihedPrm,impropPrm,nbPrm = readPRM(prm)
    #print "bondPrm";print bondPrm
    newBondPrm = updateBondPrm(bonds,atomTypes,bondPrm) 
    #print newBondPrm
    newAnglePrm = updateAnglePrm(atomTypes,anglePrm)
    #print "Dihed"
    newDihedPrm = updateDihedPrm(atomTypes,dihedPrm)	
    #print "Improp"		
    newImpPrm = updateDihedPrm(atomTypes,impropPrm)
    newNBPrm = updateNonBond(atomTypes,nbPrm)
    #print newNBPrm
    ljPath = sys.argv[0].split("GAAMP-UpdateVdw")[0]
    ljOptPrm = readOptLJ(ljPath+"LJPARAMETERS.dat")	
    writeUpdatedPrm(newBondPrm,newAnglePrm,newDihedPrm,newImpPrm,newNBPrm,ljOptPrm)

    os.system("mv tmp.rtf mol.rtf") ##lets swap new rtf & prm for further uses
    os.system("mv tmp.prm mol.prm")

    os.chdir(baseDir)    ### back to $GAAMPRUN


