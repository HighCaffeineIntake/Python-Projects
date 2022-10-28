import json
import copy

def ImportEverything():
    GNRA = Import("C:\\Users\\admin\Downloads\\HL_43074.16.json")
    GNRA1= Import("C:\\Users\\admin\\Downloads\\HL_02361.6 (1).json")
    UNCG = Import("C:\\Users\\admin\\Downloads\\HL_79299.10.json")
    TLoops1 = Import("C:\\Users\\admin\\Downloads\\HL_68493.6.json")
    TLoops2 = Import("C:\\Users\\admin\\Downloads\\HL_33239.15.json")
    return GNRA,GNRA1,UNCG,TLoops1,TLoops2

def Import(x):
    with open(x,"r") as file:
        JSON=json.load(file)
    return JSON

def Reverse(x):
    x2=copy.copy(x)
    x2.reverse()
    return x2

def CheckBPContinuityGNRA(x): #GNRA Continuity can also be used for T-Loops as they also have the same continuity criteria
    DiffList=[]
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    DiffSet=set(DiffList)
    try: 
        DiffSet.remove(1)
        if len(DiffSet) == 0:
            return True
        if len(DiffSet) != 0:
            return False
    except:
        return False

def CheckBPContinuityUNCG(x):
    DiffList=[]
    Checklist=[1,2,1,1]
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    Comparison = [i for i,j in zip(Checklist,DiffList) if i == j]
    if len(Comparison) == len(Checklist):
        return True
    else:
        return False

def CheckBPContinuityGNRA1(x):
    DiffList=[]
    Checklist=[1,1,1,2,1]
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    Comparison = [i for i,j in zip(Checklist,DiffList) if i == j]
    if len(Comparison) == len(Checklist):
        return True
    else:
        return False

def MakeGNRASeqs(x):
    GNRASeqs=[]
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        Continuity=False
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])
        Continuity=CheckBPContinuityGNRA(Basepositions)
        if Continuity == True:
            Seq=Seq[1:-1]
            RevSeq = copy.copy(Seq)
            RevSeq.reverse()
            Seq="".join(Seq)
            RevSeq="".join(RevSeq)
            GNRASeqs.append(Seq)
            GNRASeqs.append(RevSeq)
    return GNRASeqs

def MakeUNCGSeqs(x): #Both SeqDicts functions create their respective sequence dictionary by assiging IDs to their Sequence, this is going to be the GNRA Version
    UNCGSeqs=[]
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])
        Continuity=CheckBPContinuityUNCG(Basepositions)
        if Continuity == True:
            Seq=Seq[1:-1]
            SeqList=Insert_N(Seq,1,"N")
            for entry in SeqList:
                UNCGSeqs.append(entry)
    return UNCGSeqs

def Insert_N(temp,Pos_N,Type_N): #Full Insert Function, differentiates between all types of Nucleotide and possibility, Insert @Pos-1
    SeqList=[]
    if Type_N == "Y":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "R":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
    elif Type_N == "W":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
    elif Type_N == "S":
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "K":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
    elif Type_N == "M":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "D":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
    elif Type_N == "V":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "H":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "B":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "A":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
    elif Type_N == "U":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
    elif Type_N == "C":
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
    elif Type_N == "G":
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
    elif Type_N =="N":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempUR=copy.copy(tempU)
        tempUR.reverse()
        tempU2="".join(tempU)
        tempUR2="".join(tempUR)
        SeqList.append(tempU2)
        SeqList.append(tempUR2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempGR=copy.copy(tempG)
        tempGR.reverse()
        tempG2="".join(tempG)
        tempGR2="".join(tempGR)
        SeqList.append(tempG2)
        SeqList.append(tempGR2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempCR=copy.copy(tempC)
        tempCR.reverse()
        tempC2="".join(tempC)
        tempCR2="".join(tempCR)
        SeqList.append(tempC2)
        SeqList.append(tempCR2)
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempAR=copy.copy(tempA)
        tempAR.reverse()
        tempA2="".join(tempA)
        tempAR2="".join(tempAR)
        SeqList.append(tempA2)
        SeqList.append(tempAR2)
    else:print("No nucleotide type given")
    return SeqList

def MakeGNRA1Seqs(x):
    GNRA1Seqs=[]
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])
        Continuity=CheckBPContinuityGNRA1(Basepositions)
        if Continuity == True:
            Seq=Seq[1:-1]
            SeqList=Insert_N(Seq,3,"N")
            for entry in SeqList:
                GNRA1Seqs.append(entry)
    return GNRA1Seqs

def MakeTLoops1Seqs(x): #TLoops with one buldged base
    TLoops1Seqs=[]
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])
        Continuity=CheckBPContinuityTLoops1(Basepositions)
        if Continuity == True:
            Seq=Seq[1:-1]
            SeqList=Insert_N(Seq,6,"N")
            for entry in SeqList:
                TLoops1Seqs.append(entry)
    return TLoops1Seqs
             
def CheckBPContinuityTLoops1(x):
    DiffList=[]
    Checklist=[1,1,1,1,1,2]
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    Comparison = [i for i,j in zip(Checklist,DiffList) if i == j]
    if len(Comparison) == len(Checklist):
        return True
    else:
        return False

def MakeTLoops2Seqs(x):
    TLoops2Seqs=[]
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])
        Continuity=CheckBPContinuityTLoops2(Basepositions)
        if Continuity == True:
            Seq=Seq[1:-1]
            RevSeq = copy.copy(Seq)
            RevSeq.reverse()
            Seq="".join(Seq)
            RevSeq="".join(RevSeq)
            TLoops2Seqs.append(Seq)
            TLoops2Seqs.append(RevSeq)
    return TLoops2Seqs 

def CheckBPContinuityTLoops2(x):
    DiffList=[]
    Checklist=[1,1,1,1,1,1,1,1]
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    Comparison = [i for i,j in zip(Checklist,DiffList) if i == j]
    if len(Comparison) == len(Checklist):
        return True
    else:
        return False

if __name__ == "__main__":
    Biglist=[]
    GNRA,GNRA1,UNCG,TLoops1,TLoops2 = ImportEverything()
    GNRASeqs=MakeGNRASeqs(GNRA)
    GNRA1Seqs=MakeGNRA1Seqs(GNRA1)
    UNCGSeqs=MakeUNCGSeqs(UNCG)
    TLoops1Seqs = MakeTLoops1Seqs(TLoops1)
    TLoops2Seqs = MakeTLoops2Seqs(TLoops2)
    with open("HairpinMotifs.csv","w") as file:
        for entry in GNRASeqs:
            file.write(entry+"+4\n")
        for entry in GNRA1Seqs:
            file.write(entry+"+4\n")
        for entry in UNCGSeqs:
            file.write(entry+"+3\n")
        for entry in TLoops1Seqs:
            file.write(entry+"+2\n")
        for entry in TLoops2Seqs:
            file.write(entry+"+2\n")