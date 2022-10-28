import json
import copy
import glob

def ImportEverything():
    KTurns = Import("C:\\Users\\admin\\Downloads\\IL_29549.8.json")
    CLoops = Import("C:\\Users\\admin\\Downloads\\IL_63596.3.json")
    SR_A= Import("C:\\Users\\admin\\Downloads\\IL_02349.4.json")
    return KTurns, CLoops,SR_A

def Import(x):
    with open(x,"r") as file:
        JSON=json.load(file)
    return JSON

def MakeKTurnsSeqs(x):
    Seqs=[]
    FullSeqList=[]
    Break = int(x["chainbreak"])
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])#Split the sequence and the basepositions here and analyze each one
        SeqA = Seq[0:Break]
        SeqB = Seq[Break:]
        BasepositionsA = Basepositions[0:Break]
        BasepositionsB = Basepositions[Break:]
        ContinuityA = CheckBPContinuityKTurnsA(BasepositionsA)
        if ContinuityA == True:
            ContinuityB = CheckBPContinuityKTurnsB(BasepositionsB)
            if ContinuityB == True:
                SeqA2=SeqA[1:-1]
                SeqB2=SeqB[1:-1]
                SeqList=Insert_N(SeqB2,2,"N")
                SeqA2="".join(SeqA2)
                for entry2 in SeqList:
                    FullSeq=SeqA2+"$"+entry2
                    FullSeqList.append(FullSeq)
                for entry3 in FullSeqList:
                    RevSeqList=list(entry3)
                    RevSeq=Reverse(RevSeqList)
                    RevSeq="".join(RevSeq)
                    Seqs.append(entry3)
                    Seqs.append(RevSeq)
    return Seqs

def CheckBPContinuityKTurnsA(x):
    DiffList=[]
    Checklist=[1,1,1,1]
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

def CheckBPContinuityKTurnsB(x):
    DiffList=[]
    Checklist=[1,1,2,1,1,1]
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

def Insert_N(temp,Pos_N,Type_N): #Full Insert Function, differentiates between all types of Nucleotide and possibility, Insert @Pos-1,No Reverses!
    SeqList=[]
    if Type_N == "Y":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "R":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
    elif Type_N == "W":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
    elif Type_N == "S":
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "K":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
    elif Type_N == "M":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "D":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
    elif Type_N == "V":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "H":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "B":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "A":
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
    elif Type_N == "U":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
    elif Type_N == "C":
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
    elif Type_N == "G":
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
    elif Type_N =="N":
        tempU=copy.copy(temp)
        tempU.insert(Pos_N,"U")
        tempU2="".join(tempU)
        SeqList.append(tempU2)
        tempG=copy.copy(temp)
        tempG.insert(Pos_N,"G")
        tempG2="".join(tempG)
        SeqList.append(tempG2)
        tempC=copy.copy(temp)
        tempC.insert(Pos_N,"C")
        tempC2="".join(tempC)
        SeqList.append(tempC2)
        tempA=copy.copy(temp)
        tempA.insert(Pos_N,"A")
        tempA2="".join(tempA)
        SeqList.append(tempA2)
    else:print("No nucleotide type given")
    return SeqList

def CheckBPContinuityCLoopsA(x):
    DiffList=[]
    Checklist=[1,1,1,1]
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

def CheckBPContinuityCLoopsB(x):
    DiffList=[]
    Checklist=[2]
    Checklist2=[3]
    UnnecessaryVAL=0
    ChecklistVAL=2
    Checklist2VAL=3
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    Comparison = [i for i,j in zip(Checklist,DiffList) if i == j]
    if len(Comparison) == len(Checklist):
        return True,ChecklistVAL
    else:
        Comparison2 = [i for i,j in zip(Checklist2,DiffList) if i == j]
        if len(Comparison2) == len(Checklist2):
            return True,Checklist2VAL
        else:
            return False,UnnecessaryVAL

def Reverse(x):
    x2=copy.copy(x)
    x2.reverse()
    return x2

def MakeCLoopsSeqs(x):
    Seqs=[]
    FullSeqList=[]
    Break = int(x["chainbreak"])
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])#Split the sequence and the basepositions here and analyze each one
        SeqA = Seq[0:Break]
        SeqB = Seq[Break:]
        BasepositionsA = Basepositions[0:Break]
        BasepositionsB = Basepositions[Break:]
        ContinuityA = CheckBPContinuityCLoopsA(BasepositionsA)
        if ContinuityA == True:
            ContinuityB,BVAL = CheckBPContinuityCLoopsB(BasepositionsB)
            if ContinuityB == True:
                SeqA2=SeqA[1:-1]
                SeqB2=SeqB[1:-1]
                if BVAL == 2:
                    SeqList=Insert_N(SeqB2,0,"N")
                if BVAL == 3:
                    SeqList2=Insert_N(SeqB2,0,"N")
                    for entry in SeqList2:
                        entry2=list(entry)
                        SeqList=Insert_N(entry2,1,"N")
                SeqA2="".join(SeqA2)
                for entry2 in SeqList:
                    FullSeq=SeqA2+"$"+entry2
                    FullSeqList.append(FullSeq)
                for entry3 in FullSeqList:
                    RevSeqList=list(entry3)
                    RevSeq=Reverse(RevSeqList)
                    RevSeq="".join(RevSeq)
                    Seqs.append(entry3)
                    Seqs.append(RevSeq)
    return Seqs

def ImportSarcinRicinMotifs():
    SRList=[]
    for name in glob.glob("C:/Master/moRNAO/MotifFolder/*"):
        f = open(name,"r")
        file= json.load(f)
        SRList.append(file)
    return SRList

def MakeSRSeqs(x):
    Seqs=[]
    for entry in x:
        FullSeqList=[]
        Break = int(entry["chainbreak"])
        for entry2 in entry["alignment"]:
            Seq=[]
            Basepositions=[]
            for Nucleotide in entry["alignment"][entry2]:
                Datalist=Nucleotide.split("|")
                Seq.append(Datalist[3])
                Basepositions.append(Datalist[4])#Split the sequence and the basepositions here and analyze each one
            SeqA = Seq[0:Break]
            SeqB = Seq[Break:]
            BasepositionsA = Basepositions[0:Break]
            BasepositionsB = Basepositions[Break:]
            ContinuityA = CheckBPContinuitySR(BasepositionsA)
            if ContinuityA == True:
                ContinuityB = CheckBPContinuitySR(BasepositionsB)
                if ContinuityB == True:
                    SeqA2=SeqA[1:-1]
                    SeqB2=SeqB[1:-1]
                    SeqA2="".join(SeqA2)
                    SeqB2="".join(SeqB2)
                    FullSeq=SeqA2+"$"+SeqB2
                    FullSeqList.append(FullSeq)
                for FullSeq in FullSeqList:
                    FullSeq2=list(FullSeq)
                    RevSeq=Reverse(FullSeq2)
                    RevSeq="".join(RevSeq)
                    Seqs.append(FullSeq)
                    Seqs.append(RevSeq)
    return Seqs

def CheckBPContinuitySR(x):
    DiffList=[]
    for i in range(len(x)-1):
        Base_Position=int(x[i])
        Next_Base=int(x[i+1])
        Diff = Next_Base - Base_Position
        DiffList.append(Diff)
    DiffSet=set(DiffList)
    DiffSet.remove(1)
    if len(DiffSet) == 0:
        return True
    else:
        return False

def MakeSR_ASeqs(x):
    Seqs=[]
    FullSeqList=[]
    Break = int(x["chainbreak"])
    for entry in x["alignment"]:
        Seq=[]
        Basepositions=[]
        for Nucleotide in x["alignment"][entry]:
            Datalist=Nucleotide.split("|")
            Seq.append(Datalist[3])
            Basepositions.append(Datalist[4])#Split the sequence and the basepositions here and analyze each one
        SeqA = Seq[0:Break]
        SeqB = Seq[Break:]
        BasepositionsA = Basepositions[0:Break]
        BasepositionsB = Basepositions[Break:]
        ContinuityA = CheckBPContinuitySR_AA(BasepositionsA)
        if ContinuityA == True:
            ContinuityB = CheckBPContinuitySR(BasepositionsB)
            if ContinuityB == True:
                SeqA2=SeqA[1:-1]
                SeqB2=SeqB[1:-1]
                SeqList=Insert_N(SeqA2,2,"A")
                SeqB2="".join(SeqB2)
                for entry2 in SeqList:
                    FullSeq=entry2+"$"+SeqB2
                    FullSeqList.append(FullSeq)
                for entry3 in FullSeqList:
                    RevSeqList=list(entry3)
                    RevSeq=Reverse(RevSeqList)
                    RevSeq="".join(RevSeq)
                    Seqs.append(entry3)
                    Seqs.append(RevSeq)
    return Seqs

def CheckBPContinuitySR_AA(x):
    DiffList=[]
    Checklist=[1,1,2,1,1,1]
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
    KTurns,CLoops,SR_A = ImportEverything() #For now only contains motif IL_2949.8 (KTurns)
    SRs=ImportSarcinRicinMotifs()
    KTurnSeqs = MakeKTurnsSeqs(KTurns)
    CLoopSeqs = MakeCLoopsSeqs(CLoops)
    SRSeqs=MakeSRSeqs(SRs)
    SR_ASeqs=MakeSR_ASeqs(SR_A)
    with open("InternalMotifs.csv","w") as file:
        for entry in KTurnSeqs:
            file.write(entry+"+K\n")
        for entry in CLoopSeqs:
            file.write(entry+"+C\n")
        for entry in SRSeqs:
            file.write(entry+"+S\n")
        for entry in SR_ASeqs:
            file.write(entry+"+S\n")
