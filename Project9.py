from collections import Counter

rmdict={"RM00003":"C","RM00008":"4","RM00010":"K","RM00011":"K","RM00018":"S","RM00019":"S","RM00024":"2","RM00029":"3"}
mrdict={"C":"RM00003","4":"RM00008","K":"RM00010","K":"RM00011","S":"RM00018","S":"RM00019","2":"RM00024","3":"RM00029"}
Üdict={"RM00003":"C-Tetraloop","RM00008":"GNRA-Tetraloop","RM00010":"Kink-Turn","RM00011":"Kink-Turn","RM00018":"Sarcin-Ricin Motif","RM00019":"Sarcin-Ricin Motif","RM00024":"T-Loop","RM00029":"UNCG-Tetraloop"}

def CreateNamelistAndRFAMdict(x):
    Namelist=[]
    RFAMdict={}
    Foundlist=[]
    for entry in x:
        entry2=entry.split("\t")
        name=entry2[0]
        rfam=entry2[1]
        found=entry2[2]
        Namelist.append(name)
        RFAMdict[name]=list(rfam.split(","))
        Foundlist.append(found)
    return Namelist,RFAMdict,Foundlist

def CreateFoundDict(x,y,z):
    SortedSet=[]
    InDict={}
    Founddict={}
    for entry in y:
        if entry not in SortedSet:
            SortedSet.append(entry)
    i=0
    j=-1
    for entry in SortedSet:
        j+=z[entry]
        InDict[entry]=[i,j]
        Templist=Foundlist[i:j+1]
        i=j+1
        Founddict[entry]=Templist
    return Founddict

def IndvCounter(x,y,z,):
    RFAMList=[]
    for entry in x:
        RFAMS=y[entry]
        FOUNDS=z[entry]
        for rfam in RFAMS:
            if rmdict[rfam] in FOUNDS:
                RFAMList.append(rfam)
    RFAMListCounted=Counter(RFAMList)
    return RFAMListCounted

def FoundDictSets(x):
    FoundDict2={}
    for entry in x:
        templist=[]
        trulist=[]
        NewFoundSet=set(x[entry])
        for Found in NewFoundSet:
            temp=list(Found)
            templist.append(temp)
        for entry2 in templist:
            for entry3 in entry2:
                if entry3 not in trulist:
                    trulist.append(entry3)
        FoundDict2[entry]=trulist
    return FoundDict2

def RFAMDictCounter(x):
    Counterlist=[]
    for entry in x:
        for rfam in x[entry]:
            Counterlist.append(rfam)
    Countedlist=Counter(Counterlist)
    return Countedlist

def MathAndOutput(x,y,z):
    for entry in x: 
        Gesamt=x[entry]
        try:
            Gezählt=y[entry]
        except:
            Gezählt=0
        Verhältnis=Gezählt/Gesamt*100
        Verhältnis=round(Verhältnis,2)
        Name=z[entry]
        Output="{Name} wurde insgesamt {Gesamt} gesucht und {Gezählt} mal gefunden. Genauigkeit für {Name} liegt bei {Verhältnis}%".format(Name=Name,Gesamt=Gesamt,Gezählt=Gezählt,Verhältnis=Verhältnis)
        print(Output)
    return 0

if __name__ == "__main__":
    with open("C:\\Master\\Project6\\Out,0-150Updated.txt","r") as file:
        file=file.readlines()#File needs to be truncated (no empty lines at the end)
        Namelist,RFAMdict,Foundlist=CreateNamelistAndRFAMdict(file)
    Counterlist=Counter(Namelist)
    FoundDict=CreateFoundDict(Foundlist,Namelist,Counterlist)
    FoundDict2=FoundDictSets(FoundDict)
    Nameset=set(Namelist)
    RFAMListCounted=IndvCounter(Nameset,RFAMdict,FoundDict2)
    Countedlist=RFAMDictCounter(RFAMdict)
    MathAndOutput(Countedlist,RFAMListCounted,Üdict)