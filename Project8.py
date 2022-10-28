import pandas as pd
import re
from datetime import date
import os
slist=["C","4","K","S","2","3","X"]
rmdict={"RM00003":"C","RM00008":"4","RM00010":"K","RM00011":"K","RM00018":"S","RM00019":"S","RM00024":"2","RM00029":"3"}
mrdict={"C":"RM00003","4":"RM00008","K":"RM00010","K":"RM00011","S":"RM00018","S":"RM00019","2":"RM00024","3":"RM00029"}
SingleRFAMList=[]
MultiRFAMList=[]
i=0
j=0
LowEndFilter=-30 #Gewünschter Filterwert x 10 muss hier eingetragen werden, Änderung wurde vorgenommen um für Berechnungen Integer benutzen zu können statt Floats
Steps=1
today=date.today()

def ReadAndDrop(x):
    data=pd.read_csv(x,sep=",")
    data.drop(columns="IDX",inplace=True)
    datalink=x
    return data,datalink

def SeqSetAndSort(x):
    SeqSet=set(x)
    SingleDataList=[]
    for element in SeqSet:
        SingleData=data.loc[data["SeqName"]==element]
        SingleDataList.append(SingleData)
    del SeqSet,SingleData
    return SingleDataList
        
def RFAMBack2List(x):
    if len(x) > 5:
        x=x.split(",")
    else:
        x=list(x)
    return x

def RFAMBack2List2(x):
    if len(x) == 11:
        x=x[2:-2]
        x=x.split()
    else:
        x=x[1:-1]
        x=re.sub("'","",x)
        x=re.sub(" ","",x)
        x=x.split(",")
    return x

def FoundBack2List(x):
    if len(x) > 1:
        x=x.split(",")
    else:
        x=list(x)
    return x

def FoundBack2List2(x):
    if len(x) == 5:
        x=x[2:-2]
        x=x.split()
    else:
        x=x[1:-1]
        x=re.sub("'","",x)
        x=re.sub(" ","",x)
        x=x.split(",")
    return x

def MFESort(x):
    for element in x:
        element=element.sort_values("MFE")
    return x

def Lensorting(x):
    Biglist=[SingleRFAMList,MultiRFAMList]
    for element in x:
        bruh=element.iloc[:1]
        bruh2=bruh["RFAMMotif"]
        for object in bruh2:
            if len(object)==1:
                SingleRFAMList.append(element)
            else:
                MultiRFAMList.append(element)
    return Biglist

def SingleRFAMListCounting(x,y,u,w):
    for element in x:
        bruhlist=[]
        if len(element.index)==2:
            bruh1=element.iloc[:1]
            bruh2=element.iloc[1:2]
            for obj in bruh1["MFE"]:
                bruhlist.append(obj)
            for obj2 in bruh2["MFE"]:
                bruhlist.append(obj2)
            p=float(bruhlist[0])-float(bruhlist[1])
            if p in range(u,1,w):
                y+=1
        else:
            pass
    return y

def MultiRFAMListCouting(x,y,u,w):
    for element in x:
        SList=[]
        bruh=element.iloc[:1]
        bruh2=bruh["RFAMMotif"]
        element=element.iloc[1:]
        for rm in bruh2:
            for rm2 in rm:
                SList.append(rmdict[rm2])
        for idx,dub in enumerate(element["RFAMMotif"]):
            if len(dub)==len(SList): #Wenn eins der RFAMMotifen die gleiche Länge hat wie die beste Sequenz, muss sie  alle Motife enthalten da nur die Strukturen hier landen die genau die Motife enthalten
                a=float(bruh["MFE"])
                element2=element.iloc[idx]
                b=float(element2["MFE"])
                p=float(a-b)
                if p in range(u,1,w):
                    y+=1
                else:pass
            else:pass
    return y

def Output(a,b,c,e,f,g,h):
    with open("AnalysisFile.txt","a") as file:
        file.write("Analysis Done on: {day}.{month}.{year}\n".format(year=today.year,month=today.month,day=today.day))
        file.write("Input used: {}\n".format(os.path.basename(h)))
        file.write("MFE Filter: {} kcal/mol\n".format(e/10))
        file.write("Gesamte Datenpunkte: {}\n".format(len(a)))
        file.write("Datenpuntke mit einem RFAMMotif: {}\n".format(len(b)))
        file.write("Datenpunkte mit mehreren RFAMMotifen: {}\n".format(len(c)))
        file.write("Singles Found: {}\n".format(f))
        file.write("Multis Found: {}\n".format(g))
        k=i+j
        file.write("Added Founds: {}\n".format(k))
        file.write("________________________\n")
        
if __name__ == "__main__":
    data,datalink=ReadAndDrop("C:\\Master\\moRNAO\\TopOut,0-150Updated.csv")
    data["RFAMMotif"]=data["RFAMMotif"].apply(RFAMBack2List2)
    data["FoundMotif"]=data["FoundMotif"].apply(FoundBack2List2)
    SingleDataList=SeqSetAndSort(data["SeqName"])
    MFESort(SingleDataList)
    Biglist=Lensorting(SingleDataList)
    i=SingleRFAMListCounting(SingleRFAMList,i,LowEndFilter,Steps)
    j=MultiRFAMListCouting(MultiRFAMList,j,LowEndFilter,Steps)
    Output(SingleDataList,SingleRFAMList,MultiRFAMList,LowEndFilter,i,j,datalink)