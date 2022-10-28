import re
import pandas as pd
from operator import attrgetter
#Anlegen des rmdicts, welches die RM Annotation von RFAM als Key und die von mir definierten Bennennungen für die Motife als Values enthält
rmdict={"RM00003":"C","RM00008":"4","RM00010":"K","RM00011":"K","RM00018":"S","RM00019":"S","RM00024":"2","RM00029":"3"}
mrdict={"C":"RM00003","4":"RM00008","K":"RM00010","K":"RM00011","S":"RM00018","S":"RM00019","2":"RM00024","3":"RM00029"}
slist=["C","4","K","S","2","3","X"]
namelist=[]
rfamlist=[]
foundlist=[]
mfelist=[]
prettylist=[]
namelist2=[]
rfamlist2=[]
foundlist2=[]
mfelist2=[]
prettylist2=[]
rowlist=[]
namelist3=[]
rfamlist3=[]
foundlist3=[]
mfelist3=[]
prettylist3=[]
tl=[]

def kcal(x):
    x=x/10
    return x

def splitrm(x):
    x=re.split(",",x)
    return x

def splitfound(x):
    x=re.split("",x)
    for element in x:
        if element not in slist:
            x.remove(element)
    return x

def makerow(x):
    x=row(x[1][0],x[1][1],x[1][2],x[1][3],x[1][4])
    return x

def makedict(a,b,c,d,e):
    Newdict={"SeqName":a,"RFAMMotif":b,"FoundMotif":c,"MFE":d,"PrettySeq":e}
    return Newdict

def RFAM2Found(x):
    y=[]
    for entry in x:
        entry=rmdict[entry]
        y.append(entry)
    return y

def top(x,y,z):
    List2=[]
    List3=[]
    for entry in x:              #Erster Teil ist dafür da von jedem einzelnen Namen erstmal das beste Ergebnis in Liste2 zu packen
        zz=z.loc[z["SeqName"]==entry]
        zzz=zz.nsmallest(1,"MFE",keep="first")
        zzzz=zzz.iloc[0]
        List2.append(zzzz)
    print(len(List2),"len List2")
    q=y["RFAMMotif"].values                     
    p=y["FoundMotif"].values
    for i in range(len(q)):    #Ab hier wird für die Multirfam Sequenzen gearbeitet: Zuerst werden alle Ergebnisse die alle zu den Rfam passenden FoundMotifs haben in Liste3 eingetragen
        n=[]
        for entry in q[i]:
            entry=rmdict[entry]
            n.append(entry)
        t=0
        for element in n:
            if element in p[i]:
                t+=1
        if t==(len(n)) and t!=1:
            List3.append(y.iloc[i])
    Tempdf=pd.DataFrame(List3)#Tempdf enthält für alle Multiseqs die Einträge die, alle Rfams auch als FoundMotifs haben
    for element in x:
        Name=Tempdf.loc[Tempdf["SeqName"]==element]
        try:
            Best=Name.nsmallest(1,"MFE",keep="first")
            Best2=Best.iloc[0]
            List2.append(Best2)#Aus Tempdf/List3 wird jetzt jeweils für jede Multiseq das beste Ergebnis rausgesucht -> Bestes Ergebnis das für alle Motife zutrifft
        except:pass
    for element in x:
        t=y.loc[y["SeqName"]==element] #t enthält für jeden Namen alle Einträge aus data, die mindestens eins der Rfam Motife enthält
        c=t.iloc[0]["RFAMMotif"] #c enthält von jedem Eintrag in t nur die Liste der RM Motife
        if len(c) == 1:
            u=t.nsmallest(1,"MFE",keep="first") #Falls für einen Namen nur ein einziges RM Motif gesucht wird, wird mit u das Object mit der besten MFE in t ausgesucht. 
            v=u.iloc[0]
            List2.append(v) #Das passende Object wird in List2 gepackt --> Für die Single Sequenzen sind jetzt sowohl das Beste insgesamt als auch das beste für das gesuchte RMMotif drin
        if len(c) != 1: #Ab hier werden für die Multisequenzen jeweils die besten Einträge für die einzelnen RM Motife rausgesucht
            for element in c:
                List=[] #Leere Liste um darin die Einträge zu speichern welche für die jeweilige Multiseq die richtigen einzelnen Motife tragen
                celement=rmdict[element] #Celement ist die SingleLetter Version von jeweils einem RMMotifs welches in Element steckt
                t.loc[:,"RFAMMotif"] = element #Es wird ein Datenframe erstellt bei dem das RFAMMotif für die komplette Zeile gesetzt wird, print(t) zeigt es gut
                q=t["FoundMotif"] #Nur die Column aus t die die ganzen FoundMotifs enthält um im Anschluss zu vergleichen ob das jeweilige Celement enthalten ist
                i=0 #Anfangspunkt für Indices die durchgegangen werden mit t.iloc[i]
                for element in q: #Elemente sind die einzelnen FoundMotif Listen aus t, also alle FoundMotif Listen aus Data für den Namen die mindestens eins der Motife tragen
                    if celement in element:
                        List.append(t.iloc[i]) #Alle in denen eins der beiden drin ist wird abgespeichert in List, da t sortiert ist nach MFE kann im Anschluss List[0] als bestes für das 
                                                #jeweilige RMMotif verwendet
                    i+=1
                try:
                    List2.append(List[0])
                except:
                    pass
    del List3,Best,Best2,n,q,p,t,u,v
    Topdf=pd.DataFrame(List2)#List2 enhält jetzt die besten Teile, für Multis das Beste mit allen und das beste für die einzelnen sowie für Single das Beste welches auch die richtigen Founds hat
    return Topdf

def Str2List(x):
    if type(x) == str:
        x = x.split()
    return x

def SortAndResetIndex(x):
    x=x.sort_values(["SeqName","MFE"])
    x=x.reset_index()
    x=x.drop("index",axis=1)
    return x

class row:

    def __init__(self,name,rfam,found,mfe,pretty):
        self.name=name
        self.rfam=rfam
        self.found=found
        self.mfe=mfe
        self.pretty=pretty

    def get_name(self):
        return self.name

    def get_rfam(self):
        return self.rfam
    
    def get_found(self):
        return self.found

    def get_mfe(self):
        return self.mfe

    def get_pretty(self):
        return self.pretty

    def makelists(self):
        for element in self.rfam:
            if rmdict[element] in self.found:
                x=",".join(self.rfam)
                y=",".join(self.found)
                namelist.append(self.name)
                rfamlist.append(x)
                foundlist.append(y)
                mfelist.append(self.mfe)
                prettylist.append(self.pretty)
        return namelist,rfamlist,foundlist,mfelist,prettylist

    def createrowlist(self):
        if self.pretty in prettylist:
            rowlist.append(self)
        return rowlist
    
    def sortnm(self):
        x=sorted(self, key=attrgetter("name","mfe"))
        return x

if __name__=="__main__":
    pd.options.mode.chained_assignment = None
    data= pd.read_csv("C:\\Master\\Project6\\Out,0-150Updated.txt",sep="\t",header=None)
    data.columns=["SeqName","RFAMMotif","FoundMotif","MFE","PrettySeq"]
    data["MFE"]=data["MFE"].apply(kcal)
    data["RFAMMotif"]=data["RFAMMotif"].apply(splitrm)
    data["FoundMotif"]=data["FoundMotif"].apply(splitfound)
    for rows in data.iterrows():#Makelists filtert mit for element in self.rfam+if rmdict[element] in self.found alle Sequenzen raus, in dem bei keiner Struktur die gesuchte Found gefunden wurde
        rows=makerow(rows)
        rows.makelists()
        rows.createrowlist()
    sortedrowlist= row.sortnm(rowlist)
    for element in sortedrowlist:
        namelist2.append(row.get_name(element))
        rfamlist2.append(row.get_rfam(element))
        foundlist2.append(row.get_found(element))
        mfelist2.append(row.get_mfe(element))
        prettylist2.append(row.get_pretty(element))
    nameset=set(namelist2)
    rowdict2=makedict(namelist2,rfamlist2,foundlist2,mfelist2,prettylist2)#rowdict enthält SORTIERTE Listen!
    data3=pd.DataFrame(rowdict2)
    Topdf=top(nameset,data3,data)
    Topdf["RFAMMotif"]=Topdf["RFAMMotif"].apply(Str2List)
    Topdf=SortAndResetIndex(Topdf)
    with open("C:\Master\Project6\Rework\TopOutTest.csv","w") as file:
        Topdf.to_csv(file)
    file.close