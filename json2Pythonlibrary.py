import json
import sys
data = []
longdict={}
seqdict={}
iddict={}
numdict={}
with open(sys.argv[1]) as json_file: #Datei wird dem Programm als Commandline Argument übergeben.
    data = json.load(json_file) #Laden des JSON Files und erstellen von später verwendeten Dictionaries sowie Abspeicherung des JSON Files in "data" als Liste

 #Erstellen des Longdict welches als Key die ID der Sequenz und als Values die Informationseinträge dafür enthält
for i in range(len(data)):
    for element in data[i]["alignment"]:
        temp=[]
        for long in data[i]["alignment"][element]:
            temp.append (long)
        longdict[element]=temp
#Erstellen des Seqdict mit den Sequenz IDs als Keys und den dazugehörigen Sequenzen
for element in longdict:
    x=element
    y=longdict[element]
    temp2=[]
    for element in y:
        temp=[]
        a = element.split("|")
        for idx, val in enumerate(a):
            if idx == 3:
                temp2.append(val)
    temp2="".join(temp2)
    seqdict[x]=temp2
#Erstellen des Numdict mit den Sequenz IDs als Keys und den dazugehörigen Basenpaarpositionen in ihrer RNA
for element in longdict:
    x=element
    y=longdict[element]
    temp2=[]
    for element in y:
        temp=[]
        a = element.split("|")
        for idx, val in enumerate(a):
            if idx == 4:
                temp2.append(val)
    numdict[x]=temp2
#Erstellen des iddicts mit den Motif IDs als Keys und den Alignment Loop IDs als Values
for i in range(len(data)):
    alignment_keys=[]
    alignment_keys=list(data[i]["alignment"].keys())
    iddict[data[i]["motif_id"]]=alignment_keys
#Loopf ID Eingabe als try/except Schleife
try:
    Eingabe = sys.argv[2]
    for value in iddict[Eingabe]:
        print(value,",",seqdict[value],",",numdict[value])
except:
    print("Please enter a motif ID after your JSON file")

json_file.close