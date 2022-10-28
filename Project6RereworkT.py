from asyncio.subprocess import PIPE
from Bio import AlignIO
import re
import subprocess
import os
import multiprocessing
import re
import logging
from memory_profiler import memory_usage
import logging.handlers
from collections import Counter
import sys,traceback

Seqdict={}
tlist=[]
tdict={}
Startdict={}
Enddict={}
RMDict={}
CMDict={}
molist=["RM00003","RM00008","RM00010","RM00011","RM00018","RM00019","RM00024","RM00029"]


def listener_configurer():
    root=logging.getLogger()
    h=logging.handlers.RotatingFileHandler("/home/ubuntu/Master/Project6/Project6Rerework.log")
    f=logging.Formatter("%(message)s")
    h.setFormatter(f)
    root.addHandler(h)

def listener_process(queue,configurer):
    configurer()
    while True:
        try:
            record=queue.get()
            if record is None:
                break
            logger=logging.getLogger(record.name)
            logger.handle(record)
        except Exception:
            print("Problem:", file =sys.stderr)
            traceback.print_exc(file=sys.stderr)

def worker_configurer(queue):
    h = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.INFO)

def Errconv(x):
    time=re.search("\nuser\t[\S]+\n",x)
    time=time.group()
    time=time[6:-2]
    time=time.split("m")
    if float(time[0])!=0:
        Add=60*float(time[0])
        time[1]=float(time[1])+Add
    else:
        pass
    time=float(time[1])
    return time

def Subprocess(x,return_dict,queue,configurer,y):
    configurer(queue)
    if len(x) != 0:
        logger=logging.getLogger(__name__)
        level=logging.INFO
        cp=subprocess.Popen(x,text=True,stdout=PIPE,stderr=PIPE,shell=True)
        memusage=memory_usage(proc=cp.pid,max_usage=True,interval=0.1)
        out,err=cp.communicate()
        time=Errconv(err)
        sequencelength=len(x)-36
        return_dict[y]=out
        logmsg=(y,sequencelength,memusage,time)
        logger.log(level,logmsg)
    return return_dict

def Multi(x):
    manager=multiprocessing.Manager()
    return_dict=manager.dict()
    processlist=[]
    queue=multiprocessing.Queue()
    listener=multiprocessing.Process(target=listener_process,args=(queue,listener_configurer))
    listener.start()
    for element in x:
        process=multiprocessing.Process(target=Subprocess,args=(x[element],return_dict,queue,worker_configurer,element))
        process.start()
        processlist.append(process)
    for process in processlist:
        process.join()
    queue.put_nowait(None)
    listener.join()
    return return_dict

def A2S(x):
    x=str(x)
    y=str((re.sub("-*","",x)))
    return y

def Tlist2(x):
    namelist=[]
    rowlist=[]
    ConDict={}
    Fdict={}
    with open(x,"r") as file:
        c=file.readlines()
        for row in c:
            row = row.split()
            if row[0] in molist:
                rowlist.append(row)
                name=row[2]+"/"+row[3]+"-"+row[4]
                namelist.append(name)
    CList=Counter(namelist)
    for element in CList:
        if CList[element] == 1:
            c=namelist.index(element)
            ConDict[element]=c
        if CList[element] != 1:
            c=[x for x,e in enumerate(namelist) if e == element]
            ConDict[element]=c
    for entry in ConDict:
        if type(ConDict[entry]) == int:
            temp=rowlist[ConDict[entry]][0]
            templist=[]
            templist.append(temp)
            Fdict[entry]=templist
        if type(ConDict[entry]) == list:
            tlist=[]
            for entry2 in ConDict[entry]:
                tlist.append(rowlist[entry2][0])
            Fdict[entry]=tlist
    return Fdict

def Alignment(x):
    with open ("/home/ubuntu/Master/Rfam2.seed","r") as file:
        format ="stockholm"
        for alignment in AlignIO.parse(file,format):
            for record in alignment:
                if record.description in x.keys():
                    Seqdict[record.description]=A2S(record.seq)
#                if record.description not in x.keys():
#                    print(record.description)
    file.close
    return Seqdict

def MCMDict(x):
    for entry in x:
        command="time /home/ubuntu/Master/Project5/out "+x[entry]
        CMDict[entry]=command
    return CMDict

def LF(x):
    dellist=[]
    for entry in x:
        if len(x[entry]) >=100:
            dellist.append(entry)
    for entry in dellist:
        del x[entry]
    return x  

def MakeOut(x):
    for entry in x:
        d=x[entry]
        e=re.sub("Answer: \n","",d)
        z=re.split("\n",e)
        for res in z:
            res2=re.split(",",res)
            if len(res2) == 3:
                for idx,part in enumerate(res2):
                    if idx == 0:
                        u=part[4:-1]
                        if u=="":
                            u="X"
                    if idx == 1:
                        t=part[1:-3]
                    if idx == 2:
                        s=part[1:-2]
                n=Fdict[entry]
                n=",".join(n)
                with open("/home/ubuntu/Master/Project6/Out.txt","a") as file:
                    file.write(entry+"\t"+n+"\t"+u+"\t"+t+"\t"+s+"\n")
                file.close

if __name__=="__main__":
    Fdict=Tlist2("/home/ubuntu/Master/motif_matches.txt")
    Seqdict=Alignment(Fdict)
    Seqdict2=LF(Seqdict)
    CMDict=MCMDict(Seqdict2)
    return_dict=Multi(CMDict)
    MakeOut(return_dict)
    with open("/home/ubuntu/Master/Project6/Out.txt","a") as file:
        file.write("\n")
    file.close