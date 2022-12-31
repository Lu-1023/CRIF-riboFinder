from scipy import stats
import re
import getopt
import sys
import re
import os
import numpy
import argparse
import random

# get the python script input method2
parser=argparse.ArgumentParser(description="the prediction method for CRIF.")
parser.add_argument("-f","--frame",type=int,metavar="frame",required=True,help="please input 1 or 2 means +1 or -1 frame")
parser.add_argument("-s","--sp",type=str,metavar="organism",required=True,help="the sample is from which species.")
parser.add_argument("-d","--doc",type=str,metavar="sample",required=True,help="the location where your files are.pattern like this :/media/hp/disk4/shuang/Ribo_seq/XXX")
args=parser.parse_args()
frame=args.frame
sample=args.doc
organism=args.sp

#get the CDS_DNA ,and input it to a directory
gene={}
cds=open("/media/hp/disk1/lu/demo/CDS_DNA.fa","r")
for line in cds:
    line=re.sub(r'\s+$','',line)
    if (re.match(r'^>',line)):
        gene_id=re.sub(r'^>','',line)
        continue
    gene.update({gene_id:line})
cds.close()

#get the size of the cds_gene
global size
size={}
for key,value in gene.items() :
    size.update({key:len(value)})

#read the position of hidden stop codon
p1={}
p2={}
# hsc2=open("/media/hp/disk1/song/Genomes/hiden_stop_codon/"+organism+"_ALL_hsc_position_frame2.txt","r")
# for l in hsc1:
#     if "id" in l:
#         continue
#     id,motif,start,stop,same=l.strip("\n").split("\t")
#     p1.update({id+"_"+start:int(same)})
# hsc2.close()
hsc=open(organism+"_ALL_hsc_position.txt","r")
for l in hsc:
    if "id" in l:
        continue
    id,motif,start,stop,same1,same2=l.strip("\n").split("\t")
    p1.update({id+"_"+start:int(same1)})
    p2.update({id+"_"+start:int(same2)})
hsc.close()


def bin(S,E,F,ID):
    read_count=[]
    #S:start E:end F:frame
    bin_num_box=10
    bin_num=10
    bin_size = int(round((E - F - S) // (3 * bin_num_box), 0))
    bin_size_set=3
    if S-100 < bin_size_set*3*bin_num:
        return
    if size[ID]-S-100 < bin_size_set*3*bin_num:
        return
    if E-S < 30:
        return
    #upstream
    for i in range(int(S-bin_size_set*bin_num*3+F),S,int(bin_size_set*3)):
        winup=[]
        for j in range(i,int(i+bin_size_set*3),3):
            winup.append(rpf.get(ID+"_"+str(j),0))
        read_count.append(sum(winup)/bin_size_set)
    #CRAF box
    for i in range(S+F,int(S+bin_size*(bin_num_box-1)*3),int(bin_size*3)):
        box=[]
        for j in range(i,int(i+bin_size*3)+3,3):
            box.append(rpf.get(ID+"_"+str(j),0))
        read_count.append(sum(box)/bin_size)
    dot=0
    box=[]
    for i in range(int(S+bin_size*(bin_num_box-1)*3+F),E,3):
        dot+=1
        box.append(rpf.get(ID+"_"+str(i),0))
    read_count.append(sum(box)/dot)
    #downstream
    for i in range(E,int(E+bin_size_set*bin_num*3),int(bin_size_set*3)):
        windown=[]
        for j in range(i,i+bin_size_set*3,3):
            windown.append(rpf.get(ID+"_"+str(j),0))
        read_count.append(sum(windown)/bin_size_set)
    read_count=list(map(str,read_count))
    return read_count

def judge(x,y,z):
    # x: your number y:frame z:gene id
    for i in range(x+y,size[z],3):
        if gene[z][i:i+3]=="TAA" or gene[z][i:i+3]=="TGA" or gene[z][i:i+3]=="TAG":
            return i+3
    return 0

#the
def bin_fs(x,y):
    #x:sample y:position
    bg=open(y+"/b_align/noscreen/"+x+"_A.bg",'r')
    global rpf
    rpf={}
    ts=[]
    for l in bg:
        id,s,e,v=l.strip("\n").split("\t")
        ts.append(id)
        for i in range(int(s),int(e)):
            rpf.update({id+"_"+str(i):int(v)})
    bg.close()
    ts=list(set(ts))
    crif=open(y+"/b_align/noscreen/FS/"+x+"/"+x+"frame"+str(frame)+"fs_fc1.txt",'r')
    plot1 = open(y + "/b_align/noscreen/FS_plot/" + x + "/" + x + "frame" + str(frame) + "_1.txt", 'w')
    plot2 = open(y + "/b_align/noscreen/FS_plot/" + x + "/" + x + "frame" + str(frame) + "_2.txt", 'w')
    plot3 = open(y + "/b_align/noscreen/FS_plot/" + x + "/" + x + "frame" + str(frame) + "_nc.txt", 'w')
    for l in crif:
        id,motif,start,stop,same=l.strip("\n").split("\t")[:5]
        if id == "id":
            continue
        same1 = p1[id + "_" + start]
        same2 = p2[id + "_" + start]
        box1 = bin(int(start), same1, 1, id)
        box2 = bin(int(start), same2, 2, id)
        if box2:
            plot2.writelines("\t".join([id,motif,start,stop,str(same2)]+box2)+"\n")
        if box1:
            plot1.writelines("\t".join([id,motif,start,stop,str(same1)]+box1)+"\n")
    for i in range(0,5000):
        # negative control
        nc = random.sample(ts, 1)[0]  # it means we extract one dot from ts list each time
        if nc in size.keys():
            a = random.sample(range(100, size[nc] - 100), 1)[0]
        #b = random.sample(range(a, size[nc] - 100), 1)[0]
            if (a - 100) % 3 != 0:
                a = a + (3 - (a - 100) % 3)
            b = judge(a, frame, nc)
            box3 = bin(a, b, frame, nc)
            if box3:
                plot3.writelines("\t".join([nc, "nc", str(a), str(frame), str(b)] + box3) + "\n")
        else:
            continue
    plot1.close()
    plot2.close()
    plot3.close()
    crif.close()

file=open(sample+"/b_align/deadapter.txt",'r').read().splitlines()
for j in file:
    try:
        os.mkdir(sample+"/b_align/noscreen/FS_plot/"+j)
    except:
        print("file exists!")
    bin_fs(j,sample)
    print("===================finish "+j+" !===================")
