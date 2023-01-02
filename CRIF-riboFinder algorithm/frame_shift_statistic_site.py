"""
This script is used to define the frame-shift sites in repeat motif.I use the independent t.text
to test whether the baseline(-1 or +1 frame) is different from the repeat box (from the motif_stop
site to the hidden stop codon). And then compared with the average of the baseline and repeat box ,
if baseline > repeat box , we see the motif site as the frame-shifting site.
intensity = (box_avg/all0_avg)-(allfs_avg/all0_avg)
"""
from scipy import stats
import re
import getopt
import sys
import re
import os
import numpy
import argparse

numpy.seterr(divide='ignore',invalid='ignore')
# get the python script input method1
# def usage():
#     print("python frame_shift_statistic_site.py --sp <species> --doc <PATH> -f <frame>\n")
#     print("<species> : the sample from which species\n<doc> : /media/hp/disk4/shuang/Ribo_seq/XXX\n<frame> : 1 or 2 represents +1 or -1 frame")
# if len(sys.argv[1:])<1 :
#     usage()
#     sys.exit("WRONG,you must type the options firstly!")
# try:
#     opt,arg=getopt.getopt(sys.argv[1:],"hs:d:f:",["help","sp=","doc=","frame="])
# except:
#     usage()
#     sys.exit("WRONG,non-existent options you type!")
# for name,value in opt:
#     if name in ("-h","--help"):
#         usage()
#     if name in ("-s","--sp"):
#         organism=value
#     if name in ("-d","--doc"):
#         sample=value
#     if name in ("-f","--frame"):
#         frame=int(value)

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
cds=open("/media/hp/disk1/song/Genomes/"+organism+"/"+organism+"_ref/CDS_DNA.fa","r")
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

#function:get the repeat box
def R_box(x,y,z): #x : motif start site .y: +1 or -1 hiden stop codon site .z: frame.
    r_box=[]
    global num
    num=0
    for j in range(x+z,y,3): #?
        if (id+"_"+str(j) in rpf.keys()):
           r_box.append(round(rpf[id+"_"+str(j)],1))
           num+=1
        else:
            r_box.append(0)
    r_box=numpy.array(r_box)
    #set the dot : numpy.set_printoptions(precision=4) or line 74
    global r_Len
    r_Len=(y-x-z)/3 # codon length
    #take the codon counts in 0 frame at repeat box, and the the medium to normalize repeat box
    inframe_box=[]
    global rbox_in_ecodon
    rbox_in_ecodon=0
    for j in range(x,y,3):
        if (id+"_"+str(j) in rpf.keys()):
           inframe_box.append(round(rpf[id+"_"+str(j)],1))
           rbox_in_ecodon+=1
        else:
            inframe_box.append(0)
    #we want to get the framshiting expression
    ratio=numpy.sum(r_box)/r_Len
    inframe_box_m=numpy.mean(inframe_box)
    r_box=(r_box*num/r_Len)/inframe_box_m
    r_box=numpy.round(r_box,6)
    # remove 0 in r_box to aviod affecting the p-value
    list0=[]
    for i in range(0,len(r_box)):
        if (r_box[i]==0):
            list0.append(i)
    r_box=numpy.delete(r_box,list0)
    rinframe=numpy.sum(inframe_box)/r_Len
    return r_box,ratio,rinframe

try:
    os.mkdir(sample+"/b_align/noscreen/FS")
except:
    print("file exists!")

#function : get the fs sites using the same hiden stop codon
def find_seq(x):
    hsc=open("/media/hp/disk1/song/Genomes/hiden_stop_codon/"+organism +"_ALL_hsc_position_frame"+str(x)+".txt","r")
    h={}
    for l in hsc:
        if (re.match(r'^id',l)):
            continue
        ID,motif,start,stop,same=l.strip("\n").split("\t")
        if (ID+"_"+same in h.keys()):
                h[ID + "_" + same].extend([start, stop])
        else:
            h.update({ID + "_" + same: [start, stop]})
    hsc.close()
    return h

#function : get the down/up reads count as the 2th compare
def down_up(s,e,f):
    up = []
    down = []
    #jump the whole repeat
    for a in range(s - 33+f, s, 3):
        if (id + "_" + str(a) in rpf.keys()):
            up.append(round(rpf[id + "_" + str(a)], 1))
        else:
            up.append(0)
    # jump the whole repeat
    for a in range(e+f, e + 30+f, 3):
        if (id + "_" + str(a) in rpf.keys()):
            down.append(round(rpf[id + "_" + str(a)], 1))
        else:
            down.append(0)
    return up,down


#function : to get the baseline and repeat box
def pro_bg (x,y):
    print("===================now processing "+x+" !===================")
    #get the _A.bg file and make it to a diretoryrpf={}
    global rpf
    rpf={}
    gene_exist=set() #storage the gene id in the sample
    bg_file=open(y+"/b_align/noscreen/"+x+"_A.bedgraph",'r')
    for L in bg_file:
        a=L.strip('\n').split('\t')
        for i in range(int(a[1]),int(a[2])):
           if int(a[3])>0 :
            rpf.update({a[0]+"_"+str(i):float(a[3])})
            gene_exist.add(a[0])
    bg_file.close()
    H=find_seq(frame)
    fs=open(y+"/b_align/noscreen/FS/"+x+"/"+x+"frame"+str(frame)+".txt",'w')
    fs.write("id"+"\t"+"motif"+"\t"+"start"+"\t"+"stop"+"\t"+"p-value"+"\t"+"foldchange"+"\n")
    box=open(y+"/b_align/noscreen/FS/"+x+"/"+x+"frame"+str(frame)+"_repeatbox.txt",'w')
    inf=open(y+"/b_align/noscreen/FS/"+x+"/"+x+"frame"+str(frame)+"_information.txt",'w')
    inf.write("id"+"\t"+"motif"+"\t"+"start"+"\t"+"stop"+"\t"+"hsc"+"\t"+"e_codon_inframe"+"\t"+"e_codon_baseline"+"\t"+"e_codon_rbox"+"\t"+"e_codon_rbox_inframe"+"\t"+"rbox_length(codon)"+"\t"+"cds_length(codon)"+"\t"+"rbox_inframe_avg"+"\t"+"p-value"+"\t"+"foldchange"+"\t"+"intensity(%)"+"\n")
    second=open(y+"/b_align/noscreen/FS/"+x+"/"+x+"frame"+str(frame)+"_compare.txt",'w')
    second.write("id"+"\t"+"motif"+"\t"+"start"+"\t"+"stop"+"\t"+"hsc"+"\t"+"p-value"+"\t"+"foldchange"+"\t"+"down/up"+"\n")
    # get the information for the hiden stop codon
    hsc = open("/media/hp/disk1/song/Genomes/hiden_stop_codon/" + organism + "_ALL_hsc_position_frame"+str(frame)+".txt", 'r')
    for l in hsc:
        l = l.strip('\n')
        if (re.match(r'^id', l)):
            continue
        global id
        id, motif, start, stop, same= l.split('\t')
        if (id in gene_exist):
            baseline=[]
            global e_count
            e_count=0
            for i in range(100+30+frame,size[id]-100-30,3):
                if (id+"_"+str(i) in rpf.keys()) :
                    baseline.append(round(rpf[id+"_"+str(i)],1))
                    e_count+=1
                else:
                    baseline.append(0)
            #value avg(all -1)
            length = (size[id] - 200) / 3  # codon length
            baseline = numpy.array(baseline)
            value_allfs=numpy.sum(baseline)/length

            # take the codon counts in 0 frame , and the the medium to normalize baseline
            inframe=[]
            in_codon = 0
            for j in range(100+30,size[id]-100-30,3):
                if (id+"_"+str(j) in rpf.keys()) :
                    inframe.append(round(rpf[id+"_"+str(j)],1))
                    in_codon+=1
                else:
                    inframe.append(0)
            inframe_m=numpy.mean(inframe)
            # value avg(all0)
            value_all0 = numpy.sum(inframe) / length

            #array baseline
            baseline=(baseline*e_count/length)/inframe_m
            baseline = list(numpy.round(baseline, 6))

            # remove 0 in baseline to aviod affecting the p-value
            list0 = []
            for i in range(0,len(baseline)):
                if (baseline[i] == 0):
                    list0.append(i)
            baseline = numpy.delete(baseline, list0)
            # if the -1 signal all is 0 ,skip it
            if (sum(baseline)==0):
                continue
            repeat_box,value_b,value_bIn= R_box(int(start), int(same), frame)
            repeat_box=list(repeat_box)
            intensity=(value_b*100/value_all0)-(value_allfs*100/value_all0)

            #length judgenment : the length of repeat box and baseline must be more 1
            if (len(baseline)<3 or  len(repeat_box)<3):
                continue
            if (stats.levene(baseline,repeat_box)[1]>0.05):
                p=stats.ttest_ind(baseline,repeat_box)[1]
            else:
                p = stats.ttest_ind(baseline, repeat_box,equal_var=False)[1]
            foldchange=round(numpy.mean(repeat_box)/numpy.mean(baseline),3)
            for i in range(0,len(repeat_box)):
                repeat_box[i]=str(repeat_box[i])
            inf.writelines("\t".join([id,motif,start,stop,same,str(in_codon),str(e_count),str(num),str(rbox_in_ecodon),str(r_Len),str(length),str(value_bIn),str(p),str(foldchange),str(intensity)])+"\n")
            box.writelines("\t".join([id,motif,start,stop]+repeat_box)+"\n")
            fs.writelines("\t".join([id,motif,start,stop,str(p),str(foldchange)])+"\n")
            if id+"_"+same in H.keys() and len(H[id+"_"+same])>2:
                up_inframe,down_inframe=down_up(int(start),int(stop),0)
                up_outframe, down_outframe = down_up(int(start),int(stop), frame)
                #normalize up and down respectively
                up_outframe=numpy.round(up_outframe/numpy.mean(up_inframe),3)
                down_outframe = numpy.round(down_outframe / numpy.mean(down_inframe),3)
                up_outframe=numpy.nan_to_num(up_outframe)
                down_outframe=numpy.nan_to_num(down_outframe)
                foldchange_2th=round((numpy.sum(down_outframe)+1)/(numpy.sum(up_outframe)+1),3)
                second.write("\t".join([id,motif,start,stop,same,str(p),str(foldchange),str(foldchange_2th)])+"\n")
    box.close()
    fs.close()
    hsc.close()
    inf.close()
    second.close()
#run
file=open(sample+"/b_align/deadapter.txt",'r').read().splitlines()
for j in file:
    try:
        os.mkdir(sample+"/b_align/noscreen/FS/"+j)
    except:
        print("file exists!")
    pro_bg(j,sample)
    print("===================finish "+j+" !===================")

