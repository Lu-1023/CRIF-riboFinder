import re
import os
import sys
import getopt
import numpy
def usage():
    print("python hsc.py -sp <species> \n")
    print(
        "<species> : the sample from which species\n")

if len(sys.argv[1:]) < 1:
    usage()
    sys.exit("WRONG,you must type the options firstly!")
try:
    opt, arg = getopt.getopt(sys.argv[1:], "hs:", ["help", "sp="])
except:
    usage()
    sys.exit("WRONG,non-existent options you type!")
for name, value in opt:
    if name in ("-h", "--help"):
        usage()
    if name in ("-s", "--sp"):
        organism = value
#function to create dict1:{id+hsc:[start1,stop1,start2,stop2```]} or {id+hsc:[start,stop]} and respectively dict2 contain motif information
def find_seq(x):
    hsc=open("/media/hp/disk1/song/Genomes/hiden_stop_codon/"+organism +"_ALL_hsc_position.txt","r")
    h={}
    m={}
    for l in hsc:
        if (re.match(r'^id',l)):
            continue
        id,motif,start,stop,stop_1,stop_2=l.strip("\n").split("\t")
        if (i == 1):
            if (id+"_"+stop_1 in h.keys()):
                h[id + "_" + stop_1].extend([start, stop])
                m[id + "_" + stop_1].append(motif)
            else:
                h.update({id + "_" + stop_1: [start, stop]})
                m.update({id + "_" + stop_1: [motif]})
        else:
            if (id+"_"+stop_2 in h.keys()):
                h[id + "_" + stop_2].extend([start, stop])
                m[id + "_" + stop_2].append(motif)
            else:
                h.update({id + "_" + stop_2: [start, stop]})
                m.update({id + "_" + stop_2: [motif]})
    hsc.close()
    return h,m

for i in range(1,3):
    H,M=find_seq(i)
    hsc_frame = open("/media/hp/disk1/song/Genomes/hiden_stop_codon/" + organism + "_ALL_hsc_position_frame" + str(i) + ".txt", "w")
    hsc_frame.write("id"+"\t"+"motif"+"\t"+"start"+"\t"+"stop"+"\t"+"stop_"+str(i)+"\n")
    for k in H.keys():
        id = re.match(r'(.*)_\d+', k).group(1)
        same = re.match(r'.*_(\d+)', k).group(1)
        #fs sites which use same hsc
        if len(H[k])>2 :
            H[k] = list(map(int, H[k]))
            dist=[]
            record = []
            for j in range(2,len(H[k]),2):
                if H[k][j]-H[k][j-1]<30:
                    dist.extend([H[k][j-2],H[k][j-1],H[k][j],H[k][j+1]])
                    record.extend([j/2-1,j/2])
                    record=list(map(int,record))
                    if j+2==len(H[k]):
                        hsc_frame.write("\t".join([id, ",".join(M[k][min(record):max(record)+1]), str(min(dist)), str(max(dist)), same]) + "\n")
                else:
                    if j+2==len(H[k]):
                        if len(dist)>2:
                            hsc_frame.write("\t".join([id,",".join(M[k][min(record):max(record)+1]),str(min(dist)),str(max(dist)),same])+"\n")
                        elif len(dist)==0:
                            hsc_frame.write("\t".join([id, M[k][int(j/2)-1], str(H[k][j-2]), str(H[k][j-1]),same]) + "\n")
                        hsc_frame.write("\t".join([id, M[k][int(j/2)],str(H[k][j]), str(H[k][j+1]), same]) + "\n")
                    else:
                        if len(dist)>2:
                            hsc_frame.write("\t".join([id,",".join(M[k][min(record):max(record)+1]),str(min(dist)),str(max(dist)),same])+"\n")
                            dist.clear()
                            record.clear()
                        elif len(dist)==0:
                            hsc_frame.write("\t".join([id, M[k][int(j/2)-1], str(H[k][j-2]), str(H[k][j-1]), same]) + "\n")
                        continue
        #fs sites has unique hsc
        else :
            hsc_frame.write("\t".join([id,M[k][0],H[k][0],H[k][1],same])+"\n")
    hsc_frame.close()
