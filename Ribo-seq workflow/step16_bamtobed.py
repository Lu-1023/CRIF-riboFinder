import getopt
import sys
import re
import os
# get the python script input
def usage():
    print("python bamtobed.py --doc <PATH> --periodicity <int>\n" )
    print("<PATH> :/media/hp/disk4/shuang/Ribo_seq/XXX\n<int> : the periodicity cutoff you set")
if len(sys.argv[1:])<1  :
    usage()
    sys.exit("WRONG,you must type the options firstly!")
try:
    opt,arg=getopt.getopt(sys.argv[1:],"hd:p:",["help","doc=","periodicity="])
except getopt.GetoptError :
    usage()
    sys.exit("WRONG,non-existent options you type!")

for name,value in opt:
    if name in ("-h","--help"):
        usage()
    if name in ("-d","--doc"):
        sample=value
    if name in ("-p", "--periodicity"):
        threshold = value

#function : read the .bam file an the name of read its length and CIGAR (mapping result) ,then write to single end .bed file
def pro_bam(x,y):
    #print to the screen
    print("proccesing", x," convertion start")
    # judge plus 15 , 14 or 16 to reads in this sample
    frame=open(y+"/b_align/codon_stat/"+x+".txt",'r')
    trim={}
    for f in frame:
        f=f.strip('\n')
        if re.match(r'^\d\d\t\d+\t',f):
            length,num,in0,plus,minus=f.split('\t')
            if in0 >= threshold :
                trim[int(length)]=15
            elif plus >= threshold:
                trim[int(length)]=14
            elif minus >= threshold :
                trim[int(length)]=16
            else:
                continue
        else:
            continue
    #show the length of reads we record
    print(trim)
    bam = os.popen("samtools view " + y + "/align2/" + x +".bam").readlines()
    for F in bam:
        F = F.strip('\n')
        id = F.split('\t')[0]
        flag=F.split('\t')[1]
        chr = F.split('\t')[2]
        bamstart = int(F.split('\t')[3]) - 1  # start : .bam-1=.bed
        quality = F.split('\t')[4]
        cigar = F.split('\t')[5]
        if flag != '4': # 4 : un-mapped mark
            # bedtools bamtobed -split only to split reads contained CIGAR "N" ,so we calculate the mapping length of reads we must add the "D" number
            length = re.findall(r'(\d+)[MD]', cigar)
            for i in range(0, len(length)):
                length[i] = int(length[i])
            length=sum(length)
            if length in trim.keys():
                bed = open(y + "/align2/" + x + ".bed", 'a+')
                if flag in "0 256": # + strand case
                    #record the CIGAR 'M','N','D' and judge have the splicing site or not
                    strand="+"
                    total=re.findall(r'(\d+)([MND])',cigar)
                    s = 0
                    splice=0
                    for i in range(0, len(total)):
                        s +=int(total[i][0])
                        if 'M' in total[i][1] or 'D' in total[i][1]:
                            splice +=int(total[i][0])
                        if splice > 15 :
                           pos = i
                           skip = trim[length] - (splice - int(total[pos][0]))
                           start = bamstart + (s - int(total[pos][0]))  + skip
                           end = start + 1
                           start = str(start)
                           end = str(end)
                           row = "\t".join([chr, start, end, id, quality, strand])
                           bed.writelines(row + "\n")
                           break
                        else :
                            continue

                elif flag in "16 272" : # - strand case
                    strand="-"
                    total = re.findall(r'(\d+)([MND])', cigar)
                    s = 0
                    splice = 0
                    bamend=bamstart
                    for a in range(0,len(total)):
                        bamend += int(total[a][0])
                    for i in range(len(total)-1,-1,-1):
                        s += int(total[i][0])
                        if 'M' in total[i][1] or 'D' in total[i][1]:
                            splice += int(total[i][0])
                        if splice > 15:
                            pos = i
                            skip = trim[length] - (splice - int(total[pos][0]))
                            end=bamend - (s - int(total[pos][0])) - skip
                            start = end - 1
                            start = str(start)
                            end = str(end)
                            row = "\t".join([chr, start, end, id, quality, strand])
                            bed.writelines(row + "\n")
                            break
                        else:
                            continue

                else : # exceptio+-n
                    exception=open(y+"/b_align/codon_stat/"+x+"_other.txt",'a+')
                    exception.writelines(F+"\n")
                    exception.close()

            else:
                continue

        else:
            continue
    frame.close()
    bed.close()



#get the unique sample names
file=open(sample+"/b_align/deadapter.txt",'r').read().splitlines()
for j in file:
    pro_bam(j,sample)
    print("===================finish "+j+"===================")

