import sys
import pysam 

# %prog bamFile outF 
if len(sys.argv) < 3:
  exit(1)

inbamFN = sys.argv[1]
pileupFN = sys.argv[2]
outbamFN = sys.argv[3]
cutoff = int(sys.argv[4])

inbamF = pysam.Samfile(inbamFN,'rb')
outbamF = pysam.AlignmentFile(outbamFN, "wb", template=inbamF)

# initialize list. 
pileupF = open(pileupFN, 'r')
read_dict = dict()
for line in pileupF:
  items = line.strip().split('\t')
  if int(items[3]) > cutoff:
    read_dict[(items[0],int(items[1]),items[2]=="True")] = int(items[3])
pileupF.close()
print("read_dict generated:",len(read_dict)," items")

curr_pileup,count = [],0
#the bam file should be coordinate sorted. 
for idx,read in enumerate(inbamF.fetch()):
#  if idx == 100000:
#    break
  if idx% 100000 == 0:
    print( "%s lines processed"%idx)
  read_info = (read.reference_name,read.reference_start,read.is_reverse)
  # if read in pileup, skip this read. 
  if read_info == curr_pileup and count!=1:
    count -= 1
#    print(read_info)
    #continue
  elif read_info in read_dict:
    curr_pileup = read_info
    count = read_dict.pop(curr_pileup)
  else:
    outbamF.write(read)
    
