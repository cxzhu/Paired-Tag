import sys
import pysam 

# %prog bamFile outF 
if len(sys.argv) < 3:
  exit(1)

bamFN = sys.argv[1]
outFN = sys.argv[2]


bamF = pysam.Samfile(bamFN,'rb')
# initialize list. 
read_dict = dict()
for idx,read in enumerate(bamF.fetch(until_eof=True)):
  if idx% 100000 == 0:
    print( "%s lines processed"%idx)
#  print(read.reference_name,read.reference_start,read.is_reverse)
  try: 
    read_dict[(read.reference_name,read.reference_start,read.is_reverse)] += 1 
  except KeyError:
    read_dict[(read.reference_name,read.reference_start,read.is_reverse)] = 1 

# sort the items. 
sorted_list = sorted(read_dict.items(),key=lambda x: x[1],reverse=True)

outF = open(outFN,'w')
for item in sorted_list:
  outF.write(str(item[0][0])+"\t"+str(item[0][1]) + "\t" +
        str(item[0][2])+"\t"+str(item[1])+"\n")
outF.close()
