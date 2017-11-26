# Import and Create Files
## python make_fasta.py hapsFile sampleFile
import sys
haps_file=sys.argv[1]
# sample_file=haps_file.strip(".haps")+".sample"
## Nina's adjustment for filename differences between haps and sample file 
sample_file=sys.argv[2]
# changed .haps to .gen_haps 
fasta_file=haps_file.strip(".gen_haps")+".fasta"

# Make Dictionary from Haps file with RS Numbers as Keys
rs_dict={}
snps=[]
haps=open(haps_file)
haps=haps.readlines()
for i in range(len(haps)):
	col=haps[i].strip().split()
	if col[3] in ["A","C","G","T"] and col[4] in ["A","C","G","T"]:
		for j in range(5,len(col)):
			if col[j]=="0":
				col[j]=col[3]
			elif col[j]=="1":
				col[j]=col[4]
		rs_dict[col[1]]=col[5:]
		snps.append(col[1])

# Make Dictionary with Sample ID as Keys
id_dict={}
sample=open(sample_file)
sample=sample.readlines()

ids=[]
for i in range(2,len(sample)):
	col=sample[i].strip().split()
	ids.append(col[1]+"_A")
	ids.append(col[1]+"_B")

for id in ids:
	id_dict[id]=[]

for snp in snps:
	for id in ids:
		id_dict[id].append(rs_dict[snp][ids.index(id)])

# Output FASTA file
fasta=open(fasta_file,"w")

for id in ids:
	print >> fasta, ">"+id
	print >> fasta, "".join(id_dict[id])

fasta.close()