#Trying to pull fasta nucleotide entries based on a text list, with optional aa start and end points

#Takes "entry" the entry gene and subfolder and "list" a filename containing the cds of interest.  As of 9/18, refers to blast_align_tree output folder
#Example: python extract_seq.py gene_ID list.txt
#Example: python extract_seq.py AT5G20480.1 genelist.txt
## python extract_seq.py 10.csv C:\science\blast_align_tree\Zm00001d014134_T001\Zm00001d014134_T001.nt.merged.fa
## or python C:\science\blast_align_tree\extract_seq.py 10.csv Zm0001d014134_T001.nt.merged.fa 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import argparse
import io
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("entry", help="name of gene given to original blast scripts")
parser.add_argument("list", help="name of the genes list, a text file with one gene on each line.  As of 9/8/18, tree.R can generate this")
parser.add_argument("aa_start", type=int, nargs='?', default=0)
parser.add_argument("aa_end", type=int, nargs='?')
args = parser.parse_args()
cwd = os.getcwd()

print('\n')
print("Extract_seqs is using {} as entry".format(args.entry))
orthogroup_seqs = cwd + "/" + str(args.entry) + "/output/Orthogroup_Sequences"
seqlist = cwd + "/" + str(args.entry) + "/output/" + str(args.list)
output = cwd + "/" + str(args.entry) + "/output/" + str(args.list) + ".orthogroups.txt"
print("Looking for all sequences in Orthofinder output Orthogroup_Sequences: {}".format(orthogroup_seqs))

print(seqlist)

#quit()

directory = os.fsencode(orthogroup_seqs)

save_file = open(output, 'a')
save_file.write('taxa'+'\t'+'OG'+'\t'+'OG_type') #line break
save_file.write('\n')
save_file.close()

with open(seqlist) as f:
    for line2 in f:
        line = line2.rstrip('\n') #Strip line break
        for file in os.listdir(directory):
            filename = cwd + "/" + str(args.entry) + "/output/Orthogroup_Sequences/"+os.fsdecode(file)
            for seq_record in SeqIO.parse(filename, "fasta"):
                save_file = open(output, 'a+')
                if line == seq_record.id:
                        save_file = open(output, 'a')
                        save_file.write(seq_record.id)
                        save_file.write('\t') #tab break
                        temp = str((file))
                        save_file.write(str(os.fsdecode(file)))
                        save_file.write('\t') #tab break
                        save_file.write('nonunique-OG') #tab break
                        save_file.write('\n')
                        break
            else:
                continue
            break
    save_file.close()

