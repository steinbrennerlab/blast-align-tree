from Bio import SeqIO
import sys
from Bio import AlignIO

input_file = sys.argv[1]
output_file = sys.argv[2]

for record in AlignIO.read(input_file, "fasta"):
    print(record.description)


def modified(records):
    for record in records:
        #Clear the description
        record.description=""
        yield record

records = modified(SeqIO.parse(input_file,"fasta"))
count = SeqIO.write(records, output_file , "fasta")
print("Converted %i records" % count)

with open(output_file, "w") as o:
    for record in AlignIO.read(input_file, "fasta"):
       print(record.name)
       record.description=""
       SeqIO.write(record, o, "fasta")