# Adam Steinbrenner
# astein10@uw.edu
# UW Seattle, Dept of Biology

#

# See readme for usage. You can run this script from the terminal by adding the first two standard arguments entry and entrydb, lists of arguments -n -h, -dbs, with optional arguments -dbs -add -add_db -aa
### pre-generated example sh blast-align-tree.sh AT4G33430.1 TAIR10cds.fa -n 10 10 10 -dbs TAIR10cds.fa Pvul218cds.fa Vung469cds.fa -hdr gene: polypeptide= locus=
#### example sh blast-align-tree.sh AT1G19220.1 TAIR10cds.fa -n 30 30 -dbs TAIR10cds.fa Pvul218cds.fa -hdr gene: polypeptide=

# Argparse turns arguments into variables
# argparse-bash https://github.com/nhoffman/argparse-bash
source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('entry')
parser.add_argument('entrydb')
parser.add_argument('-n', '--seqs', nargs='+',
                    help='input for blast -max_target_seqs; multiple values allowed')
parser.add_argument('-hdr', '--header', nargs='+',
                    help='input for which parts of the fasta header to pull; do these in the order of dbs. For example, you want 10 for Sbic3cds.fa to pull SOBIC IDs')
parser.add_argument('-dbs', '--database', nargs='+',
                    help='input for which blastDBs to search; (e.g. Vung469.cds.fa TAIR10cds.fa); multiple values allowed')
parser.add_argument('-add', '--add_seqs', nargs='*',
					help='input additional sequences to get from the add_database (see below) to add to the final file')
parser.add_argument('-add_db', '--add_dbs', nargs='*',
					help='input additional databases to get the addseqs from (see above)')
parser.add_argument('-aa', '--slice', nargs='*',
					help='input amino acids to cutoff')
EOF

#Echo input arguments
echo Sequence to find, entry: "$ENTRY"
echo Database to search, entrydb: "$ENTRYDB"
echo "-dbs, Databases to blast against:"
for a in "${DATABASE[@]}"; do
    echo "$a"
done
echo
echo Number of subject seqs to pull from tblastn:
for a in "${SEQS[@]}"; do
    echo "$a"
done
echo dbs length is "${#DATABASE[@]}"

#make a subdirectory with the gene name $ENTRY for storing output files
mkdir $PWD/$ENTRY
mkdir $PWD/$ENTRY/output


#Finds the sequence with header ENTRY in the fasta file ENTRYDB and extracts a translation file "ENTRY.seq.fa" to the correct subfolder
if ((${#SLICE[@]}>0)); then
	python $PWD/scripts/extract_translate.py $ENTRY $ENTRYDB ${SLICE[0]} ${SLICE[1]}
else
	python $PWD/scripts/extract_translate.py $ENTRY $ENTRYDB
fi
	
#The main BLASTP function identifies blast hits, up to a maximum -max_target_seqs defined in the arguments.  We call BLASTP locally to run on blast databases built with protein fasta files, which are located in the subfolder "genomes". We then generate blast outputs with various -outfmt options to ONLY include the hits names (sseqid) in text format. This list of genes is used in future steps
#Please note that -max_hsps is set to one; this means only one BLAST hit per target sequence will be recovered
#The for loop will repeat this for each blast database specified in the command line arguments "-db" 
for ((i=0; i<"${#DATABASE[@]}"; i++)); do
	echo
	echo blastp on: "${DATABASE[i]}";
	echo and will pull this many hits, sorted by best e-value...... "${SEQS[i]}";
	#This blast outputs ONLY the hit names (sseqid)
	blastp -query $PWD/$ENTRY/$ENTRY.seq.fa -db $PWD/genomes/${DATABASE[i]} -max_target_seqs ${SEQS[i]} -max_hsps 1 -outfmt "6 sseqid"  -out $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn
	
	#outputs the full BLAST results (outfmt 1).  You can look at this file to gauge completeness based on e-values, etc
	blastp -query $PWD/$ENTRY/$ENTRY.seq.fa -db $PWD/genomes/${DATABASE[i]} -max_hsps 1 -outfmt 1  -out $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.full
	
	#Defines which genome each sequence came from and appends this information
	head -n 5 $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.full
	sed -i '1s/^/hit query_id\tsubject_id\tpct_identity\taln_length\tn_of_mismatches\tgap_openings\tq_start q_end\ts_start   s_end\te_value bit_score\n/' $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.full 
	
done

#Optional: uncomment the next line to give a pause point
#read -p "Pause and check the blast output lists, then press [Enter] key to start..."

	
#Pull the full length CDS sequences from the blast databases using blastdbcmd (part of the BLAST package) and removes their stop codons using remove_stop.py
for ((i=0; i<"${#DATABASE[@]}"; i++))
do
blastdbcmd -db $PWD/genomes/${DATABASE[i]} -entry_batch $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn -out $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.blastdb.stop.fa

python $PWD/scripts/remove_stop.py $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.blastdb.stop.fa $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.blastdb.fa
python $PWD/scripts/remove_stop.py $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.blastdb.stop.fa $PWD/$ENTRY/$ENTRY.${DATABASE[i]}.seq.tblastn.blastdb.translate.fa
done

	
#Calls python fasta header parser for each of the translations, to get only the part of the header specified in the argument -hdr
for ((i=0; i<"${#DATABASE[@]}"; i++))
do
python $PWD/scripts/pull_id_fasta.py $ENTRY ${DATABASE[i]} ${HEADER[i]} .seq.tblastn.blastdb.fa
python $PWD/scripts/pull_id_fasta.py $ENTRY ${DATABASE[i]} ${HEADER[i]} .seq.tblastn.blastdb.translate.fa
done

#Creates "coding" files, this a list of genes with the genome they originally came from
for ((i=0; i<"${#DATABASE[@]}"; i++))
do
python $PWD/scripts/pull_id_fasta_coding.py $ENTRY ${DATABASE[i]} ${HEADER[i]}
done

#Combines coding files into one merged_coding.txt file, then gives it a header "taxa \t genome" which is needed later in ggtree
ls -v $PWD/$ENTRY/*.txt | xargs cat > $PWD/$ENTRY/merged_coding.txt
sleep 1
sed -i '1s/^/taxa\tgenome\n/' $PWD/$ENTRY/merged_coding.txt

#Merges all original header translated fasta files and removes any duplicates (awk).  Also simplifies all filenames .
ls -v $PWD/$ENTRY/*.seq.tblastn.blastdb.fa | xargs cat > $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.merged.fa
awk '/^>/{f=!d[$1];d[$1]=1}f' $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.merged.fa > $PWD/$ENTRY/$ENTRY.nt.merged.fa

ls -v $PWD/$ENTRY/*.seq.tblastn.blastdb.fa.parse.fa | xargs cat > $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.nt.parse.merged.fa

awk '/^>/{f=!d[$1];d[$1]=1}f' $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.nt.parse.merged.fa > $PWD/$ENTRY/$ENTRY.nt.parse.merged.fa

ls -v $PWD/$ENTRY/*.seq.tblastn.blastdb.translate.fa | xargs cat > $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.translate.merged.fa
awk '/^>/{f=!d[$1];d[$1]=1}f' $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.translate.merged.fa > $PWD/$ENTRY/$ENTRY.merged.fa

ls -v $PWD/$ENTRY/*.seq.tblastn.blastdb.translate.fa.parse.fa | xargs cat > $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.translate.fa.parse.merged.fa
awk '/^>/{f=!d[$1];d[$1]=1}f' $PWD/$ENTRY/$ENTRY.seq.tblastn.blastdb.translate.fa.parse.merged.fa > $PWD/$ENTRY/$ENTRY.parse.merged.fa

#To add an outgroup or other special sequence, use the -add and -add_dbs options. If a sequence is specified in the -add and -add_dbs option, this section adds all add_seq files to the translation file
echo "(${#ADD_SEQS[0]})"
if ((${#ADD_SEQS[0]}>0)); then
	for ((i=0; i<"${#ADD_SEQS[@]}"; i++))
	do
	python $PWD/scripts/add_seq.py $ENTRY ${ADD_DBS[i]} ${ADD_SEQS[i]}
	done
fi

#Runs clustal omega on the parsed header fasta file
clustalo -i $PWD/$ENTRY/$ENTRY.parse.merged.fa -o $PWD/$ENTRY/$ENTRY.parse.merged.clustal.fa

#Runs FastTree and then makes a copy "combinedtree.nwk" for feeding into ggtree R script
echo "Running FastTree"
FastTree -out $PWD/$ENTRY/$ENTRY.parse.merged.clustal.fa.nwk $PWD/$ENTRY/$ENTRY.parse.merged.clustal.fa
cp $PWD/$ENTRY/$ENTRY.parse.merged.clustal.fa.nwk $PWD/$ENTRY/combinedtree.nwk

#creates variables blasted and font that scales with how many seqeuences are in the output.  These variables are then fed to ggtree
blasted="" 
declare -i font=0
blasted+=$ENTRY
for ((i=0; i<"${#DATABASE[@]}"; i++)); do
	blasted+="_"
	blasted+=${DATABASE[i]}
	blasted+="_"
	blasted+=${SEQS[i]}
	font+=${SEQS[i]}
done
echo $blasted
mkdir $PWD/$ENTRY/$blasted

echo $font

#Initiaties the tree visualization R script within the subfolder $ENTRY and with the number of sequences $BLASTED
Rscript visualize-tree.r --entry $ENTRY --write $blasted 

#Calls prank for codon alignment based on the aligned amino acid sequences
prank -convert -d=$PWD/$ENTRY/$ENTRY.parse.merged.clustal.fa -dna=$PWD/$ENTRY/$ENTRY.nt.parse.merged.fa -o=$PWD/$ENTRY/output/codon.fa