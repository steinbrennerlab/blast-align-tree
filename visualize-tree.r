# Adam Steinbrenner
# astein10@uw.edu
# UW Seattle, Dept of Biology

#Make sure these are installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("EBImage")
#BiocManager::install("treeio")
#BiocManager::install("ggtree")
#quiBiocManager::install("Biostrings")
#install.packages("phytools")
#install.packages("optparse")
#install.packages("tidyselect")
#install.packages("tidyselect")
#install.packages("labeling")
#install.packages("tidyverse")

#Certain ggtree functions may require the development version of treeio; make sure it is current from github
#install.packages("devtools")
#devtools::install_github("GuangchuangYu/treeio")
library("ggplot2")
library("treeio")
library("phytools") # for sims and ASRs
library("EBImage") # for images
library("ggtree")
library("optparse")
library("tidytree")
library("ape")

sessionInfo()

### Options. Example usage below:
### Rscript visualize-tree.R -e AT4G33430.1 -b output 
option_list <- list( 
    make_option(c("-e", "--entry"), action="store", type="character",
        help="entry"),
	make_option(c("-b", "--write"), action="store", type="character", 
        help="parameters term for file writing"),
	make_option(c("-a", "--reroot"), action="store", type="character", 
        help="node to reroot"),
	make_option(c("-n", "--node"), action="store", default=0, type="integer", 
        help="feed script a node and it will subset one node higher"),
	make_option(c("-f", "--width"), action="store", default=6, type="numeric", 
        help="width for pdf, optional"),
	make_option(c("-g", "--height"), action="store", default=3, type="numeric", 
        help="height for pdf, optional"), #code below replaces this with a calculated value
	make_option(c("-i", "--line"), action="store", default=0.1, type="numeric", 
        help="line_width"),
	make_option(c("-x", "--push"), action="store", default=0.3, type="numeric", 
        help="push right side to make room"),
	make_option(c("-r", "--symbol_offset"), action="store", default=1, type="numeric", 
        help="font size for gene symbols"),
	make_option(c("-y", "--label_offset"), action="store", default=0.05, type="numeric", 
        help="font size for gene symbols"),
	make_option(c("-m", "--symbol_size"), action="store", default=1, type="numeric", 
        help="offset for gene symbols"),
	make_option(c("-w", "--bootstrap_offset"), action="store", default=0.05, type="numeric", 
        help="bootstrap_offset"),
	make_option(c("-q", "--heatmap_width"), action="store", default=0.3, type="numeric", 
        help="heatmap_width"),
	make_option(c("-p", "--heatmap_offset"), action="store", default=5.9, type="numeric", 
        help="heatmap_offset"),
	make_option(c("-z", "--size"), action="store", default=1, type="numeric", 
        help="size of font"),
	make_option(c("-l", "--labels_boolean"), action="store", default=1, type="numeric", 
        help="whether to include node labels"),
	make_option(c("-k", "--bootstrap_boolean"), action="store", default=0, type="numeric", 
        help="whether to include bootstrap labels")
		)
message(option_list)
opt <- parse_args(OptionParser(option_list=option_list))

#setWD -- replace this with your own or use getwd()
working_dir <- getwd()
setwd(working_dir)

#Takes in the output from FastTree fed by blast_align_tree bash script
tree_newick <- paste(opt$entry,"/","combinedtree.nwk", sep='')
message(tree_newick)
tree <- read.tree(tree_newick)

#Takes in merged_coding file from blast_align_tree
gene_species_list <- paste(opt$entry,"/","merged_coding.txt", sep='')
message(gene_species_list)

#define output filenames
file <- paste(opt$entry,"/output/",opt$write,".pdf", sep='')
message(file)
file_csv <- paste(opt$entry,"/output/",opt$write,".csv", sep='')
message(file_csv)
file_nwk <- paste(opt$entry,"/output/",opt$write,".nwk", sep='')
message(file_nwk)


#optional reroot
if (length(opt$reroot)>0) {
	tree <- reroot(tree, which(tree$tip.label == opt$reroot))
}

###If option -n is provided, this will take a subset of the original tree!
nodenum <- opt$node
if (opt$node>0) {
message("Node input detected.  Subtree based on node:")
message(opt$node)
tree <- tree_subset(tree, nodenum, levels_back = 0)
}
write.tree(tree,file_nwk)

p <- ggtree(tree, size=opt$line) #size specifies line size thickness

#add species list to the tree object
dd <- read.table(gene_species_list, sep="\t", header = TRUE, stringsAsFactor=F)
p <- p %<+% dd 

#How many tips does the tree have
a <- as.integer(length(tree$tip.label))
print(a)

#set dimensions for output based on number of tips
node_count <- length(tree$tip.label)
print("The number of nodes in tree p is:")
print(node_count)
opt$height <- ((node_count/11.27)+0.2)
xmax <- (max(p$data$x) + 0.07 + opt$symbol_offset)

opt$width <- max(p$data$x) + 5.5

#Define tip and node label sizes
size <- opt$size
size2 <- (size)
print("The tip label font size is")
message(size)
print("The node label font size is")
message(size2)

###OPTIONAL: adds hmm coding to ggtree object p
###OPTIONAL: takes hmm_coding and adds to a separate dataframe dd2
#dir3 <- paste(getwd(),"/",opt$entry,"/","hmm_coding.txt", sep='')
#message(dir3)
#dd2 <- read.table(dir3, sep="\t", header = TRUE, stringsAsFactor=F) #optional, reads a second merged file!
#p <- p %<+% dd2 + geom_tippoint(aes(size=size2/2,color=type,shape=type)) + scale_size_identity() #weird, you need scale_size_identity! https://groups.google.com/forum/#!topic/bioc-ggtree/XHaq9Sk3b00

#Takes the tree object and converts it to a dataframe using fortify and data.frame. then reorders it according to the graphical position
#Apparently fortify might deprecate and switch to the "broom" package for tidying data. In the future it would be good to do this on the ggtree object "p", not "tree", so that flip functions will be reflected in the output
tips <- fortify(tree)
tips <- data.frame(tips$label,tips$y,tips$isTip)
tips <- tips[order(tips[,3],tips[,2],decreasing=T),]

#Writes the graphically ordered tips to a csv file
for (i in 1:node_count) {
write(as.matrix(tips)[,1][i],file=file_csv,sep=",",append=T)
}

#Using the output csv, this line calls python scripts to get the original nucleotide sequences from the parsed, merged, fasta file.  This output fa is in the same order as the tree
system(paste("python scripts/extract_seq.py ",opt$entry," ",opt$write,".csv",sep=""))
system(paste("python scripts/extract_seq_alignment.py ",opt$entry," ",opt$write,".csv.",sep=""))
system(paste("python scripts/extract_seq_aa.py ",opt$entry," ",opt$write,".csv",sep=""))

#Read datasets for tree visualizatioon
print("Reading labels file")
dir3 <- paste("datasets/gene_symbols.txt", sep='')
dd3 <- read.table(dir3, sep="\t", header = TRUE, stringsAsFactor=F, quote="")
p <- p %<+% dd3

print("Reading DEG list file")
dir4 <- paste("datasets/Vu_DEGs.txt", sep='')
dd4 <- read.table(dir4, sep="\t", header = TRUE, stringsAsFactor=F, quote="")
p <- p %<+% dd4

#Choose color scheme, Either a pre-defined one (values = col) or just enough to include all species	

#create a character vector to define the color schemes for each aes color factor
standard_colors <- c("black","blue","darkgoldenrod","purple","orange","darkgreen","black","blue","purple","darkgreen","black","blue","orange","purple","darkgreen","cadetblue","deeppink","darkgoldenrod","brown4","olivedrab2","cyan","magenta","#008080","lavender","#bcf60c","#aaffc3","#ffd8b1","#fabebe","#fffac8","green")


#Reads RNAseq data for cowpea genes
counts_file <- read.table("datasets/Vu_inceptin_1hr.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactor=F)

print("Reading counts file")
counts_file2 <- read.delim2("datasets/Vu_counts.txt", sep="\t", header = TRUE, stringsAsFactor=F)
p <- p %<+% counts_file2

counts_file3 <- read.delim2("datasets/Pv_counts.txt", sep="\t", header = TRUE, stringsAsFactor=F)
p <- p %<+% counts_file3

#Add categories from Steinbrenner 2022, Plant Journal -- cowpea transcriptomic responses to wounding and In11 peptide
p <- p + 
	geom_treescale(fontsize = opt$size, linesize = opt$line,x=0, y=max(p$data$y), width=0.1) + 
	#cowpeaDEGs
	geom_tiplab(aes(label=onehr_dmg),color="black", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+0.4)) +
	geom_tiplab(aes(label=in_1_subset),color="darkgoldenrod1", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+0.6)) +
	geom_tiplab(aes(label=sixhr_dmg),color="black", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+1.4)) +
	geom_tiplab(aes(label=in_6_subset),fontface="bold",color="red1", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+1.6)) +
	
	geom_tiplab(aes(label=D_onehr_dmg), color="black", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+0.4)) +
	geom_tiplab(aes(label=D_in_1_subset), color="steelblue1", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+0.6)) +
	geom_tiplab(aes(label=D_sixhr_dmg), color="black", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+1.4)) +
	geom_tiplab(aes(label=D_in_6_subset), color="slateblue1", size=1,align=T, linetype=NA, offset=(opt$symbol_offset+1.6)) 

xlabs_positions <- c(xmax+0.4,xmax+0.75,xmax+1.4,xmax+1.7,xmax+2.4,xmax+2.75,xmax+3.0,xmax+3.3,xmax+3.6,xmax+3.8,xmax+4.0,xmax+4.2)

#Add data from Steinbrenner 2022, Plant Journal -- cowpea transcriptomic responses to wounding and In11 peptide
p <- p +
	geom_tiplab(size=size,offset=opt$label_offset,aes(color=hit,fontface="bold")) + #tip labels (gene names)
	geom_tiplab(aes(label=symbol,color=hit), size=opt$symbol_size,align=T, linetype=NA, offset=opt$symbol_offset) + #gene symbol names
	scale_colour_manual(values=standard_colors) + 
	#cowpeaDEGs
	#geom_tiplab(aes(label=onehr_dmg), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+1)) +
	#geom_tiplab(aes(label=in_1_subset), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+1.1)) +
	#geom_tiplab(aes(label=sixhr_dmg), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+2.1)) +
	#geom_tiplab(aes(label=in_6_subset), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+2.2)) +

	#cowpea counts data
	annotate("text",size=size,,x=c(xmax+0.4),y=max(p$data$y)+1.6,hjust=0,label=c("Steinbrenner et al 2022 DOI:10.1111/tpj.15732 3'RNASeq DEG categories and log2(counts)"))+
	geom_tiplab(aes(label=U1_avg), color="gray50",size=1,align=T, linetype=NA, offset=(opt$symbol_offset+2.4)) +
	geom_tiplab(aes(label=H1_avg,fontface=onehr_dmg_font), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+2.7)) +
	geom_tiplab(aes(label=I1_avg,fontface=onehr_in_font), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+3.0)) +
	geom_tiplab(aes(label=H6_avg,fontface=sixhr_dmg_font), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+3.3)) +
	geom_tiplab(aes(label=I6_avg,fontface=sixhr_in_font), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+3.6)) +
	
	#bean l2FCdata
	annotate("text",size=size,,x=c(xmax+3.6),y=max(p$data$y)+2.1,hjust=0,label=c("Unpublished Phaseolus vulgaris"))+
	annotate("text",size=size,,x=c(xmax+3.6),y=max(p$data$y)+1.6,hjust=0,label=c("PE-RNAseq log2(counts)"))+
	geom_tiplab(aes(label=Pv_Undmg), size=1, align=T, linetype=NA, offset=(opt$symbol_offset+3.8)) +
	geom_tiplab(aes(label=Dmg1hr), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+4.0)) +
	geom_tiplab(aes(label=Inc1hr), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+4.2)) +
	theme(legend.position = "none") +
	annotate("text",size=size,x=xlabs_positions,y=max(p$data$y)+0.7,label=c("Dmg1","Inceptin_effect","Dmg6","Inceptin_effect","Undmg","Dmg1","In1","Dmg6","In6","Undmg","Dmg1hr","In1hr"))


#Add unpublished data from Natalia Guayazan-Palacios. Email us for more info astein10@uw.edu -- uncomment if needed
DEG_file <- read.table("datasets/additional_FC_data.txt", sep="\t", header = TRUE, stringsAsFactor=F)
head(DEG_file)
names(DEG_file)
p <- p %<+% DEG_file
p <- p + geom_tiplab(aes(label=NGP_WH_1h_l2FC), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+4.4)) + 
	geom_tiplab(aes(label=NGP_InH_1h_l2FC), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+4.6)) + 
	geom_tiplab(aes(label=NGP_flg22H_1h_l2FC), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+4.8)) + 
	geom_tiplab(aes(label=NGP_WH_6h_l2FC), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.0)) + 
	geom_tiplab(aes(label=NGP_InH_6h_l2FC), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.2)) + 
	geom_tiplab(aes(label=NGP_flg22H_6h_l2FC), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.4))+
	annotate("text",size=size,x=c(xmax+4.4,xmax+4.6,xmax+4.8,xmax+5.0,xmax+5.2,xmax+5.4),y=max(p$data$y)+0.7,label=c("WH_1h","InH","flg22H","WH_6h","InH","flg22H"))+
	annotate("text",size=size,,x=c(xmax+4.4),y=max(p$data$y)+2.8,hjust=0,label=c("Guayazan-Palacios unpubl. data log2FC"))

#Add data from Bjornsen et al 2021. 
PAMP_degs <- read.delim2("datasets/bjornsen_pamps.txt", sep="\t", header = TRUE, stringsAsFactor=F)
p <- p %<+% PAMP_degs
p <- p + 
	annotate("text",size=size,,x=c(xmax+5.6),y=max(p$data$y)+1.6,hjust=0,label=c("Bjornsen_PAMP-max-log2(FC)"))+
	geom_tiplab(aes(label=flg22), size=1, align=T, linetype=NA, offset=(opt$symbol_offset+5.6)) +
	geom_tiplab(aes(label=elf18), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.8)) +
	geom_tiplab(aes(label=nlp20), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+6.0)) +
	geom_tiplab(aes(label=og), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+6.2)) +
	geom_tiplab(aes(label=chitin), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+6.4)) +
	annotate("text",size=size,x=c(xmax+5.6,xmax+5.8,xmax+6.0,xmax+6.2,xmax+6.4),y=max(p$data$y)+0.7,label=c("flg22","elf18","nlp20","OG","chitin"))

#print("Reading ciselements file")
#cis_elements <- read.delim2("C:/science/blast_align_tree/datasets/cis_element_counts.csv", sep=",", header = TRUE, stringsAsFactor=F)
#p <- p %<+% cis_elements
#p <- p + geom_tiplab(aes(label=LuxT), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.4)) +
#	geom_tiplab(aes(label=LuxA), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.6)) +
#	geom_tiplab(aes(label=Gbox), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+5.8)) +
#	geom_tiplab(aes(label=EE), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+6.0)) +
#	geom_tiplab(aes(label=CBS_A), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+6.2)) +
#	geom_tiplab(aes(label=CBS_B), size=1,align=T, linetype=NA, offset=(opt$symbol_offset+6.4)) +
#		annotate("text",size=size,x=c(xmax+5.4,xmax+5.6,xmax+5.8,xmax+6.0,xmax+6.2,xmax+6.4),y=max(p$data$y)+0.7,label=c("LuxT","LuxA","Gbox","EE","CBS-A","CBS-B"))
#}

#If option -l is present, include node numbering
if (opt$labels_boolean > 0) {
	p <- p + geom_text2(aes(subset=!isTip, label=node), color="red", hjust=-.3, size=size2) #node labels
}

#converts "labels" bootstrap to a 1-100 integer
d <- p$data
d <- d[!d$isTip,]
d$label <- as.integer(100*(as.numeric(d$label)))
message(d$label)

#p <- viewClade(p,595)
#d <- d[d$label < 75,] #option to only show some threshold bootstrap
if (opt$bootstrap_boolean > 0) {
	p <- p + geom_text(data=d, aes(label=label), size=size, nudge_x=opt$bootstrap_offset)  #bootstraps
}

#Create a visualization of the multiple sequence alignment with any amino acid indicated as black, and any gap indicated as grey

##msa_colors <- c("gray85","red","orange","green",rep(c("black"),each=20)) ##version of color scale for domains
msa_colors <- c("gray85",rep(c("black"),each=30))
aa <- paste(opt$entry,"/output/",opt$write,".csv.aa.fa", sep='')
msa_pre <- paste(opt$entry,"/output/",opt$write,".csv.aa.ungapped.fa", sep='')
msa <- paste(opt$entry,"/output/",opt$write,".csv.aa.ungapped.headers.fa", sep='')
system(paste("trimal -in ",aa," -out ",msa_pre," -noallgaps",sep=""))
system(paste("python scripts/remove_header.py ",msa_pre," ",msa,sep=""))

#Set visualization parameters for heatmaps
##heatmap colors and limits
upper <- 7 #these are the upper and lower limits used to set colors. If outside the limits cell is gray
lower <- -7
low_color<-"#000099" #blue
high_color<-"#FF0000" #red
#high_color<-"#ffa500" #orange
#high_color<-"#FFCC33" #yellow


#open pdf to write
if (opt$width>0) {
pdf(file, height=opt$height, width=opt$width)
} else {
pdf(file)
}

##Choose one: heatmap or just the tree
p
#gheatmap(p,counts_file, offset = opt$heatmap_offset, width=opt$heatmap_width+3, font.size=size, colnames_angle=-20, hjust=0, color="black") + scale_fill_gradient2(low=low_color,high=high_color,mid="white",limits=c(lower,upper)) #heatmap
dev.off()

if (opt$width>0) {
pdf(paste(file,".MSA.pdf",sep=""), height=opt$height, width=opt$width)
} else {
pdf(paste(file,"_MSA",sep=""))
}
print("Using this MSA sequence")
print(msa)
msaplot(p,msa,offset=opt$heatmap_offset,color=msa_colors)

dev.off()


