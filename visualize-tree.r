# Adam Steinbrenner
# astein10@uw.edu
# UW Seattle, Dept of Biology

#Make sure these are installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("EBImage")
#BiocManager::install("treeio")
#BiocManager::install("ggtree")
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
	make_option(c("-w", "--bootstrap_offset"), action="store", default=0.07, type="numeric", 
        help="bootstrap_offset"),
	make_option(c("-q", "--heatmap_width"), action="store", default=0.3, type="numeric", 
        help="heatmap_width"),
	make_option(c("-p", "--heatmap_offset"), action="store", default=5.9, type="numeric", 
        help="heatmap_offset"),
	make_option(c("-z", "--size"), action="store", default=1, type="numeric", 
        help="size of font"),
	make_option(c("-l", "--labels_boolean"), action="store", default=1, type="numeric", 
        help="whether to include node labels, cis elements")
		)
message(option_list)
opt <- parse_args(OptionParser(option_list=option_list))

#setWD -- replace this with your own or use getwd()
working_dir <- getwd()
setwd(working_dir)

#Takes in the output from FastTree fed by blast_align_tree
dir1 <- paste(opt$entry,"/","combinedtree.nwk", sep='')
message(dir1)
tree <- read.tree(dir1)

#Takes in merged_coding file from blast_align_tree
dir2 <- paste(opt$entry,"/","merged_coding.txt", sep='')
message(dir2)

file <- paste(opt$entry,"/output/",opt$write,".pdf", sep='')
message(file)

file_csv <- paste(opt$entry,"/output/",opt$write,".csv", sep='')
message(file_csv)

file_nwk <- paste(opt$entry,"/output/",opt$write,".nwk", sep='')
message(file_nwk)




###Several heatmaps to choose from
counts_file <- read.table("datasets/both_beans_I1H1.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactor=F)


###If option -n is provided, this will take a subset of the original tree; use a node 1 deeper than you want!
nodenum <- opt$node
if (opt$node>0) {
message("Node input detected.  Subtree based on node:")
message(opt$node)
tree <- tree_subset(tree, nodenum, levels_back = 1) #Right now a bug with levels_back=0 is preventing me from specifying the node ITSELF
}
write.tree(tree,file_nwk)

#optional reroot
#tree <- root(tree, which(tree$tip.label == "AT5G03620"))

q <- ggtree(tree, size=opt$line) #size specifies line size thicknessz

dd <- read.table(dir2, sep="\t", header = TRUE, stringsAsFactor=F)
q <- q %<+% dd 

a <- as.integer(length(tree$tip.label))
print(a)

#set dimensions
node_count <- length(tree$tip.label)
print("The number of nodes in tree q is:")
print(node_count)
opt$height <- ((node_count/11.27)+0.2)

opt$width <- max(q$data$x) + 5.5

size <- opt$size
size2 <- (size)
print("The tip label font size is")
message(size)
print("The node label font size is")
message(size2)

###OPTIONAL: adds hmm coding to ggtree object q
###OPTIONAL: takes hmm_coding and adds to a separate dataframe dd2
#dir3 <- paste(getwd(),"/",opt$entry,"/","hmm_coding.txt", sep='')
#message(dir3)
#dd2 <- read.table(dir3, sep="\t", header = TRUE, stringsAsFactor=F) #optional, reads a second merged file!
#q <- q %<+% dd2 + geom_tippoint(aes(size=size2/2,color=type,shape=type)) + scale_size_identity() #weird, you need scale_size_identity! https://groups.google.com/forum/#!topic/bioc-ggtree/XHaq9Sk3b00

#Takes the tree object and converts it to a dataframe using fortify.  Simplifies it down using data.frame and then reorders it according to the graphical position!
#Apparently fortify might deprecate and switch to the "broom" package for tidying data. In the future it would be good to do this on the ggtree object "q", not "tree", so that flip functions will be reflected in the output
tips <- fortify(tree)
tips <- data.frame(tips$label,tips$y,tips$isTip)
tips <- tips[order(tips[,3],tips[,2],decreasing=T),]

#Writes the tips to a csv file.  Name is based on the option -b specified when the script is called
for (i in 1:node_count) {
write(as.matrix(tips)[,1][i],file=file_csv,sep=",",append=T)
}

#Using the output csv, this line calls python scripts to get the original nucleotide sequences from the parsed, merged, fasta file!  This output fa is in the same order as the tree!
system(paste("python scripts/extract_seq.py ",opt$entry," ",opt$write,".csv",sep=""))
system(paste("python scripts/extract_seq_alignment.py ",opt$entry," ",opt$write,".csv.",sep=""))
system(paste("python scripts/extract_seq_aa.py ",opt$entry," ",opt$write,".csv",sep=""))

#Take in datasets
print("Reading labels file")
dir3 <- paste("datasets/labels2.txt", sep='')
dd3 <- read.table(dir3, sep="\t", header = TRUE, stringsAsFactor=F, quote="")
q <- q %<+% dd3

print("Reading DEG list file")
dir4 <- paste("datasets/cowpea_venn3.txt", sep='')
dd4 <- read.table(dir4, sep="\t", header = TRUE, stringsAsFactor=F, quote="")
q <- q %<+% dd4

print("Reading counts file")
counts_file2 <- read.delim2("datasets/cowpea_counts_2.txt", sep="\t", header = TRUE, stringsAsFactor=F)
q <- q %<+% counts_file2

counts_file3 <- read.delim2("datasets/2021_04_26_pvul_1hr.txt", sep="\t", header = TRUE, stringsAsFactor=F)
q <- q %<+% counts_file3

xmax <- (max(q$data$x) + 0.07 + opt$symbol_offset)

#If option -l is present, include cis_elements
if (opt$labels_boolean > 0) {
	print("Reading ciselements file")
	cis_elements <- read.delim2("C:/science/blast_align_tree/datasets/cis_element_counts.csv", sep=",", header = TRUE, stringsAsFactor=F)
	q <- q %<+% cis_elements
	q <- q + geom_tiplab(aes(label=LuxT), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+4)) +
		geom_tiplab(aes(label=LuxA), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+4.2)) +
		geom_tiplab(aes(label=Gbox), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+4.4)) +
		geom_tiplab(aes(label=EE), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+4.6)) +
		geom_tiplab(aes(label=CBS_A), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+4.8)) +
		geom_tiplab(aes(label=CBS_B), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+5.0)) +
		annotate("text",size=size,x=c(xmax+4,xmax+4.2,xmax+4.4,xmax+4.6,xmax+4.8,xmax+5.0),y=max(q$data$y)+0.7,label=c("LuxT","LuxA","Gbox","EE","CBS-A","CBS-B"))
}

#If option -l is present, include node numbering
if (opt$labels_boolean > 0) {
	q <- q + geom_text2(aes(subset=!isTip, label=node), color="red", hjust=-.3, size=size2) #node labels
}

#converts "labels" bootstrap to a 1-100 integer
d <- q$data
d <- d[!d$isTip,]
d$label <- as.integer(100*(as.numeric(d$label)))
message(d$label)

#q <- viewClade(q,595)
d <- d[d$label < 75,] #option to only show some threshold bootstrap


#Choose color scheme, Either a pre-defined one (values = col) or just enough to include all species	

#create a character vector to define the color schemes for each aes color factor
col1 <- c("black","darkgoldenrod1","gray50","red1","steelblue1","slateblue1","black","green4","blue","green")
col2 <- c("black","darkgoldenrod1","gray50","red1","steelblue1","slateblue1","green4","black","blue","green","black")
names(col2) <- c("black","darkgoldenrod1","gray50","red1","steelblue1","slateblue1","TAIR10cds.fa","Vung469cds.fa","Pvul218cds.fa","Slyc250cds.fa","Pkinase")
copb <- c("black","#4a6630","#0014c2","#ffa519")
standard_colors <- c("black","blue","blue","darkgoldenrod","purple","orange","darkgreen","black","blue","purple","darkgreen","black","blue","orange","purple","darkgreen","cadetblue","deeppink","darkgoldenrod","brown4","olivedrab2","cyan","magenta","#008080","lavender","#bcf60c","#aaffc3","#ffd8b1","#fabebe","#fffac8","green")

	



q <- q + 
	geom_treescale(fontsize = opt$size, linesize = opt$line,x=0, y=max(q$data$y), width=0.1) + 
	#cowpeaDEGs
	geom_tiplab(aes(label=onehr_dmg),color="black", size=(opt$symbol_size),align=T, linetype=NA, offset=(opt$symbol_offset+0.4)) +
	geom_tiplab(aes(label=in_1_subset),color="darkgoldenrod1", size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+0.6)) +
	geom_tiplab(aes(label=sixhr_dmg),color="black", size=(opt$symbol_size),align=T, linetype=NA, offset=(opt$symbol_offset+1.4)) +
	geom_tiplab(aes(label=in_6_subset),fontface="bold",color="red1", size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+1.6)) +
	
	geom_tiplab(aes(label=D_onehr_dmg), color="black", size=(opt$symbol_size),align=T, linetype=NA, offset=(opt$symbol_offset+0.4)) +
	geom_tiplab(aes(label=D_in_1_subset), color="steelblue1", size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+0.6)) +
	geom_tiplab(aes(label=D_sixhr_dmg), color="black", size=(opt$symbol_size),align=T, linetype=NA, offset=(opt$symbol_offset+1.4)) +
	geom_tiplab(aes(label=D_in_6_subset), color="slateblue1", size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+1.6)) 

xlabs_positions <- c(xmax+0.4,xmax+0.75,xmax+1.4,xmax+1.7,xmax+2.4,xmax+2.75,xmax+3.0,xmax+3.3,xmax+3.6)


q <- q +
	geom_tiplab(size=size,offset=opt$label_offset,aes(color=hit,fontface="bold")) + #tip labels (gene names)
	geom_tiplab(aes(label=symbol,color=hit), size=opt$symbol_size,align=T, linetype=NA, offset=opt$symbol_offset) + #gene symbol names
	
	#cowpeaDEGs
	#geom_tiplab(aes(label=onehr_dmg), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+1)) +
	#geom_tiplab(aes(label=in_1_subset), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+1.1)) +
	#geom_tiplab(aes(label=sixhr_dmg), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+2.1)) +
	#sgeom_tiplab(aes(label=in_6_subset), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+2.2)) +

	#cowpea counts data
	geom_tiplab(aes(label=U1_avg), color="gray50",size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+2.4)) +
	geom_tiplab(aes(label=H1_avg,color=onehr_dmg_color,fontface=onehr_dmg_font), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+2.7)) +
	geom_tiplab(aes(label=I1_avg,color=onehr_in_color,fontface=onehr_in_font), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+3.0)) +
	geom_tiplab(aes(label=H6_avg,color=sixhr_dmg_color,fontface=sixhr_dmg_font), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+3.3)) +
	geom_tiplab(aes(label=I6_avg,color=sixhr_in_color,fontface=sixhr_in_font), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+3.6)) +

	geom_tiplab(aes(label=log1hr), size=opt$symbol_size,align=T, linetype=NA, offset=(opt$symbol_offset+3.8)) + #common bean
	theme(legend.position = "none") +
	
	#add column labels at the top, customize this for your dataset
	annotate("text",size=size,x=xlabs_positions,y=max(q$data$y)+0.7,label=c("Dmg1","Inceptin_effect","Dmg6","Inceptin_effect","Undmg","Dmg1","In1","Dmg6","In6")) +
	scale_colour_manual(values=standard_colors) #std
	



#version of color scale for domains
#msa_colors <- c("gray85","red","orange","green",rep(c("black"),each=20))
msa_colors <- c("gray85",rep(c("black"),each=30))
aa <- paste(opt$entry,"/output/",opt$write,".csv.aa.fa", sep='')
msa_pre <- paste(opt$entry,"/output/",opt$write,".csv.aa.ungapped.fa", sep='')
msa <- paste(opt$entry,"/output/",opt$write,".csv.aa.ungapped.headers.fa", sep='')
system(paste("trimal -in ",aa," -out ",msa_pre," -noallgaps",sep=""))
system(paste("python scripts/remove_header.py ",msa_pre," ",msa,sep=""))

#heatmap colors and limits
upper <- 7 #these are the upper and lower limits, used to set colors!  Use 8 for lipids_deg, less for both beans
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
#q

gheatmap(q,counts_file, offset = opt$heatmap_offset, width=opt$heatmap_width+3, font.size=size, colnames_angle=-20, hjust=0, color="black") + scale_fill_gradient2(low=low_color,high=high_color,mid="white",limits=c(lower,upper)) #heatmap
dev.off()

if (opt$width>0) {
pdf(paste(file,".MSA.pdf",sep=""), height=opt$height, width=opt$width)
} else {
pdf(paste(file,"_MSA",sep=""))
}
print("Using this MSA sequence")
print(msa)
msaplot(q,msa,offset=opt$heatmap_offset,color=msa_colors)

dev.off()


