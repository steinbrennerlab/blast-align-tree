# Adam Steinbrenner
# astein10@uw.edu
# UW Seattle, Dept of Biology

#Make sure these are installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("EBImage")
#BiocManager::install("treeio")
#BiocManager::install("ggtree")
#BiocManager::install("Biostrings")
#install.packages("phytools")
#install.packages("optparse")
#install.packages("tidyselect")
#install.packages("tidyselect")
#install.packages("labeling")
#install.packages("tidyverse")
#install.packages("broom")

#Certain ggtree functions may require the development version of treeio; make sure it is current from github
#install.packages("devtools")
#devtools::install_github("GuangchuangYu/treeio")
library("ggplot2")
library("treeio")
library("phytools") # for sims and ASRs
#("EBImage") # for images
library("ggtree")
library("optparse")
library("tidytree")
library("ape")
library("broom")

sessionInfo()

### Options. Example usage below:
### Rscript visualize-tree.R -e AT5G45250.1 -b rps4 -k 1 -l 0 -m 1 -n 852
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
	make_option(c("-r", "--symbol_offset"), action="store", default=0.5, type="numeric", 
        help="font size for gene symbols"),
	make_option(c("-y", "--label_offset"), action="store", default=0.05, type="numeric", 
        help="font size for gene symbols"),
	make_option(c("-m", "--symbol_size"), action="store", default=1, type="numeric", 
        help="offset for gene symbols"),
	make_option(c("-w", "--bootstrap_offset"), action="store", default=0.05, type="numeric", 
        help="bootstrap_offset"),
	make_option(c("-q", "--heatmap_width"), action="store", default=0.3, type="numeric", 
        help="heatmap_width"),
	make_option(c("-p", "--heatmap_offset"), action="store", default=2.0, type="numeric", 
        help="heatmap_offset"),
	make_option(c("-z", "--size"), action="store", default=1, type="numeric", 
        help="size of font"),
	make_option(c("-l", "--labels_boolean"), action="store", default=1, type="numeric", 
        help="whether to include node labels"),
	make_option(c("-k", "--bootstrap_boolean"), action="store", default=0, type="numeric", 
        help="whether to include bootstrap labels"),
	make_option("--genomeLabel_boolean", action="store", default=0, type="numeric", 
        help="whether to include genome labels"),
	make_option(c("-j", "--tree_type"), action="store", default=1, type="numeric", 
        help="Specify a (1) plain tree, (2) heatmap (2), or (3) multiple sequence alignment"),
	make_option(c("--features"),type = "character",default = NULL,
		help = "Optional path to features.txt. If omitted, will try <entry>/output/features.txt. If missing/unreadable, feature overlay is skipped."),
	make_option(c("--subdir"),type = "character",default = NULL,
		help = "Relative subdirectory to append under both --entry and --write (e.g., 'runs/20250903_1557'). Leading slashes will be stripped."),
	make_option(c("--dist_type"), action="store", type="character", default="patristic", 
		help="Distance type to compute: 'patristic' or 'phenetic' (case-insensitive)"),
	make_option(c("--dist_digits"), action="store", type="integer", default=3, help="Number of digits to round distances when printed")
	)
message(option_list)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- path resolution (does not mutate opt$entry/opt$write) ----
append_rel <- function(base, subdir) {
  if (is.null(subdir) || !nzchar(subdir)) return(base)
  file.path(base, gsub("^[/\\\\]+", "", subdir))
}
ENTRY_DIR <- append_rel(opt$entry, opt$subdir)
WRITE_DIR <- append_rel(opt$subdir, opt$write)
message(sprintf("ENTRY_DIR=%s", ENTRY_DIR))
message(sprintf("WRITE_DIR=%s", WRITE_DIR))

###
# Check for presence of features.txt in the output folder
###

`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

feature_file_default <- file.path(ENTRY_DIR, "output", "features.txt")
feature_file <- opt$features %||% feature_file_default

read_features_safe <- function(path) {
  if (!is.character(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  if (isTRUE(file.info(path)$size == 0)) return(NULL)
  tryCatch(
    utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = "",
      strip.white = TRUE
    ),
    error = function(e) NULL
  )
}

required_cols_present <- function(df) {
  req <- c("label", "aa_start", "aa_end", "feature")
  !is.null(df) && all(req %in% names(df))
}

features_df <- read_features_safe(feature_file)
has_features <- required_cols_present(features_df)
if (!has_features) {
  message(sprintf("No usable features file found at '%s'. Skipping feature overlay.", feature_file))
}

###
# End feature file check
###

#setWD
working_dir <- getwd()
setwd(working_dir)

#Takes in the output from FastTree fed by blast_align_tree bash script
tree_newick <- paste(ENTRY_DIR,"/","combinedtree.nwk", sep='')
message(tree_newick)
tree <- read.tree(tree_newick)

#Takes in merged_genome_mapping file from blast_align_tree
gene_species_list <- paste(ENTRY_DIR,"/","merged_genome_mapping.txt", sep='')
message(gene_species_list)

#define output filenames
file <- paste(ENTRY_DIR,"/output/",opt$write,".pdf", sep='')
message(file)
file_csv <- paste(ENTRY_DIR,"/output/",opt$write,".csv", sep='')
message(file_csv)
file_nwk <- paste(ENTRY_DIR,"/output/",opt$write,".nwk", sep='')
message(file_nwk)



#optional reroot if option -a is specified
if (length(opt$reroot)>0) {
	tree <- reroot(tree, which(tree$tip.label == opt$reroot))
}

###If option -n is provided, this will take a subset of the original tree
nodenum <- opt$node
if (opt$node>0) {
message("Node input detected.  Subtree based on node:")
message(opt$node)
tree <- tree_subset(tree, nodenum, levels_back = 0)
}
write.tree(tree,file_nwk)

p <- ggtree(tree, size=opt$line) #size specifies line size thickness

p <- p + theme(legend.title = element_text(size = 3), 
               legend.text = element_text(size = 3))

#add species list to the tree object; this is the "merged coding" file from the blast-align-tree bash script, a tab separated text file of taxa and genome
dd <- read.table(gene_species_list, sep="\t", header = TRUE, stringsAsFactor=F)
p <- p %<+% dd 

#How many tips does the tree have
a <- as.integer(length(tree$tip.label))

#set height dimension for output based on number of tips
node_count <- length(tree$tip.label)
print("The number of nodes in tree p is:")
print(node_count)
opt$height <- ((node_count/11.27)+0.5)

#Define tip and node label sizes
size <- opt$size
size2 <- (size)
print("The tip label font size is")
message(size)
print("The node label font size is")
message(size2)

# If option -l is specified, include node numbering
if (opt$labels_boolean > 0) {
	p <- p + geom_text2(aes(subset=!isTip, label=node), color="red", hjust=-.3, size=size2) #node labels
}

xmax <- max(p$data$x)

###########################

#Convert the tree to a dataframe then reorder according to the graphical position
get_tree_df <- function(tree, layout = "rectangular") {
  pp <- ggtree::ggtree(tree, layout = layout)
  pp$data  # tibble with x, y, parent, isTip, node, label, etc.
}
tips <- get_tree_df(tree)
tips <- data.frame(tips$label,tips$y,tips$isTip)
tips <- tips[order(tips[,3],tips[,2],decreasing=T),]

# Write the tip names to a csv file in the same order as the tree visualization
# Pull out the first column
v <- as.character(as.matrix(tips)[,1])
# Write all lines
write(v, file = file_csv, sep = "\n")

#Using the csv of gene names, extract_seq.py extracts original nucleotide sequences from the parsed, merged, fasta file.  The output fasta file is in the same order as the tree!
system(paste("python scripts/extract_seq.py ",ENTRY_DIR," ",opt$entry," ",opt$write,".csv",sep=""))
system(paste("python scripts/extract_seq_aa.py ",ENTRY_DIR," ",opt$entry," ",opt$write,".csv",sep=""))

# Read gene symbols to display next to annotated genes
  ## Read gene symbols file from root directory
print("Reading labels file")
dir3 <- paste("gene_symbols.txt", sep='')
dd3 <- read.table(dir3, sep="\t", header = TRUE, stringsAsFactor=F, quote="")
p <- p %<+% dd3

## create a character vector to define the color schemes for each aes color factor
standard_colors <- c(
  "black", "blue", "darkgoldenrod", "purple", "orange", "darkgreen",  # 1-6
  "slategray4",    #  7
  "steelblue4",    #  8
  "sienna4",       #  9
  "darkslateblue", # 10
  "mediumslateblue",#11
  "peru",          # 12
  "darkolivegreen4",#13
  "seagreen4",     # 14
  "maroon4",       # 15
  "midnightblue",  # 16
  "darkorchid4",   # 17
  "chocolate4",    # 18
  "cadetblue4",    # 19
  "rosybrown4",    # 20
  "slateblue4",    # 21
  "olivedrab4",    # 22
  "darkcyan",      # 23
  "tomato4",       # 24
  "forestgreen",   # 25
  "indianred4",    # 26
  "paleturquoise4",# 27
  "mediumpurple4", # 28
  "darkmagenta",   # 30
  "gray35"         # 31
)
standard_colors_generator <- function(n) {
  rep(standard_colors, length.out = n)
}
standard_colors <- standard_colors_generator(100)

  ## Add gene symbol labels to the tree using offsets specified in the options
p <- p +
	geom_tiplab(size=size,offset=opt$label_offset,aes(color=genome,fontface="bold")) + 
	#tip labels (gene names) colored by species
	geom_tiplab(aes(label=symbol,color=genome), size=opt$symbol_size,align=T, linetype=NA, offset=opt$symbol_offset) + 
	#gene symbol names
	scale_colour_manual(values=standard_colors)

if (opt$genomeLabel_boolean > 0) {
	p <- p +
	geom_tiplab(aes(label=genome,color=genome), size=opt$symbol_size,align=T, linetype=NA, offset=opt$symbol_offset+1.3) + 
	#genome label names
	scale_colour_manual(values=standard_colors)
}

#converts "labels" bootstrap to a 1-100 integer
d <- p$data
d <- d[!d$isTip,]
d$label <- as.integer(100*(as.numeric(d$label)))
message(d$label)
#option to only show some threshold bootstrap
#d <- d[d$label < 75,] 

#If option -k is specified, include bootstraps
if (opt$bootstrap_boolean > 0) {
	p <- p + geom_text(data=d, aes(label=label), size=size, nudge_x=opt$bootstrap_offset)  #bootstraps
}



# Read all datasets in the /datasets folder
## Function to read datasets from a folder

read_datasets <- function(folder_path) {
  print("read datasets function")
  dataset_files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)
  datasets <- list()
  
  for (file in dataset_files) {
    print(file)
	dataset_name <- tools::file_path_sans_ext(basename(file))
    datasets[[dataset_name]] <- read.delim2(file, sep = "\t", header = TRUE, stringsAsFactor = FALSE)
	print(dataset_name)
  }
  
  return(datasets)
}

## Function to add datasets to the tree, calculating the text offset for each column
add_datasets_to_tree <- function(p, datasets, base_offset) {
  cumulative_offset <- opt$symbol_offset + 0.5
  for (i in seq_along(datasets)) {
    dataset_name <- names(datasets)[i]
    dataset <- datasets[[i]]
    
    p <- p %<+% dataset
    
    # Add dataset name as text annotation
    # This should stay fixed for each dataset
    p <- p + annotate("text", size = opt$size, x = max(p$data$x) + base_offset + cumulative_offset,
                      y = max(p$data$y) + 2.4, label = dataset_name, fontface = "bold", hjust = 0)
    
    # Add columns as tip labels
    columns_to_display <- names(dataset)[-1]  # Exclude the first column
    for (j in seq_along(columns_to_display)) {
      column_name <- columns_to_display[j]
      column_offset <- (j - 1) * 0.2
      
      # The column data should use the cumulative offset plus its own offset
      p <- p + geom_tiplab(aes(label = .data[[column_name]]), size = 1, align = TRUE,
                           linetype = NA, offset = base_offset + cumulative_offset + column_offset)
      
      # Column names should align with their data
      p <- p + annotate("text", size = opt$size * 0.8, 
                        x = max(p$data$x) + base_offset + cumulative_offset + column_offset,
                        y = max(p$data$y) + 0.7, 
                        label = column_name, 
                        angle = 45, 
                        hjust = 0)
    }
    
    # Update the cumulative offset after processing each dataset
    cumulative_offset <- cumulative_offset + length(columns_to_display) * 0.2 + 0.15
  }
  total_offset <<- cumulative_offset
  return(p)
}


  # Read datasets
datasets <- read_datasets("datasets/")

  ## Create a duplicate tree for heatmap and MSA trees below -- it will not have datasets
p2 <- p

# Add datasets to the tree
p <- add_datasets_to_tree(p, datasets, base_offset = 0)
  


opt$width <- max(p$data$x) + total_offset

#Creates a visualization of the multiple sequence alignment with any amino acid indicated as black, and any gap indicated as grey
##msa_colors <- c("gray85","red","orange","green",rep(c("black"),each=20)) ##version of color scale for domains
msa_colors <- c("gray85",rep(c("black"),each=30))
aa <- paste(ENTRY_DIR,"/output/",opt$write,".csv.aa.fa", sep='')
msa_pre <- paste(ENTRY_DIR,"/output/",opt$write,".csv.aa.ungapped.fa", sep='')
msa <- paste(ENTRY_DIR,"/output/",opt$write,".csv.aa.ungapped.headers.fa", sep='')

#trimal to trim alignment to an ungapped version -- this allows a subtree to become a comprehensible alignment for the MSA visualization
system(paste("trimAl -in ",aa," -out ",msa_pre," -noallgaps",sep=""))
system(paste("python scripts/remove_header.py ",msa_pre," ",msa,sep=""))


#Prints 2 PDFs: data as text, and multiple sequence alignment cartoon

pdf(file, height=opt$height, width=opt$width)
p
dev.off()

#Optional -- read a datafile for a heatmap

#Set visualization parameters for heatmaps
##heatmap colors and limits
#upper <- 7 #these are the upper and lower limits used to set colors. If outside the limits cell is gray
#lower <- -7
#low_color<-"#000099" #blue
#high_color<-"#FF0000" #red
#high_color<-"#ffa500" #orange
#high_color<-"#FFCC33" #yellow

#print("Reading heatmap data file")
#counts_file <- read.table("datasets/In11_log2FC.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactor=F)
#p <- p %<+% counts_file

#pdf(paste(file,".heatmap.pdf",sep=""), height=opt$height, width=opt$width)
#gheatmap(p2,counts_file, offset = opt$heatmap_offset, width=opt$heatmap_width+3, font.size=size, colnames_angle=-20, hjust=0, color="black") + #scale_fill_gradient2(low=low_color,high=high_color,mid="white",limits=c(lower,upper)) 
#dev.off()


### Version 2 of the tree will have an MSA cartoon

p_msa<-msaplot(p2,msa,offset=opt$heatmap_offset,color=msa_colors,bg_line=FALSE,) + guides(fill = "none")

if (has_features) {

# --- Add feature rectangles over the existing MSA (drop-in; UNALIGNED -> ALIGNED) ----
# This version reads unaligned positions and maps them to aligned MSA columns.
# Requirements: Biostrings, ggnewscale (installed on demand below).

feature_file <- file.path(ENTRY_DIR, "output", "features.txt")

# 0) Read features (unaligned coordinates)
feat <- utils::read.delim(feature_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE, quote = "")
stopifnot(all(c("label","aa_start","aa_end","feature") %in% names(feat)))
feat$aa_start <- as.integer(feat$aa_start)
feat$aa_end   <- as.integer(feat$aa_end)
feat <- feat[!is.na(feat$label) & !is.na(feat$aa_start) & !is.na(feat$aa_end) & feat$aa_end >= feat$aa_start, , drop = FALSE]

# 1) Build alignment mapping: for each sequence, map unaligned AA index -> aligned column
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
alnAA <- Biostrings::readAAStringSet(msa)    # 'msa' path defined above
lab_order <- names(alnAA)

build_map_list <- function(aa_set) {
  ml <- vector("list", length(aa_set))
  names(ml) <- names(aa_set)
  for (nm in names(aa_set)) {
    s <- as.character(aa_set[[nm]])
    ung <- 0L
    # m[unaligned_index] = aligned_column_index
    m <- integer(0)
    for (i in seq_len(nchar(s))) {
      ch <- substr(s, i, i)
      if (ch != "-") {
        ung <- ung + 1L
        m[ung] <- i
      }
    }
    ml[[nm]] <- m
  }
  ml
}
map_list <- build_map_list(alnAA)

# Helpers to convert unaligned positions to aligned columns, clamping inside seq length
map_pos_start <- function(label, pos) {
  m <- map_list[[label]]
  if (is.null(m) || length(m) == 0L) return(NA_integer_)
  if (pos < 1L) pos <- 1L
  if (pos > length(m)) pos <- length(m)
  m[pos]
}
map_pos_end <- function(label, pos) {
  m <- map_list[[label]]
  if (is.null(m) || length(m) == 0L) return(NA_integer_)
  if (pos < 1L) pos <- 1L
  if (pos > length(m)) pos <- length(m)
  m[pos]
}

# 2) Convert feature starts/ends (unaligned) -> MSA columns (aligned)
feat$aa_start_col <- mapply(map_pos_start, feat$label, feat$aa_start)
feat$aa_end_col   <- mapply(map_pos_end,   feat$label, feat$aa_end)
feat <- feat[!is.na(feat$aa_start_col) & !is.na(feat$aa_end_col), , drop = FALSE]

# ---- Export aligned coordinates for provenance --------------------------------
aligned_out <- file.path(ENTRY_DIR, "output", paste(opt$write, "_features_aligned.txt", sep=" "))
export_cols <- c("label", "feature",
                 "aa_start", "aa_end",           # original unaligned coords
                 "aa_start_col", "aa_end_col")   # aligned MSA columns

# Add a length column (aligned span) for convenience
feat$len_aa_aligned <- feat$aa_end_col - feat$aa_start_col + 1L
export_cols <- c(export_cols, "len_aa_aligned")

# Ensure output dir exists
dir.create(dirname(aligned_out), recursive = TRUE, showWarnings = FALSE)

# Write tab-delimited, no quotes
utils::write.table(
  feat[, export_cols, drop = FALSE],
  file      = aligned_out,
  sep       = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote     = FALSE
)

message(sprintf("[features] wrote aligned coordinates → %s", aligned_out))
# ------------------------------------------------------------------------------

# 3) Extract the per-column x-bounds from p_msa’s drawn grid
pb <- ggplot2::ggplot_build(p_msa)

normalize_rects <- function(d) {
  if (!("PANEL" %in% names(d))) d$PANEL <- 1L
  if (all(c("xmin","xmax","ymin","ymax") %in% names(d))) {
    d[, c("xmin","xmax","ymin","ymax","PANEL"), drop = FALSE]
  } else if (all(c("x","y","width","height") %in% names(d))) {
    data.frame(
      xmin  = d$x - d$width/2,
      xmax  = d$x + d$width/2,
      ymin  = d$y - d$height/2,
      ymax  = d$y + d$height/2,
      PANEL = d$PANEL
    )
  } else if (all(c("x","y") %in% names(d))) {
    ux <- sort(unique(d$x)); uy <- sort(unique(d$y))
    wx <- if (length(ux) > 1) min(diff(ux)) else 1
    hy <- if (length(uy) > 1) min(diff(uy)) else 1
    data.frame(
      xmin  = d$x - wx/2,
      xmax  = d$x + wx/2,
      ymin  = d$y - hy/2,
      ymax  = d$y + hy/2,
      PANEL = d$PANEL
    )
  } else NULL
}
rect_list <- Filter(Negate(is.null), lapply(pb$data, normalize_rects))
stopifnot(length(rect_list) > 0)
msa_rect <- rect_list[[ which.max(vapply(rect_list, nrow, integer(1))) ]]

col_center <- (msa_rect$xmin + msa_rect$xmax) / 2
u_centers  <- sort(unique(round(col_center, 6)))
col_id     <- match(round(col_center, 6), u_centers)  # 1..N columns

xleft_by_col  <- tapply(msa_rect$xmin, col_id, min, na.rm = TRUE)
xright_by_col <- tapply(msa_rect$xmax, col_id, max, na.rm = TRUE)

col_edges <- data.frame(
  aa_col = seq_along(u_centers),
  xleft  = as.numeric(xleft_by_col),
  xright = as.numeric(xright_by_col)
)

# 4) Map tip labels to y positions from p_msa
tip_df <- p_msa$data
if ("isTip" %in% names(tip_df)) {
  tip_y <- tip_df[tip_df$isTip %in% TRUE & !is.na(tip_df$label), c("label","y")]
} else {
  tmp <- tip_df[!is.na(tip_df$label), c("label","y")]
  tip_y <- tmp[!duplicated(tmp$label), , drop = FALSE]
}

# 5) Join: features -> y, then to column edges using aligned columns
feat_pos <- merge(feat, tip_y, by = "label", all.x = TRUE)
feat_pos <- merge(feat_pos, setNames(col_edges, c("aa_start_col","xleft",".z1")), by = "aa_start_col", all.x = TRUE)
feat_pos <- merge(feat_pos, setNames(col_edges, c("aa_end_col",".z2","xright")), by = "aa_end_col", all.x = TRUE)
feat_pos <- feat_pos[!is.na(feat_pos$xleft) & !is.na(feat_pos$xright) & !is.na(feat_pos$y), , drop = FALSE]

# 6) Compute vertical extents from drawn grid
row_h <- stats::median(msa_rect$ymax - msa_rect$ymin, na.rm = TRUE)
feat_pos$ymin <- feat_pos$y - row_h/2
feat_pos$ymax <- feat_pos$y + row_h/2

# 7) Draw: rectangles for long spans, midpoint dots for short (<10 aa)
if (!requireNamespace("ggnewscale", quietly = TRUE)) utils::install.packages("ggnewscale")
p_msa <- p_msa + ggnewscale::new_scale_fill()

feat_pos$len_aa <- feat_pos$aa_end - feat_pos$aa_start + 1L
long_pos  <- feat_pos[feat_pos$len_aa >= 10, , drop = FALSE]
short_pos <- feat_pos[feat_pos$len_aa <  10, , drop = FALSE]
if (nrow(short_pos) > 0) {
  short_pos$x_mid <- (short_pos$xleft + short_pos$xright)/2
  short_pos$y_mid <- short_pos$y
}

feature_colors <- c(
  "#FF0000", # bright red
  "#FF7F00", # vivid orange
  "#00FF00", # neon green
  "#00CED1", # bright turquoise
  "#00BFFF", # vivid sky blue
  "#FFD700", # bright yellow-gold
  "#1E90FF", # dodger blue
  "#0000FF", # pure bright blue
  "#8A2BE2", # vivid blue violet
  "#FF00FF", # magenta
  "#FF1493", # deep pink
  "#FF69B4", # hot pink
  "#FF4500", # orange red
  "#ADFF2F", # green yellow
  "#7CFC00", # lawn green
  "#40E0D0", # turquoise
  "#00FA9A", # medium spring green
  "#DA70D6", # orchid
  "#FF6347", # tomato
  "#FFFF00"  # pure yellow
)
levs <- unique(feat_pos$feature)
fills <- setNames(rep(feature_colors, length.out = length(levs)), levs)

if (nrow(long_pos) > 0) {
  p_msa <- p_msa +
    ggplot2::geom_rect(
      data = long_pos,
      ggplot2::aes(xmin = xleft, xmax = xright, ymin = ymin, ymax = ymax, fill = feature),
      inherit.aes = FALSE,
      alpha = 0.55,
      color = "black",
      linewidth = 0.2
    )
}
if (nrow(short_pos) > 0) {
  p_msa <- p_msa +
    ggplot2::geom_point(
      data = short_pos,
      ggplot2::aes(x = x_mid, y = y_mid, fill = feature),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.2,
      stroke = 0.3,
      color = "black"
    )
}
p_msa <- p_msa + ggplot2::scale_fill_manual(values = fills, drop = FALSE, name = "Feature")
# -------------------------------------------------------------------------------------

} else {
  # Do nothing. Keep the rest of the script unchanged.
}

pdf(paste(file,".MSA.pdf",sep=""), height=opt$height, width=opt$width)
print("Using this MSA sequence")
print(msa)
p_msa
dev.off()

### ------------------- Version 3: nearest-by-genome distances -------------------
### Requires: tree (ape::phylo), dd with columns 'taxa' and 'genome',
###           p2 (your base tree with labels), msaplot inputs already prepared
### Output:   <entry>/output/<write>.nearest.pdf

# Helper: compute a full pairwise distance matrix
compute_distance_matrix <- function(dist_type = "patristic", tree, msa_path) {
  dist_type <- tolower(dist_type)
  tips <- tree$tip.label

  if (dist_type == "patristic") {
    D <- ape::cophenetic.phylo(tree)
    # ensure full square over current tips
    D <- D[tips, tips, drop = FALSE]
    return(D)
  }

  # phenetic = uncorrected p-distance on your trimmed ungapped AA alignment
  if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
  aa <- Biostrings::readAAStringSet(msa_path)
  seqs <- setNames(as.character(aa), names(aa))

  common <- intersect(tips, names(seqs))
  D <- matrix(NA_real_, length(tips), length(tips), dimnames = list(tips, tips))

  # simple p-distance over equal-length ungapped AA strings
  get_pd <- function(a, b) {
    a <- strsplit(a, "", fixed = TRUE)[[1]]
    b <- strsplit(b, "", fixed = TRUE)[[1]]
    L <- min(length(a), length(b))
    if (L == 0) return(NA_real_)
    mean(a[seq_len(L)] != b[seq_len(L)])
  }

  for (i in seq_along(common)) {
    for (j in i:length(common)) {
      li <- common[i]; lj <- common[j]
      d  <- get_pd(seqs[[li]], seqs[[lj]])
      D[li, lj] <- D[lj, li] <- d
    }
  }
  D
}

# Helper: build a wide table: first column 'label', then one column per genome
nearest_by_genome_table <- function(tree, dd, D, round_digits = 3) {
  tips <- tree$tip.label
  dd2  <- dd[dd$taxa %in% tips, c("taxa","genome")]
  rownames(dd2) <- NULL
  genomes <- sort(unique(dd2$genome))
  by_gen  <- split(dd2$taxa, dd2$genome)

  # prepare string matrix so printing looks clean
  out_mat <- matrix("", nrow = length(tips), ncol = length(genomes),
                    dimnames = list(tips, genomes))

  for (lab in tips) {
    g_self <- dd2$genome[match(lab, dd2$taxa)]
    for (g in genomes) {
      if (identical(g, g_self)) {
        out_mat[lab, g] <- ""               # blank for own genome
      } else {
        candidates <- intersect(by_gen[[g]], rownames(D))
        if (length(candidates)) {
          val <- suppressWarnings(min(D[lab, candidates], na.rm = TRUE))
          if (is.finite(val)) {
            out_mat[lab, g] <- formatC(val, digits = round_digits, format = "f")
          } else {
            out_mat[lab, g] <- ""
          }
        } else {
          out_mat[lab, g] <- ""
        }
      }
    }
  }

  data.frame(label = tips, as.data.frame(out_mat, check.names = FALSE), check.names = FALSE)
}

# 1) Build the chosen distance matrix
Dmat <- compute_distance_matrix(opt$dist_type, tree, msa)

# 2) Build the dataset (first column MUST be 'label' for %<+%)
nearest_df <- nearest_by_genome_table(tree, dd, Dmat, round_digits = opt$dist_digits)

# Reorder rows to match the plotted tree's tip order (top-to-bottom)
tip_order <- tryCatch({
  dfp <- p2$data
  dfp <- dfp[!is.na(dfp$label) & dfp$isTip %in% TRUE, c("label","y")]
  dfp <- dfp[order(dfp$y, decreasing = TRUE), , drop = FALSE]
  unique(as.character(dfp$label))
}, error = function(e) tree$tip.label)

nearest_df <- nearest_df[match(tip_order, nearest_df$label), , drop = FALSE]
row.names(nearest_df) <- NULL

# 2.5) Also write the nearest-by-genome table to CSV (wide format)
nearest_csv <- file.path(ENTRY_DIR, "output",
                         paste0(opt$write, ".nearest_", tolower(opt$dist_type), ".csv"))
dir.create(dirname(nearest_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  nearest_df,
  file      = nearest_csv,
  sep       = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote     = TRUE,
  na        = ""
)
message(sprintf("[nearest] wrote table → %s", nearest_csv))

# 3) Add to a clean copy of your base plot and render
datasets_nearest <- list()
datasets_nearest[[paste0("Nearest (", tolower(opt$dist_type), ")")]] <- nearest_df

p_nearest <- add_datasets_to_tree(p2, datasets_nearest, base_offset = 0)

width_nearest <- max(p_nearest$data$x) + total_offset
pdf(paste0(file, ".nearest.pdf"), height = opt$height, width = width_nearest)
p_nearest
dev.off()
### ---------------------------------------------------------------------------
