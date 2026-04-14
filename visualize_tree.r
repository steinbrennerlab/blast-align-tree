# Adam Steinbrenner
# astein10@uw.edu
# UW Seattle, Dept of Biology

# Required packages: ggplot2, treeio, phytools, ggtree, optparse, tidytree, ape, broom
# BiocManager packages: treeio, ggtree, Biostrings
suppressPackageStartupMessages({
  library("ggplot2")
  library("treeio")
  library("phytools")
  library("ggtree")
  library("optparse")
  library("tidytree")
  library("ape")
  library("broom")
})

# Determine the directory where this R script is located (for robust script paths)
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("--file=", "", file_arg))
    return(dirname(script_path))
  }
  return(getwd())
}
SCRIPT_BASE <- get_script_dir()

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
	make_option(c("--heatmap_font_size"), action="store", default=1.2, type="numeric",
        help="font size for heatmap column labels"),
	make_option(c("--heatmap_label_angle"), action="store", default=45, type="numeric",
        help="angle for heatmap column labels"),
	make_option(c("-d", "--datasets"), action="store", type="character", default=NULL,
        help="Comma- or semicolon-separated dataset files/directories to display and heatmap. Supports .txt, .tsv, and .csv. Directories are searched recursively for supported files. If omitted, top-level files in datasets/ are used."),
	make_option(c("--heatmap_file"), action="store", type="character", default=NULL,
        help="Deprecated alias for one dataset file. Used only when --datasets is omitted."),
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
		help = "Optional path to features.txt. If omitted, will try <entry>/features.txt. If missing/unreadable, feature overlay is skipped."),
	make_option(c("--subdir"),type = "character",default = NULL,
		help = "Relative subdirectory to append under both --entry and --write (e.g., 'runs/20250903_1557'). Leading slashes will be stripped."),
	make_option(c("--dist_type"), action="store", type="character", default="patristic",
		help="Distance type to compute: 'patristic' or 'phenetic' (case-insensitive)"),
	make_option(c("--dist_digits"), action="store", type="integer", default=3, help="Number of digits to round distances when printed"),
	make_option(c("--nearest_from"), type = "character", default = NULL,
		help = "Comma- or semicolon-separated list of genomes to include in the nearest-distance panel. Use exact names as in merged_genome_mapping.txt (column 'genome'). If omitted, all genomes (except a tip's own) are considered."),
make_option(c("--genome_colors"), action="store", type="character", default=NULL,
  help="Optional genome-specific color overrides. Format: 'GenomeA=red,GenomeB=#00FF00'")

	)
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

parse_csv_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  unique(trimws(unlist(strsplit(x, "[,;]"))))
}

parse_genome_colors <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)

  pairs <- trimws(unlist(strsplit(x, ",")))
  res <- lapply(pairs, function(p) {
    kv <- trimws(unlist(strsplit(p, "=")))
    if (length(kv) != 2) return(NULL)
    data.frame(genome = kv[1], color = kv[2], stringsAsFactors = FALSE)
  })

  do.call(rbind, Filter(Negate(is.null), res))
}

###
# Check for presence of features.txt in the run folder
###

`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

feature_file_default <- file.path(ENTRY_DIR, "features.txt")
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
tree <- read.tree(tree_newick)

#Takes in merged_genome_mapping file from blast_align_tree
gene_species_list <- paste(ENTRY_DIR,"/","merged_genome_mapping.txt", sep='')

#define output filenames
file <- file.path(ENTRY_DIR, paste0(opt$write, ".pdf"))
file_csv <- file.path(ENTRY_DIR, paste0(opt$write, ".csv"))
file_nwk <- file.path(ENTRY_DIR, paste0(opt$write, ".nwk"))



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

p <- p + theme(legend.title = element_text(size = 6),
               legend.text = element_text(size = 5),
               legend.key.size = unit(0.3, "cm"),
               legend.margin = margin(2, 2, 2, 2))

#add species list to the tree object; this is the "merged coding" file from the blast-align-tree bash script, a tab separated text file of taxa and genome
dd <- read.table(gene_species_list, sep="\t", header = TRUE, stringsAsFactor=F)
p <- p %<+% dd

#How many tips does the tree have
a <- as.integer(length(tree$tip.label))

#set height dimension for output based on number of tips
node_count <- length(tree$tip.label)
opt$height <- ((node_count/11.27)+0.5)
message(sprintf("[R] Tree loaded: %d tips", node_count))

#Define tip and node label sizes
size <- opt$size
size2 <- (size)

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

# Extract sequences in tree order: NT (for CDS) and AA (for alignment visualization)
system2("python", shQuote(c(file.path(SCRIPT_BASE, "scripts/extract_seq.py"), ENTRY_DIR, opt$entry, paste0(opt$write, ".csv"))), stdout = FALSE, stderr = FALSE)
system2("python", shQuote(c(file.path(SCRIPT_BASE, "scripts/extract_seq.py"), ENTRY_DIR, opt$entry, paste0(opt$write, ".csv"), "--aa")), stdout = FALSE, stderr = FALSE)

# Read gene symbols to display next to annotated genes
  ## Read gene symbols file from root directory
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

# ---- Genome color overrides ----
override_df <- parse_genome_colors(opt$genome_colors)

genomes_present <- sort(unique(dd$genome))

# Default mapping (existing behavior)
default_map <- setNames(standard_colors_generator(length(genomes_present)),
                        genomes_present)

if (!is.null(override_df)) {

  override_df <- override_df[override_df$genome %in% genomes_present, , drop = FALSE]

  if (nrow(override_df)) {
    message(sprintf("[colors] applying overrides: %s",
      paste(paste(override_df$genome, override_df$color, sep="="), collapse=", ")))

    default_map[override_df$genome] <- override_df$color
  } else {
    message("[colors] overrides provided but no matching genomes found")
  }
}

final_color_map <- default_map

  ## Add gene symbol labels to the tree using offsets specified in the options
p <- p +
	geom_tiplab(size=size,offset=opt$label_offset,aes(color=genome,fontface="bold")) +
	#tip labels (gene names) colored by species
	geom_tiplab(aes(label=symbol,color=genome), size=opt$symbol_size,align=T, linetype=NA, offset=opt$symbol_offset) +
	#gene symbol names
	scale_colour_manual(values = final_color_map, guide = guide_legend(ncol = 1, keywidth = 0.3, keyheight = 0.3))

if (opt$genomeLabel_boolean > 0) {
	p <- p +
	geom_tiplab(aes(label=genome,color=genome), size=opt$symbol_size,align=T, linetype=NA, offset=opt$symbol_offset+1.3) +
	#genome label names
	scale_colour_manual(values = final_color_map)
}

#converts "labels" bootstrap to a 1-100 integer
d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(as.character(d$label))
d$label <- ifelse(d$label <= 1, as.integer(d$label * 100), as.integer(d$label))
#option to only show some threshold bootstrap
#d <- d[d$label < 75,]

#If option -k is specified, include bootstraps
if (opt$bootstrap_boolean > 0) {
	p <- p + geom_text(data=d, aes(label=label), size=size, nudge_x=opt$bootstrap_offset)  #bootstraps
}


# Read datasets selected by --datasets/--heatmap_file.
## Files must contain one ID column followed by one or more data columns. The ID
## column is normalized to "taxa" so ggtree can join values to tip labels.

default_dataset_dir <- "datasets"
supported_dataset_pattern <- "\\.(txt|tsv|csv)$"

split_path_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  parts <- trimws(unlist(strsplit(x, "[,;]")))
  parts[nzchar(parts)]
}

resolve_existing_path <- function(path) {
  candidates <- c(path, file.path(SCRIPT_BASE, path))
  candidates <- candidates[!duplicated(candidates)]
  for (candidate in candidates) {
    if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }
  stop(sprintf("Dataset path not found: %s", path), call. = FALSE)
}

expand_dataset_inputs <- function(inputs) {
  files <- character(0)

  for (input in inputs) {
    expanded <- Sys.glob(input)
    if (length(expanded) == 0) expanded <- Sys.glob(file.path(SCRIPT_BASE, input))
    if (length(expanded) == 0) expanded <- input

    for (item in expanded) {
      path <- resolve_existing_path(item)

      if (dir.exists(path)) {
        found <- list.files(
          path,
          pattern = supported_dataset_pattern,
          full.names = TRUE,
          recursive = TRUE,
          ignore.case = TRUE
        )
        files <- c(files, found)
      } else if (grepl(supported_dataset_pattern, path, ignore.case = TRUE)) {
        files <- c(files, path)
      } else {
        stop(sprintf("Unsupported dataset extension for '%s'. Use .txt, .tsv, or .csv.", path), call. = FALSE)
      }
    }
  }

  unique(normalizePath(files, winslash = "/", mustWork = TRUE))
}

default_dataset_files <- function() {
  dataset_dir <- resolve_existing_path(default_dataset_dir)
  files <- list.files(
    dataset_dir,
    pattern = supported_dataset_pattern,
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  unique(normalizePath(files, winslash = "/", mustWork = TRUE))
}

read_dataset_file <- function(path, dataset_name) {
  ext <- tolower(tools::file_ext(path))

  df <- tryCatch(
    if (ext == "csv") {
      utils::read.csv(
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = ""
      )
    } else {
      utils::read.delim(
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "",
        comment.char = ""
      )
    },
    error = function(e) {
      stop(sprintf("Could not read dataset '%s': %s", path, e$message), call. = FALSE)
    }
  )

  if (ncol(df) < 2) {
    stop(sprintf("Dataset '%s' must contain an ID column plus at least one data column.", path), call. = FALSE)
  }

  names(df)[1] <- "taxa"
  df$taxa <- trimws(as.character(df$taxa))
  df <- df[!is.na(df$taxa) & nzchar(df$taxa), , drop = FALSE]

  if (nrow(df) == 0) {
    stop(sprintf("Dataset '%s' has no usable taxa IDs in the first column.", path), call. = FALSE)
  }

  display_cols <- trimws(names(df)[-1])
  empty_cols <- !nzchar(display_cols)
  if (any(empty_cols)) {
    display_cols[empty_cols] <- paste0("column_", which(empty_cols))
  }

  safe_prefix <- make.names(dataset_name)
  internal_cols <- make.unique(make.names(paste(safe_prefix, display_cols, sep = "__")), sep = "_")
  names(df)[-1] <- internal_cols
  raw_df <- df

  collapse_values <- function(x) {
    values <- trimws(as.character(x))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 0) return("")
    paste(unique(values), collapse = ";")
  }

  display_df <- raw_df
  for (col in internal_cols) {
    display_df[[col]] <- ifelse(is.na(display_df[[col]]), "", as.character(display_df[[col]]))
  }
  if (anyDuplicated(display_df$taxa)) {
    message(sprintf("[R] Dataset '%s' has duplicate taxa IDs; collapsing duplicate values for direct tree labels.", dataset_name))
    display_df <- stats::aggregate(
      display_df[, internal_cols, drop = FALSE],
      by = list(taxa = display_df$taxa),
      FUN = collapse_values
    )
  }

  list(
    name = dataset_name,
    path = path,
    data = display_df,
    raw_data = raw_df,
    columns = data.frame(
      internal = internal_cols,
      label = display_cols,
      stringsAsFactors = FALSE
    )
  )
}

read_datasets <- function(files) {
  dataset_names <- tools::file_path_sans_ext(basename(files))
  dataset_names <- make.unique(dataset_names, sep = "_")
  datasets <- vector("list", length(files))

  for (i in seq_along(files)) {
    datasets[[i]] <- read_dataset_file(files[i], dataset_names[i])
  }

  names(datasets) <- dataset_names
  datasets
}

## Function to add datasets to the tree, calculating the text offset for each column
add_datasets_to_tree <- function(p, datasets, base_offset) {
  cumulative_offset <- opt$symbol_offset + 0.5

  if (length(datasets) == 0) {
    return(list(plot = p, total_offset = cumulative_offset))
  }

  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    dataset_name <- dataset$name
    dataset_df <- dataset$data

    p <- p %<+% dataset_df

    # Add dataset name as text annotation
    p <- p + annotate("text", size = opt$size, x = max(p$data$x) + base_offset + cumulative_offset,
                      y = max(p$data$y) + 2.4, label = dataset_name, fontface = "bold", hjust = 0)

    # Add columns as tip labels
    columns_to_display <- dataset$columns
    for (j in seq_len(nrow(columns_to_display))) {
      column_name <- columns_to_display$internal[j]
      column_label <- columns_to_display$label[j]
      column_offset <- (j - 1) * 0.2

      # The column data should use the cumulative offset plus its own offset.
      p <- local({
        column_name_local <- column_name
        p + geom_tiplab(aes(label = .data[[column_name_local]]), size = 1, align = TRUE,
                        linetype = NA, offset = base_offset + cumulative_offset + column_offset)
      })

      # Column names should align with their data.
      p <- p + annotate("text", size = opt$size * 0.8,
                        x = max(p$data$x) + base_offset + cumulative_offset + column_offset,
                        y = max(p$data$y) + 0.7,
                        label = column_label,
                        angle = 45,
                        hjust = 0)
    }

    # Update the cumulative offset after processing each dataset.
    cumulative_offset <- cumulative_offset + nrow(columns_to_display) * 0.2 + 0.15
  }

  list(plot = p, total_offset = cumulative_offset)
}


dataset_inputs <- split_path_list(opt$datasets)
if (length(dataset_inputs) > 0 && !is.null(opt$heatmap_file) && nzchar(opt$heatmap_file)) {
  message("[R] --datasets provided; ignoring deprecated --heatmap_file.")
}
if (length(dataset_inputs) == 0) {
  dataset_inputs <- split_path_list(opt$heatmap_file)
}
if (length(dataset_inputs) == 0) {
  dataset_files <- default_dataset_files()
} else {
  dataset_files <- expand_dataset_inputs(dataset_inputs)
}

if (length(dataset_files) == 0) {
  stop("No supported dataset files found. Use .txt, .tsv, or .csv files.", call. = FALSE)
}
message(sprintf("[R] Datasets selected: %s", paste(basename(dataset_files), collapse = ", ")))
datasets <- read_datasets(dataset_files)

  ## Create a duplicate tree for heatmap and MSA trees below -- it will not have datasets
p2 <- p

# Add datasets to the tree
dataset_tree <- add_datasets_to_tree(p, datasets, base_offset = 0)
p <- dataset_tree$plot
total_offset <- dataset_tree$total_offset

heatmap_width <- opt$width
opt$width <- max(opt$width, max(p$data$x) + total_offset)

#Creates a visualization of the multiple sequence alignment with any amino acid indicated as black, and any gap indicated as grey
##msa_colors <- c("gray85","red","orange","green",rep(c("black"),each=20)) ##version of color scale for domains
msa_colors <- c("gray85",rep(c("black"),each=30))
aa <- file.path(ENTRY_DIR, paste0(opt$write, ".csv.aa.fa"))
msa_pre <- file.path(ENTRY_DIR, paste0(opt$write, ".csv.aa.no_all_gap_columns.fa"))
msa <- file.path(ENTRY_DIR, paste0(opt$write, ".csv.aa.no_all_gap_columns.ids_only.fa"))

#trimal removes columns where every sequence is a gap, leaving within-sequence gaps intact for the MSA visualization
system2("trimal", shQuote(c("-in", aa, "-out", msa_pre, "-noallgaps")), stdout = FALSE, stderr = FALSE)
# Keep only FASTA IDs in the cleaned alignment so tree labels and MSA labels match exactly
system2("python", shQuote(c(file.path(SCRIPT_BASE, "scripts/remove_header.py"), msa_pre, msa)), stdout = FALSE, stderr = FALSE)

# Ensure MSA and tree tips match exactly (avoids msaplot crash from mismatches)
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
msa_seqs <- Biostrings::readAAStringSet(msa)

# Remove MSA sequences not in the tree
msa_keep <- names(msa_seqs) %in% tree$tip.label
if (sum(msa_keep) < length(msa_seqs)) {
  message(sprintf("[R] Filtering MSA: keeping %d of %d sequences to match tree tips",
                  sum(msa_keep), length(msa_seqs)))
  msa_seqs <- msa_seqs[msa_keep]
}

# Add gap-only sequences for tree tips missing from the MSA
missing_tips <- setdiff(tree$tip.label, names(msa_seqs))
if (length(missing_tips) > 0) {
  aln_width <- unique(Biostrings::width(msa_seqs))[1]
  gap_str <- paste(rep("-", aln_width), collapse = "")
  gap_seqs <- Biostrings::AAStringSet(setNames(rep(gap_str, length(missing_tips)), missing_tips))
  msa_seqs <- c(msa_seqs, gap_seqs)
  message(sprintf("[R] Added %d gap-only sequences for tree tips missing from MSA",
                  length(missing_tips)))
}

Biostrings::writeXStringSet(msa_seqs, msa)


###
#   Tree Version 1 -- Display all data in /datasets. txt files should contain column "taxa" followed by relevant columns
###

message(sprintf("[R] Writing tree PDF: %s", file))
pdf(file, height=opt$height, width=opt$width)
print(p)
dev.off()

###
#   Tree Version 2 -- Display a specific file as a heatmap
###

#Set visualization parameters for heatmaps
##heatmap colors and limits
# The scale type is auto-detected from the data:
#   - If any negative values exist   -> diverging scale (blue-white-red, symmetric around 0)
#                                       suits fold-change data (e.g., Bjornsen PAMP L2FC)
#   - If all values are non-negative -> sequential scale (white-red, 0 to max)
#                                       suits absolute expression (e.g., Klepikova log2 counts)
# No user configuration required. Pass any supported BAT-format dataset files and the heatmap will adapt.

coerce_numeric_for_heatmap <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))

  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "N/A", "NaN", "na", "n/a", "null", "NULL",
                     "#N/A", "#NUM!", "#DIV/0!", "#VALUE!", "#REF!")] <- NA_character_
  converted <- suppressWarnings(as.numeric(x_chr))

  non_missing <- !is.na(x_chr)
  if (!any(!is.na(converted))) return(NULL)
  if (any(non_missing & is.na(converted))) return(NULL)

  converted
}

mean_or_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

build_heatmap_table <- function(datasets) {
  tables <- list()
  skipped <- character(0)

  for (dataset in datasets) {
    df <- dataset$raw_data
    numeric_df <- data.frame(taxa = df$taxa, stringsAsFactors = FALSE, check.names = FALSE)
    heatmap_col_names <- character(0)

    for (j in seq_len(nrow(dataset$columns))) {
      internal <- dataset$columns$internal[j]
      display <- dataset$columns$label[j]
      numeric_values <- coerce_numeric_for_heatmap(df[[internal]])

      if (is.null(numeric_values)) {
        skipped <- c(skipped, sprintf("%s:%s", dataset$name, display))
        next
      }

      heatmap_name <- if (length(datasets) > 1) {
        paste(dataset$name, display, sep = ": ")
      } else {
        display
      }
      heatmap_col_names <- c(heatmap_col_names, heatmap_name)
      numeric_df[[heatmap_name]] <- numeric_values
    }

    if (length(heatmap_col_names) == 0) next

    if (anyDuplicated(numeric_df$taxa)) {
      message(sprintf("[R] Heatmap dataset '%s' has duplicate taxa IDs; averaging duplicate numeric values.", dataset$name))
    }

    agg <- stats::aggregate(
      numeric_df[, heatmap_col_names, drop = FALSE],
      by = list(taxa = numeric_df$taxa),
      FUN = mean_or_na
    )
    rownames(agg) <- agg$taxa
    agg$taxa <- NULL
    names(agg) <- make.unique(names(agg), sep = "_")
    tables[[dataset$name]] <- agg
  }

  if (length(skipped) > 0) {
    message(sprintf("[R] Heatmap skipped nonnumeric columns: %s", paste(skipped, collapse = ", ")))
  }
  if (length(tables) == 0) return(NULL)

  all_taxa <- unique(unlist(lapply(tables, rownames), use.names = FALSE))
  combined <- data.frame(row.names = all_taxa, check.names = FALSE)

  for (tbl in tables) {
    for (col in names(tbl)) {
      combined_col <- col
      while (combined_col %in% names(combined)) {
        combined_col <- make.unique(c(names(combined), col), sep = "_")[length(names(combined)) + 1]
      }
      combined[[combined_col]] <- NA_real_
      combined[rownames(tbl), combined_col] <- tbl[[col]]
    }
  }

  combined
}

heatmap_file <- build_heatmap_table(datasets)

if (!is.null(heatmap_file) && ncol(heatmap_file) > 0) {
  # Inspect data to choose scale type
  .numeric_vals <- suppressWarnings(as.numeric(unlist(heatmap_file, use.names = FALSE)))
  .numeric_vals <- .numeric_vals[is.finite(.numeric_vals)]
  data_min <- if (length(.numeric_vals)) min(.numeric_vals) else 0
  data_max <- if (length(.numeric_vals)) max(.numeric_vals) else 1
  is_diverging <- data_min < 0

  if (is_diverging) {
    # Values crossing 0 get a diverging scale centered at 0.
    limit <- max(abs(data_min), abs(data_max))
    heatmap_scale <- scale_fill_gradient2(
      low = "#000099", mid = "white", high = "#FF0000",
      midpoint = 0, limits = c(-limit, limit),
      na.value = "grey85"
    )
    scale_label <- "value"
    message(sprintf("[R] Heatmap auto-detect: diverging scale, limits +/-%.2f (data range %.2f to %.2f)",
                    limit, data_min, data_max))
  } else {
    # Non-negative data: render 0 as missing/undetected so it is visually distinct.
    heatmap_file[heatmap_file == 0] <- NA
    sequential_max <- if (data_max > 0) data_max else 1
    heatmap_scale <- scale_fill_gradient(
      low = "#FFFFFF", high = "#FF0000",
      limits = c(0, sequential_max),
      na.value = "grey85"
    )
    scale_label <- "value"
    message(sprintf("[R] Heatmap auto-detect: sequential scale, limits 0 to %.2f", sequential_max))
  }

  p_heatmap <- p2 + labs(fill = scale_label, size=1) + guides(
    fill = guide_colorbar(
      title.theme  = element_text(size = 6),
      label.theme  = element_text(size = 6),
      barwidth     = unit(0.3, "cm"),
      barheight    = unit(1, "cm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  )

  p_heatmap <- p_heatmap +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10))

  heatmap_width <- max(
    heatmap_width,
    max(p2$data$x, na.rm = TRUE) + opt$heatmap_offset + ncol(heatmap_file) * 0.15 + 2
  )

  message(sprintf("[R] Writing heatmap PDF: %s.heatmapL2Counts.pdf", file))
  pdf(paste(file,".heatmapL2Counts.pdf",sep=""), height=opt$height, width=heatmap_width)

  heatmap_plot <- gheatmap(
    p_heatmap,
    heatmap_file,
    offset = opt$heatmap_offset,
    width = opt$heatmap_width,
    colnames_offset_y = 0.5,
    colnames_angle = opt$heatmap_label_angle,
    font.size = opt$heatmap_font_size,
    hjust = 0,
    color="black",
    colnames_position = "top"
  ) + heatmap_scale
  print(heatmap_plot)

  dev.off()
} else {
  message("[R] No numeric dataset columns available for the heatmap tree. Skipping heatmap PDF.")
}

###
#   Tree Version 3 -- Display a multiple sequence alignment cartoon
###

p_msa<-msaplot(p2,msa,offset=opt$heatmap_offset,color=msa_colors,bg_line=FALSE,) + guides(fill = "none")

if (has_features) {

# --- Add feature rectangles over the existing MSA (drop-in; UNALIGNED -> ALIGNED) ----
# This version reads unaligned positions and maps them to aligned MSA columns.
# Requirements: Biostrings, ggnewscale (installed on demand below).

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
feat$aa_start_col <- as.integer(mapply(map_pos_start, feat$label, feat$aa_start))
feat$aa_end_col   <- as.integer(mapply(map_pos_end,   feat$label, feat$aa_end))
feat <- feat[!is.na(feat$aa_start_col) & !is.na(feat$aa_end_col), , drop = FALSE]

if (nrow(feat) == 0L) {
  message("[features] No features mapped to aligned columns. Skipping feature overlay.")
  has_features <- FALSE
} else {

# ---- Export aligned coordinates for provenance --------------------------------
aligned_out <- file.path(ENTRY_DIR, paste0(opt$write, ".features_aligned.txt"))
export_cols <- c("label", "feature",
                 "aa_start", "aa_end",           # original unaligned coords
                 "aa_start_col", "aa_end_col")   # aligned MSA columns

# Add a length column (aligned span) for convenience
feat$len_aa_aligned <- feat$aa_end_col - feat$aa_start_col + 1L
export_cols <- c(export_cols, "len_aa_aligned")

# Ensure run dir exists
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

# 3) Extract the per-column x-bounds from p_msa's drawn grid
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
#    Use annotate() to avoid ggnewscale/msaplot fill-scale conflicts.

feat_pos$len_aa <- feat_pos$aa_end - feat_pos$aa_start + 1L
long_pos  <- feat_pos[feat_pos$len_aa >= 10, , drop = FALSE]
short_pos <- feat_pos[feat_pos$len_aa <  10, , drop = FALSE]
if (nrow(short_pos) > 0) {
  short_pos$x_mid <- (short_pos$xleft + short_pos$xright)/2
  short_pos$y_mid <- short_pos$y
}

feature_colors <- c(
  "#FF0000", "#FF7F00", "#00FF00", "#00CED1", "#00BFFF",
  "#FFD700", "#1E90FF", "#0000FF", "#8A2BE2", "#FF00FF",
  "#FF1493", "#FF69B4", "#FF4500", "#ADFF2F", "#7CFC00",
  "#40E0D0", "#00FA9A", "#DA70D6", "#FF6347", "#FFFF00"
)
levs <- unique(feat_pos$feature)
fills <- setNames(rep(feature_colors, length.out = length(levs)), levs)

for (i in seq_len(nrow(long_pos))) {
  row <- long_pos[i, ]
  p_msa <- p_msa +
    annotate("rect",
      xmin = row$xleft, xmax = row$xright,
      ymin = row$ymin, ymax = row$ymax,
      fill = fills[row$feature],
      alpha = 0.55, color = "black", linewidth = 0.2)
}
for (i in seq_len(nrow(short_pos))) {
  row <- short_pos[i, ]
  p_msa <- p_msa +
    annotate("point",
      x = row$x_mid, y = row$y_mid,
      fill = fills[row$feature],
      shape = 21, size = 2.2, stroke = 0.3, color = "black")
}
# ---- Manual legend for feature overlays (annotate has no legend) ----
msa_left <- min(col_edges$xleft)
msa_w    <- max(col_edges$xright) - msa_left
legend_y_base <- max(p_msa$data$y, na.rm = TRUE) + 1.8

for (fi in seq_along(levs)) {
  ly <- legend_y_base + (length(levs) - fi) * 1.0
  p_msa <- p_msa +
    annotate("point",
      x = msa_left, y = ly,
      fill = fills[levs[fi]], color = "black",
      shape = 22, size = 3, stroke = 0.4) +
    annotate("text",
      x = msa_left + msa_w * 0.03, y = ly,
      label = levs[fi], hjust = 0, size = 2.5)
}

p_msa <- p_msa +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 25, r = 5, b = 5, l = 5))
# -------------------------------------------------------------------------------------

} # end inner else (nrow(feat) > 0)

} else {
  # Do nothing. Keep the rest of the script unchanged.
}

msa_pdf_height <- if (has_features) opt$height + 0.5 else opt$height
message(sprintf("[R] Writing MSA PDF: %s.MSA.pdf", file))
pdf(paste(file,".MSA.pdf",sep=""), height=msa_pdf_height, width=opt$width)
print(p_msa)
dev.off()

