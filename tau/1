library("DESeq2")

## read the command line argumnts
# 1: file path for sample sheet and count table, 
# 2: output directory to write the file
# 3: output directory base

## before running this scrpit we need to have sample sheets and count tables prepared
## These files should be in csv format. sample sheets include the information about the sample 
## name and the treatment. Count table has the information about the gene expression. columns are genes and rows are samples


args <- commandArgs(trailingOnly = TRUE)

input_file_path = args[1]
output_file_path = args[2]
outBase = args[3]

# set working directory 
setwd(input_file_path)

# Create the name for the output base directory
outName=paste(outBase,"_deseq_results",sep = "")

# Create the location to save the files
output_dir <- file.path(output_file_path,outBase)

# Check if the directory is already there and if it exists stop! 
if (!dir.exists(output_dir)){
dir.create(output_dir)
} else {
    stop("Output directory already exists!")
}




# get the list of samples sheets
# sample sheets need to have the pattern "sammple_sheet" in the file name
sample_sheet <- list.files(path=input_file_path,pattern="*sample_sheet*")

# get the list of count table files
# count table files need to have pattern "count" in the file name
count_files <- list.files(path=input_file_path,pattern="*count*")



for (i in 1:length(sample_sheet))
{
    print (c(count_files[i],sample_sheet[i]))
    sampleInfo <- read.csv(sample_sheet[i],check.names = FALSE)
    rownames(sampleInfo) <- sampleInfo$sample
    sampleInfo$sample <- NULL
    print(c(rownames(sampleInfo)))

    print(nrow(sampleInfo))

    countInfo <- read.csv(count_files[i],check.names = FALSE)
    rownames(countInfo) <- countInfo$index
    countInfo$index <- NULL
    print(c(colnames(countInfo)))
    print(ncol(countInfo))

    dds <- DESeqDataSetFromMatrix(countData = round(countInfo),colData = sampleInfo,design = ~ Treatment)
    # dds <- DESeqDataSetFromMatrix(countData = round(countInfo+1),colData = sampleInfo,design = ~ Treatment)
    # res <- lfcShrink(dds, coef=2, res=res) 
    # resNorm <- lfcShrink(dds, coef=2, type="normal")
    
    dds <- DESeq(dds)

    res <- results(dds, contrast=c("Treatment","T","C"))
    #res_shrink <- lfcShrink(dds, coef=2, res=res)
    resNorm <- lfcShrink(dds, coef=2, type="normal")
    print(c(head(res)))

    
    write.csv(as.data.frame(resNorm),file=paste(sample_sheet[i], "deseq2_res.csv", sep = "."))
    normalized_data <- counts(dds, normalized=TRUE)

    write.csv(as.data.frame(normalized_data),file=paste(sample_sheet[i], "deseq2_norm_counts.csv", sep = "."))
}
