#!/usr/bin/env Rscript
options(warn = -1)

suppressMessages(library(edgeR))
suppressMessages(library(argparser)) #https://github.com/cran/argparser


# 检查R包
allkgs <- installed.packages()[,"Package"]
need <- c("argparser","edgeR")
nedd_install <- c("conda install -c conda-forge r-argparser","conda install -c bioconda bioconductor-edger")
n=0
for (x in 1:length(need)) {
  if(!need[x] %in% allkgs){
    message("Please run: ",nedd_install[x])
    n=n+1
  }
}
if(n != 0)stop("Please install the R packages first!")
# conda install -c conda-forge r-argparser
# conda install -c bioconda bioconductor-edger

# 参数设置
p <- argparser::arg_parser("Merge featureCounts(linux version) and calculate FPKM/TPM")

p <- argparser::add_argument(p, "--input_path", help="input: a directory containing the counts matrix. Recommended named: '<sample>.count'", type="character",default = "./")
p <- argparser::add_argument(p, "--feature_length", help="input: the file of gene length. The file has two column,<Feautreid> and <Length>. 'default' for Length of result in featureCounts!", type="character",default = "default")
p <- argparser::add_argument(p, "--pattern", help="limit pattern of input files using regular expression in R language",type="character",default = "*count$")
p <- argparser::add_argument(p, "--output_path", help="output: an existent directory", type="character",default = "./")
p <- argparser::add_argument(p, "--prefix", help="give the file of output matrix a prefix like '<output_prefix>_.*'", type="character",default = "result_gene",short = "-f")
p <- argparser::add_argument(p, "--feature_name", help="feature name of output.", type="character",default = "gene_id")

# 参数解析与定义
argv <- argparser::parse_args(p)

input_path <- argv$input_path
feature_length <- argv$feature_length
pattern <- argv$pattern
output_path <- argv$output_path
output_prefix <- argv$prefix
feature_name <- argv$feature_name

# 测试用

if(F){
  input_path <- "./Quant_featureCounts/known"
  feature_length <- "./GRCh38.110.trans_length.txt" # "default" 或者ID文件
  pattern <- "*\\.count$"
  output_path <- "../"
  output_prefix <- "result_transcript"
  feature_name <- "transcript_id"
}

# 代码主体
file_name <- dir(path = input_path,pattern = pattern, full.names = T)
if(length(file_name) == 0){stop("\nCould not find files with pattern of '",pattern,"'!")}
# file <- file_name
# merge all count matrix
if(feature_length != "default"){
  if(file.exists(feature_length)){
    gene_length <- read.table(feature_length, header = T,sep = "\t")
    stopifnot(ncol(gene_length) == 2)
    colnames(gene_length) <- c("Geneid","Length")
  }else{
    stop("Could not find the feature length file: ",feature_length," !")
  }
}else{
  gene_length <- read.table(file_name[1], header = T,comment.char = "#") 
  gene_length <- subset(gene_length, select = c("Geneid","Length"))
}

count_df <- data.frame(Geneid = matrix(ncol = 1,nrow = nrow(gene_length)) )
count_df$Geneid <- gene_length$Geneid
for (i in 1:length(file_name)) {
  df_tmp <- read.table(file_name[i], header = T,comment.char = "#")
  df_tmp <- subset(df_tmp, select = c(1,7))
  colnames(df_tmp)[2] <- sub(pattern = pattern,replacement = "",x = basename(file_name[i]))
  # colnames(df_tmp)[2] <- basename(file_name[i]) %>% str_remove("\\.\\w+$")
  # count_df <- count_df %>% full_join(df_tmp,by="Geneid")
  count_df <- merge(count_df,df_tmp, by = "Geneid")
  cat(i," count matrixs have merged!\n")
}
cat("Congratulations! All count matrixs have merged!\n")

# count_mat <- count_df %>% tibble::column_to_rownames(var = "Geneid")
rownames(count_df) <- count_df$Geneid # 关联gene_length与count_df的"Geneid"，相同顺序
count_mat <- count_df[,2:ncol(count_df)]
# count_mat <- count_mat[match(gene_length$Geneid, rownames(count_mat)),]
gene_length <- gene_length[match(rownames(count_mat), gene_length$Geneid),]
# gene_length <- df %>% select(1,2) %>% tibble::column_to_rownames(var = "Geneid")

cat("There are ",nrow(count_mat)," genes in the annotation!\n")
# 仅作查看非0值基因数
filter_count_mat <- count_mat[rowSums(count_mat)>0,]
cat(nrow(filter_count_mat),"(",format( round(nrow(filter_count_mat)/nrow(count_mat)*100, digits = 2), nsmall = 2),"%) genes expressed in all samples!\n")

stopifnot(identical(rownames(count_mat), gene_length$Geneid))
# caculate FPKM
cat("Calculating the FPKM…\n")
fpkm <- edgeR::rpkm(count_mat,gene.length = gene_length$Length)
fpkm <- as.data.frame(fpkm)

# caculate TPM
cat("Calculating the TPM…\n")
fpkm2tpm <- function(fpkm){
  if((is.matrix(fpkm) | is.data.frame(fpkm)) & all(fpkm>=0)){ #fpkm所有值为非负且为矩阵或者数据框
    tpm <- t(t(fpkm)/colSums(fpkm))*10^6
  }
  return(tpm)
}
tpm <- as.data.frame(fpkm2tpm(fpkm = fpkm)) 

# write out 
cat("writing out raw counts matrix\n")
# out_count_mat <- count_mat %>% rownames_to_column(var = feature_name)
out_count_mat <- data.frame(id = rownames(count_mat), count_mat)
colnames(out_count_mat) <- c(feature_name, colnames(count_mat))
count_file <- paste0(output_path,"/",output_prefix,".counts")
write.table(out_count_mat, file = count_file,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("writing out FPKM\n")
# out_fpkm <- fpkm %>% rownames_to_column(var = feature_name)
out_fpkm <- data.frame(id = rownames(fpkm), fpkm)
colnames(out_fpkm) <- c(feature_name, colnames(fpkm))
fpkm_file <- paste0(output_path,"/",output_prefix,".fpkm")
write.table(out_fpkm, file = fpkm_file,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("writing out TPM\n")
# out_tpm <- tpm %>% rownames_to_column(var = feature_name)
out_tpm <- data.frame(id = rownames(fpkm), tpm)
colnames(out_tpm) <- c(feature_name, colnames(tpm))
tpm_file <- paste0(output_path,"/",output_prefix,".tpm")
write.table(out_tpm, file = tpm_file,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("Congratulations! All of the missions have been completed!\n")
