#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(edgeR))
suppressMessages(library(argparser)) #https://github.com/cran/argparser


# 参数设置
p <- arg_parser("Merge featureCounts(linux version) and calculate FPKM/TPM")

p <- add_argument(p, "--input_path", help="input: a directory containing the counts matrix named with '<sample>.count'", type="character",default = "./")
p <- add_argument(p, "--pattern", help="limit pattern of input files using regular expression in R language",type="character",default = "*count$")
p <- add_argument(p, "--output_path", help="output: an existent directory", type="character",default = "./")
p <- add_argument(p, "--prefix", help="give the file of output matrix a prefix like '<output_prefix>_genes.*'", type="character",default = "my",short = "-f")

# 参数解析与定义
argv <- parse_args(p)

path <- argv$input_path
pattern <- argv$pattern
output_path <- argv$output_path
output_prefix <- argv$prefix

# 代码主体
file_name <- dir(path = path,pattern = pattern)
file <- paste0(path,"/",file_name)
# merge all count matrix
df <- read.table(file[1], header = T,comment.char = "#") %>% select(c(1,6,7))
colnames(df)[3] <- basename(file[1]) %>% str_remove("\\.\\w+$")
cat("1 count matrix has merged!\n")
for (i in 2:length(file)) {
  df_tmp <- read.table(file[i], header = T,comment.char = "#") %>% select(c(1,7))
  colnames(df_tmp)[2] <- basename(file[i]) %>% str_remove("\\.\\w+$")
  df <- df %>% full_join(df_tmp,by="Geneid")
  cat(i," count matrixs have merged!\n")
}
cat("Congratulations! All count matrixs have merged!\n")

count_mat <- df %>% select(-2) %>% tibble::column_to_rownames(var = "Geneid")
gene_length <- df %>% select(1,2) %>% tibble::column_to_rownames(var = "Geneid")
filter_count_mat <- count_mat[rowSums(count_mat)>0,]
cat(nrow(count_mat)-nrow(filter_count_mat)," genes were not expressed in all samples!\n")

# caculate FPKM
cat("Calculating the FPKM…\n")
fpkm <- rpkm(count_mat,gene.length = gene_length$Length) %>% as.data.frame()

# caculate TPM
cat("Calculating the TPM…\n")
fpkm2tpm <- function(fpkm){
  if((is.matrix(fpkm) | is.data.frame(fpkm)) & all(fpkm>=0)){ #fpkm所有值为非负且为矩阵或者数据框
    tpm <- t(t(fpkm)/colSums(fpkm))*10^6
  }
  return(tpm)
}
tpm <- fpkm2tpm(fpkm = fpkm) %>% as.data.frame() 

# write out 
cat("writing out raw counts matrix\n")
out_count_mat <- count_mat %>% rownames_to_column(var = "gene_id")
count_file <- paste0(output_path,"/",output_prefix,"_","genes.counts")
write.table(out_count_mat, file = count_file,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("writing out FPKM\n")
out_fpkm <- fpkm %>% rownames_to_column(var = "gene_id")
fpkm_file <- paste0(output_path,"/",output_prefix,"_","genes.fpkm")
write.table(out_fpkm, file = fpkm_file,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("writing out TPM\n")
out_tpm <- tpm %>% rownames_to_column(var = "gene_id")
tpm_file <- paste0(output_path,"/",output_prefix,"_","genes.tpm")
write.table(tpm, file = tpm_file,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("Congratulations! All of the missions have been completed!\n")
