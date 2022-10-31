#!/usr/bin/env Rscript

Sys.setenv(LANGUAGE = "en")
#mamba  install -c conda-forge -c bioconda r-argparse r-jsonlite r-tidyr  bioconductor-annotationforge bioconductor-go.db r-pkgbuild r-vroom r-dplyr
#library(jsonlite)
#library(readr)
#library(dplyr)
#library(vroom)
#library(tidyr)
#library(stats)
#library(utils)
#library(AnnotationForge)
#library(pkgbuild)
#library(argparse)
start_time <- Sys.time() #记录开始时间
suppressPackageStartupMessages(library(magrittr, quietly = T))

# 0. 参数

parser <- argparse::ArgumentParser(description = "Creating ORGDB package and get some files from emapper2 and 'ko00001.json' of KEGG for non-model organisms!")
parser$add_argument("-i","--emapper2_annotations",action = "store",type = "character", 
                    help = "Annotations file from emapper2.")
# type: 'logical', 'integer', 'double' or 'character'
parser$add_argument("-k", "--ko_json", action = "store",type = "character", default = NULL, required = FALSE,
                    help = "Specify the 'ko00001.json'. If not given, it will be downloaded from 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir='.")
parser$add_argument("-g","--genus", type = "character", default = "Actinidia",
                    help = "The genus of organism.[default: Actinidia]")
parser$add_argument("-s","--species", type = "character", default = "'chinensis var. chinensis'",
                    help = "The species of organism.[default: 'chinensis var. chinensis']")
parser$add_argument("-t", "--tax_id", type = "integer", default = 1590841, #这个参数integer或者character应该都没问题。
                    help = "Taxonomy ID  that represents your organism. Taxonomy's address: 'https://www.ncbi.nlm.nih.gov/taxonomy'[default: 1590841]")
parser$add_argument("-o", "--out_dir", type = "character", default = NULL, 
                    help = "Specify the directory of output. If not existed, create it." )
parser$add_argument("-v", "--version", type = "character", default = "0.1",
                    help = "The version of ORGDB package.format: \'x.y.z\'[default: 0.1]")
parser$add_argument("-w", "--write_out", action='store_true', 
                    help = "Write out or don't write out some files from emapper2 and 'ko00001.json' to '--out_dir'. If the option '-w' exists, all files will be written out."  )
# 没参数提示
# 默认选项default设置好像不能自动添加到帮助文档中。
# help 内不知道如何换行
# help 内禁止手动回车换行！！！否则win系统写完的脚本，在linux下无法识别，同时也无法用dos2unix转换。
args <- parser$parse_args()
#parser$print_help() #打印帮助文档,用于调试

# 尝试模拟解析给出参数后的现有赋值。
#parser$parse_args()
#parser$parse_args(c('-omy','-w'))

emapper2_annotations_file <- args$emapper2_annotations
ko_json <- args$ko_json
genus <- args$genus
species <- args$species
tax_id <- args$tax_id
out_dir <- args$out_dir
version <- args$version
write_out <- args$write_out

# 参数检查
if(is.null(emapper2_annotations_file)){
  parser$print_help() 
  stop("Option '-i' is requited!")
  }
if(is.null(out_dir)){
  parser$print_help()
  stop("Option '-o' is requited!")
  }
if(!file.exists(emapper2_annotations_file)){
  parser$print_help()
  stop(emapper2_annotations_file,": No such file!\n")
  }
if(!is.null(ko_json)){
  if(!file.exists(ko_json)) {
    parser$print_help()
    stop(ko_json,": No such file!\n")
  }else{
      message(ko_json,": ready!")
    }
}

# 参数示例，写脚本前期用于调试
if(F){
  emapper2_annotations_file <- "./Red5.emapper.annotations"
  ko_json <- NULL
  genus = "Actinidia"
  species = "chinensis var. chinensis"
  tax_id = "1590841"
  out_dir <- "./my_restlts"
  version = "0.1"
  write_out <- TRUE
}


# 1.1 prepare the name of ord.db package----
# 放到后面
#genus <- genus %>% str_sub(1,1)
#species <- species %>% str_sub(1,2)

# 1.4 提前检查是否有同名文件，有的话必须提前终止。
dbname <- paste0(stringr::str_sub(genus,1,1),stringr::str_sub(species,1,2))
orgdb <- paste0("org.",dbname,".eg.db")
orgdb_sqlite <- paste0("org.",dbname,".eg.sqlite")
if(dir.exists(orgdb)){stop("Directory of ",orgdb," exists, delete it or give anthor '--genus' or '--species'.")}
if(file.exists(orgdb_sqlite)){stop("Directory of ",orgdb_sqlite," exists, delete it or give anthor '--genus' or '--species'.\n")}

# 2.1 从annotation获取gene2go----
get_go_info <- function(emapper2_annotations_file,out_dir, write_out){
  # 1.2 check output directory!
  if(!file.exists(out_dir)){
    dir.create(out_dir)
    message(out_dir,": no such directory!\n",
            out_dir,": created!\n")
  }else{
    message(out_dir,": ready!\n")
  }
  message("Note:",emapper2_annotations_file,": must be an annotation file from emapper!\n",
          "    And its formats refer to emapper 2.1.7 and eggnog 5.0.2\n")
  suppressWarnings(
    emmaper <- readr::read_delim(emapper2_annotations_file, delim = "\t",
                                 col_names = TRUE, comment = "##",
                                 trim_ws = TRUE,
                                 escape_double = FALSE,
                                 na = "-",
                                 show_col_types = FALSE))
  
  colnames(emmaper)[1] <- "query"
  
  extract_colname <- c("query","COG_category","Description",
                       "Preferred_name","GOs","KEGG_ko","KEGG_Pathway","PFAMs")
  
  new_colname <- c("GID","COG","GENENAME","SYMBOL","GO","KO","Pathway","Pfam")
  
  emmaper2 <- emmaper %>% dplyr::select(match(extract_colname,colnames(emmaper)))
  colnames(emmaper2) <- new_colname
  if(write_out){
    vroom::vroom_write(emmaper2, file = paste0(out_dir,"/emmaper2_eff_info.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  }
  
  ## 简单统计
  ### 总基因条数
  num_gene <- emmaper2 %>% nrow() #32561
  message("Total nummber of genes: ",num_gene)
  ### 有SYMBOL
  num_SYMBOL <- emmaper2 %>% dplyr::filter(!is.na(SYMBOL)) %>% nrow() #2580
  message("Number of genes containing the SYMBOL: ",num_SYMBOL)
  ### 有COG
  num_gene_cog <- emmaper2 %>% dplyr::filter(!is.na(COG)) %>% nrow() #30242
  message("Number of genes containing the COG: ",num_gene_cog)
  ### 有Pfam
  num_pfam <- emmaper2 %>% dplyr::filter(!is.na(Pfam)) %>% nrow() #29174
  message("Number of genes detecting in Pfam: ",num_pfam)
  ### 有GO
  num_go <- emmaper2 %>% dplyr::filter(!is.na(GO)) %>% nrow() #17706
  message("Number of genes detecting in GO: ",num_go)
  ### 有KO
  num_KEGG_KO <- emmaper2 %>% dplyr::filter(!is.na(KO)) %>% nrow() #16263
  message("Number of genes containing KEGG KO: ",num_KEGG_KO)
  
  ### 有Pathway
  num_KEGG_Pathway <- emmaper2 %>% dplyr::filter(!is.na(Pathway)) %>% nrow() #10170
  message("Number of genes containing KEGG Pathway: ",num_KEGG_Pathway)
  
  ## 提取GO信息
  ### gene ID 与推测的 Gene Symbol
  gene_info <- emmaper2 %>% dplyr::select(GID,GENENAME ) %>% dplyr::filter(!is.na(GENENAME)) #差点选反了
  gene2symbol <- emmaper2 %>% dplyr::select(GID,SYMBOL ) %>% dplyr::filter(!is.na(SYMBOL))
  
  ### 提取gene ID与GO ID
  gene2go <- emmaper2 %>% dplyr::select(GID,GO) %>% 
    dplyr::filter(!is.na(GO)) %>% #去掉空值，就是留下判断为真的内容。
    tidyr::separate_rows(GO,sep = ",", convert = TRUE) %>%  #拆分GO列，使其one2one;convert 有利于列的类型准确显示 numeric, integer, or logical.
    dplyr::mutate(EVIDENCE = 'IEA') %>%  #添加列
    dplyr::filter(!is.na(GO))
  ## 提取KEGG Pathway信息
  gene2pathway <- emmaper2 %>% dplyr::select(GID, Pathway) %>% 
    dplyr::filter(!is.na(Pathway)) %>% 
    tidyr::separate_rows(Pathway, sep = ",", convert = TRUE) %>% 
    dplyr::filter(stringr::str_detect(Pathway,"ko")) %>%   #过滤掉map编号。map和ko的数字编号对应
    dplyr::filter(!is.na(Pathway))
  
  ## 提取KEGG KO信息
  gene2ko <- emmaper2 %>% dplyr::select(GID,KO) %>% 
    dplyr::filter(!is.na(KO)) %>% 
    tidyr::separate_rows(KO, sep = ",",convert = TRUE)  
  gene2ko$KO <- gene2ko$KO %>% stringr::str_remove("ko:")
  gene2ko <- gene2ko %>% dplyr::distinct() #无重复呀
  
  #str(gene2pathway)
  #summary(gene2pathway)
  #tibble(gene2pathway)
  results <- list()
  results$gene_info <- gene_info
  results$gene2symbol <- gene2symbol
  results$gene2go <- gene2go
  results$gene2pathway <- gene2pathway
  results$gene2ko <- gene2ko
  return(results)
}


# 3.1 定义函数----
get_kegg_info <- function(ko_json){
  # 1.3 prepare ko00001.json file! 
  
  if(!is.null(ko_json)){
    ko_json <- ko_json
    if(!file.exists(ko_json)){stop(ko_json,": No such file!")}
  }else {
    
    if(!dir.exists("./ko00001")){dir.create("ko00001")}
    
    if(!file.exists("./ko00001/ko00001.json")){ #必须force = TRUE 放在()中。| (force = TRUE)
      mess
      utils::download.file("https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=",
                           destfile = "./ko00001/ko00001.json")
    }
    
    if(file.exists("./ko00001/ko00001.json")) {message("\nThe ko00001.json exists !\n")}
    ko_json <- "./ko00001/ko00001.json"
    if(!file.exists(ko_json)){stop("Failed to download ",ko_json,". ")}
    
  }
  
  # 提取
  results <- list()
  kegg <- jsonlite::fromJSON(ko_json)
  pathway2name <- tibble::tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble::tibble(KO = character(), Pathway = character())
  gene_KO_Symbol_Name <- tibble::tibble(Entry = character(), Symbol = character(), Name = character())
  for (a in seq_along(kegg[["children"]][["children"]])) { #索引1-8
    A <- kegg[["children"]][["name"]][[a]] #A的name
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) { #A的children
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]  #B的name
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) { #B的children
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]] #C的name
        pathway_id <- stringr::str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- stringr::str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% stringr::str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble::tibble(Pathway = pathway_id, Name = pathway_name)) #需要值
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos <- stringr::str_match(kos_info, "K[0-9]*")[,1] #[,1]转为向量
        ko_symbol <- stringr::str_match(kos_info,"K[0-9]*.*;") %>% stringr::str_remove("K[0-9]*  ") %>% stringr::str_remove(";")
        ko_name <- stringr::str_match(kos_info,"; .*") %>% stringr::str_remove("; ")
        ko2pathway <- rbind(ko2pathway, tibble::tibble(KO = kos, Pathway = rep(pathway_id, length(kos)))) #需要值
        gene_KO_Symbol_Name <- rbind(gene_KO_Symbol_Name,tibble::tibble(Entry = kos, Symbol = ko_symbol, Name = ko_name)) #需要值，有用值
        
      }}}
  results$pathway2name <- pathway2name %>% stats::na.omit()  
  results$ko2pathway <- ko2pathway %>% dplyr::distinct() %>% stats::na.omit() 
  results$gene_KO_Symbol_Name <- gene_KO_Symbol_Name
  return(results)
}
# 描述大KO：kos_info
# gene_KO_Symbol_Name 提取了基因大KO的ID，Symbo和Name。
# pathway2name 通路ID 小ko与通路名称 #小ko在kegg官网以map开头。
# ko2pathway  基因大Ko 与通路小ko

get_go_kegg_info <- function(emapper2_annotations_file, ko_json, out_dir,write_out){
  # 2.2 运行----
  results_go <- get_go_info(emapper2_annotations_file = emapper2_annotations_file,out_dir=out_dir,write_out = write_out)
  # 2.2 运行
  results_kegg <- get_kegg_info(ko_json = ko_json)
  
  results <- list()
  results$results_go <- results_go
  results$results_kegg <- results_kegg
  return(results)
}

message("Extracting some effective information!\n")
results <- get_go_kegg_info(emapper2_annotations_file = emapper2_annotations_file, ko_json = ko_json, out_dir = out_dir, write_out = write_out)
invisible(gc())



## 物种信息taixd 1590841
## genus = "Actinidia"
## species = "chinensis var. chinensis"

# 3.1 构建OrgDb----
## 构建OrgDb
#gene_info <- results$results_go$gene_info
#gene2symbol <- results$results_go$gene2symbol
#gene2go <- results$results_go$gene2go 
####gene2pathway <- results_go$gene2pathway#gene2pathway 来自eggnog注释，40955行
gene2ko <- results$results_go$gene2ko
ko2pathway <- results$results_kegg$ko2pathway 
new_gene2pathway <- gene2ko %>% dplyr::left_join(ko2pathway,by = "KO") %>% dplyr::select(1,3) %>% dplyr::distinct() %>% na.omit() #eggnog的gene2ko与KEGG的ko2pathway
message("'gene2pathway' comes from gene2ko of eggnog and ko2pathway of eggnog!\n")

create_ORGDB <- function(gene_info = results$results_go$gene_info, symbol = results$results_go$gene2symbol, ko = new_gene2pathway, go = results$results_go$gene2go,
                         genus = genus, species = species,tax_id = tax_id, version = version, out_dir = out_dir){
  dbname <- paste0(stringr::str_sub(genus,1,1),stringr::str_sub(species,1,2))
  orgdb <- paste0("org.",dbname,".eg.db")
  
  message("Creating ORGDB package: ",orgdb,"\n")
  AnnotationForge::makeOrgPackage(gene_info = results$results_go$gene_info,
                                  symbol = results$results_go$gene2symbol,
                                  ko = new_gene2pathway, #gene2pathway
                                  go = results$results_go$gene2go, #gene2go
                                  version = version,
                                  maintainer = "Shimao Zheng <zhengshimao007@163.com>",
                                  author = "Shimao Zheng  <zhengshimao007@163.com>",
                                  outputDir = ".",
                                  tax_id = tax_id,
                                  genus = stringr::str_sub(genus,1,1),
                                  species = stringr::str_sub(species,1,2),
                                  goTable = "go",
                                  verbose = TRUE)
  # makeOrgPackage() 函数注意事项
  #unlink = TRUE, The 1st column must always be the gene ID 'GID'
  # 函数内最后一个参数不能以“,”结尾。
  
  #system("R CMD build org.Ach.eg.db") #这种运行方式也可以
  if(file.exists(orgdb)){
    pkgbuild::build(orgdb,dest_path = out_dir)
  }
}

# 构建
create_ORGDB(gene_info = results$results_go$gene_info, symbol = results$results_go$gene2symbol, ko = new_gene2pathway, go = results$results_go$gene2go,
  genus = genus, species = species,tax_id = tax_id, version = version, out_dir = out_dir)
invisible(gc())


# write out some files 
if(write_out){
  
  message("Writing out some file from emapper2 and 'ko00001.json'.")
  vroom::vroom_write(results$results_go$gene_info, file = paste0(out_dir,"/emapper2_gene_info.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  vroom::vroom_write(results$results_go$gene2symbol, file = paste0(out_dir,"/emapper2_gene2symbol.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  vroom::vroom_write(results$results_go$gene2go, file = paste0(out_dir,"/emapper2_gene2go.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  vroom::vroom_write(results$results_go$gene2ko, file = paste0(out_dir,"/emapper2_gene2ko.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  vroom::vroom_write(results$results_kegg$pathway2name, file = paste0(out_dir,"/kegg_pathway2name.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  vroom::vroom_write(results$results_kegg$ko2pathway, file = paste0(out_dir,"/kegg_ko2pathway.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
  vroom::vroom_write(results$results_kegg$gene_KO_Symbol_Name, file = paste0(out_dir,"/kegg_gene_KO_Symbol_Name.tsv"), quote = "none", delim = "\t", na = "", col_names = T)
}

invisible(gc())

message("\n\nDONE!\n")

runtime <- difftime(Sys.time(), start_time, units = "mins")
message("Time consuming: ",round(runtime,2)," mins\n")

