#1. 载入需要的工具包
rm(list = ls())#清除内存变量
gc()
library(cytofkit)
library(dplyr)
library(gplots)
library(ggrepel)
library(Rmisc)
library(colorRamps)
library(RColorBrewer)
library(car)
library(ggpubr)
library(cytofexplorer)
library(reshape2)
library(ggdendro)
library(cowplot)

#2.数据预处理（读取->Downsample->合并->Transform）

#2.1 读取->Downsample->合并
projectdir<-"/media/zhujing/DATA2/1-Increased IL-10-producing regulatory T cells are"
mergeMethod<-"all"
fixedNum<-5000

#目录设置
wdir<-"Pipeline"
raw_fcs_dir<-"PG_output_panel_all"

#读取文件 downsample，merge
wdir<-paste0(projectdir,"/",wdir)
raw_fcs_dir<-paste0(wdir,"/",raw_fcs_dir)
metadata_dir<-paste0(wdir,"/metadata")
#out_Files<-paste0(projectdir,"/transfer_station")
setwd(wdir)


inFile_1<-list.files(raw_fcs_dir,pattern = '.fcs$',full.names = TRUE)
combined_data_raw<-cytof_exprsMerge(fcsFiles=inFile_1,
                                    transformMethod = "none",
                                    mergeMethod = mergeMethod,
                                    fixedNum = fixedNum)
#cytofAsinh
#arcsinh

#简化列表
paraname<-colnames(combined_data_raw)
paraname<-sub(".*<","",paraname)
paraname<-sub(">.*","",paraname)
paraname<-sub("-","_",paraname)
#paraname<-sub("+","_",paraname)
#as.data.frame(combined_data_raw,row.names = NULL,optional = FALSE)
colnames(combined_data_raw)<-paraname

#增加File_ID
File_ID<-sub("_[0-9]*$","",row.names(combined_data_raw))
File_ID
combined_data_raw<-data.frame(combined_data_raw,File_ID)

#Chekpoint
head(combined_data_raw)

#transform
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
transform_id=which(all_markers$transform==1)
transform_id
combined_data_transformed<-combined_data_raw
combined_data_transformed[,transform_id]<-apply(combined_data_transformed[,transform_id,drop=FALSE],
                                                2,cytofAsinh)
#checkpoint
head(combined_data_transformed)

#2.2 marker选项
#打开all_marker.csv，增加expr_para，heatmap两列分别用来指定用来分析差异表达和出现在heatmap中的marker
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"),header = TRUE)
all_markers$markers<-colnames(combined_data_raw)[1:nrow(all_markers)]

#checkpoint
all_markers

#2.3 确定分组
groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors = FALSE)
groups$File_ID<-unique(File_ID)

#checkpoint
groups

#3 针对各个cluster进行统计分析
#source("./backup/Cluster_Expr.R")
#source("./backup/Cluster_Abundance.R")

#3.1 数据初步整理,全局参数设定
cluster_stat<-stat_by_cluster(combined_data_transformed,
                              all_markers,
                              cluster_name="PhenoGraph",
                              summerise_method="median",
                              groups=groups,
                              major_cond="Patient_Type",
                              group_seq=c("HC","W","ICU_W","ICU"), #设置major_cod在各个统计中的顺序
                              stat.paired = T,
                              stat.method = "t-test")
#head(cluster_stat)
#3.2亚群丰度 ***
#绘制boxplot
#source("backup/Cluster_Abundance.R")
draw_abundance_report(cluster_stat,
                      heatmap_ctrl = "W") #heatmap_ctrl选择样本数量多的组

#绘制亚群丰度火山图
draw_abundance_volcano(cluster_stat,
                       cond1="ICU",
                       cond2="W",
                       #stat.method = "t-test"
                       )

#3.3 marker差异表达分析 ***
#需要all_markers.csv中含有expr_para一列
#生成每个cluster的expr_para的heatmap，boxplot 
cluster_expr_report(cluster_stat,
                    #cluster_id=c(3,39), #根据细胞类型选择cluster_id
                    heatmap_ctrl=c("ICU_W","ICU","W","HC"),
                    subgroup_cond = "Patient_Type",
                    #heatmap1_trans_method = "0_to_Max"
                    )

#绘制差异表达火山图
draw_expr_volcano(cluster_stat,
                  cond1="W",
                  cond2="ICU",
                  cluster_id = 7,
                  #cluster_id = c(23,30),
                  dif.level = 1.5,
                  stat.paired = FALSE,
                  #stat.method = "auto",
                  conf.level = 0.95,
                  )
#3.4 绘制cluster的heatmap ***
#需要all_markers.csv中含有heatmap一列
#source("./backup/Heatmap_Output.R")
cluster_expr_matrix<-cluster_stat$cluster_expr_matrix
#'到metadata目录中找到samplename.csv,打开后，
#'在后面添加一列，列名为heatmap，所有需要在热图上显示的marker上标记1
#'多余的列可以保留，只要列名不重复，就不会影响程序运行。
#'保存后退出。

all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"),header = TRUE)
heatmap_ID<-as.character(dplyr::filter(all_markers,heatmap==1)$markers)
draw_expr_heatmap(xdata=cluster_expr_matrix[,heatmap_ID],
                  Rowv=T,
                  Colv=T,
                  trans_method="simpleAsinh_Min_to_Max",
                  color_style=3
                  )
"""
PhenoGraph_ID<-as.character(dplyr::filter(all_markers,PhenoGraph==1)$markers)
cluster_marker_preview(cluster_expr_matrix[,PhenoGraph_ID],
                       groups = groups,
                       all_markers = all_markers,
                       cluster_name = "PhenoGraph",
                       major_cond = "HC",
                       cluster_id=c(1,2,3,4)
)
"""