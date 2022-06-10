#1、载入需要的工具包
rm(list=ls()) #清除内存变量
library(flowCore)
library(Rcpp)
library(cytofkit)
library(igraph)
library(ggplot2)
library(ggthemes)
library(Rtsne)
library(dplyr)
library(cytofexplorer)

#2 数据预处理（读取->Downsample->合并->Transform）
#2.1 读取FCS文件，Downsample, Merge等
# 输入
projectdir<-"/media/zhujing/DATA2/5-A Streamlined CyTOF Workflow To Facilitate Standardized Multi-Site Immune Profiling of COVID-19 Patients"
mergeMethod<-"ceil" #合并方法："ceil",'all','fixed','min'
fixedNum<-5000 #设置从每个文件抽取的细胞数

#目录设置
wdir<-"Pipeline"
raw_fcs_dir<-"Singelets"
wdir<-paste0(projectdir,"/",wdir)
raw_fcs_dir<-paste0(projectdir,"/",raw_fcs_dir)
metadata_dir<-paste0(wdir,"/metadata")
out_Files<-paste0(wdir,"/transfer_station")
setwd(wdir)

"""
#读取文件,Downsample,Merge
inFile<-list.files(raw_fcs_dir,pattern = '.fcs$',full.names = TRUE)
foo<-read.flowSet(inFile, emptyValue = FALSE)
#foo<-read.FCS(inFile,emptyValue = FALSE)
outFile <- file.path(out_Files)
write.flowSet(foo,outFile)
#write.FCS(foo,outFile)
#运行后注释掉以上代码
"""
inFile_1<-list.files(out_Files,pattern = '.fcs$',full.names = TRUE)
combined_data_raw<-cytof_exprsMerge(fcsFiles=inFile_1,
                                    transformMethod = "none",
                                    mergeMethod = mergeMethod,
                                    fixedNum = fixedNum)

#简化列表
paraname<-colnames(combined_data_raw)
#paraname
paraname<-sub(".*<","",paraname)
paraname<-sub(">.*","",paraname)
paraname<-sub("-","_",paraname)
#paraname<-sub("+","_",paraname)
#as.data.frame(combined_data_raw,row.names = NULL,optional = FALSE)
colnames(combined_data_raw)<-paraname

#增加File_ID
File_ID<-sub("_[0-9]*$","",row.names(combined_data_raw))
#File_ID
combined_data_raw<-data.frame(combined_data_raw,File_ID)

#Chekpoint
head(combined_data_raw)

#2.2 Transform 数据转换:针对部分通道，方法 cytofAsinh
# input:
#重要：通道选择，仅仅首次需要input
#write.csv(colnames(combined_data_raw),paste0(metadata_dir,"/all_markers.csv"), row.names=FALSE) 
#仅第一次运行，通道一旦运行后注释掉上一行代码。
#到工作目录中找到all_markers.csv,打开后，首行A列标明“markers”，B列标明“transform”
#B列中在要transform的标签行标记1，保存退出
#后面可以把write.csv......这行注释，就不会再运行此行。

#确定分组,文件分组，仅仅首次需要input
#write.csv(unique(File_ID),paste0(metadata_dir,"/all_samples.csv"),row.names = FALSE)

#到metadata目录中找到samplename.csv，打开后
#A列名设为File_ID，根据样本分组情况设置后面几列；
#B列名设为Short_name,每个样本的简称，可以方便的出现在输出的图片中；
#C列以后根据实验设计，输入实际的分组情况，例如：Timepoint,Tissue_Type,Patient_Type ect.
#第一次运行，通道一旦选好可以将此行注释

#数据转换
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
transform_id=which(all_markers$transform==1)
transform_id
combined_data_transformed<-combined_data_raw
combined_data_transformed[,transform_id]<-apply(combined_data_transformed[,transform_id,drop=FALSE],
                                                2,cytofAsinh)

#checkpoint
head(combined_data_transformed)

#3. 运行phenograph聚类（仅有一个参数K，默认30）
#输入
#phenograph参数设置：
k=30 #计算knn网络的邻居个数
#通道选择：到metadata目录中找到all_markers.csv,打开后，在右边选择一个空列输入“PhenoGraph”，
#并在要进行phenogtraph聚类的标签行标记1，保存退出。

#模块主体部分
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
PhenoGraph_id=which(all_markers$PhenoGraph==1)  #给出向量、数组或类似列表的对象中被认为是真的值的索引，在类似数组的对象中允许使用数组索引。 
PhenoGraph_input_data=combined_data_transformed[,PhenoGraph_id]

#phenograph的elbow测试，由于要进行多伦聚类测试，耗时较长，建议在total event较少的情况下使用
#PG_elbow(PhenoGraph_input_data)

PhenoGraph_result<-as.numeric(membership(Rphenograph(data=PhenoGraph_input_data,k=k)))

#Checkpoint
hist(PhenoGraph_result)


#4. 将结果导出成FCS文件
##input：整合待导出数据
combined_data_output<-data.frame(PhenoGraph=PhenoGraph_result)
head(combined_data_output)

#根据File_ID将合并数据还原成单个样本数据，导出FCS文件
row.names(combined_data_output)<-row.names(combined_data_raw)
head(combined_data_output)
cytof_addToFCS_modified(data=combined_data_output,
                        rawFCSdir = out_Files,
                        analyzedFCSdir = paste(wdir,"/PG_output_panel_all",sep=""),
                        newHeader = "PG_")

#运行结束后注释掉以上代码。

#5 整合数据，从各个meta cluster或文件中抽取定量细胞进行可视化（downsample）
#source("./backup/FlowSOM_tSNE_Output.R")
combined_data_analysed<-data.frame(combined_data_transformed,
                                   PhenoGraph=PhenoGraph_result)
head(combined_data_analysed)
combined_data_sampled<-equal_sample(x=combined_data_analysed,
                                    sample_method="ceil", #两种取值：“ceil”，“all”
                                    sample_type="by_file",
                                    sample_num=1000,
                                    cluster_name="PhenoGraph")
head(combined_data_sampled)
str(combined_data_sampled)

#6. 运行tSNE(BH-sne)降维（two options，二选一，本步骤需要时间较长）
#Input：tSNE参数设置：
max_iter=10000  #迭代次数
perplexity=30  #困惑度
#seed=123  #随机种子
theta=0.5  #权衡速度与准确度，越小越精确，越大速度越快
dims=2  #降维输出维度（默认=2）

#通道选择：找到all_markers.csv,打开后，在新一列首行输入“tSNE”,
#并在要进行tSNE降维的标签行标记1，保存退出。

#tsne分析
if (exists('seed')) set.seed(seed)
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
tSNE_para_id=which(all_markers$tSNE==1)
tSNE_input_data=combined_data_sampled[,tSNE_para_id]
#head(tSNE_input_data)

tsne_result<-Rtsne(tSNE_input_data,
                   initial_dims=ncol(tSNE_input_data),
                   pca=FALSE,
                   dims=dims,
                   check_duplicates=FALSE,
                   perplexity=perplexity,
                   max_iter=max_iter,
                   theta=theta)$Y
row.names(tsne_result)<-row.names(combined_data_sampled)
colnames(tsne_result)<-c("tsne_1","tsne_2")
head(tsne_result)

#checkpoint:一下语句可以绘制出降维分析的草图
plot(tsne_result)

#7 将结果导出成FCS文件
#source("./backup/modified_cytofkit_functions.R")

#Input:整合待导出数据
combined_data_output<-data.frame(tsne_result,
                                 PhenoGraph=combined_data_sampled$PhenoGraph)
head(combined_data_output)

#根据File_ID将合并数据还原为单个样本数据，导出FCS文件
#row.names(combined_data_output)<-row.names(combined_data_raw)
cytof_addToFCS_modified(data = combined_data_output,
                        rawFCSdir = out_Files,
                        analyzedFCSdir = paste(wdir,"/PGtsne_output",sep = ""),
                        newHeader = "PG_tsne_"
)

#8.降维和聚类结果可视化
#确定分组,文件分组，仅仅首次需要input
#write.csv(unique(File_ID),paste0(metadata_dir,"/all_samples.csv"),row.names = FALSE)

#到metadata目录中找到samplename.csv，打开后
#A列名设为File_ID，根据样本分组情况设置后面几列；
#B列名设为Short_name,每个样本的简称，可以方便的出现在输出的图片中；
#C列以后根据实验设计，输入实际的分组情况，例如：Timepoint,Tissue_Type,Patient_Type ect.
#第一次运行，通道一旦选好可以将此行注释

groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors = FALSE)
groups$File_ID<-unique(File_ID)

#checkpoint
groups

#整合聚类和降维数据
#source("./backup/FlowSOM_tSNE_Output.R")
combined_data_plot<-cbind(combined_data_sampled,
                           tsne_result)
#checkpoint
head(combined_data_plot)

#8.1.降维和聚类结果作图（R语言脚本直接生成）
#生成图片
combined_data_plot$PhenoGraph<-as.factor(combined_data_plot$PhenoGraph)
centers<-combined_data_plot%>%
  group_by(PhenoGraph) %>%
  summarise(tsne_1=median(tsne_1),tsne_2=median(tsne_2))

mytheme<-theme(panel.background = element_rect(fill = "white",colour = "black",size=0.2),
               legend.key = element_rect(fill = "white",colour = "white"), #图标
               legend.background = (element_rect(colour = "white",fill="white")))

#output1
##生成合并tSNE-Phenograph图像
ggplot(combined_data_plot)+
  geom_point(aes(x=tsne_1,y=tsne_2,colour=PhenoGraph))+
  mytheme+
  geom_text(data=centers,aes(x=tsne_1,y=tsne_2),label=centers$PhenoGraph,colour="black",size=5)

#output2
##生成单个文件的tSNE-Phenograph图像
ggplot(combined_data_plot)+
  geom_point(aes(x=tsne_1,y=tsne_2,colour=PhenoGraph))+
  mytheme+
  facet_wrap(~File_ID)

#8.2生成tsne-Phenograph系列图片 ***
#source("./backup/FlowSOM_tSNE Analysis backup (Pipeline2-2).R")
draw_tsne_figs(combined_data_plot = combined_data_plot,
               groups = groups,
               cluster_color = dif_seq_rainbow,
               cluster_name = "PhenoGraph",
               major_cond = "Patient_Type",
               reduction_dm1 = "tsne_1",
               reduction_dm2 = "tsne_2",
               #dot_size = 5,
               #text_size = 20
               )

#生成marker-Cluster 系列density图片 all_markers.csv中标注density_plot
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
draw_density_plots(combined_data_plot,
                   groups = groups,
                   all_markers = all_markers,
                   cluster_color = dif_seq_rainbow, #选择颜色，其他选择：rainbow, dif_seq_rainbow,brewer_color_sets
                   cluster_name = "PhenoGraph",     #cluster通道的名称
                   cluster_id=c(1,2,3,4)
                   )
#运行后注释掉以上代码


#8.3生成所有marker的tsne热图 ***
#打开all_markers.csv，在后面的一个空列首行标记“tsne_heatmap”,
#然后在要生成热图的所在行标记1，保存退出
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
all_markers$markers<-colnames(combined_data_transformed)[1:nrow(all_markers)]
heatmap_tsne_id<-which(all_markers$tsne_heatmap==1)
heatmap_tsne_markers=colnames(combined_data_transformed[,heatmap_tsne_id])
heatmap_tsne_markers
draw_tsne_heatmaps(combined_data_plot=combined_data_plot,
                   heatmap_tsne_markers=heatmap_tsne_markers,
                   groups=groups,
                   trans_method = "Min_to_Max",
                   single_file=TRUE,
                   major_cond="Patient_Type",
                   #output_format = "pdf"
                   )

