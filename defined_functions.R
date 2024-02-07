library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(magrittr)
library(DoubletFinder)
library(future)
# this is linux version

check_gene = function(membrane.data){
  gene = c()
  gene_lst = unique(membrane.data$features.plot)
  for (current_gene in gene_lst){
    scope = subset(membrane.data,subset = features.plot ==current_gene)
    avg.exp = scope$avg.exp
    pct.exp = scope$pct.exp
    # cat(avg.exp)
    flag.a = 0
    flag.b = 0
    flag.c = 0

    # universal low expression
    for (individual in avg.exp){
      #
      # cat(str_interp('${individual}\n'))
      if (individual >1){
        flag.a = flag.a + 1
      }
    }
    #
    for (individual in pct.exp){
      # universal low percentage
      if (individual >5){
        flag.b = flag.b + 1
      }
      # universal high percentage
      if (individual >75){
        flag.c = flag.c + 1
      }
    }

    cnt = 1
    if ((flag.a == 0) | (flag.b == 0) | (flag.c == length(pct.exp))){
      # cat(str_interp('current gene ${current_gene} out!\n'))
      # cat(str_interp('${flag.a} ${flag.b};'))
      flag.a
    }else{
      cat(str_interp("'${current_gene}',"))
      gene = c(gene,current_gene)
      cnt = cnt + 1
    }
    rm(scope)
  }
  return (gene)
}


plan("multisession",workers=20) # 72 cores, 29mins, for seurat.object_1 # 16 cores, 20mins?
# availableCores()
# nbrOfWorkers()
options(future.globals.maxSize = 1 * 1024^4) # 1T
# use S3 to manage plots, for simplicity, flexibility
## constructor
create.Plots = function(){
  structure(
    list(count=0,
         plot_list=list()),
    class="Plots"
  )
}
## methods
add.Plots = function(plots,plot,name=NULL,file_path=NULL){
  plots$count = plots$count+1
  plot_name = if (is.null(name)) paste0("p",plots$count) else name
  plots$plot_list[[plot_name]] = list(plot=plot,file_path=file_path)
  return(plots)
}

get.Plots = function(plots,name){
  if(name%in% names(plots$plot_list)){
    plots$plot_list[[name]]
  }else{
    stop("Plot not found")
  }
}

draw.Plots = function(plots,plot_sizes= list(
                        list(width=8,height=8),
                        list(width=8,height=8),
                        list(width=8,height=8),
                        list(width=8,height=8),
                        list(width=30,height=100),
                        list(width=30,height=100),
                        list(width=15,height=10),
                        list(width=15,height=10),
                        list(width=15,height=10),
                        list(width=8,height=8)
                      ),
                      path='./',suffix='pdf')
{
  print(paste("Saving",plots$cnt,"plots"))
  for(i in seq_along(plots$plot_list)){
    plot_info = plots$plot_list[[i]]
    plot = plot_info$plot
    func = plot_info$file_path
    date_stamp = format(Sys.Date(),"%y%m%d")
    plot_cnt_padded = sprintf("%02d",as.numeric(sub("p","",names(plots$plot_list)[i])))
    filename = paste0(path,plot_cnt_padded,"_",func,"_",date_stamp,".",suffix)
    print(filename)
    width = plot_sizes[[i]]$width
    height = plot_sizes[[i]]$height
    if (suffix %in% c("png","jpg")){
      ggsave(filename,plot,width=width,height=height,dpi=300)
    }else{
      pdf(filename,width = width,height = height);print(plot);dev.off()
    }
    print(paste("Plotting:",filename,"\tWidth:",width,"\tHeight:",height))
  }
}

plot_sizes = list(
  list(width=8,height=8),
  list(width=8,height=8),
  list(width=8,height=8),
  list(width=8,height=8),
  list(width=30,height=100),
  list(width=30,height=100),
  list(width=15,height=10),
  list(width=15,height=10),
  list(width=15,height=10),
  list(width=8,height=8)
)

cellranger_qc = function(dir_path,project,output_dir="D:/fzy/10x_DMH/0_qc/"){
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  system(paste("tree -L 2", output_dir))
  print(paste(output_dir,project,sep = ""))
  setwd(paste(output_dir,project,sep = ""))
  start_time_cellranger_qc = Sys.time()
  #initialize the Seurat object with the raw (non-normalized data)
  seurat.object <- CreateSeuratObject(counts = Read10X(dir_path),
                                      min.cells = 3, min.features =  200,
                                      project = project)
  seurat.object = seurat.object %>%
    PercentageFeatureSet(object = ., pattern = "^mt-", col.name = "percent.mt") %>%
    PercentageFeatureSet(object = ., pattern = "^Rp", col.name = "percent.ribo") %>%
    PercentageFeatureSet(object = ., pattern = "^Hb[^(p)]", col.name = "percent.hb") %>%
    PercentageFeatureSet(object = ., pattern = "Pecam1|Pf4", col.name = "percent.plat") %>%
    PercentageFeatureSet(object = ., pattern = "Fos", col.name = "percent.cFos")

  feats = c('nFeature_RNA','nCount_RNA','percent.mt','percent.ribo','percent.hb','percent.plat','percent.cFos')
  seurat.object.plots = create.Plots()
  seurat.object.plots = add.Plots(plots=seurat.object.plots,
                                  plot={VlnPlot(seurat.object, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +NoLegend()},
                                  file_path = "VlnPlot")
  seurat.object.plots =  add.Plots(plots=seurat.object.plots,
                                  plot=FeatureScatter(seurat.object, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5),
                                  file_path = "FeatureScatter")

  seurat.object = seurat.object %>% subset(.,subset = nFeature_RNA > 500 & percent.mt < 5 & percent.hb <1) %>%
    `[`(.,!grepl("Malat1", rownames(.)),) %>%
    NormalizeData(., normalization.method = 'LogNormalize', scale.factor = 10000,verbose=F) %>%
    FindVariableFeatures(., selection.method = 'vst', nfeatures = 4000,verbose=F) %>%
    ScaleData(., features = rownames(.),verbose=F) %>%
    RunPCA(., features = VariableFeatures(.),verbose=F) %>%
    FindNeighbors(.,dims=1:30,verbose=F) %>%
    FindClusters(.,resolution = 0.6,verbose=F) # granularity, larger datasets with larger values

  #Step2: 确定参数，寻找最优pK值
  sweep.res.list <- paramSweep(seurat.object, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  ## (3) Homotypic Doublet Proportion Estimate -------------------------------------
  annotations <- seurat.object$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  (DoubletRate = ncol(seurat.object)*8*1e-6) #按每增加1000个细胞，双细胞比率增加千分之8来计算
  #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞
  (nExp_poi <- round(DoubletRate*length(seurat.object$seurat_clusters))) #最好提供celltype，而不是seurat_clusters。
  # 计算双细胞比例
  (nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)))

  ## (4)最后，使用确定好的参数鉴定Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat.object <- doubletFinder(seurat.object, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  seurat.object <- doubletFinder(seurat.object, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  # 使用nExp = nExp_poi和nExp = nExp_poi.adj,分别进行doublets鉴定，以便后续确定哪些细胞是Doublet-High Confidience
  ## Plot results ---------------------------------------------------------------------------
  tmp = seurat.object@meta.data
  seurat.object@meta.data[,"DF_hi.lo"] <- seurat.object[[colnames(tmp)[length(colnames(tmp))-2]]]
  seurat.object@meta.data$DF_hi.lo[which(seurat.object@meta.data$DF_hi.lo == "Doublet" & seurat.object[[colnames(tmp)[length(colnames(tmp))]]] == "Singlet")] <- "Doublet-Low Confidience"
  seurat.object@meta.data$DF_hi.lo[which(seurat.object@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
  ## 结果展示，分类结果在pbmc@meta.data中
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
                                plot=DimPlot(seurat.object, group.by ="DF_hi.lo",cols =c("black","red","gold")),
                                file_path = "DimplotDoublet")

  numDoublet = dim(seurat.object[,seurat.object$DF_hi.lo == "Doublet-High Confidience"])[2]
  numSinglet = dim(seurat.object)[2]
  print(paste("Disacard ",numDoublet,"high confident doublets."))
  print(paste("Remaining ",numSinglet,"cells."))
  seurat.object = seurat.object[,seurat.object$DF_hi.lo == "Singlet"]
  table(seurat.object@meta.data$DF_hi.lo)
  # Doublet-High Confidence Doublet-Low Confidence Singlet
  # 198 24  5041


  seurat.object= RunUMAP(seurat.object,dims = 1:30,verbose=F)

  seurat.object.plots=add.Plots(plots=seurat.object.plots,
            plot=DimPlot(seurat.object,label=T,pt.size=0.1),
            file_path = "Dimplot")
  # save(seurat.object,file=str_interp('${project}_qc.Rdata'))
  seurat.object.markers = FindAllMarkers(seurat.object, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25,test.use='roc',verbose = F)
  write.table(seurat.object.markers,file = str_interp("${project}_DEGs.txt"),sep = "\t",col.names = NA)
  seurat.object.markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC) -> top10
  markers = c("Snap25","Slc17a6","Gad1","Gad2","Hdc","Gpc5","Mobp","Dnah12","Col23a1","Inpp5d","Pdgfra","Dlc1","Adgrl4")
  gene_use = union(markers,top10$gene)
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
            plot={DoHeatmap(seurat.object,features = gene_use,slot = 'data')+ scale_fill_gradient2('legend name', low='white', high='red')},
            file_path = "Heatmap_raw")
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
            plot={DoHeatmap(seurat.object, features = gene_use, slot = "scale.data")+
              scale_fill_gradient2('legend name',low ='yellow',high='red',mid ='white')},
            file_path = "Heatmap_scaled")
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
            plot={DotPlot(seurat.object,features = markers)+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust =0.5))+
                ggtitle(str_interp('${project} neurons (before integration)'))} %T>% print(),
            file_path = "Dotplot")
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
            plot=VlnPlot(seurat.object,feats,group.by = 'seurat_clusters'),
            file_path = "VlnPlot_annotated")

  seurat.object=seurat.object%>%
  {
    cluster_annotations <- annotate_clusters(DotPlot(.,features = markers)$data)
    cluster_annotations$cluster <- as.character(cluster_annotations$cluster)
    identities_vector <- setNames(cluster_annotations$annotation, cluster_annotations$cluster)
    RenameIdents(., identities_vector)
  } %>%
  {
    # tryCatch(
    #   {
    #     excluded = subset(seurat.object,seurat_clusters%in%c("Debris","Doublet"))
    #     unique(excluded$RNA_snn_res.0.6)
    #     print(paste(sprintf("%.2f",100*as.numeric(dim(excluded)[2]/dim(seurat.object)[2])),"%"))
    #   }
    # )
    .$seurat_clusters = .@active.ident
    .[,!.$seurat_clusters%in%c("Debris","Doublet")]
  }

  saveRDS(seurat.object,str_interp('${project}.${format(Sys.Date(),format="%y%m%d")}.rds'))
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
                                plot={DotPlot(seurat.object,features = markers)+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust =0.5))+
                                    ggtitle(str_interp('${project} neurons (after integration)'))} %T>% print(),
                                file_path = "Dotplot")
  seurat.object.plots=add.Plots(plots=seurat.object.plots,
                plot=DimPlot(seurat.object,reduction = 'umap',label=T,pt.size=0.5),
                file_path = "Dimplot_annotated")
  draw.Plots(seurat.object.plots)
  end_time_cellranger_qc = Sys.time()
  print(end_time_cellranger_qc-start_time_cellranger_qc)
  gc()
}



annotate_clusters <- function(data_p, markers=
                                list(
                                    Neuron = c("Snap25"),
                                    Astrocyte = c("Gpc5"),
                                    Ependymal = c("Dnah12"),
                                    OPC = c("Gpc5","Pdgfra"),
                                    Oligo = c("Mobp"),
                                    Tanycyte = c("Col23a1"),
                                    Microglia = c("Inpp5d"),
                                    Pericyte = c("Dlc1"),
                                    VEC = c("Adgrl4")
                                    # ... 添加更多细胞类型和它们的标记 ...
                                    )
                            ) {
  # 创建一个空的数据框来存储簇注释
  cluster_annotations <- data.frame(cluster=unique(data_p$id), annotation=as.character(""), stringsAsFactors=FALSE)

  # 遍历每个簇
  for (cluster_id in cluster_annotations$cluster) {
    # 获取当前簇的数据子集
    cluster_data <- data_p[data_p$id == cluster_id, ]

    # 初始化标记记录器
    marker_counts <- setNames(rep(0, length(markers)), names(markers))

    # 检查每个细胞类型的标记基因
    for (cell_type_name in names(markers)) {
      marker_genes <- markers[[cell_type_name]] #markers[['Neuron']] = c('Snap25')
      if (length(marker_genes)==1){ # 对于只有一个标记基因的细胞类型，如果 出现两个
        gene_data <- cluster_data[cluster_data$features.plot == marker_genes, ]
        if (gene_data$pct.exp > 70) {
          marker_counts[cell_type_name] <- marker_counts[cell_type_name] + 1
        }
      } else{
        flag = 0
        for (gene in marker_genes){
          gene_data <- cluster_data[cluster_data$features.plot == gene, ]
          if (gene_data$pct.exp > 70) {
            flag = flag + 1
          }
        }
        if (flag == length(marker_genes)){
          marker_counts[cell_type_name] <- marker_counts[cell_type_name] + 1
          marker_counts["Astrocyte"] = 0 # because OPC and Astrocyte are in competition
        }
      }
    }
    print(marker_counts)
    # 应用注释规则
    high_markers <- names(marker_counts)[marker_counts > 0]
    print(high_markers)
    if (length(high_markers) == 1 ) { #&& all(cluster_data$pct_exp[!cluster_data$features.plot %in% unlist(markers[high_markers])] < 30)
      cluster_annotations$annotation[cluster_annotations$cluster == cluster_id] <- high_markers
    } else if (length(high_markers) > 1) {
      cluster_annotations$annotation[cluster_annotations$cluster == cluster_id] <- "Doublet"
    } else {
      cluster_annotations$annotation[cluster_annotations$cluster == cluster_id] <- "Debris"
      warning("[Incomplete classification] There is unclassified clusters. Check markers list!")
    }
  }
  # summary_df =data.frame(
  #   Cluster = unique(cluster_data$cluster),
  #   Annotation = cluster_annotations
  # )
  return(cluster_annotations)
}




# functions used in the following steps -----------------------------------
`qc` = function(data){
  return(data%>%
           subset(x = ., subset = percent.mt < 5& nFeature_RNA>200) %>%
           NormalizeData(., normalization.method = 'LogNormalize', scale.factor = 10000) %>%
           FindVariableFeatures(., selection.method = 'vst', nfeatures = 4000) %>%
           ScaleData(., features = rownames(.)) %>%
           RunPCA(., features = VariableFeatures(.)) %>%
           FindNeighbors(.,dims=1:30) %>%
           FindClusters(.,resolution = 0.6) %>%
           RunUMAP(.,dims = 1:30))
}


`visualization_before_annotation` = function(data){
  (proj = data@project.name)
  # p1 ----------------------------------------------------------------------
  p1 <- UMAPPlot(object = data, label = TRUE)
  pdf(file = str_interp('1_dim_reduction_cluster.${format(Sys.Date(),format="%y%m%d")}.pdf'),
      height = 8,width = 9);p1;dev.off()

  # p2 ----------------------------------------------------------------------
  markerGene <- c("Gpc5", #Astrocyte
                  "Dnah12", "Aqp4", #Ependymocyte
                  "Inpp5d", "Cx3cr1", #Microglia
                  "Rbfox3", "Snap25", #Neuron
                  "Slc17a6",
                  "Slc32a1", "Gad1", "Gad2",
                  "Hdc",
                  "Mog", "Mobp", #Oligo
                  "Pdgfra", #OPC
                  "Dlc1", #pericytes
                  "Col23a1", #Tancyte
                  "Adgrl4", "Adgrl") #VEC

  p2 <- DoHeatmap(object = data, features = markerGene,slot = "scale.data")+
    scale_fill_gradient2('legend name', low = 'yellow', high = 'red', mid = 'white')

  pdf(str_interp('2_heatmap_${proj}_cluster.${format(Sys.Date(),format="%y%m%d")}.pdf'),
      width = 10,height = 5);p2;dev.off()


  # p3 ----------------------------------------------------------------------
  p3 <- DoHeatmap(object = data, features = markerGene,
                  slot = "data") +
    scale_fill_gradient2('legend name', low = 'yellow', high = 'red', mid = 'white')
  pdf(str_interp('3_heatmap_${proj}_cluster_raw.${format(Sys.Date(),format="%y%m%d")}.pdf'),
      width = 10,height = 5);p3;dev.off()


  # t1 ----------------------------------------------------------------------
  data.markers <- FindAllMarkers(object = data, test.use = "wilcox",
                                  min.pct = 0.2,only.pos = T)

  write.table(data.markers,file = str_interp('${proj}.ClusterMarkers.txt'),sep = "\t",
              col.names = NA)

  top50 <- data.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

  # p4 ----------------------------------------------------------------------
  p4 <- DoHeatmap(object = data, features = c(markerGene, top50$gene),
                  slot = "scale.data")  +
    scale_fill_gradient2('legend name', low = 'yellow', high = 'red', mid = 'white')
  pdf(str_interp('4_heatmap_data_Cluster_top50.${format(Sys.Date(),format="%y%m%d")}.pdf'),
      width=20,height=80);p4;dev.off()

  # p5 ----------------------------------------------------------------------
  p5 <- DoHeatmap(object = data, features = c(markerGene, top50$gene),
                  slot = "data") +
    scale_fill_gradient2('legend name', low = 'yellow', high = 'red', mid = 'white')
  pdf(str_interp('5_heatmap_data_Cluster_top50_raw.${format(Sys.Date(),format="%y%m%d")}.pdf'),
      width=20,height=80);p5;dev.off()

  # t2 ----------------------------------------------------------------------
  data.markers <- read.table(str_interp('${proj}.ClusterMarkers.txt'), header = T)
  top50 <- data.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

  (table(data@active.ident))

  p_return = DotPlot(data,markerGene)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust =0.5))+
    ggtitle(str_interp('snRNA-seq ${proj} neurons (after integration)'))

  return(list(p_return,data.markers))
}
