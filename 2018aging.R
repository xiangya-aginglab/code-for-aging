library(data.table)
library(Seurat)
raw.data <- read.csv ('./2018aging/Skin_6Control_rawUMI.csv',header=T, row.names = 1)
raw.data = aging_2018
Sys.time()
dim(raw.data)
raw.data[1:4,1:4]
head(colnames(raw.data)) 
# Load metadata 
metadata <- read.csv("./2018aging/Skin_6Control_Metadata.csv", row.names=1, header=T)
head(metadata) 

# Find ERCC's, compute the percent ERCC, and drop them from the raw data.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
fivenum(percent.ercc)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,]
dim(raw.data) 


main_tiss <- CreateSeuratObject(counts = raw.data)
# add rownames to metadta 
row.names(metadata) <- metadata$cell_id
# add metadata to Seurat object 
main_tiss <- AddMetaData(object = main_tiss, metadata = metadata)
main_tiss <- AddMetaData(object = main_tiss, percent.ercc, col.name = "percent.ercc")
# Head to check
head(main_tiss@meta.data)

# Calculate percent ribosomal genes and add to metadata
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = main_tiss@assays$RNA@data), value = TRUE)
percent.ribo <- Matrix::colSums(main_tiss@assays$RNA@counts[ribo.genes, ])/Matrix::colSums(main_tiss@assays$RNA@data)
fivenum(percent.ribo)
main_tiss <- AddMetaData(object = main_tiss, metadata = percent.ribo, col.name = "percent.ribo")
main_tiss

VlnPlot(object = main_tiss, 
        features = "HOXC10", 
        group.by = "orig.ident")

VlnPlot(aging_2018, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
plot1 <- FeatureScatter(main_tiss, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(main_tiss, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

main_tiss <- subset(main_tiss, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mito < 20)
main_tiss <- NormalizeData(main_tiss)
main_tiss <- FindVariableFeatures(main_tiss, selection.method = "vst", nfeatures = 2000)
main_tiss <- ScaleData(main_tiss)
main_tiss <- RunPCA(main_tiss, features = VariableFeatures(object = main_tiss))
DimHeatmap(main_tiss, dims = 1:20, cells = 500, balanced = TRUE)

main_tiss <- JackStraw(main_tiss, num.replicate = 100)
main_tiss <- ScoreJackStraw(main_tiss, dims = 1:30)
JackStrawPlot(main_tiss, dims = 1:20)
ElbowPlot(main_tiss)
main_tiss <- FindNeighbors(main_tiss, dims = 1:19)
main_tiss <- FindClusters(main_tiss, resolution = 0.2)

library(harmony)
library(clustree)
seuratObj <- RunHarmony(main_tiss, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:19, 
                     reduction = "harmony")

DimPlot(seuratObj,reduction = "umap",label=T, group.by = 'orig.ident') 


sce=seuratObj
main_tiss <- FindNeighbors(main_tiss,
                          dims = 1:19) 

#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.6,0.8,1)) {
  main_tiss=FindClusters(main_tiss, #graph.name = "CCA_snn", 
                        resolution = res, algorithm = 1)
}
colnames(main_tiss@meta.data)
apply(sce_all1@meta.data[,grep("RNA_snn",colnames(sce_all1@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce_all1, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce_all1, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce_all1, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)

p1_dim=plot_grid(ncol = 3, DimPlot(sce_all1, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)


p2_tree=clustree(main_tiss@meta.data, prefix = "RNA_snn_res.")
p2_tree
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")

sce_all1 <- FindNeighbors(sce_all1, dims = 1:30)
sce_all1 <- FindClusters(sce_all1, resolution = 0.2)
main_tiss <- RunUMAP(main_tiss, dims = 1:19)

main_tiss@meta.data = main_tiss@meta.data %>%
  mutate(group = case_when(
    orig.ident == 'SC1control' | orig.ident == 'SC4control' |orig.ident == 'SC18control' | orig.ident == 'SC33control'~ "old",
    orig.ident == 'SC32control' | orig.ident == 'SC34control' ~ "young"
   )
  )

plot4=DimPlot(main_tiss, reduction = "umap",group.by="group")
plot3=DimPlot(main_tiss, reduction = "umap")
plot3+plot4
main_tiss = RunTSNE(main_tiss, dims = 1:19)
embed_tsne <- Embeddings(main_tiss, 'tsne')

plot1 = DimPlot(main_tiss, reduction = "tsne", label=T)
plot2 = DimPlot(main_tiss, reduction = "tsne", group.by='group') 
plot1 + plot2

saveRDS(main_tiss,file = 'aging_2018.rds')


plot1 = FeaturePlot(sce, features = c("CPE",'MT2A','FMOD','FST'))
plot2 = DimPlot(sce, reduction = "umap", group.by='group')
Plot3 = DimPlot(main_tiss, reduction = "umap", group.by = "celltype",label = T)
 plot1 + plot2

head(main_tiss@meta.data)

##MARKERS
library(ggplot2) 
genes_to_check = c('AIF1','HLA-DRA',  # macrophage 
                   'MPZ',  # schwann cells
                   'PMEL',# melanocyte
                   'CD3D','CXCR4',## lymphocyte
                   'CLDN5','VWF', # endothelial cell
                   'DES','ACTG2', # smooth muscle cells 
                   'COL1A1','APOD','POSTN', 'IGFBP3', 'SFRP2', 'COMP', 'CXCL12', ## fibo 
                   'KRT14' ,'KRT1','KRT2','KRT16', ## keratinocytes 
                   'LOR', #cornified envelope cells
                   'RGS5','ACTA2', 'TPM2', ## pericyte 
                  'IGJ', #B cell
                   'HMGB2',#proliferating keratinocyte
                   'SCGB1B2P','MUCL1')#secretory cells/epithelial cells

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
p <- DotPlot(main_tiss, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p 
ggsave('check_last_markers.pdf',height = 11,width = 11)

# 需要自行看图，定细胞亚群： 
celltype=data.frame(ClusterID=0:14,
                    celltype= 0:14 ) 
celltype[celltype$ClusterID %in% c(2,5,6),2]='keratinocytes' # CD19, CD79A, MS4A1
celltype[celltype$ClusterID %in% c(3),2]='endothelial cell' # TOP2A, MKI67 
celltype[celltype$ClusterID %in% c(4),2]='pericyte' # CD8
celltype[celltype$ClusterID %in% c(7),2]='macrophage' # CD4
celltype[celltype$ClusterID %in% c(8),2]='lymphocyte' # 'C1QA', 'C1QB'
celltype[celltype$ClusterID %in% c(0,1),2]='FB' # 'S100A9', 'S100A8'
celltype[celltype$ClusterID %in% c(9),2]='secretory cells/epithelial cells' # 'KLRB1','NCR1'
celltype[celltype$ClusterID %in% c(10),2]='proliferating keratinocyte' # 'CD1E','CD1C'
celltype[celltype$ClusterID %in% c(11),2]='smooth muscle cells' #'JCHAIN','LILRA4'
celltype[celltype$ClusterID %in% c(12),2]='melanocyte' 
celltype[celltype$ClusterID %in% c(13),2]='cornified envelope cells'
celltype[celltype$ClusterID %in% c(14),2]='B cell'

head(celltype)
celltype
table(celltype$celltype)
main_tiss@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  main_tiss@meta.data[which(main_tiss@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(main_tiss@meta.data$celltype)

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 

p <- DotPlot(main_tiss, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()  +th

p
ggsave(plot=p, filename="./human_aging/qc/check_marker_by_celltype.pdf",
       width = 7 ,height = 8)

DotPlot(main_tiss, features = genes_to_check,
        assay='RNA' ,group.by = 'celltype' ) + 
  coord_flip()+ scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  NULL

table(sce_all1@meta.data$celltype,sce_all1@meta.data$RNA_snn_res.0.3)

library(patchwork)
p_all_markers=DotPlot(main_tiss, features = genes_to_check,
                      assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_umap=DimPlot(main_tiss, reduction = "umap", group.by = "celltype",label = T)
p_all_markers+p_umap
ggsave('./human_aging/qc/markers_umap_by_celltype.pdf',width = 12,height = 6)
phe=sce_all1@meta.data
save(phe,file = './human_aging/qc/phe-by-markers.Rdata')


sce_all1
table(Idents(main_tiss))  
Idents(main_tiss)=main_tiss$celltype
sce.markers <- FindAllMarkers(object = main_tiss, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
pro='markers'
write.csv(sce.markers,file=('./2018aging/markers_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce_all1,top10$gene,size=3)
ggsave(filename='./human_aging/qc/markers_sce.markers_heatmap.pdf',height = 15)

library(dplyr) 
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(main_tiss,top3$gene,size=3)
ggsave(filename='./2018aging/markers_sce.markers_heatmap.pdf',height = 8,width = 10)
p <- DotPlot(main_tiss, features = unique(top3$gene),
             assay='RNA'  )  + coord_flip()+th

p
ggsave(filename='./2018aging/DotPlot_check_top3_markers_by_clusters.pdf',height = 8,width = 10)
save(sce.markers,file = './2018aging/sce.markers.Rdata')


# 
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce_all1, features = feats, pt.size = 0.01, ncol = 2) + 
  NoLegend()
p1
library(ggplot2) 
ggsave(filename="./human_aging/qc/Vlnplot1.pdf",plot=p1)

feats <- c("percent.mt")
p2=VlnPlot(sce_all1,  features = feats, pt.size = 0.01, ncol = 3, same.y.lims=T) + 
  NoLegend()
p2	
ggsave(filename="Vlnplot2.pdf",plot=p2)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", 
                  pt.size = 0.5)
ggsave(filename="Scatterplot.pdf",plot=p3)

saveRDS(sce_all1,file="sce.all.final.rds")



library(Seurat)

ct = names(table(main_tiss$celltype)) 

table(main_tiss@meta.data$group)
head(main_tiss@meta.data)
table(Idents(main_tiss))
#统计各个亚群PD-DLB vs HC 差异基因
library(DESeq2)
all_sig_markers = lapply(ct, function(x){
  # x = ct[1]
  print(x)
  markers <- FindMarkers(main_tiss, subset.ident = x,
                         ident.1 = 'old',ident.2 = 'young',
                         group.by = 'group', logfc.threshold = 0.15, test.use = 'wilcox')
  markers_sig <- subset(markers, p_val_adj < 0.05)
  return(markers_sig)
})

markers_fb <- FindMarkers(main_tiss, subset.ident = 'FB',
                       ident.1 = 'old',ident.2 = 'young',
                       group.by = 'group', logfc.threshold = 0.15, test.use = 'wilcox')

write.csv(markers_fb,'fb_diff_markers.csv')
VlnPlot(main_tiss, features = c("CPE",'MT2A','FMOD','FST'), split.by = 'group')

+?FindMarkers
names(all_sig_markers)=ct
save(all_sig_markers,file="DEGs_in_all_celltype.Rdata")

## UpSet图
# 制作UpSet输入数据inputlist列表
library(UpSetR)
inputlist=lapply(ct, function(x){
  y=rownames(all_sig_markers[[x]])
  return(y)
})

names(inputlist)=ct

#按照基因数量从大到小排序，方便后面upset图横向从大到小排列。
#mylist=inputlist[order(sapply(inputlist, length),decreasing = T)]

pdf('UpSet.pdf',width = 24,height = 16)
p<-upset(fromList(inputlist),nsets=10,order.by = 'freq',
         #matrix.color = "blue", #点的颜色
         #main.bar.color = rainbow(59), #竖bar颜色
         #sets.bar.color = rainbow(9), #横bar颜色
         shade.alpha = 0.4, #点图中阴影深浅
         mainbar.y.label = "Distinct DEGs",
         sets.x.label = "Total DEGs"
)
p
dev.off()

#markers <- subset(FindMarkers(scRNA, subset.ident = 'CD4 T',ident.1 = 'PD-DLB',ident.2 = 'HC', group.by="csfgroup"),p_val_adj < 0.05)

markers=all_sig_markers[['FB']]


markers_fb$change = ifelse(markers_fb$avg_log2FC>0,'Up', 
                        ifelse(markers_fb$avg_log2FC<0,'Down','Stable'))
markers_fb$symbol=rownames(markers_fb)

library(ggplot2)
p <- ggplot(data = markers, 
            aes(x = avg_log2FC, 
                y = -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.005),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

cg=c("HOXC10",'FZD6')
for_label=markers[cg,]

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave("CD4_Volcano.pdf")

head(aging_2018@meta.data)
sce_fb=(aging_2018[,(aging_2018$celltype %in% c( 'FB' ))])

table(sce_fb@meta.data$celltype)
sce_fb=CreateSeuratObject(counts = sce_fb@assays$RNA@counts,
                          meta.data = sce_fb@meta.data) 
sce_fb <- NormalizeData(sce_fb, normalization.method =  "LogNormalize",  
                        scale.factor = 1e4)
GetAssay(sce_fb,assay = "RNA")
sce_fb <- FindVariableFeatures(sce_fb, 
                               selection.method = "vst", nfeatures = 2000)  
sce_fb <- ScaleData(sce_fb) 
sce_fb <- RunPCA(object = sce_fb, pc.genes = VariableFeatures(sce_fb)) 

DimHeatmap(sce_fb, dims = 1:20, cells = 100, balanced = TRUE)
ElbowPlot(sce_fb) 
sce_fb <- JackStraw(sce_fb, num.replicate = 100)
sce_fb <- ScoreJackStraw(sce_fb, dims = 1:20)
JackStrawPlot(sce_fb, dims = 1:20)


sce_fb <- FindNeighbors(sce_fb, dims = 1:20)
sce_fb <- FindClusters(sce_fb, resolution = 0.2)
head(sce_fb@meta.data)
table(sce_fb@meta.data$RNA_snn_res.0.2)  

set.seed(123)
sce_fb <- RunTSNE(object = sce_fb, dims = 1:20, do.fast = TRUE)
DimPlot(sce_fb,reduction = "tsne",label=T)


sce_fb <- RunUMAP(object = sce_fb, dims = 1:20, do.fast = TRUE)
DimPlot(sce_fb,reduction = "umap",label=T) 
head(sce@meta.data)
DimPlot(sce,reduction = "umap",label=T,group.by = 'orig.ident')

head(sce_fb[["RNA"]]@counts)
library(Seurat)
library(harmony)
seuratObj <- RunHarmony(sce_fb, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:20, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=T)  


sce_fb=seuratObj
sce_fb <- FindNeighbors(sce_fb, reduction = "harmony",
                     dims = 1:20)
sce.all =sce
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4,0.5,0.8,1)) {
  sce_fb=FindClusters(sce_fb, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
colnames(sce_fb@meta.data)
apply(sce.all@meta.data[,grep("RNA_snn",colnames(sce.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)

library(clustree)
p2_tree=clustree(sce_fb_1@meta.data, prefix = "RNA_snn_res.")
p2_tree
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")
sel.clust = "RNA_snn_res.0.2"
sce_fb <- SetIdent(sce_fb, value = sel.clust)
table(sce_fb@active.ident) 
saveRDS(sce.all, "sce.all_int.rds")
DimPlot(sce.all, reduction = "umap", group.by = "seurat_clusters",label = T) 
DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T) 
VlnPlot(sce.all, features = c('FMO1'), split.by = 'group')

sce.all.marker <- FindAllMarkers(sce_fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sce.all.marker %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

library(ggplot2) 
genes_to_check = c('AIF1','HLA-DRA',  # macrophage 
                   'MPZ',  # schwann cells
                   'PMEL',# melanocyte
                   'CD3D','CXCR4',## lymphocyte
                   'CLDN5','VWF', # endothelial cell
                   'DES','ACTG2', # smooth muscle cells 
                   'COL1A1','APOD','POSTN', 'IGFBP3', 'SFRP2', 'COMP', 'CXCL12', ## fibo 
                   'KRT14' ,'KRT1','KRT2','KRT16', ## keratinocytes 
                   'LOR', #cornified envelope cells
                   'RGS5','ACTA2', 'TPM2', ## pericyte 
                   'IGJ', #B cell
                   'HMGB2',#proliferating keratinocyte
                   'SCGB1B2P','MUCL1')#secretory cells/epithelial cells

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
p <- DotPlot(main_tiss, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p 
ggsave('check_last_markers.pdf',height = 11,width = 11)

# 需要自行看图，定细胞亚群： 
celltype=data.frame(ClusterID=0:3,
                    celltype= 0:3 ) 
celltype[celltype$ClusterID %in% c(0),2]='A' # CD19, CD79A, MS4A1
celltype[celltype$ClusterID %in% c(1),2]='B' # TOP2A, MKI67 
celltype[celltype$ClusterID %in% c(2),2]='C' # CD8
celltype[celltype$ClusterID %in% c(3),2]='other' # CD8

head(celltype)
celltype
table(celltype$celltype)
main_tiss@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 

p <- DotPlot(sce.all, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()  +th

p
ggsave(plot=p, filename="./human_aging/qc/check_marker_by_celltype.pdf",
       width = 7 ,height = 8)

DotPlot(main_tiss, features = genes_to_check,
        assay='RNA' ,group.by = 'celltype' ) + 
  coord_flip()+ scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  NULL

table(sce_all1@meta.data$celltype,sce_all1@meta.data$RNA_snn_res.0.3)

library(patchwork)
p_all_markers=DotPlot(main_tiss, features = genes_to_check,
                      assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_umap=DimPlot(main_tiss, reduction = "umap", group.by = "celltype",label = T)
p_all_markers+p_umap
ggsave('./human_aging/qc/markers_umap_by_celltype.pdf',width = 12,height = 6)
phe=sce_all1@meta.data
save(phe,file = './human_aging/qc/phe-by-markers.Rdata')


sce_all1
table(Idents(sce.all))  
Idents(sce.all)=sce.all$celltype
sce.markers <- FindAllMarkers(object =sce.all, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
pro='markers'
write.csv(sce.markers,file=('./2018aging/markers_sce.markers.csv'))
library(dplyr) 
top10 <- sce.all.marker %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce_fb,top10$gene,size=3)
ggsave(filename='./human_aging/qc/markers_sce.markers_heatmap.pdf',height = 15)

library(dplyr) 
top3 <- sce.all.marker %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(sce.all,top3$gene,size=5)
ggsave(filename='./2018aging/markers_sce_fb.markers_heatmap.pdf',height = 8,width = 10)
p <- DotPlot(sce_fb, features = unique(top3$gene),
             assay='RNA'  )  + coord_flip()+th

p
ggsave(filename='./2018aging/DotPlot_check_top5_markers_by_fb_clusters.pdf',height = 8,width = 10)
save(sce.all.marker,file = './2018aging/sce.fb_marker.Rdata')



sce <- FindClusters(sce, resolution = 0.6)
table(sce@meta.data$seurat_clusters)
DimPlot(sce,reduction = "umap",label=T) 
ggsave(filename = 'harmony_umap_recluster_by_0.1.pdf') 
DimPlot(sce.all,reduction = "umap",label=T) 
ggsave(filename = 'harmony_umap_sce_recluster_by_orig.ident.pdf') 

saveRDS(sce_fb,file = './2018aging/sce_fb_1.rds')

sce_fb = sce_fb_1

p1 = DimPlot(sce.all,reduction = "umap",label=T)
p2 = DimPlot(sce.all,reduction = "umap",label=T,group.by = 'group')
p1 + p2 +p3
p3 = FeaturePlot(sce.all, features = c("CPE",'MT2A','FMOD','FST'))

VlnPlot(sce.all,features = c("CPE",'MT2A','FMOD','FST'), ncol = 2, split.by = 'group')

sce_fb_1
DimPlot(sce_fb,reduction = "umap",label=T, group.by = 'group')


FeaturePlot(sce_fb, features = c("CPE",'MT2A','FMOD','FST'))


library(Seurat)

ct = names(table(sce.all$celltype)) 
ct = c('A','B','C')

table(main_tiss@meta.data$group)
head(main_tiss@meta.data)
table(Idents(main_tiss))
#统计各个亚群PD-DLB vs HC 差异基因
library(DESeq2)
all_sig_markers = lapply(ct, function(x){
  # x = ct[1]
  print(x)
  markers <- FindMarkers(sce.all, subset.ident = x,
                         ident.1 = 'old',ident.2 = 'young',
                         group.by = 'group', logfc.threshold = 0.15, test.use = 'wilcox')
  markers_sig <- subset(markers, p_val_adj < 0.05)
  return(markers_sig)
})

markers_A = all_sig_markers [['A']]
markers_C = all_sig_markers [['C']]
markers_B = all_sig_markers [['B']]

write.csv(markers_fb,'fb_diff_markers.csv')
VlnPlot(main_tiss, features = c("CPE",'MT2A','FMOD','FST'), split.by = 'group')

+?FindMarkers
names(all_sig_markers)=ct
save(all_sig_markers,file="DEGs_in_all_celltype.Rdata")

## UpSet图
# 制作UpSet输入数据inputlist列表
library(UpSetR)
inputlist=lapply(ct, function(x){
  y=rownames(all_sig_markers[[x]])
  return(y)
})

names(inputlist)=ct

#按照基因数量从大到小排序，方便后面upset图横向从大到小排列。
#mylist=inputlist[order(sapply(inputlist, length),decreasing = T)]

pdf('UpSet.pdf',width = 24,height = 16)
p<-upset(fromList(inputlist),nsets=10,order.by = 'freq',
         #matrix.color = "blue", #点的颜色
         #main.bar.color = rainbow(59), #竖bar颜色
         #sets.bar.color = rainbow(9), #横bar颜色
         shade.alpha = 0.4, #点图中阴影深浅
         mainbar.y.label = "Distinct DEGs",
         sets.x.label = "Total DEGs"
)
p
dev.off()

#markers <- subset(FindMarkers(scRNA, subset.ident = 'CD4 T',ident.1 = 'PD-DLB',ident.2 = 'HC', group.by="csfgroup"),p_val_adj < 0.05)

markers=all_sig_markers[['FB']]


markers_C$change = ifelse(markers_C$avg_log2FC>0,'Up', 
                           ifelse(markers_C$avg_log2FC<0,'Down','Stable'))
markers_C$symbol=rownames(markers_C)
write.csv(markers_C, './2018aging/markers_C.csv')



library(ggplot2)
p <- ggplot(data = markers, 
            aes(x = avg_log2FC, 
                y = -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.005),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

cg=c("HOXC10",'FZD6')
for_label=markers[cg,]

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave("CD4_Volcano.pdf")

plot1 <- FeatureScatter(sce.all, feature1 = "CPE", feature2 = "MMP2", group.by = 'group')
plot2 <- FeatureScatter(sce.all, feature1 = "CPE", feature2 = "MMP11")
plot1 + plot2

head(yanda@meta.data,30)
samplename <-  as.vector(lapply(rownames(yanda@meta.data), function(x){
  as.character(strsplit(x,'-')[[1]][2])
}))
samplename <- as.data.frame(unlist(samplename))
rownames(samplename) <- colnames(yanda@meta.data)
colnames(samplename)[1] <- 'sample'
head(samplename)

ElbowPlot(sce_fb_1)

head(sce_fb_1@meta.data)
table(sce_fb_1@meta.data$RNA_snn_res.0.2)

VlnPlot(sce_fb_1, features = c("APOD",'SFRP2'), ncol = 3)
DimPlot(sce_fb_1)
                                  sce_fb_1 <- JackStraw(sce_fb_1, num.replicate = 100,dims = 40)
sce_fb_1 <- ScoreJackStraw(sce_fb_1, dims = 1:40)
JackStrawPlot(sce_fb_1, dims = 1:20)
