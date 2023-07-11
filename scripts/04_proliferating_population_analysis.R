###############################
##### Proliferating Cells #####
###############################
library(Seurat)
library(remotes)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(monocle3)
library(dplyr)
library(ggpubr)
library(clusterProfiler)

# Locate the Proliferaitng Cells on UMAP #
# Figure 4 C
DimPlot(kidneys,
		  cols = c("#999999", "#299CA8", "#999999",
		  			"#999999", "#999999","#999999",
		  			"#999999", "#999999", "#999999", 
		  			"#999999",
		  			"#999999"))


# Analyze this cluster separately
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(kidneys) <- "RNA"
PT_Prolif <- subset(kidneys, subset = cell_types == "PT_Prolif")
# cell cycle scoring of the proliferating cells
PT_Prolif <- CellCycleScoring(PT_Prolif, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
PT_Prolif <- NormalizeData(PT_Prolif)
PT_Prolif <- FindVariableFeatures(PT_Prolif, selection.method = "vst")
PT_Prolif <- ScaleData(PT_Prolif, features = rownames(PT_Prolif))
PT_Prolif <- RunPCA(PT_Prolif)
PT_Prolif <- RunUMAP(PT_Prolif, dims = 1:30)
DimPlot(PT_Prolif)

# Look at how many cells are in each condition
table(PT_Prolif$Condition)
# kidney_1st_untreated_rt kidney_1st_treated kidney_2nd_untreated_rt      kidney_2nd_treated 
# 2                     122                      86                     323 

DefaultAssay(PT_Prolif) <- "SCT"

# Figure 4 D
FeaturePlot(PT_Prolif,
				reduction = "umap",
				features = c(
					"STMN1",
					"SLC34A1",
					"CD3E"),
				split.by = "Condition")

# Test how many cells are T cells
PT_Prolif_T <- subset(PT_Prolif, CD3G > 1 | CD3E > 1)
table(PT_Prolif_T$Condition)
# kidney_1st_untreated_rt      kidney_1st_treated kidney_2nd_untreated_rt 
# 0                       0                      68 
# kidney_2nd_treated 
# 0 

# Look at cell cycle scoring
# Figure 4 G
DimPlot(PT_Prolif,
		  reduction = "umap") + 
			theme(axis.title = element_blank()) 

# Look at organism assignment
# Figure 4 H
DimPlot(PT_Prolif,
		  group.by = "origin") + theme(axis.title = element_blank(),
		  													plot.title = element_blank())

# Visualize top markers across cell cycle identities
# Figure 4 I
RidgePlot(PT_Prolif, 
			 features = c("STMN1", "PCNA", "HMGN2", "HMGB2"), 
			 ncol = 2)
