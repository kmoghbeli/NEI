# ## Marker Genes
# 
marker_genes_L1 <- c("Epcam", "Ptprc")

# Source: https://www.cellsignal.com/pathways/immune-cell-markers-mouse
# https://www.thermofisher.com/us/en/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work/macrophage-cell-overview.html
marker_genes_L2 <- c("Epcam", "Ptprc", 
                     "Itgam", # Cd11b - Myeloid marker (mouse)
                     "Adgre1", # F4/80, 
                     "Cd74", # MHC-II mouse marker (used by Renthal 2022 to identify immune cells in TG)
                     "Pax5",  # B 
                     "Ighd", "Cd27", # Naive (mouse) B cell markers (IgD+, CD27-)
                     "Cd83",  # DCs
                     "Cd14", "Cd68",  # Macs - note that Cd16 never comes up 
                     "Ptgs2", "Irf5", "Nos2",  # Mouse M1 Mac Markers 
                     "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a",   # T
                     "Cd44", "Sell", # Naive (mouse) T cell markers (CD44 -, CD62L +)
                     "Fcgr3", "Il17a", "Ifng",    # delta-gamma T
                     "Klrb1c", "Klrk1", "Gzma", "Gzmb", "Prf1", # NK
                     "Itga2", "Ncam1",  #NK-T
                     "Pdgfra",  # Fibroblasts
                     "Vwf")     # Endothelial

marker_genes_tg <- 
  c(#"Itgam", # Cd11b - Myeloid marker (mouse)
    #"Adgre1", # F4/80, 
    "Pax5",  # B 
    "Ighd", "Cd27", # Naive (mouse) B cell markers (IgD+, CD27-)
    "Cd3d", # T
    "Klrb1c", "Prf1", "Klrk1", "Gzma", "Gzmb",  # NK 
    #"Itga2", "Ncam1",  #NK-T
    "Cd83",  # DCs
    "Cd14", "Cd68",  # Macs - note that Cd16 never comes up 
    #"Itgax", # DCs
    "Ly6c1", 
    #"Cd74", # MHC-II mouse marker (used by Renthal 2022 to identify immune cells in TG)
    "Ptgs2", "Irf5", "Nos2",  # Mouse M1 Mac Markers 
    # "Stat1", "Retnla",  # Mouse M1 Mac Markers (less helpful)
    #"Il12a", "Il23a", "Cd163",  # M1 vs M2 (M1: IL-12 and IL23 high with CD163 neg and M2 the opposite)
    "Cd163",  # M2
    #"Arg1", # M2a
    "Socs3", "Cd86", # M2b
    "Ccr2", "Slamf6",   #M2c
    # "Tlr1", "Tlr8", "Scarb1", #M2c (less helpful)
    "Vegfa",    # M2d, 
    "Cx3cr1"  # Tissue-res Mac
  )

# CD64 - mouse cornea macrophages (maybe activation marker, i.e., low in M0)
# CCR2 maybe also an activation marker
# Ly6c - inflammatory monocytes
# Cd68
# Nature Communications: Qie et al ("integrated proteomic... macrophages") paper and supplemental data: macrophage.mousprotein.cn

# 
# Macrophage markers:
#   - A total of 13 surface markers, including CD45, F4/80, CD11b, CD117, Siglec-F, CD11c, Cx3cr1, MHCII, CD64, CD115, CD24, B220, and Ly6g, were used for cell sorting according to the published literature (Sup- plementary Fig. 1a and Supplementary Tables 1 and 2)1,8,29
# - F4/80 and CD11b bright and dim phenotypes were used to distinguish tissue- resident and recruited populations in the lung, liver, and spleen
# - Proteins that are critical for macrophage develop- ment or serve as identity markers, such as Sfpi1 (PU.1), Itgam (CD11b), Adgre1 (F4/80), Myd88, and Mertk, were found to be highly or mod- erately abundant (Supplementary Fig. 3).
# - Proteins with high mRNA correlation: TFs: Stat5a, Stat6, Runx3 // Adaptor: Myd88, Fadd // Receptor: Fcrls, Tlr5, Tlr12, Clec4f, Nod1, Nlrp6, Ccr2, Itga6, Vcam1 // Cytokine: Cxcl2
# 
# (References: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9940311/; 
#   Macrophages:
#     - F4/80, MHCII, CD68
#   
#   M2: Arg-1
#   iNOS: M1
#   Gal-3: M2
#   TREM2: M2 (but also myeloid lineage cells)
#   Mature Tissue M2: CD163
#   CD169: secondary lymphoid
#   CD206: DCs + Macs (M2)
#   CD9: anti-inflamm macs
#   CX3CR1: tissue-res macs
#   Lyve1: macrophages in vasculature
#   
#   https://www.rndsystems.com/resources/cell-markers/immune-cells/macrophages/m2d-macrophage-activation-state-markers
#   
#   https://www.thermofisher.com/us/en/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work/macrophage-cell-overview.html
#   
#   https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
#     
#     https://www.bio-rad-antibodies.com/macrophage-polarization-minireview.html?JSESSIONID_STERLING=C4A0854D67915175B7810AAE3BD9BFA9.ecommerce2&evCntryLang=US-en&cntry=US&thirdPartyCookieEnabled=true
#   
#   M1: 
#     COX2, iNOS, IRF5, STAT1, CD80, CD86, CD163 (neg), CD32, CD16, MHCII high, IL-12 high, IL-23 high, IL-10 low
#   - IFNy, IL-1, IL-6, IL-12, IL-23, TNFa, CD14, CD16, CD32, CD64, CD68, CD80, CD86, CD204, CD369, Ly-6C, Mer, MHCII, IRF5
#   
#   M2: IL-12 and IL-23 low, IL-10 and IL-1RA high, IL-1 low, IL-6 low, TNFa low
#   - CD206, CD163, CD209, FIZZ1, Ym1/2 
#   - IDO, Arginase, IL-10, TGF-b, YM1, CD14, CD115, CD163, CD204, CD206, CD209, CSF1R, FceR1, IRF4, RELM-a, STAT6, 
#   
#   M2a: Arginase1, IRF4, PPARy, STAT6, CD163+, CD301, CXCR1/2, Dectin-1+, MHCII low, Fizz1, Arg-1, Ym1/2
#   
#   M2b: SOCS3, CD86, SPHK1/2
#   M2c: TLR8, CCR2, SLAM, TLR1, SR-BI, CD163
#   M2d: VEGF
#   