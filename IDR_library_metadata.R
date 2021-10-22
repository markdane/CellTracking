library(tidyverse)
library(readxl)

data_path <- "../../../images/"
data_path <- "../images/"

plateIDs <- c("AU00601","AU00602","AU00701","AU00702","AU00801","AU00802","AU00901","AU00902","AU01001","AU01002","AU01101","AU01102")
data_paths <- paste0(data_path, plateIDs)
library_filepath <- "../metadata/AU565_library.csv"

if(!file.exists(library_filepath)){
  #read file names into a dataframe
  df <- map(data_paths, dir, pattern = "m.tif", recursive = TRUE) %>%
    unlist() %>%
    bind_cols() %>%
    rename("filename" = "...1") %>%
    mutate(filename = str_remove(filename, ".*/"),
           plate_name = str_remove_all(filename, "_.*"),
           channel_name = str_extract(filename, "_[RPG]_"),
           channel_name = str_remove_all(channel_name,"_"),
           well = str_extract(filename, "_[ABCD123456]*_"),
           well = str_remove_all(well,"_"),
           row = str_extract(well, "[[:alpha:]]"),
           column = str_extract(well, "[[:digit:]*]"),
           column = sprintf("%02d", as.integer(column)),
           well = paste0(row,column),
           field = str_extract(filename, "_[1234]_"),
           field = str_remove_all(field,"_"),
           timepoint = str_extract(filename, "[[:alnum:]]*.tif"),
           timepoint = str_remove_all(timepoint,".tif"))
  write_csv(df, "../metadata/IDR/AU565_library.csv")
} else {
  df <- read_csv(library_filepath)
}

#Lapatinib
# 1.7 InChkey
# BCFGMOOMADDAQU-UHFFFAOYSA-N
# 
# 1.8 Canonical Smiles
# CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl
# 
# 1.9 Isomers Smiles
# CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl
# Doxorubicin
# 1.8 InChkey
# MWWSFMDVAYGXBV-BXPPNZEESA-N
# 
# 1.9 Canonical Smiles
# CC1C(C(CC(O1)OC2CC(CC3=C(C4=C(C(=C23)O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O.Cl
# 
# 1.10 Isomers Smiles
# C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C(C4=C(C(=C23)O)C(=O)C5=
#                                                     C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O
# 
#Palbociclib
# 1.8 InChkey
# HJQCAEDIUJXGCQ-UHFFFAOYSA-N
# 
# 1.9 Canonical Smiles
# CC1=C(C(=O)N(C2=NC(=NC=C12)Cl)C3CCCC3)Br

#Paclitaxol
# 1.8 InChkey
# RCINICONZNJXQF-HLPYKUKBSA-N
# 
# 1.9 Canonical Smiles
# CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C
# 
# 1.10 Isomers Smiles
# CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C
#                                                                       [C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=
#                                              CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C

#Gemcitabine
# 1.8 InChkey
# OKKDEIYWILRZIA-OSZBKLCCSA-N
# 
# 1.9 Canonical Smiles
# C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)(F)F
# 
# 1.10 Isomers Smiles
# C1=CN(C(=O)N=C1N)[C@H]2C([C@@H]([C@H](O2)CO)O)(F)F


drug_metadata <- tibble(name = c("Laptinib", "Doxorubicin", "Palbociclib", "Paclitaxel", "Gemcitabine"),
                        InChlkey = c("BCFGMOOMADDAQU-UHFFFAOYSA-N", "MWWSFMDVAYGXBV-BXPPNZEESA-N", "HJQCAEDIUJXGCQ-UHFFFAOYSA-N", "RCINICONZNJXQF-HLPYKUKBSA-N", "OKKDEIYWILRZIA-OSZBKLCCSA-N"),
                        SMILES = c("CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
                                   "CC1C(C(CC(O1)OC2CC(CC3=C(C4=C(C(=C23)O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O.Cl",
                                   "CC1=C(C(=O)N(C2=NC(=NC=C12)Cl)C3CCCC3)Br",
                                   "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
                                   "C1=CN(C(=O)N=C1N)[C@H]2C([C@@H]([C@H](O2)CO)O)(F)F"),
                        "PubChem CID" = c("208908", "31703","5330286","36314","60750"),
                        "PubChem URL" = c("https://pubchem.ncbi.nlm.nih.gov/compound/208908","https://pubchem.ncbi.nlm.nih.gov/compound/31703","https://pubchem.ncbi.nlm.nih.gov/compound/5330286", "https://pubchem.ncbi.nlm.nih.gov/compound/36314","https://pubchem.ncbi.nlm.nih.gov/compound/60750"))

#add treatment metadata
path_to_metadata <- "../../../../images/"
df_data_metadata <- lapply(plateIDs, function(plateID){
  metadata <- read_excel(paste0(path_to_metadata, plateID,"/metadata/",plateID,".xlsx")) %>%
    mutate(plateID = plateID,
           well = str_remove(Well, "0"),
           "Organism" = "Homo sapiens",
           "Term Source 1 REF" = "NCBITaxon",
           "Term Source 1 Accession" = "",
           "Control Type" = case_when(Drug1 == "control" ~"negative control"))
}) %>%
  bind_rows() %>%
  right_join(df, by = c("plateID"="plate_name", "well" = "well")) %>%
  left_join(drug_metadata, by = c("Drug1" ="name")) %>%
  rename("Drug1_InChlkey" = "InChlkey",
         "Drug1_SMILES" = "SMILES",
         "Compound1 PubChem CID" = "PubChem CID",
         "Compound1 PubChem URL" = "PubChem URL") %>%
  left_join(drug_metadata, by = c("Drug2" ="name")) %>%
  rename("Drug2_InChlkey" = "InChlkey",
         "Drug2_SMILES" = "SMILES",
         "Compound2 PubChem CID" = "PubChem CID",
         "Compound2 PubChem URL" = "PubChem URL") %>%
  mutate(channel_label = case_when(channel_name == "P" ~"Phase contrast",
                                   channel_name == "G" ~"Cell cycle reporter",
                                   channel_name == "R" ~"Nuclear reporter"),
         "Tissue Type" = "Breast; Mammary gland",
         Drug1 = case_when(Drug1 == "control" ~"",
                           TRUE ~ Drug1)) %>%
  select(filename, plateID, well, Organism, "Term Source 1 REF",
         "Term Source 1 Accession", Cellline, "Tissue Type", field, channel_name, channel_label, timepoint,  "Control Type",Drug1, Drug1Concentration, Drug1ConcentrationUnits, Drug1_InChlkey, Drug1_SMILES, "Compound1 PubChem CID", "Compound1 PubChem URL", Drug2, Drug2Concentration, Drug2ConcentrationUnits, Drug2_InChlkey, Drug2_SMILES,"Compound2 PubChem CID", "Compound2 PubChem URL") %>%
  rename("Plate" = "plateID",
         "Well" = "well",
         "Cell Line" = "Cellline",
         "Channels" = "channel_name")

result <- write_tsv(df_data_metadata, "../metadata/AU565_library_all_metadata_test.tsv")

IDR_well_level_library <- df_data_metadata %>%
  select(Plate, Well, Organism, "Term Source 1 REF",
         "Term Source 1 Accession", "Cell Line", "Tissue Type","Control Type",Drug1, Drug1Concentration, Drug1ConcentrationUnits, Drug1_InChlkey, Drug1_SMILES, "Compound1 PubChem CID", "Compound1 PubChem URL", Drug2, Drug2Concentration, Drug2ConcentrationUnits, Drug2_InChlkey, Drug2_SMILES,"Compound2 PubChem CID", "Compound2 PubChem URL") %>%
  distinct() %>%
  mutate(Channels = "P:Phase contrast, R:Nuclear reporter, G:Cell cycle reporter")


result <- write_tsv(IDR_well_level_library, "../metadata/AU565_library_metadata.tsv")
