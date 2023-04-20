library(tidyverse)
library(readxl)

#cellLine <- "HCC1143"
cellLine <- "21MT1"
#cellLine <- "MDAMB157"

data_path <- "../../../../images/"
#data_path <- "../images/"

#plateIDs <- c("AU00601","AU00602","AU00701","AU00702","AU00801","AU00802","AU00901","AU00902","AU01001","AU01002","AU01101","AU01102")
#plateIDs = c("HC00701","HC00801","HC00901","HC01001","HC01301","HC01401")
plateIDs = c("2101001","2101201","2101301","2101401","2101501","2101601","2101701","2101801")
#plateIDs = c("MD00301","MD00401","MD00501","MD00601","MD00701","MD00801")

data_paths <- paste0(data_path, plateIDs)

#The library file has a row for every raw image in this cell lines dataset
#it decodes the file name to create a dataframe
library_filepath <- paste0("../metadata/",cellLine,"_file_level_library.csv")
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
  write_csv(df, library_filepath)
} else {
  df <- read_csv(library_filepath)
}

drug_metadata <- tibble(name = c("Lapatinib", "Doxorubicin", "Palbociclib", "Paclitaxel", "Gemcitabine", "BEZ235","Trametinib" ),
                        InChlkey = c("BCFGMOOMADDAQU-UHFFFAOYSA-N", "MWWSFMDVAYGXBV-BXPPNZEESA-N", "HJQCAEDIUJXGCQ-UHFFFAOYSA-N", "RCINICONZNJXQF-HLPYKUKBSA-N", "OKKDEIYWILRZIA-OSZBKLCCSA-N","JOGKUKXHTYWRGZ-UHFFFAOYSA-N","LIRYPHYGHXZJBZ-UHFFFAOYSA-N"),
                        SMILES = c("CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
                                   "CC1C(C(CC(O1)OC2CC(CC3=C(C4=C(C(=C23)O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O.Cl",
                                   "CC1=C(C(=O)N(C2=NC(=NC=C12)Cl)C3CCCC3)Br",
                                   "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
                                   "C1=CN(C(=O)N=C1N)[C@H]2C([C@@H]([C@H](O2)CO)O)(F)F","N#CC(c1ccc(cc1)n1c(=O)n(c2c1c1cc(ccc1nc2)c1cnc2c(c1)cccc2)C)(C)C",
                                   "CC(=O)Nc1cccc(c1)n1c(=O)n(C2CC2)c(=O)c2c1c(C)c(=O)n(c2Nc1ccc(cc1F)I)C "),
                        "PubChem CID" = c("208908", "31703","5330286","36314","60750","11977753","11707110"),
                        "PubChem URL" = c("https://pubchem.ncbi.nlm.nih.gov/compound/208908","https://pubchem.ncbi.nlm.nih.gov/compound/31703","https://pubchem.ncbi.nlm.nih.gov/compound/5330286", "https://pubchem.ncbi.nlm.nih.gov/compound/36314","https://pubchem.ncbi.nlm.nih.gov/compound/60750","https://pubchem.ncbi.nlm.nih.gov/compound/11977753", "https://pubchem.ncbi.nlm.nih.gov/compound/11707110"))

#add treatment metadata
path_to_metadata <- "../../../../images/"
df_data_metadata <- lapply(plateIDs, function(plateID){
  metadata <- read_excel(paste0(path_to_metadata, plateID,"/metadata/",plateID,".xlsx")) %>%
    mutate(plateID = plateID,
           well = str_remove(Well, "0"),
           "Organism" = "Homo sapiens",
           "Term Source 1 REF" = "NCBITaxon",
           "Term Source 1 Accession" = "NCBITaxon_9606",
           "Term Source 2 REF" = "EFO",
           #Need case when: AU565 EFO_0001087, MDAMB157 EFO:0000305
           "Term Source 2 Accession" = case_when(Cellline == "AU565" ~"EFO_0001087",
                                                 Cellline == "HCC1143" ~"EFO_0001169",
                                                 Cellline == "21MT1" ~"EFO_000xxxx",
                                                 Cellline == "MDAMB157" ~"EFO_0000305"),
           "Term Source 3 REF" = "UBERON",
           "Term Source 3 Accession" = "UBERON_0000310",
           "Term Source 4 REF" = "UBERON",
           "Term Source 4 Accession" = "UBERON_0001911",
           "Control Type" = case_when(Drug1 == "control" ~"negative control"))
}) %>%
  bind_rows() %>%
  right_join(df, by = c("plateID"="plate_name", "Well" = "well")) %>%
  left_join(drug_metadata, by = c("Drug1" ="name")) %>%
  rename("Compound 1 InChlkey" = "InChlkey",
         "Compound 1 SMILES" = "SMILES",
         "Compound 1 PubChem CID" = "PubChem CID",
         "Compound 1 PubChem URL" = "PubChem URL") %>%
  left_join(drug_metadata, by = c("Drug2" ="name")) %>%
  rename("Compound 2 InChlkey" = "InChlkey",
         "Compound 2 SMILES" = "SMILES",
         "Compound 2 PubChem CID" = "PubChem CID",
         "Compound 2 PubChem URL" = "PubChem URL") %>%
  mutate(channel_label = case_when(channel_name == "P" ~"Phase contrast",
                                   channel_name == "G" ~"Cell cycle reporter",
                                   channel_name == "R" ~"Nuclear reporter"),
         "Tissue Type" = "Breast",
         "Tissue Type 2" = "Mammary gland",
         "Control Type" = case_when(Drug1 == "Untreated" ~"Negative Control",
                                    TRUE ~""),
         Drug1 = case_when(Drug1 == "control" ~"",
                           TRUE ~ Drug1)
         ) %>%
  select(filename, plateID, well, Organism, "Term Source 1 REF",
         "Term Source 1 Accession", Cellline, "Term Source 2 REF",
         "Term Source 2 Accession","Tissue Type", "Term Source 3 REF",
         "Term Source 3 Accession","Tissue Type 2", "Term Source 4 REF",
         "Term Source 4 Accession",field, channel_name, channel_label, timepoint,  "Control Type",
         Drug1, Drug1Concentration, "Compound 1 InChlkey", "Compound 1 SMILES" , "Compound 1 PubChem CID", "Compound 1 PubChem URL",
         Drug2, Drug2Concentration, Drug2ConcentrationUnits, "Compound 2 InChlkey", "Compound 2 SMILES" ,"Compound 2 PubChem CID", "Compound 2 PubChem URL", "Control Type") %>%
  rename("Plate" = "plateID",
         "Well" = "well",
         "Cell Line" = "Cellline",
         "Channels" = "channel_name") 

result <- write_tsv(df_data_metadata, paste0("../metadata/",cellLine,"_file_level_library_metadata.tsv"))

IDR_well_level_library <- df_data_metadata %>%
  select(Plate, Well, 
         Organism, "Term Source 1 REF", "Term Source 1 Accession",
         "Cell Line",  "Term Source 2 REF",
         "Term Source 2 Accession",
         "Tissue Type", "Term Source 3 REF","Term Source 3 Accession",
         "Tissue Type 2", "Term Source 4 REF","Term Source 4 Accession",

         Drug1, Drug1Concentration, "Compound 1 InChlkey", "Compound 1 SMILES" , "Compound 1 PubChem CID", "Compound 1 PubChem URL",
         Drug2, Drug2Concentration, "Compound 2 InChlkey", "Compound 2 SMILES" , "Compound 2 PubChem CID", "Compound 2 PubChem URL", "Control Type") %>%
  distinct() %>%
  mutate(Drug1 = case_when(Drug1 == "Untreated" ~"",
                           TRUE ~Drug1)) %>%
  rename("Compound 1 Name" = Drug1,
         "Compound 1 Concentration (microMolar)" = Drug1Concentration,
         "Compound 2 Name" = Drug2,
         "Compound 2 Concentration (microMolar)" = Drug2Concentration) %>%
  mutate(Channels = "P:Phase contrast, R:Nuclear reporter, G:Cell cycle reporter")

result <- write_csv(IDR_well_level_library, paste0(data_path,"IDR_",cellLine,"/",cellLine,"_library_metadata.csv"))
