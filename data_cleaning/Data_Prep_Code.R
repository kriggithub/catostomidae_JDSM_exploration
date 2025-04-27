#### Data Preparation for JSDM Application to O/E Predictive Modeling
#### K. Voss April 11, 2025
#### This code produces three files from raw 2018-2019 NRSA data ( )
#### (1) fish_table = All fish data (will need to subset out the catostomids)
#### (2) environment_table = The environmental variables unaffected by stress (see sources for description)
#### (3) reference_table = Determines reference condition of sites and indicates 
#### whether site should be used in model constuction ("C") or validation ("V")
####
#### NOTES: Using UID and SITE_ID as unique sample keys depending on data set. 
#### Only including first visit to a site.
#### Reminder: We are only modeling the reference sites in the calibration set

####Step 1: Load Packages
library(tidyverse)

####Step 2: Prepare Fish Taxa Data
fish_data <- read.csv(file="NRSA1819_Fish Count.csv", header = T)

fish_table <- fish_data %>% 
  filter(VISIT_NO == 1 & AG_ECO9 %in% c("SAP","NAP")) %>%
  pivot_wider(id_cols = SITE_ID, names_from = FINAL_NAME, values_from = TOTAL, values_fill = 0)

####Step 3: Abiotic Variables Unaffected by Stress (See Meador & Carlisle (2009))  
#### https://doi-org.dml.regis.edu/10.1577/T08-132.1)
landscape_data <- read.csv(file="NRSA1819_Landscape Data.csv", header = T)
habitat_data <- read.csv(file="NRSA1819_Physical Habitat Metrics.csv", header = T)
site_data <- read.csv(file="NRSA1819_Site Information.csv", header = T)

environmental_table <- habitat_data %>%
  inner_join(site_data, by =c("UID")) %>%
  filter(VISIT_NO.x == 1 & AG_ECO9.x %in% c("SAP","NAP")) %>%
  inner_join(landscape_data, by = c("SITE_ID.x" = "SITE_ID")) %>%
  select(SITE_ID = SITE_ID.x, lat = LAT_DD83, long = LON_DD83, ecoregion_3 = US_L3CODE, elevation = ELEV_PT, 
         slope = XSLOPE_use , ws_area = WSAREASQKM, ws_precip = PRECIP8110WS, ws_temp = TMEAN8110WS, 
         ws_runoff = RUNOFFWS, bfi_ws = BFIWS, perm_ws = PERMWS, erod_ws = KFFACTWS) 

#Notes: lat & long in decimal degrees, ecoregion = Omernik Level III ecoregion 
# Elevation in cm; Watershed Area in km2; Slope in % at site
# Precipitation in mm (30 year watershed average);Temperature in deg C (30 year watershed average); 
# Runoff in ??? (watershed average); Base Flow Index in % (watershed average); 
# Permeability in in/h (watershed average); Erodibility in unitless K factors (watershed average)


####Step 4: Variables needed for Reference Determination (See USEPA (2024) pp. 34)
#### https://www.epa.gov/system/files/documents/2024-12/nrsa-2018-19-tsd-final-11252024.pdf

set.seed(seed=447)
chemical_data <- read.csv(file="NRSA1819_Water Chemistry_Chlorophyll a.csv", header = T)
chemical_data$VISIT_NO <- as.character(chemical_data$VISIT_NO) # To join properly below

reference_table <- habitat_data %>% 
  inner_join(chemical_data, by =c("UID")) %>%
  filter(VISIT_NO.x == "1" & AG_ECO9.x %in% c("SAP", "NAP")) %>%
  inner_join(environmental_table, by = c("SITE_ID.x" = "SITE_ID")) %>%
  select(SITE_ID = SITE_ID.x, AG_ECO9 = AG_ECO9.x, ecoregion_3, total_N = NTL_RESULT, total_P = PTL_RESULT, 
         chloride = CHLORIDE_RESULT, sulfate = SULFATE_RESULT, ANC = ANC_RESULT, DOC = DOC_RESULT, 
         turbidity = TURB_RESULT, RDI = W1_HALL, pct_fine = PCT_FN) %>%
  mutate(total_N_ug = total_N*1000, chloride_ueq = chloride*1000/35.45, sulfate_ueq = sulfate*1000/96.06*2) %>%
  mutate(pass_nap = (total_P <= 20) +  (total_N_ug <= 750) + (chloride_ueq <= 250 | ecoregion_3 %in% c(1,59)) + 
           (sulfate_ueq <= 250) +  (ANC >= 50 | DOC >= 5) + (turbidity <=5) + (RDI <= 2) + (pct_fine <= 25),
         pass_sap = (total_P <= 20) +  (total_N_ug <= 750) + (chloride_ueq <= 200) + (sulfate_ueq <= 400) +  
           (ANC >= 50 | DOC >= 5) + (turbidity <=5) + (RDI <= 2) + (pct_fine <= 25)) %>%
  mutate(ref_cond = case_when(
    AG_ECO9 == "NAP" & pass_nap == 8 ~ "REF",
    AG_ECO9 == "SAP" & pass_sap == 8 ~ "REF",
    .default = "NON-REF"), val_cond = sample(x=c("C","V"),size=n(), replace="TRUE",prob=c(0.70,0.30)))

#Notes: Total N in mg/L; Total P in ug/L; chloride in mg/L; sulfate in mg/L; ANC in ueq/L;
# DOC in mg/L; turbidity in NTU; RDI as described in Ch. 7 above/

####Step 5: Check the distribution of reference and validation sites by aggregated ecoregion 
table(reference_table$ref_cond, reference_table$val_cond, reference_table$AG_ECO9)



