## Frequency calculation using GnomAD api
#
# Author: Samuel Jorquera
#
# Estimates the absolute frequency (AF) of a set of variants using the GnomAD api
#
#

## Libraries
library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)



## Variants of interest:
vars <- c(
  "4-23795377-C-T",   #PPARGC1A rs6821591
  "2-177265309-T-G",  #NFE2L2 rs6721961
  "7-95316772-A-T",   #PON1 rs854560
  "7-87509329-A-G",   #ABCB1 rs1045642
  "14-20456995-T-G",  #APEX1 rs1130409
  "3-9757089-C-G",    #OGG1 rs1052133
  "3-165773492-C-T",  #BCHE rs1803274
  "12-117256215-C-T", #NOS1 rs12829185
  "12-117236534-C-T", #NOS1 rs10774910
  "12-117247465-G-A", #NOS1 rs1047735
  "12-117215033-G-A", #NOS1 rs2682826
  "12-117267969-C-T", #NOS1 rs3741480
  "12-117246975-G-T", #NOS1 rs816353
  "11-27658369-C-T",  #BDNF rs6265
  "6-32441753-G-A",   #HLA-DRA rs3129882
  "5-134169184-T-C",  #SKP1 rs2284312
  "5-1447745-C-T",    #SLC6A3 (Upstream) rs2652510
  "5-1447726-G-A"    #SLC6A3 (Upstream) rs2550956 
)


gnomAD <- function(id) {
  
  #Query for GnomAD GraphiQL api
  query <- paste0('{
    variant(variantId: "', id, '", dataset: gnomad_r4) {
      genome {
        populations {
          id
          ac
          an
        }
      }
    }
  }')
  
  # GnomAD api
  res <- POST(
    url = "https://gnomad.broadinstitute.org/api",
    body = list(query = query),
    encode = "json"
  )
  
  # Read
  parsed <- fromJSON(content(res, "text", encoding = "UTF-8"))
  
  # AF calculation
  pops <- parsed$data$variant$genome$populations
  
  df <- data.frame(
    variant = id,
    population = pops$id,
    AC = pops$ac,
    AN = pops$an
  )
  
  df$AF <- df$AC / df$AN
  
  return(df)
}

## Apply gnomAD() over variants
results <- lapply(vars, gnomAD)
df <- bind_rows(results)


## filter for gnomAD ancestries.
AFs <- df %>% filter(population %in% c(
  "amr", #American Admixed
  "fin", #European (finish)
  "ami", #Amish
  "eas", #East Asian
  "mid", #Middle Eastern
  "afr", #African/African American
  "sas", #South Asian
  "asj", #Ashkenazi Jewish
  "nfe"  #European (non-finish)
))

## Rename ancestries
AFs$population <- factor(AFs$population, 
                         levels = c("amr", "fin","ami","eas","mid","afr", "sas", "asj", "nfe"),
                         labels = c("AMR", "FIN", "AMI", "EAS", "MID", "AFR","SAS","ASJ","EUR"))

## Rename variants
AFs <- AFs %>%
  mutate(variant = recode(variant,
                          "4-23795377-C-T"="PPARGC1A rs6821591",
                          "7-95316772-A-T"="PON1 rs854560",
                          "14-20456995-T-G"="APEX1 rs1130409",
                          "3-9757089-C-G"="OGG1 rs1052133",
                          "3-165773492-C-T"="BCHE rs1803274",
                          "12-117256215-C-T"="NOS1 rs12829185",
                          "12-117236534-C-T"="NOS1 rs10774910",
                          "12-117247465-G-A"="NOS1 rs1047735",
                          "12-117215033-G-A"="NOS1 rs2682826",
                          "12-117267969-C-T"="NOS1 rs3741480",
                          "12-117246975-G-T"="NOS1 rs816353",
                          "11-27658369-C-T"="BDNF rs6265",
                          "2-177265309-T-G"="NFE2L2 rs6721961",
                          "7-87509329-A-G"="ABCB1 rs1045642",
                          "6-32441753-G-A"="HLA-DRA rs3129882",
                          "5-134169184-T-C"="SKP1 rs2284312",
                          "5-1447745-C-T"="SLC6A3 (Upstream) rs2652510",
                          "5-1447726-G-A"="SLC6A3 (Upstream) rs2550956"))

## Data rearrangement
Pops <- data.frame(variant = levels(as.factor(AFs$variant)))
for (ancestry in levels(AFs$population)) {
  aux <- AFs %>% filter(population == ancestry)
  aux <- data.frame(variant = aux$variant,
                    AF = aux$AF)
  colnames(aux)[2] <- as.character(ancestry)
  Pops <- merge(Pops, aux, by.x ="variant", by.y = "variant")
}

write.csv(Pops, file = "GPEI_freq_ancestries.csv")
