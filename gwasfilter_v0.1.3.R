# gwasfilter
# Functions: Filtering genome-wide association studies & Obtaining related SNP-trait associations database from GWAS Catalog
# Date of release: 2020-07-10
# Authors:  Songchun Yang (Peking University) & Chongyang Li (Dalian Maritime University) 
# Version: 0.1.3
#==========================================================================================================
message("==================================================================================================")
message("# Welcome to use gwasfilter !")
message("# Version: 0.1.3 (for test)")
message("# Date: 2020-07-10")
message("# Example (total cholesterol):")
message("source(\"gwasfilter.R\")                       #source R script")
message("get_gwasdata(path=YOUR_PATH, checknew=TRUE)  #step1: Downloading the GWAS Catalog database")
message("get_efo(trait=\"Cholesterol\")                 #step2: Querying the EFO ID of a trait")
message("obtain_trait(efoindex=2, append=T)           #step3: Obtaining all reported traits of an aimed trait")
message("store_trait(traitindex=c(1,2,7,8,12))        #step4: Storing the right reported traits of an aimed trait")
message("gwasfilter(efofile=\"TMP_EFO_LIST.csv\", traitfile=\"TMP_EFO_TRAITS.csv\", exclude0=F, replication=T, sample=0, ancestry.eas=T, eas.only=F, journal=c(1), association=T)")
message("                                             #step5: Starting filtering GWASs on the basis of self-defined criterions")
message("==================================================================================================")

require_packages <- function(pkg) {
  ERR <- try( library(pkg, character.only=TRUE), silent=TRUE )
  if (class(ERR)=="try-error") {
    message(paste("installing package",pkg))
    install.packages(pkg)
    library(pkg, character.only=TRUE)
  }
}

require_packages(pkg="readr")
require_packages(pkg="dplyr")
require_packages(pkg="stringr")

check_allele <- function(allele_vector){
  a <- gsub("A","1",allele_vector)
  a <- gsub("T","2",a)
  a <- gsub("C","3",a)
  a <- gsub("G","4",a)
  a <- ifelse(is.na(a),"9",a)
  a <- ifelse(a=="","9",a)
  a <- as.numeric(a)
  return(is.na(a))
}

## Function 1: Downloading or refreshing the GWAS Catalog database
get_gwasdata <- function(path, checknew=TRUE) {
  FURL <- "https://www.ebi.ac.uk/gwas/api/search/downloads/"
  DFS <- c("studies_alternative", "alternative", "ancestry")
  ck.path <- try (is.null(path), silent=TRUE)
  if (class(ck.path)=="try-error") {
    path <- getwd()
    message(paste("Downloading GWAS Catalog Databases to", path, "as default."))
  } else {
    dir.create(path)
  }
  ck <- file.access(paste0(path, "/", DFS[1]), mode = 0)
  if (ck==-1) {
    download.file(url=paste0(FURL, DFS[1]), destfile=paste0(path, "/", DFS[1]))
    download.file(url=paste0(FURL, DFS[2]), destfile=paste0(path, "/", DFS[2]))
    download.file(url=paste0(FURL, DFS[3]), destfile=paste0(path, "/", DFS[3]))
  } else {
    if (checknew==TRUE) {
      download.file(url=paste0(FURL, DFS[1]), destfile=paste0(path, "/", DFS[1], ".tmp"))
      OLD <- read_tsv(file = paste0(path, "/", DFS[1]))
      NEW <- read_tsv(file = paste0(path, "/", DFS[1], ".tmp"))
      if ( nrow(NEW) > nrow(OLD) ) {
        file.rename(paste0(path, "/", DFS[1], ".tmp"), paste0(path, "/", DFS[1]))
        message("==================================")
        message("Updating GWAS Catalog Databases...")
        message("==================================")
        download.file(url=paste0(FURL, DFS[2]), destfile=paste0(path, "/", DFS[2]))
        download.file(url=paste0(FURL, DFS[3]), destfile=paste0(path, "/", DFS[3]))
      } else {
        message("================================================")
        message("GWAS Catalog Databases do not need to be updated")
        message("================================================")
        file.remove(paste0(path, "/", DFS[1], ".tmp"))
      }
    }
  }
  catalog_df <- list()
  catalog_df[[1]] <- read_tsv(file = paste0(path, "/", DFS[1]), progress=FALSE)
  catalog_study <- catalog_df[[1]]
  save(catalog_study, file="catalog_study.RData")
  catalog_df[[2]] <- read_tsv(file = paste0(path, "/", DFS[2]), progress=FALSE)
  catalog_df[[3]] <- read_tsv(file = paste0(path, "/", DFS[3]), progress=FALSE)
  
  df <- subset(catalog_df[[1]], select=c("MAPPED_TRAIT","MAPPED_TRAIT_URI","DISEASE/TRAIT"))
  df$EFO <- gsub("http://www.ebi.ac.uk/efo/", "", df$MAPPED_TRAIT_URI)
  df$EFO <- gsub("http://www.orpha.net/ORDO/", "", df$EFO)
  df$EFO <- gsub("http://purl.obolibrary.org/obo/", "", df$EFO)
  EFOs <- as.data.frame(str_split_fixed(df$EFO, ",", 2), stringsAsFactors=F)
  TRAIT_EFO <- subset(df, EFOs$V2=="" & !duplicated(EFO), select=c("MAPPED_TRAIT","EFO"))
  save(TRAIT_EFO, file="TRAIT_EFO.RData")
  EFO_TRAIT_INFO <- subset(df, select=c("EFO","DISEASE/TRAIT","MAPPED_TRAIT"))
  save(EFO_TRAIT_INFO, file="EFO_TRAIT_INFO.RData")
  
  ances <- subset(catalog_df[[3]], STAGE == "initial", select=c(1,8,9) ) 
  eas <- grep("East Asian", ances$`BROAD ANCESTRAL CATEGORY`)
  ances[eas, "is_EAS"] <- 1
  ances[is.na(ances$is_EAS), "is_EAS"] <- 0
  tra <- grep(",", ances$`BROAD ANCESTRAL CATEGORY`)
  ances[tra, "is_TRA"] <- 1
  ances[is.na(ances$is_TRA), "is_TRA"] <- 0
  ances$ances_count <- ances$is_TRA + 1
  
  sample_size <- aggregate(ances$`NUMBER OF INDIVDUALS`, by=list(ances$`STUDY ACCESSION`), FUN=sum)
  colnames(sample_size) <- c("STUDY ACCESSION", "NUMBER OF INDIVDUALS")
  is_EAS <- aggregate(ances$is_EAS, by=list(ances$`STUDY ACCESSION`), FUN=sum)
  colnames(is_EAS) <- c("STUDY ACCESSION", "is_EAS")
  is_TRA <- aggregate(ances$ances_count, by=list(ances$`STUDY ACCESSION`), FUN=sum)
  colnames(is_TRA) <- c("STUDY ACCESSION", "ances_count")
  catalog_ances <- merge(sample_size, is_EAS, by="STUDY ACCESSION")
  catalog_ances <- merge(catalog_ances, is_TRA, by="STUDY ACCESSION")
  catalog_ances$is_TRA <- ifelse(catalog_ances$ances_count>1, 1, 0)
  save(catalog_ances, file="catalog_ancestry.RData")
  
  asso <- as.data.frame(catalog_df[[2]][, c(37, 11, 12, 13, 15, 21, 27, 22, 25, 31, 32, 28, 30)])
  EA <- as.data.frame(str_split_fixed(asso$`STRONGEST SNP-RISK ALLELE`, "-", 2), stringsAsFactors=F)
  asso$EA <- as.character(EA[,2])
  asso[grep("increase", asso$`95% CI (TEXT)`), "beta_symbol"] <- 1
  asso[grep("decrease", asso$`95% CI (TEXT)`), "beta_symbol"] <- -1
  asso[is.na(asso$beta_symbol), "beta_symbol"] <- 9
  asso$beta <- ifelse(abs(asso$beta_symbol)==1, asso$`OR or BETA` * asso$beta_symbol, NA)
  asso$OR <- ifelse(asso$beta_symbol==9, asso$`OR or BETA`, NA)
  CI <- as.data.frame(str_split_fixed(asso$`95% CI (TEXT)`, "]", 2))
  CI2 <- gsub("[[]", "", CI[,1])
  CI3 <- as.data.frame(str_split_fixed(CI2, "-", 2))
  asso$LCI <- as.numeric(as.character(CI3[,1]))
  asso$UCI <- as.numeric(as.character(CI3[,2]))
  asso.org <- subset(asso, select=c(1,2,3,4,5,8,14,7,9,16:19,12,13))
  colnames(asso.org) <- c("STUDY ACCESSION","Region","Chr","Pos(hg38)","Locus","rsid","Effect_allele","Effect_allele_freq","Function","beta","OR","LCI","UCI","P-VALUE","P-VALUE(TEXT)")
  asso.org$Effect_allele <- ifelse(asso.org$Effect_allele=="?", NA, asso.org$Effect_allele)
  asso.org$Effect_allele_freq <- ifelse(asso.org$Effect_allele=="NR", NA, asso.org$Effect_allele_freq)
  asso.org$Effect_allele_freq <- as.numeric(asso.org$Effect_allele_freq)
  asso.org$Chr.missing <- ifelse(is.na(asso.org$Chr), 1, 0)
  asso.org$Pos.missing <- ifelse(is.na(asso.org$`Pos(hg38)`), 1, 0)
  asso.org[-grep(x = asso.org$rsid,"^rs"), "snpname.invalid"] <- 1
  asso.org[is.na(asso.org$snpname.invalid),"snpname.invalid"] <- 0
  asso.org$EA.missing <- ifelse(is.na(asso.org$Effect_allele), 1, 0)
  EA_invalid <- check_allele(asso.org$Effect_allele)
  asso.org$EA.invalid <- ifelse(EA_invalid==TRUE, 1, 0)
  asso.org$EAF.missing <- ifelse(is.na(asso.org$Effect_allele_freq), 1, 0)
  asso.org$betaOR.missing <- ifelse(is.na(asso.org$beta) & is.na(asso.org$OR), 1, 0)
  asso.org$P.missing <- ifelse(is.na(asso.org$`P-VALUE`), 1, 0)
  save(asso.org, file="catalog_associations.RData")
  return(catalog_df)
}

## Function 2: Querying the EFO ID of a trait 
get_efo <- function(trait){
  load("TRAIT_EFO.RData")
  INDEX <- grep(pattern=trait, x=TRAIT_EFO$MAPPED_TRAIT, ignore.case=TRUE)
  TRAIT.EFO <- TRAIT_EFO[INDEX, c("MAPPED_TRAIT","EFO")]
  #print(TRAIT.EFO)
  save(TRAIT.EFO, file="TMP_TRAIT_EFO.RData")
  write.csv(TRAIT.EFO, file="tmp_get_efo.csv")
  message("============================================")
  message("For details please view file tmp_get_efo.csv")
  message("============================================")
  return(TRAIT.EFO)
}

## Function 3: Obtaining all reported traits of an aimed trait 
obtain_trait <- function(efolist, efoindex, append = T){
  ck1 <- try(is.null(efolist), silent=TRUE)
  if ( class(ck1)=="try-error" ) {
    load("TMP_TRAIT_EFO.RData")
    efolist <- TRAIT.EFO
    #efolist <- read.csv("tmp_get_efo.csv", as.is=T, header=T)
    efolist$N_reported.traits <- NA
  }
  AIM_EFO <- as.character(efolist[efoindex, "EFO"])
  ck2 <- file.access("TMP_EFO_LIST.csv", mode = 0)
  
  if (ck2 == -1) {
    TMP_EFO_LIST <- efolist[efoindex, ]
  } else {
    tmp <- read.csv("TMP_EFO_LIST.csv", header=TRUE, as.is=TRUE)		
    if ( AIM_EFO %in% tmp$EFO == FALSE ) {
      TMP_EFO_LIST <- rbind(tmp, efolist[efoindex, ])
    } else {
      TMP_EFO_LIST <- tmp
    }
  }
  
  ##add in v0.1.2
  if(!append){
    TMP_EFO_LIST = TMP_EFO_LIST[dim(TMP_EFO_LIST)[1],]
  }
  ##
  
  write.csv(TMP_EFO_LIST, file="TMP_EFO_LIST.csv", row.names=FALSE)
  load("EFO_TRAIT_INFO.RData")
  INDEX <- grep(pattern=AIM_EFO, x=EFO_TRAIT_INFO$EFO)
  SUB <- EFO_TRAIT_INFO[INDEX, ]
  SUB2 <- SUB[!duplicated(SUB$`DISEASE/TRAIT`), ]
  #print(SUB2)
  save(SUB2, file="TMP_REPORTED_TRAIT.RData")
  write.csv(SUB2, file="tmp_obtain_trait.csv")
  message("=================================================")
  message("For details please view file tmp_obtain_trait.csv")
  message("=================================================")
  save(AIM_EFO, file="TMP_AIM_EFO.RData")
  return(SUB2)
}

## Function 4: Storing the right reported traits of an aimed trait 
store_trait <- function(aim.efo, traitlist, traitindex){
  ck1 <- try(is.null(aim.efo), silent=TRUE)
  if ( class(ck1)=="try-error") {
    load("TMP_AIM_EFO.RData")
  } else {
    AIM_EFO <- aim.efo
  }	
  ck2 <- try(is.null(traitlist), silent=TRUE)
  if ( class(ck2)=="try-error") {
    load("TMP_REPORTED_TRAIT.RData")
    traitlist <- SUB2
  }
  KEEPLIST <- traitlist[traitindex, c("EFO", "DISEASE/TRAIT")]
  colnames(KEEPLIST) <- c("related.EFO", "reported.traits")
  KEEPLIST$EFO <- AIM_EFO
  TMP_EFO_LIST <- read.csv("TMP_EFO_LIST.csv", header=TRUE, as.is=TRUE)
  TMP_EFO_LIST[TMP_EFO_LIST$EFO==AIM_EFO, "N_reported.traits"] <- nrow(KEEPLIST)
  write.csv(TMP_EFO_LIST, file="TMP_EFO_LIST.csv", row.names=FALSE)
  ck <- file.access("TMP_EFO_TRAITS.csv", mode = 0)
  if (ck == -1) {
    TMP_EFO_TRAITS <- KEEPLIST
  } else {
    tmp <- read.csv("TMP_EFO_TRAITS.csv", header=TRUE, as.is=TRUE)
    tmp2 <- rbind(tmp, KEEPLIST)
    TMP_EFO_TRAITS <- tmp2[!duplicated(tmp2),]
  }
  write.csv(TMP_EFO_TRAITS, file="TMP_EFO_TRAITS.csv", row.names=FALSE)
  #print(TMP_EFO_LIST)
  return(TMP_EFO_LIST)
}

## Function 5: Starting filtering GWASs on the basis of self-defined criterions
gwasfilter <- function(efofile="TMP_EFO_LIST.csv", traitfile="TMP_EFO_TRAITS.csv", exclude0=F, replication=F, sample=0, ancestry.eas=F, eas.only=F, journal=c(1,2,3,4,9), association=F) {
  load("catalog_study.RData")
  load("catalog_ancestry.RData")
  load("catalog_associations.RData")
  journal_rank <- read.csv(file="Catalog_JCR_2018.csv", header=TRUE, as.is=TRUE)

  qctab <- read.csv(efofile, header=TRUE, as.is=TRUE)
  kept_traits <- read.csv(traitfile, header=TRUE, as.is=TRUE)
  
  ## Function Settings
  if (replication==TRUE) {
    REPLICATION <- 1
  } else {
    REPLICATION <- c(0,1)
  }
  
  SAMPLE <- sample
  
  if (ancestry.eas==TRUE) {
    ANCESTRY <- 1
  } else {
    ANCESTRY <- c(0,1)
  }
  
  if (eas.only==TRUE) {
    EAS.ONLY <- 0
  } else {
    EAS.ONLY <- c(0,1)
  }
  
  RANK <- journal

  ASSOCIATION <- association

  qctab$status <- NA 
  qctab$EFO_record <- NA
  qctab$merged_record <- NA
  qctab$nonzero_merged_record <- NA
  qctab$is_replicated <- NA
  qctab$Sample.size.median <- NA
  qctab$is_EAS <- NA
  qctab$JCR.Q1 <- NA
  qctab$total_kept <- NA
  final_kept_study <- list()
  final_snp <- list()
  
  for (i in 1:nrow(qctab) ) {
    TraitID <- qctab[i, "MAPPED_TRAIT"]
    EFO <- qctab[i, "EFO"]
    message(paste("Dealing with Trait",i ,":", TraitID))
    if (is.na(qctab[i,"status"])) {
      keep1 = c() 
      for (j in 1:nrow(catalog_study) ) {
        uris = catalog_study$MAPPED_TRAIT_URI[j]
        if( grepl(pattern = EFO, x = uris) ){
          keep1 = append(keep1, j)
        }
      }
      keep_study1 <- catalog_study[keep1, ]
      qctab[i, "EFO_record"] <- nrow(keep_study1)
      if ( nrow(keep_study1)==0 ) {
        message("Warning message:")
        message(paste("Trait", i, ":", TraitID, "EFO has no merged record"))
        qctab[i, "status"] <- 0
      } else {
        right_traits_name <- kept_traits[kept_traits$EFO==EFO, "reported.traits"]
        keep_study2 <- keep_study1[keep_study1$`DISEASE/TRAIT` %in% right_traits_name, ]
        qctab[i,"merged_record"] <- nrow(keep_study2)
        if ( nrow(keep_study2)==0 ) {
          message("Warning message:")
          message(paste("Trait", i, ":", TraitID, "reported traits has no merged record"))
          qctab[i, "status"] <- 0
        } else {
          for (x in 1:length(right_traits_name)) {
            NAME = right_traits_name[x]
            ck1 <- nrow(keep_study2[keep_study2$`DISEASE/TRAIT`==NAME,])
            kept_traits[kept_traits$EFO==EFO & kept_traits$reported.traits==NAME, "N_merged_record"] <- ck1
            if (ck1==0) {
              message("Warning message:")
              message(paste("Trait", i, ":", TraitID, "reported trait", NAME, "has no merged record"))
            }
          }
          
          keep_study2_nonzero <- subset(keep_study2, `ASSOCIATION COUNT`>0)
          qctab[i,"nonzero_merged_record"] <- nrow(keep_study2_nonzero)

          if (exclude0==TRUE) {
            keep_study3 <- as.data.frame(keep_study2_nonzero)
          } else if (exclude0==FALSE) {
            keep_study3 <- as.data.frame(keep_study2)
          }
          
          if ( nrow(keep_study3)>0 ) {
            keep_study3$replicated <- ifelse(!is.na(keep_study3$`REPLICATION SAMPLE SIZE`), 1, 0)
            keep_study3[, "Sample.size"] <- NA
            keep_study3[, "is_EAS"] <- NA
            keep_study3[, "is_TRA"] <- NA
            keep_study3[, "JCRrank"] <- NA
            for (x2 in 1:nrow(keep_study3)) {
              GCST <- as.character(keep_study3[x2, "STUDY ACCESSION"])
              ck.ances <- catalog_ances[catalog_ances$`STUDY ACCESSION`==GCST, ]
              if (nrow(ck.ances)>0) {
                keep_study3[x2, "Sample.size"] <- catalog_ances[catalog_ances$`STUDY ACCESSION`==GCST, "NUMBER OF INDIVDUALS"]
                keep_study3[x2, "is_EAS"] <- catalog_ances[catalog_ances$`STUDY ACCESSION`==GCST, "is_EAS"]
                keep_study3[x2, "is_TRA"] <- catalog_ances[catalog_ances$`STUDY ACCESSION`==GCST, "is_TRA"]
              }
              JOUR = as.character(keep_study3[x2, "JOURNAL"])
              ck.jour <- journal_rank[journal_rank$JOURNAL==JOUR, ]
              if (nrow(ck.jour)>0) {
                keep_study3[x2, "JCRrank"] <- journal_rank[journal_rank$JOURNAL==JOUR, "JCRrank"]
              }
            }
            qctab[i, "is_replicated"] <- sum(keep_study3$replicated)
            qctab[i, "Sample.size.median"] <- median(keep_study3$Sample.size, na.rm=T)
            qctab[i, "is_EAS"] <- sum(keep_study3$is_EAS, na.rm=T)
            JCR.Q1 <- keep_study3[keep_study3$JCRrank==1, ]
            qctab[i, "JCR.Q1"] <- nrow(JCR.Q1)
            
            keep_study4 <- subset(keep_study3, replicated%in%REPLICATION & (Sample.size>=SAMPLE | is.na(Sample.size))
                                  & (is_EAS%in%ANCESTRY | is.na(is_EAS) ) & (is_TRA%in%EAS.ONLY | is.na(is_TRA)) & JCRrank%in%RANK )
            qctab[i, "total_kept"] <- nrow(keep_study4)
            
            if (nrow(keep_study4)>0) {
              keep_study4$EFO <- EFO
              if (ASSOCIATION==TRUE) {
                snp_database <- c()
                for (st in 1:nrow(keep_study4)) {
                  GCST <- as.character(keep_study4[st, "STUDY ACCESSION"])
                  snp <- subset(asso.org, `STUDY ACCESSION`==GCST)
                  if (nrow(snp)>0) {
                    snp$EFO <- EFO
                    keep_study4[st, "Chr.missing"] <- sum(snp$Chr.missing, na.rm=T)
                    keep_study4[st, "Pos.missing"] <- sum(snp$Pos.missing, na.rm=T)
                    keep_study4[st, "snpname.invalid"] <- sum(snp$snpname.invalid, na.rm=T)
                    keep_study4[st, "EA.missing"] <- sum(snp$EA.missing, na.rm=T)
                    keep_study4[st, "EA.invalid"] <- sum(snp$EA.invalid, na.rm=T)
                    keep_study4[st, "betaOR.missing"] <- sum(snp$betaOR.missing, na.rm=T)
                    keep_study4[st, "P.missing"] <- sum(snp$P.missing, na.rm=T)
                    snp_database <- rbind(snp_database, snp)
                  }
                }
                final_snp[[i]] <- snp_database
              }
            } else {
              if (ASSOCIATION==TRUE) {
                final_snp[[i]] <- asso.org[0, ]
              }
            }
            
            final_kept_study[[i]] <- keep_study4
            qctab[i, "status"] <- 1
            ns <- nrow(keep_study4)
            if (ns>1) {
              str = "studies have"
            } else {
              str = "study has"
            }
            message(paste0("(",i,"/",nrow(qctab),") Trait <",TraitID,"> finished! ", ns, " ", str, " been filtered."))
          } else {
            qctab[i, "status"] <- 0
          }
        }
      }     
    } else {
      message(paste0("(",i,"/",nrow(qctab),") Trait <",TraitID,"> has finished before!"))
    }
  }
  output_check <- subset(qctab, total_kept>0)
  if (nrow(output_check)>0) {
    STUDY_DATABASE <- c()
    for (i in 1:nrow(qctab) ) {
      STUDY_DATABASE <- rbind(STUDY_DATABASE, final_kept_study[[i]])
    }
    if (ASSOCIATION==TRUE) {
      STUDY_DATABASE <- STUDY_DATABASE[,c(22,1:21,23:29)]
    } else {
      STUDY_DATABASE <- STUDY_DATABASE[,c(22,1:21)]
    }
    file1 = paste0(Sys.Date(), "_study_database.csv")
    write.csv(STUDY_DATABASE, file=file1, row.names=F)
    message(paste0("Study database has been exported to ", getwd(), "/", file1) )
    
    if (ASSOCIATION==TRUE) {
      SNP_DATABASE <- c()
      for (i in 1:nrow(qctab) ) {
        SNP_DATABASE <- rbind(SNP_DATABASE, final_snp[[i]])
      }
      if (nrow(SNP_DATABASE)>0) {
        SNP_DATABASE <- SNP_DATABASE[,c(24,1:23)]
        file2 = paste0(Sys.Date(), "_snp_database.csv")
        write.csv(SNP_DATABASE, file=file2, row.names=F)
        message(paste0("SNP database has been exported to ", getwd(), "/", file2) )
      }		
    }
  } else {
    message("[Warning] No study has been kept after filtering. Maybe try another filtering strategy.")
  }
  
  file3 = paste0(Sys.Date(), "_qc_table1.csv")
  write.csv(qctab, file=file3, row.names=F)
  message(paste0("Quality control table 1 has been exported to ", getwd(), "/", file3) )
  file4 = paste0(Sys.Date(), "_qc_table2.csv")
  write.csv(kept_traits, file=file4, row.names=F)
  message(paste0("Quality control table 2 has been exported to ", getwd(), "/", file4) )
  
}