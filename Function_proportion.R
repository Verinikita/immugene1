#' @title fun_prop
#' @description This function estimate the proportion of CDR3 regions.
#' @param data_list is a objet with a list of dataframes with the informations on the samples.
#' @param T.name this parameter can be  "TRBGD"
#' @return a data frame with two columns, a first with Samples id and a second with a proportion index
#' @export dataframe with proportion region TRB and TRGD
#' @examples
#'
#'
#'
fun_prop<- function(data_list, T.name){
  if(T.name == "TRBGD") {
    #
    data_TCR<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    #
    df_TCR <- lapply(data_TCR, function(x) x[which(is.na(x[ ,("cloneIG")])==TRUE), ])
    #
    data_TCRA<- purrr::map(df_TCR, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRA"), sep = "TRA"))
    #Borro las filas de cloneTRA con TRA
    df_TCRA <- lapply(data_TCRA, function(x) x[which(is.na(x[ ,("cloneTRA")])==TRUE), ])
    #
    data2<- purrr::map(df_TCRA, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRG"), sep = "TRG"))
    data3<- purrr::map(data2,   ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRD"), sep = "TRD"))
    data4<- purrr::map(data3,   ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRB"), sep = "TRB"))
    #  
    df_CD  <- lapply(data4, function(x) x[which(is.na(x[ ,("cloneTRD")])==TRUE), ]) 
    df_CD1 <- lapply(data4, function(x) x[which(is.na(x[ ,("cloneTRG")])==TRUE), ]) 
    df_CD2 <- lapply(data4, function(x) x[which(is.na(x[ ,("cloneTRB")])==F), ])
    
    data22<- NULL
    for(i in names(data4)){
      df1 <- nrow(data4[[i]][complete.cases(data4[[i]]["cloneTRG"]), ])
      df2 <- nrow(data4[[i]][complete.cases(data4[[i]]["cloneTRD"]), ])
      df3 <- nrow(data4[[i]][complete.cases(data4[[i]]["cloneTRB"]), ])
      if (lengths(data4[[i]]["Proportion"]) >= 1 & df1 > df2) {
        data22[[i]] <- df_CD[[i]][,c("cloneTRB","cloneTRG")]
        data22[[i]]["Proportions_BGyD"]<-prop.table(df_CD[[i]]["Clones"])
        colnames(data22[[i]]) <-c("cloneTRB","cloneTRGD","Proportions_BGyD")
      }else if (lengths(data4[[i]]["Proportion"]) >= 1  & df2 >=df1) {
        data22[[i]] <- df_CD1[[i]][,c("cloneTRB","cloneTRD")]
        data22[[i]]["Proportions_BGyD"]<-prop.table(df_CD1[[i]]["Clones"])
        colnames(data22[[i]]) <-c("cloneTRB","cloneTRGD","Proportions_BGyD")
      }else if (lengths(data4[[i]]["Proportion"]) >= 1 & df2==0 & df1==0 & df3>0){
        data22[[i]] <- df_CD2[[i]][,"cloneTRB"]
        data22[[i]]["Proportions_BGyD"] <- prop.table(df_CD2[[i]]["Clones"])
        colnames(data22[[i]]) <-c("cloneTRB","cloneTRGD","Proportions_BGyD")
      }else {
        data22[[i]]$cloneTRB <- NA
        data22[[i]]$cloneTRGD <- NA
        data22[[i]]$Proportions_BGyD <- NA
        data22[[i]]<- as.data.frame(data22[[i]])
        # colnames(data22[[i]]) <-c("cloneTRB","cloneTRGD","Proportions_BGyD")
      }
    }
    #
    df_CD33 <- lapply(data22, function(x) x[which(is.na(x[ ,("cloneTRGD")])==FALSE), ]) 
    df_CD55 <- lapply(data22, function(x) x[which(is.na(x[ ,("cloneTRB")])==FALSE), ])
    
    data_fitB <- lapply(1:length(df_CD55), function(i){
      if (lengths(df_CD55[[i]]["Proportions_BGyD"]) >= 1){
        sum(df_CD55[[i]]["Proportions_BGyD"])
      }else{
        return(0)
      }
    })
    
    data_fitGD <- lapply(1:length(df_CD33), function(i){
      if (lengths(df_CD33[[i]]["Proportions_BGyD"]) >= 1){
        sum(df_CD33[[i]]["Proportions_BGyD"])
        
      }else{
        return(0)
      }
    })
    
    
    #######################################
    #Editing dataframe TRB
    data_TCRB <- t(as.data.frame(data_fitB))
    Sample<- as.list(names(data4))
    colnames(data_TCRB)<- data.frame("Proportions_BGyD")
    rownames(data_TCRB)<- c(Sample)
    data_TCRBf<- cbind(data_TCRB, "Sample"= Sample)
    data_TCRBff<- data.frame(data_TCRBf) #Como vector
    data_TCRB_f5<- data.frame((data_TCRBff$Sample))
    data_TCRB_f5<- data.frame(t(data_TCRB_f5))
    colnames(data_TCRB_f5)<- "Sample"
    data_TCRB_f4<- data.frame((data_TCRBff$Proportions_BGyD))
    data_TCRB_f4<- data.frame(t(data_TCRB_f4))
    colnames(data_TCRB_f4)<- "Proportions_TRB_GD"
    data_TCR_B<- cbind(data_TCRB_f4, "Sample"= data_TCRB_f5$Sample)
    
    #Editing dataframe TRGD
    data_TCRGD <- t(as.data.frame(data_fitGD))
    Sample<- as.list(names(data4))
    colnames(data_TCRGD)<- data.frame("Proportion_M_TRGyD")
    rownames(data_TCRGD)<- c(Sample)
    data_TCRGDf<- cbind(data_TCRGD, "Sample"= Sample)
    data_TCRGDff<- data.frame(data_TCRGDf) #Como vector
    data_TCRGD_f5<- data.frame((data_TCRGDff$Sample))
    data_TCRGD_f5<- data.frame(t(data_TCRGD_f5))
    colnames(data_TCRGD_f5)<- "Sample"
    data_TCRGD_f4<- data.frame((data_TCRGDff$Proportion_M_TRGyD))
    data_TCRGD_f4<- data.frame(t(data_TCRGD_f4))
    colnames(data_TCRGD_f4)<- "Proportion_TRGD_B"
    data_TCR_GD<- cbind(data_TCRGD_f4, "Sample"= data_TCRGD_f5$Sample)
    dataTCRBGD<-merge(data_TCR_B, data_TCR_GD, by = "Sample")
    return(dataTCRBGD)
  }else if(T.name == "TRX"){
    #Filter IG data
    data_TCR<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    #Delete files of IG
    df_TCR <- lapply(data_TCR, function(x) x[which(is.na(x[ ,("cloneIG")])==TRUE), ])
    
    data1<- purrr::map(df_TCR, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRA"), sep = "TRA"))
    data2<- purrr::map(data1,  ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRG"), sep = "TRG"))
    data3<- purrr::map(data2,  ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRD"), sep = "TRD"))
    data4<- purrr::map(data3,  ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRB"), sep = "TRB"))

    data22<- NULL
    for(i in names(data4)){
      if (lengths(data4[[i]]["Clones"]) >= 1) {
        data22[[i]] <- data4[[i]][,c("Clones", "cloneTRA","cloneTRB","cloneTRG","cloneTRD")]
        data22[[i]]["Proportions_TRX"]<-prop.table(data4[[i]]["Clones"])
      }else {
        data22[[i]]$Proportions_TRX <- NA
        data22[[i]]$cloneTRA        <- NA
        data22[[i]]$cloneTRB        <- NA
        data22[[i]]$cloneTRG        <- NA
        data22[[i]]$cloneTRD        <- NA
        data22[[i]]$Clones          <- NA
        data22[[i]]<- as.data.frame(data22[[i]])
        # colnames(data22[[i]]) <-c("cloneTRB","cloneTRGD","Proportions_BGyD")
      }
    }
   
    df_C    <- lapply(data22, function(x) x[which(is.na(x[ ,("cloneTRA")])==F), ])
    df_CD   <- lapply(data22, function(x) x[which(is.na(x[ ,("cloneTRB")])==F), ])
    df_CD1  <- lapply(data22, function(x) x[which(is.na(x[ ,("cloneTRG")])==F), ]) 
    df_CD2  <- lapply(data22, function(x) x[which(is.na(x[ ,("cloneTRD")])==F), ]) 
    
    
    data_fitA <- lapply(1:length(df_C), function(i){
      if (lengths(df_C[[i]]["Proportions_TRX"]) >= 1){
        sum(df_C[[i]]["Proportions_TRX"])
      }else{
        return(0)
      }
    })
    
    data_fitB <- lapply(1:length(df_CD), function(i){
      if (lengths(df_CD[[i]]["Proportions_TRX"]) >= 1){
        sum(df_CD[[i]]["Proportions_TRX"])
        
      }else{
        return(0)
      }
    })
    
    data_fitG <- lapply(1:length(df_CD1), function(i){
      if (lengths(df_CD1[[i]]["Proportions_TRX"]) >= 1){
        sum(df_CD1[[i]]["Proportions_TRX"])
        
      }else{
        return(0)
      }
    })
    
    data_fitD <- lapply(1:length(df_CD2), function(i){
      if (lengths(df_CD2[[i]]["Proportions_TRX"]) >= 1){
        sum(df_CD2[[i]]["Proportions_TRX"])
        
      }else{
        return(0)
      }
    })
    
    #######################################
    #Editing dataframe TRA
    data_TCRA <- t(as.data.frame(data_fitA))
    Sample<- as.list(names(data4))
    colnames(data_TCRA)<- data.frame("Proportions_TRX")
    rownames(data_TCRA)<- c(Sample)
    data_TCRAf<- cbind(data_TCRA, "Sample"= Sample)
    data_TCRAff<- data.frame(data_TCRAf) #vector
    data_TCRA_f5<- data.frame((data_TCRAff$Sample))
    data_TCRA_f5<- data.frame(t(data_TCRA_f5))
    colnames(data_TCRA_f5)<- "Sample"
    data_TCRA_f4<- data.frame((data_TCRAff$Proportions_TRX))
    data_TCRA_f4<- data.frame(t(data_TCRA_f4))
    colnames(data_TCRA_f4)<- "Proportions_TRA"
    data_TCR_A<- cbind(data_TCRA_f4, "Sample"= data_TCRA_f5$Sample)
    #######################################
    #Editing dataframe TRB
    data_TCRB <- t(as.data.frame(data_fitB))
    Sample<- as.list(names(data4))
    colnames(data_TCRB)<- data.frame("Proportions_TRX")
    rownames(data_TCRB)<- c(Sample)
    data_TCRBf<- cbind(data_TCRB, "Sample"= Sample)
    data_TCRBff<- data.frame(data_TCRBf) #vector
    data_TCRB_f5<- data.frame((data_TCRBff$Sample))
    data_TCRB_f5<- data.frame(t(data_TCRB_f5))
    colnames(data_TCRB_f5)<- "Sample"
    data_TCRB_f4<- data.frame((data_TCRBff$Proportions_TRX))
    data_TCRB_f4<- data.frame(t(data_TCRB_f4))
    colnames(data_TCRB_f4)<- "Proportions_TRB"
    data_TCR_B<- cbind(data_TCRB_f4, "Sample"= data_TCRB_f5$Sample)
    #######################################
    #Editing dataframe TRG
    data_TCRG <- t(as.data.frame(data_fitG))
    Sample<- as.list(names(data4))
    colnames(data_TCRG)<- data.frame("Proportions_TRX")
    rownames(data_TCRG)<- c(Sample)
    data_TCRGf<- cbind(data_TCRG, "Sample"= Sample)
    data_TCRGff<- data.frame(data_TCRGf) #vector
    data_TCRG_f5<- data.frame((data_TCRGff$Sample))
    data_TCRG_f5<- data.frame(t(data_TCRG_f5))
    colnames(data_TCRG_f5)<- "Sample"
    data_TCRG_f4<- data.frame((data_TCRGff$Proportions_TRX))
    data_TCRG_f4<- data.frame(t(data_TCRG_f4))
    colnames(data_TCRG_f4)<- "Proportions_TRG"
    data_TCR_G<- cbind(data_TCRG_f4, "Sample"= data_TCRG_f5$Sample)
    #Editing dataframe TRD
    data_TCRD <- t(as.data.frame(data_fitD))
    Sample<- as.list(names(data4))
    colnames(data_TCRD)<- data.frame("Proportions_TRX")
    rownames(data_TCRD)<- c(Sample)
    data_TCRDf<- cbind(data_TCRD, "Sample"= Sample)
    data_TCRDff<- data.frame(data_TCRDf) #vector
    data_TCRD_f5<- data.frame((data_TCRDff$Sample))
    data_TCRD_f5<- data.frame(t(data_TCRD_f5))
    colnames(data_TCRD_f5)<- "Sample"
    data_TCRD_f4<- data.frame((data_TCRDff$Proportions_TRX))
    data_TCRD_f4<- data.frame(t(data_TCRD_f4))
    colnames(data_TCRD_f4)<- "Proportions_TRD"
    data_TCR_D<- cbind(data_TCRD_f4, "Sample"= data_TCRD_f5$Sample)
    
    dataTCRAB<-merge(data_TCR_A,data_TCR_B, by = "Sample")
    dataTCRGD<-merge(data_TCR_G,data_TCR_D, by = "Sample")
    dataTCRABGD<-merge(dataTCRAB, dataTCRGD, by = "Sample")
    return(dataTCRABGD)
    
  }
}
