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
    #1FiltroIG de los datos completos
    data_TCR<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    #Borro las filas de IG
    df_TCR <- lapply(data_TCR, function(x) x[which(is.na(x[ ,("cloneIG")])==TRUE), ])
    #1FiltroTRA de los datos sin IG
    data_TCRA<- purrr::map(df_TCR, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRA"), sep = "TRA"))
    #Borro las filas de cloneTRA con TRA
    df_TCRA <- lapply(data_TCRA, function(x) x[which(is.na(x[ ,("cloneTRA")])==TRUE), ])
    #1 Esta si me filtra los TCRG y D en otra columna
    data2<- purrr::map(df_TCRA, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRG"), sep = "TRG"))
    data3<- purrr::map(data2,   ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRD"), sep = "TRD"))
    data4<- purrr::map(data3,   ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRB"), sep = "TRB"))
    #2 Filtro por los NA de la columna cloneTRBGyD Dejando todo
    df_CD  <- lapply(data4, function(x) x[which(is.na(x[ ,("cloneTRD")])==TRUE), ]) #Esta queda con TRB y TRG
    df_CD1 <- lapply(data4, function(x) x[which(is.na(x[ ,("cloneTRG")])==TRUE), ]) #Esta queda con TRB y TRD
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
    #Acomodo el dataframe para TRB
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

    #Acomodo el dataframe para TRGD
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
    }
}



