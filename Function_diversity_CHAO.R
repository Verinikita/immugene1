#' @title divCHAO
#' @description This function analyzes the input data and gets the estimate
#' on the diversity CHAO of CDR3 regions.
#' @param data_list is a objet with a list of dataframes with the informations on the samples.
#' @param column this value can be "ntIG+TCR", "TRX.nt", "TRX.aa", "TRA.nt", "TRA.aa", "TRB.nt",
#' "TRB.aa", "TRG.nt", "TRG.aa", "TRD.nt", "TRD.aa". "ntIG+TCR"  for the complete use of the information,
#' without filterin the isolation residues for example IG data, then you can filter according
#' your preferences .nt only the nucleotide CDR3 sequences and .aa only CDR3 aminocidic regions,
#' of the different chains TRA TRB TRG TRD.
#' @return a data frame with two columns, a first with Samples id and a second with a diversity index
#' @export normaldivCHAO
#' @examples
#'
#'
divCHAO<- function(data_list, column){
  calcular_diversidad_chao1 <- function(x) {
    # Crear una matriz con los datos de abundancia de cada clon en cada muestra
    matriz_abundancia <- as.matrix(plyr::ldply(x))
    so <- length(x[x > 0])
    s1 <- length(x[x == 1])
    s2 <- length(x[x == 2])
    if ((s1 - s2)^2 == (s1 + s2)^2)
      chao1<- (so + s1 * (s1 - 1)/((s2 + 1) * 2))
    else
      chao1<-(so + s1^2/(s2 * 2))}

  if(column == "TRX.aa") {
    m1<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    m2 <- lapply(m1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4<- purrr::map(m3, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    da_CDR3aa_filterim<- purrr::map(m5,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::count(CDR3.aa, wt = Clones))
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){

        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_aa_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_aa_sinIG = 0)
        return(df_diversidad2)
      }

    })
   purrr::map_dfr(list(df6), dplyr::bind_rows)
  }else if(column == "ntIG+TCR") {
    #ingreso la lista de dataframes
    muestra1<- data_list
    df1 = data.frame()
    for (i in 1:length(muestra1)) {
      # Calcular la diversidad con Chao1 para la muestra i
      diversidad_chao1 <- calcular_diversidad_chao1(muestra1[[i]]$Clones)
      # Crear un data frame con los valores de diversidad
      df_diversidad <- data.frame(Sample = names(muestra1[i]), diversity_ntIGTCR = diversidad_chao1)
      df1 = rbind(df1, df_diversidad)
    }
    purrr::map_dfr(list(df1), dplyr::bind_rows)
  }else if(column == "TRX.nt"){
    m1<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    m2 <- lapply(m1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos

    #ingreso la lista de dataframes
    muestra2<- m3
    df11 = data.frame()
    for (i in 1:length(muestra2)) {
      # Calcular la diversidad con Chao1 para la muestra i
      diversidad_chao1 <- calcular_diversidad_chao1(muestra2[[i]]$Clones)
      # Crear un data frame con los valores de diversidad
      df_diversidad <- data.frame(Sample = names(muestra2[i]), diversity_nt = diversidad_chao1)
      df11 = rbind(df11, df_diversidad)
    }
    purrr::map_dfr(list(df11), dplyr::bind_rows)

  }else if(column == "TRA.nt"){
    m1<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    m2 <- lapply(m1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4<- purrr::map(m3, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRA"), sep = "TRA"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[ ,("cloneTRA")])==FALSE), ])

    #ingreso la lista de dataframes
    muestra5<- m5
    div_TRA_nt = data.frame()
    for (i in 1:length(muestra5)) {
      # Calcular la diversidad con Chao1 para la muestra i
      diversidad_chao1 <- calcular_diversidad_chao1(muestra5[[i]]$Clones)
      # Crear un data frame con los valores de diversidad
      df_diversidad5 <- data.frame(Sample = names(muestra5[i]), diversity_ntTRA_sIG = diversidad_chao1)
      div_TRA_nt <- rbind(div_TRA_nt, df_diversidad5)
    }
    return(div_TRA_nt)
  }else if(column == "TRA.aa"){
    m4<- purrr::map(data_list, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6<- purrr::map(m5, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRA"), sep = "TRA"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[ ,("cloneTRA")])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::count(CDR3.aa, wt = Clones))
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRA_aa_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRA_aa_sinIG = 0)
        return(df_diversidad2)
      }
    })

    purrr::map_dfr(list(df6), dplyr::bind_rows)

  }else if(column == "TRB.nt"){
    m1<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    m2 <- lapply(m1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4<- purrr::map(m3, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRB"), sep = "TRB"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[ ,("cloneTRB")])==FALSE), ])
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    da_CDR3aa_filterim<- m5
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRB_nt_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRB_nt_sinIG = 0)
        return(df_diversidad2)
      }
    })
    purrr::map_dfr(list(df6), dplyr::bind_rows)

  }else if(column == "TRB.aa"){
    m4<- purrr::map(data_list, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6<- purrr::map(m5, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRB"), sep = "TRB"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[ ,("cloneTRB")])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::count(CDR3.aa, wt = Clones))
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRB_aa_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRB_aa_sinIG = 0)
        return(df_diversidad2)
      }
    })

    purrr::map_dfr(list(df6), dplyr::bind_rows)
  }else if(column == "TRG.nt"){
    m1<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    m2 <- lapply(m1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4<- purrr::map(m3, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRG"), sep = "TRG"))
    da_CDR3aa_filterim <- lapply(m4, function(x) x[which(is.na(x[ ,("cloneTRG")])==FALSE), ])
    #Cambio el nombre de la columna 2 = count  por la lista de nombres

    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
        })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRG_nt_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRG_nt_sinIG = 0)
        return(df_diversidad2)
      }
    })
    purrr::map_dfr(list(df6), dplyr::bind_rows)
  }else if(column == "TRG.aa"){
    m4<- purrr::map(data_list, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6<- purrr::map(m5, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRG"), sep = "TRG"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[ ,("cloneTRG")])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::count(CDR3.aa, wt = Clones))
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRG_aa_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRG_aa_sinIG = 0)
        return(df_diversidad2)
      }
    })

    purrr::map_dfr(list(df6), dplyr::bind_rows)
  }else if(column == "TRD.nt"){
    m1<- purrr::map(data_list, ~ tidyr::separate(.x, V.name, c("V.name", "cloneIG"), sep = "IG"))
    m2 <- lapply(m1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4<- purrr::map(m3, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRD"), sep = "TRD"))
    da_CDR3aa_filterim <- lapply(m4, function(x) x[which(is.na(x[ ,("cloneTRD")])==FALSE), ])
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRD_nt_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRD_nt_sinIG = 0)
        return(df_diversidad2)
      }
    })
    purrr::map_dfr(list(df6), dplyr::bind_rows)
  }else if(column == "TRD.aa"){
    m4<- purrr::map(data_list, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6<- purrr::map(m5, ~ tidyr::separate(.x, V.name, c("V.name", "cloneTRD"), sep = "TRD"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[ ,("cloneTRD")])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::count(CDR3.aa, wt = Clones))
    #Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_fiter <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #ingreso la lista de dataframes
    data_conteoff<- data_fiter
    df6 <- lapply(1:length(data_conteoff), function(i){
      if (lengths(data_conteoff[[i]][2]) >= 1){
        # Calcular la diversidad con Chao1 para la muestra i
        diversidad_chao1 <- calcular_diversidad_chao1(data_conteoff[[i]][2])
        # Crear un data frame con los valores de diversidad
        df_diversidad1 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRD_aa_sinIG = diversidad_chao1)
        return(df_diversidad1)
      }else {
        df_diversidad2 <- data.frame(Sample = names(data_conteoff[[i]][2]), div_TRD_aa_sinIG = 0)
        return(df_diversidad2)
      }
    })

    purrr::map_dfr(list(df6), dplyr::bind_rows)
  }
}

