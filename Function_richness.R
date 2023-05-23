#' @title fun_rich
#' @description This function analyzes the input data and gets the estimate
#' on the richness of CDR3 regions.
#' @param data_list is a objet with a list of dataframes with the informations on the samples.
#' @param T.name this parameter can be  "TRX.nt", "TRX.aa", "TRA.nt", "TRA.aa", "TRB.nt",
#' "TRB.aa", "TRG.nt", "TRG.aa", "TRD.nt", "TRD.aa"
#' @return a data frame with two columns, a first with Samples id and a second with a richness index
#' @export dataframe with richness region
#' @examples
#'
#'
fun_rich<- function(data_list, T.name){
  if(T.name == "TRX.nt") {
    mmm<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneIG"), sep = "IG"))
    #2 Filtro por los NA de la columna cloneIG
    df_CDR3filter <- lapply(mmm, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    #3 Elimino la columna cloneIG
    df_CDR3TRX <- purrr::map(df_CDR3filter, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    #4 Sumo las filas para CDR3nt
    data_CDR3nt_sum<- purrr::map(df_CDR3TRX,  ~ dplyr::group_by(.x, CDR3.nt) %>% dplyr::summarize(cont=dplyr::n()))
    #5 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3nt <-
      lapply(names(data_CDR3nt_sum), function(i){
        x <- data_CDR3nt_sum[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #6 Une los dataframes CDR3
    data_CDR3ntf <- purrr::reduce(data_CDR3nt, dplyr::full_join, by = "CDR3.nt")
    #7 Sumo las cantidades
    richness_TRX_CDR3nt<- colSums(!is.na(data_CDR3ntf[,-1]))
    richness_TRX_CDR3nt<- as.data.frame(richness_TRX_CDR3nt)
    richness_TRX_CDR3nt["Sample"]<- row.names(richness_TRX_CDR3nt)
    richness_TRX_CDR3nt

  }else if(T.name == "TRA.nt") {
    ###Ahora sumo TRAnt
    #1 Esta si me filtra los TRA en otra columna
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRA"), sep = "TRA"))
    #2 Filtro por los NA de la columna cloneTRA
    df_CDR3filter <- lapply(filter1, function(x) x[which(is.na(x[,"cloneTRA"])==FALSE), ])
    #3 Sumo las filas para CDR3nt
    data_CDR3nt_sum<- purrr::map(df_CDR3filter,  ~ dplyr::group_by(.x, CDR3.nt) %>% dplyr::summarize(cont=dplyr::n()))
    #4 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3nt <-
      lapply(names(data_CDR3nt_sum), function(i){
        x <- data_CDR3nt_sum[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #5 Une los dataframes CDR3
    data_CDR3ntf <- purrr::reduce(data_CDR3nt, dplyr::full_join, by = "CDR3.nt")
    #6 Sumo las cantidades
    richness_TRA_CDR3nt<- colSums(!is.na(data_CDR3ntf[,-1]))
    richness_TRA_CDR3nt<- as.data.frame(richness_TRA_CDR3nt)
    richness_TRA_CDR3nt["Sample"]<- row.names(richness_TRA_CDR3nt)
    richness_TRA_CDR3nt
  }else if(T.name == "TRB.nt") {
    ###Ahora sumo TRBnt
    #1 Esta si me filtra los TRB en otra columna
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRB"), sep = "TRB"))
    #2 Filtro por los NA de la columna cloneTRA
    df_CDR3filter <- lapply(filter1, function(x) x[which(is.na(x[,"cloneTRB"])==FALSE), ])
    #3 Sumo las filas para CDR3nt
    data_CDR3nt_sum<- purrr::map(df_CDR3filter,  ~ dplyr::group_by(.x, CDR3.nt) %>% dplyr::summarize(cont=dplyr::n()))
    #4 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3nt <-
      lapply(names(data_CDR3nt_sum), function(i){
        x <- data_CDR3nt_sum[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #5 Une los dataframes CDR3
    data_CDR3ntf <- purrr::reduce(data_CDR3nt, dplyr::full_join, by = "CDR3.nt")
    #6 Sumo las cantidades
    richness_TRB_CDR3nt<- colSums(!is.na(data_CDR3ntf[,-1]))
    richness_TRB_CDR3nt<- as.data.frame(richness_TRB_CDR3nt)
    richness_TRB_CDR3nt["Sample"]<- row.names(richness_TRB_CDR3nt)
    richness_TRB_CDR3nt
  }else if(T.name == "TRG.nt") {
    ###Ahora sumo TRGnt
    #1 Esta si me filtra los TRG en otra columna
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRG"), sep = "TRG"))
    #2 Filtro por los NA de la columna cloneTRA
    df_CDR3filter <- lapply(filter1, function(x) x[which(is.na(x[,"cloneTRG"])==FALSE), ])
    #3 Sumo las filas para CDR3nt
    data_CDR3nt_sum<- purrr::map(df_CDR3filter,  ~ dplyr::group_by(.x, CDR3.nt) %>% dplyr::summarize(cont=dplyr::n()))
    #4 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3nt <-
      lapply(names(data_CDR3nt_sum), function(i){
        x <- data_CDR3nt_sum[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #5 Une los dataframes CDR3
    data_CDR3ntf <- purrr::reduce(data_CDR3nt, dplyr::full_join, by = "CDR3.nt")
    #6 Sumo las cantidades
    richness_TRG_CDR3nt<- colSums(!is.na(data_CDR3ntf[,-1]))
    richness_TRG_CDR3nt<- as.data.frame(richness_TRG_CDR3nt)
    richness_TRG_CDR3nt["Sample"]<- row.names(richness_TRG_CDR3nt)
    richness_TRG_CDR3nt
  }else if(T.name == "TRD.nt") {
    ###Ahora sumo TRDnt
    #1 Esta si me filtra los TRA en otra columna
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRD"), sep = "TRD"))
    #2 Filtro por los NA de la columna cloneTRD
    df_CDR3filter <- lapply(filter1, function(x) x[which(is.na(x[,"cloneTRD"])==FALSE), ])
    #3 Sumo las filas para CDR3nt
    data_CDR3nt_sum<- purrr::map(df_CDR3filter,  ~ dplyr::group_by(.x, CDR3.nt) %>% dplyr::summarize(cont=dplyr::n()))
    #4 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3nt <-
      lapply(names(data_CDR3nt_sum), function(i){
        x <- data_CDR3nt_sum[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #5 Une los dataframes CDR3
    data_CDR3ntf <- purrr::reduce(data_CDR3nt, dplyr::full_join, by = "CDR3.nt")
    #6 Sumo las cantidades
    richness_TRD_CDR3nt<- colSums(!is.na(data_CDR3ntf[,-1]))
    richness_TRD_CDR3nt<- as.data.frame(richness_TRD_CDR3nt)
    richness_TRD_CDR3nt["Sample"]<- row.names(richness_TRD_CDR3nt)
    richness_TRD_CDR3nt
  }else if(T.name == "TRX.aa") {
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneIG"), sep = "IG"))
    m2 <- lapply(filter1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4<- purrr::map(m3, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    da_CDR3aa_filterim<- purrr::map(m5,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::summarize(cont=dplyr::n())) #Tiene en cuenta los unique para cada secuencia CDR3aa si esta 2 veces n=2
    #5 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3aa <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #6 Une los dataframes CDR3
    data_CDR3aaf <- purrr::reduce(data_CDR3aa, dplyr::full_join, by = "CDR3.aa")
    #7 Sumo las cantidades
    richness_TRX_CDR3aa<- colSums(!is.na(data_CDR3aaf[,-1]))
    richness_TRX_CDR3aa<- as.data.frame(richness_TRX_CDR3aa)
    richness_TRX_CDR3aa["Sample"]<- row.names(richness_TRX_CDR3aa)
    richness_TRX_CDR3aa
  }else if(T.name == "TRA.aa") {
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneIG"), sep = "IG"))
    m2 <- lapply(filter1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4 <- purrr::map(m3, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6 <- purrr::map(m5, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRA"), sep = "TRA"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[,"cloneTRA"])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::summarize(cont=dplyr::n())) #Tiene en cuenta los unique para cada secuencia CDR3aa si esta 2 veces n=2
    #5 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3aa <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #6 Une los dataframes CDR3
    data_CDR3aaf <- purrr::reduce(data_CDR3aa, dplyr::full_join, by = "CDR3.aa")

    #7 Sumo las cantidades
    richness_TRA_CDR3aa<- colSums(!is.na(data_CDR3aaf[,-1]))
    richness_TRA_CDR3aa<- as.data.frame(richness_TRA_CDR3aa)
    richness_TRA_CDR3aa["Sample"]<- row.names(richness_TRA_CDR3aa)
    richness_TRA_CDR3aa

  }else if(T.name == "TRB.aa") {
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneIG"), sep = "IG"))
    m2 <- lapply(filter1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4 <- purrr::map(m3, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6 <- purrr::map(m5, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRB"), sep = "TRB"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[,"cloneTRB"])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::summarize(cont=dplyr::n())) #Tiene en cuenta los unique para cada secuencia CDR3aa si esta 2 veces n=2
    #5 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3aa <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #6 Une los dataframes CDR3
    data_CDR3aaf <- purrr::reduce(data_CDR3aa, dplyr::full_join, by = "CDR3.aa")
    #7 Sumo las cantidades
    richness_TRB_CDR3aa<- colSums(!is.na(data_CDR3aaf[,-1]))
    richness_TRB_CDR3aa<- as.data.frame(richness_TRB_CDR3aa)
    richness_TRB_CDR3aa["Sample"]<- row.names(richness_TRB_CDR3aa)
    richness_TRB_CDR3aa

  }else if(T.name == "TRG.aa") {
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneIG"), sep = "IG"))
    m2 <- lapply(filter1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map(m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4 <- purrr::map(m3, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply(m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6 <- purrr::map(m5, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRG"), sep = "TRG"))
    m7 <- lapply(m6, function(x) x[which(is.na(x[,"cloneTRG"])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map(m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::summarize(cont=dplyr::n())) #Tiene en cuenta los unique para cada secuencia CDR3aa si esta 2 veces n=2
    #5 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3aa <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #6 Une los dataframes CDR3
    data_CDR3aaf <- purrr::reduce(data_CDR3aa, dplyr::full_join, by = "CDR3.aa")
    #7 Sumo las cantidades
    richness_TRG_CDR3aa<- colSums(!is.na(data_CDR3aaf[,-1]))
    richness_TRG_CDR3aa<- as.data.frame(richness_TRG_CDR3aa)
    richness_TRG_CDR3aa["Sample"]<- row.names(richness_TRG_CDR3aa)
    richness_TRG_CDR3aa
  }else if(T.name == "TRD.aa") {
    filter1<- purrr::map(data_list, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneIG"), sep = "IG"))
    m2 <- lapply(filter1, function(x) x[which(is.na(x[,"cloneIG"])==TRUE), ])
    m3 <- purrr::map( m2, ~ dplyr::select(.x,-one_of(c("cloneIG"))))  #elimino la columna otros creada arriba para separar clonotypos
    m4 <- purrr::map( m3, ~ tidyr::separate(.x, CDR3.aa, c("CDR3.aa", "cl"), sep = "\\~|\\*"))
    m5 <- lapply( m4, function(x) x[which(is.na(x[,"cl"])==TRUE), ])
    m6 <- purrr::map( m5, ~ tidyr::separate(.x, allVHitsWithScore, c("allVHitsWithScore", "cloneTRD"), sep = "TRD"))
    m7 <- lapply( m6, function(x) x[which(is.na(x[,"cloneTRD"])==FALSE), ])
    da_CDR3aa_filterim<- purrr::map( m7,  ~ dplyr::group_by(.x, CDR3.aa) %>% dplyr::summarize(cont=dplyr::n())) #Tiene en cuenta los unique para cada secuencia CDR3aa si esta 2 veces n=2
    #5 Cambio el nombre de la columna 2 = count  por la lista de nombres
    data_CDR3aa <-
      lapply(names(da_CDR3aa_filterim), function(i){
        x <- da_CDR3aa_filterim[[ i ]]
        # cambiar el nombre de la segunda columna
        names(x)[2] <- i
        # devuelve
        x
      })
    #6 Une los dataframes CDR3
    data_CDR3aaf <- purrr::reduce(data_CDR3aa, dplyr::full_join, by = "CDR3.aa")
    #7 Sumo las cantidades
    richness_TRD_CDR3aa<- colSums(!is.na(data_CDR3aaf[,-1]))
    richness_TRD_CDR3aa<- as.data.frame(richness_TRD_CDR3aa)
    richness_TRD_CDR3aa["Sample"]<- row.names(richness_TRD_CDR3aa)
    richness_TRD_CDR3aa

  }
}
