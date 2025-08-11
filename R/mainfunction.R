#--------------------------------------------------------------------------------------------------##
# Import packages  
#--------------------------------------------------------------------------------------------------##

##Ensure that the version of iNEXT.3D is 1.0.7.


################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################
#--------------------------------------------------------------------------------------------------##
# Main function iNEXTmeta_beta() for meta analysis  
#--------------------------------------------------------------------------------------------------##

#' iNEXTmeta_beta() is a function that estimates the difference of standardized 3D (taxonomic, phylogenetic and
#' functional) gamma, alpha, beta diversity, and four dissimilarity measures with two treatments (e.g., enhanced vs. control), and perform meta analysis for several
#' studies/sites. 
#' @param data (a) For datatype = "abundance", data can be input as a matrix/data.frame (species by assemblages).
#' Here an assemblage refers to a combination of site and treatment.
#' The first row lists the names of sites and second row lists the names of treatment. 
#' (b) For datatype = "incidence_raw", data can be input as a list with several lists (assemblages) of
#' matrices/data.frames, each matrix represents species-by-sampling units incidence data . 
#' The first list must be a dataframe with names of site-by-treatment for each list. 
#' @param diversity selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, 
#' and 'FD' = Functional diversity. 
#' @param order.q a number specifying the diversity order. Default is \code{order.q = 0}. 
#' @param datatype data type of input data: individual-based abundance data (datatype = "abundance"), 
#' species by sampling-units incidence matrix (datatype = "incidence_raw") with all entries being 0 (non-detection)
#' or 1 (detection).
#' @param level a number specifying a particular value of sample coverage (between 0 and 1) that will be used to compute standardized 3D estimates. 
#' If level = NULL, then this function computes the standardized 3D diversity estimates and four dissimilarity measures
#' for the minimal quartile among all SC(2n) values that are greater than 0.5.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. 
#' Bootstrap replications are generally time consuming. 
#' Default is nboot = 10. If more accurate results are required, set nboot = 100 (or nboot = 200).
#' @param treatment_order a character vector for the names of treatment. The difference of standardized 3D 
#' diversity and the four dissimilarity measures will be computed as diversity of the first treatment minus the diversity of second treatment. 
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param PDtree (required only when diversity = "PD"), a phylogenetic tree in Newick format for all observed
#' species in the pooled data. 
#' @param PDreftime (required only when diversity = "PD"), a number specifying reference 
#' times for PD. Default is NULL (i.e., the age of the root of PDtree).
#' @param PDtype (required only when diversity = "PD"), select PD type: PDtype = "PD" (effective total branch length)
#' or PDtype = "meanPD" (effective number of equally divergent lineages). Default is "meanPD", where meanPD = PD/tree depth.
#' @param FDdistM (required only when diversity = "FD"), a species pairwise distance matrix for all species in
#' the pooled data.
#' @param FDtype (required only when diversity = "FD"), select FD type: FDtype = "tau_value" for FD under 
#' specified a threshold value, or FDtype = "AUC" (area under the curve of tau-profile) for an overall FD which
#' integrates all threshold values between zero and one. Default is "AUC".
#' @param FDtau (required only when diversity = "FD" and FDtype = "tau_value"), a number between 0 and 1 
#' specifying tau value (threshold levels). If NULL (default), then threshold is set to be the mean distance 
#' between any two individuals randomly selected from the pooled data (i.e., quadratic entropy).
#' @param FDcut_number (argument only for diversity = "FD" and FDtype = "AUC"), a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is FDcut_number = 30. A larger value can be set to obtain more accurate AUC value.
#'
#'
#' @import iNEXT.3D
#' @importFrom dplyr group_by group_split summarise filter select mutate ungroup rename rename_with bind_rows across everything
#' @importFrom purrr map map2_dfr map_dfr map_df imap
#' @importFrom stringr str_to_title
#' @importFrom tibble tibble
#' @importFrom stats qnorm
#' @importFrom metafor rma
#' @importFrom iNEXT.beta3D iNEXTbeta3D
#' 
#' @return a list. The list contains eight dataframes (coverage summary, gamma, alpha, beta diversity, and four dissimilarity measures) that show the results for each site.



iNEXTmeta_beta <- function(data, diversity="TD", order.q=0, datatype="incidence_raw", level = NULL, nboot = 10, treatment_order, conf=0.95, 
                           PDtree, PDreftime = NULL, PDtype = "meanPD",
                           FDdistM, FDtype = "AUC", FDtau = NULL, FDcut_number = 30){
  
  if (datatype=="incidence_raw"&diversity=="FD"){
    
    
    nT = sapply(2:length(data), function(i) ncol(data[[i]]))
    
    if (sum(nT < 8) > 0){
      
      
      message("For functional diversity, the number of sampling units in your data is fewer than eight, which may be insufficient for making reliable inferences using incidence-based models.\n 
              In this case, we treat species incidence frequencies as “abundances”, and all outputs are generated using abundance-based models.\n 
              Accordingly, the sample coverage estimation is also based on the abundance data.")
      
      
      tmp = data[[2]] %>% rowSums %>% data.frame %>% rownames_to_column
      for (i in 3:length(data)) {
        tmp = full_join(tmp, 
                        data[[i]] %>% rowSums %>% data.frame   %>% rownames_to_column,
                        by = "rowname") 
      }
      rownames(tmp) = tmp[,1]
      tmp = tmp[,-1]
      colnames(tmp) <- NULL
      
      
      data = rbind(data[[1]], as.matrix(tmp)) %>% as.matrix
      
      datatype = "abundance" 
    }
    
  }
  
  if(is.null(nboot))nboot=10
  
  if(is.null(conf))conf=0.95
  
  if (datatype=="abundance"){
    
    name <- rownames(data)[1]
    name_site <- unique(as.character(data[1,]))
    n_site <- length(name_site)
    order_site <- order(factor(data[1,], levels=name_site))
    data <- data[,order_site]
    row_nn <- row.names(data)[-c(1,2)]
    
    data_T1 = lapply(name_site, function(x) {
      
      data[-c(1,2), data[2,] == treatment_order[1] & data[1,] == x] %>% apply(2, function(w) as.numeric(w)) %>% set_rownames(row_nn)
      
    }) %>% set_names(name_site)
    
    data_T2 = lapply(name_site, function(x) {
      
      data[-c(1,2), data[2,] == treatment_order[2] & data[1,] == x] %>% apply(2, function(w) as.numeric(w)) %>% set_rownames(row_nn)
      
    }) %>% set_names(name_site)
    
    
    site_names <- unique(data[1,])
    index <- data.frame(Site = rep(site_names,each=2), Treatment = rep(unique(data[2,]),length(site_names)))
    
    
    
    
    
    C_vector <-sapply(1:nrow(index), function(i) data[, data[1,] == index[i,1] & data[2,] == index[i,2]] %>% .[-(1:2),] %>%
                        unlist %>% as.numeric %>% iNEXT.3D:::Coverage(., datatype, 2*sum(.)))
    
    
    
    
    
    
    summary_C <- quantile(C_vector, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    
    summary_C <- as.data.frame(t(summary_C))
    colnames(summary_C) <- c("C_Max", "C_25quantile", "C_50quantile", "C_75quantile", "C_100quantile")
    
    
    
    if(is.null(level)){
      
      level=min(summary_C[summary_C>0.5])
      
    }
    
    
    
    if(level<=0.5) {
      
      print("Due to the low coverage level, we switch to the default coverage level.")
      level=min(summary_C[summary_C>0.5])
    }
    
    
    
    
  }else if (datatype=="incidence_raw"){
    
    name <- rownames(data[[1]])[1]
    name_site <- unique(as.matrix(data[[1]])[1,])
    n_site <- length(name_site)
    order_site <- order(factor(data[[1]][1,], levels=name_site))
    data_remove <- data[-1]
    data_remove <- data_remove[order_site]
    data1 <- data[[1]][,order_site]
    
    data_T1 <- lapply(name_site, function(x) {
      
      data_remove[data1[2,] == treatment_order[1] & data1[1,] == x]
      
    }) %>% set_names(name_site)
    
    data_T2 <- lapply(name_site, function(x) {
      
      data_remove[data1[2,] == treatment_order[2] & data1[1,] == x]
      
    }) %>% set_names(name_site)
    
    site_names <- unique(as.vector(data1[1,]))  
    index <- data.frame(Site =  rep(site_names,each=2), Treatment = rep(unique(data1[2,]),length(site_names)))
    
    
    
    
    
    C_vector <- sapply(1:nrow(index), function(i) {
      
      data[-1][data[[1]][1,] == index[i,1] & data[[1]][2,] == index[i,2]] %>% 
        do.call(rbind,.) %>% 
        iNEXT.3D:::Coverage(., "incidence_raw", 2 * ncol(.))
      
    })
    
    
    
    
    summary_C <- quantile(C_vector, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    
    summary_C <- as.data.frame(t(summary_C))
    colnames(summary_C) <- c("C_Max", "C_25quantile", "C_50quantile", "C_75quantile", "C_100quantile")
    
    
    
    if(is.null(level)){
      
      level=min(summary_C[summary_C>0.5])
    }
    
    
    
    if(level<=0.5) {
      
      print("Due to the low coverage level, we switch to the default coverage level.")
      level=min(summary_C[summary_C>0.5])
    }
    
    
  }else{
    NULL
  }
  
  
  
  
  
  
  
  
  
  
  
  CC <- qnorm((1-conf)/2, lower.tail = F)
  
  ### each dataset
  ibeta_div_T1 <- iNEXTbeta3D(data_T1, diversity=diversity, q=order.q, datatype=datatype, level=level, nboot=nboot, PDtree=PDtree, PDreftime=PDreftime, PDtype=PDtype, FDdistM=FDdistM, FDtype=FDtype, FDtau=FDtau, FDcut_number=FDcut_number)
  ibeta_div_T2 <- iNEXTbeta3D(data_T2, diversity=diversity, q=order.q, datatype=datatype, level=level, nboot=nboot, PDtree=PDtree, PDreftime=PDreftime, PDtype=PDtype, FDdistM=FDdistM, FDtype=FDtype, FDtau=FDtau, FDcut_number=FDcut_number)
  
  
  
  
  gab_T1 = lapply(ibeta_div_T1, function(x) rbind(x$gamma %>% rename("qD" = "Gamma")          %>% mutate(Type = "Gamma"),
                                                  x$alpha %>% rename("qD" = "Alpha")          %>% mutate(Type = "Alpha"),
                                                  x$beta  %>% rename("qD" = "Beta")           %>% mutate(Type = "Beta"),
                                                  x$"1-C" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-C"),
                                                  x$"1-U" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-U"),
                                                  x$"1-V" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-V"),
                                                  x$"1-S" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-S"))) %>% do.call(rbind,.)
  
  gab_T2 = lapply(ibeta_div_T2, function(x) rbind(x$gamma %>% rename("qD" = "Gamma")          %>% mutate(Type = "Gamma"),
                                                  x$alpha %>% rename("qD" = "Alpha")          %>% mutate(Type = "Alpha"),
                                                  x$beta  %>% rename("qD" = "Beta")           %>% mutate(Type = "Beta"),
                                                  x$"1-C" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-C"),
                                                  x$"1-U" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-U"),
                                                  x$"1-V" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-V"),
                                                  x$"1-S" %>% rename("qD" = "Dissimilarity")  %>% mutate(Type = "1-S"))) %>% do.call(rbind,.)
  
  finaldata <- lapply(c("Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", "1-S"), function(y) {
    
    div_T1 = gab_T1 %>% filter(Type == y)
    div_T2 = gab_T2 %>% filter(Type == y)
    
    D1 <- div_T1 %>% .$qD
    D2 <- div_T2 %>% .$qD
    
    Diff <- D1 - D2
    Var <- div_T1$s.e.^2 + div_T2$s.e.^2
    UCL <- Diff+CC*sqrt(Var)
    LCL <- Diff-CC*sqrt(Var)
    
    ### fixed model
    den_fixed <- (sum(1/Var))
    delta_fixed <- sum(Diff/Var)/den_fixed
    
    
    Totaldata <- data.frame(Difference=c(Diff,delta_fixed), 
                            Coverage=level,
                            SE=c(sqrt(Var),sqrt(1/den_fixed)),
                            LCL=c(LCL,delta_fixed-CC*sqrt(1/den_fixed)),
                            UCL=c(UCL,delta_fixed+CC*sqrt(1/den_fixed)),
                            Site=c(name_site,"meta_analysis"),
                            Order.q=rep(order.q,(n_site+1)),
                            Diversity=rep(diversity,(n_site+1)),
                            T1=as.numeric(c(D1,"")), T2=as.numeric(c(D2,"")),
                            weight_fixed=as.numeric(c(((1/Var)/sum(1/Var))*100,"100")))
    colnames(Totaldata)[c(9,10)] <- treatment_order
    colnames(Totaldata)[6] <- name
    
    return(Totaldata)
  })
  
  
  
  names(finaldata) <- paste("Summary", c("Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", "1-S"), sep = "_")
  
  finaldata <- c(list(Summary_Coverage = summary_C), finaldata)
  return(finaldata)
}


#--------------------------------------------------------------------------------------------------##
# Main function ggiNEXTmeta() for graphical displays  
#--------------------------------------------------------------------------------------------------##

#' ggiNEXTmeta() is a function that provides forest plot for the difference of standardized 3D
#' diversity and four dissimilarity measures with two treatments.
#' 
#' @param data the outcome of the iNEXTmeta_beta function. 
#' @param range the range of the forest plot. 
#' @param num_round a number that the values show on the plot are rounded to the specified number of decimal places.
#' @param type (argument only for the output from iNEXTmeta_beta), selection of diversity type for plotting forest plot: 'Gamma' = Gamma diversity, 'Alpha' = Alpha diversity, 'Beta' = Beta diversity,
#' 1-C' = Sørensen-type non-overlap dissimilarity, '1-U' = Jaccard-type non-overlap dissimilarity, '1-V' = Sørensen-type turnover dissimilarity, and '1-S' = Jaccard-type turnover dissimilarity.
#' 
#' 
#' @importFrom dplyr filter
#' @importFrom forestplot forestplot fp_set_style fp_add_header fp_append_row fp_decorate_graph fp_set_zebra_style fp_add_lines fpTxtGp
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom grid grid.text gpar
#' @importFrom gridExtra grid.arrange
#' 
#' 
#' @return A forest plot that visualizing the output of iNEXTbeta3Dmeta. In the plot, it shows the difference of diversity between two treatments for each study/site and meta analysis (fixed-effects model).

ggiNEXTmeta <- function(data, range, num_round, type = NULL)
{
  if (!is.null(type)) data = data[paste("Summary", type, sep = "_")]
  name <- colnames(data[[1]])[6]
  name_treat <- colnames(data[[1]])[c(9,10)]
  ll <- length(data[[1]][,6])
  name_site <- data[[1]][,6][-ll]
  diversity <- data[[1]]$Diversity[1]
  order <- data[[1]]$Order.q[1]
  Coverage=unique(round(data[[1]]$Coverage[-ll],num_round))
  
  forestplot_q <- tibble::tibble(mean  = round(data[[1]]$Difference[-ll],num_round),
                                 lower = round(data[[1]]$LCL[-ll],num_round),
                                 upper = round(data[[1]]$UCL[-ll],num_round),
                                 study = name_site,
                                 q_T1 = round(data[[1]][,9][-ll],num_round),
                                 q_T2 = round(data[[1]][,10][-ll],num_round), 
                                 diff = round(data[[1]]$Difference[-ll],num_round), 
                                 LCL = round(data[[1]]$LCL[-ll],num_round),
                                 UCL= round(data[[1]]$UCL[-ll],num_round),
                                 w_fixed= paste(round(data[[1]]$weight_fixed[-ll],2),"%",sep=""))
  
  forestplot_q |>
    forestplot(labeltext = c(study, q_T1, q_T2, diff, LCL, UCL, w_fixed),title=paste(type,"with Coverage =",Coverage),
               clip = range,
               xlog = F) |>
    fp_set_style(box = "royalblue",
                 line = "darkblue",
                 summary = "royalblue",
                 txt_gp = fpTxtGp(ticks = gpar(cex=1))) |> 
    fp_add_header(study = c("", name),
                  q_T1 = c(paste(diversity,"(q=",order,")",sep=""), name_treat[1]),
                  q_T2 = c(paste(diversity,"(q=",order,")",sep=""), name_treat[2]),
                  diff=c("","Difference"),
                  LCL=c("","LCL"),
                  UCL=c("","UCL"),
                  w_fixed=c("","W(fixed)")) |>
    fp_append_row(mean  = round(data[[1]]$Difference[ll],num_round),
                  lower = round(data[[1]]$LCL[ll],num_round),
                  upper = round(data[[1]]$UCL[ll],num_round),
                  study = "Meta analysis",
                  diff = round(data[[1]]$Difference[ll],num_round),
                  LCL = round(data[[1]]$LCL[ll],num_round),
                  UCL = round(data[[1]]$UCL[ll],num_round),
                  is.summary = TRUE) |> 
    fp_decorate_graph(graph.pos = 4) |> 
    fp_set_zebra_style("#EFEFEF")
}


#' @examples
#' \donttest{
#' # Meta-analysis of incidence data for taxonomic diversity at orders q = 0, 1, and 2
#' data(Bat_incidence_data)
#' meta_inci_TD_q0 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "TD", order.q = 0, 
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50, 
#'                                  treatment_order = c("E","C"), conf = 0.95)
#'                                  
#'                                  
#' meta_inci_TD_q1 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "TD", order.q = 1, 
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50, 
#'                                  treatment_order = c("E","C"), conf = 0.95)
#'                                  
#'                                  
#' meta_inci_TD_q2 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "TD", order.q = 2, 
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50, 
#'                                  treatment_order = c("E","C"), conf = 0.95)
#' 
#' 
#' # Meta-analysis of incidence data for phylogenetic diversity at orders q = 0, 1, and 2
#' data(Bat_incidence_data)
#' data(Bat_tree)
#' 
#' meta_inci_PD_q0 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "PD", order.q = 0,
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50,
#'                                  treatment_order = c("E","C"), conf = 0.95, PDtree = Bat_tree)
#'
#'
#' meta_inci_PD_q1 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "PD", order.q = 1,
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50,
#'                                  treatment_order = c("E","C"), conf = 0.95, PDtree = Bat_tree)
#'   
#'                                 
#' meta_inci_PD_q1 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "PD", order.q = 2,
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50,
#'                                  treatment_order = c("E","C"), conf = 0.95, PDtree = Bat_tree)                                 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Meta-analysis of incidence data for functional diversity at orders q = 0, 1, and 2
#' data(Bat_incidence_data)
#' data(Bat_distM)
#' 
#' 
#' meta_inci_FD_q0 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "FD", order.q = 0, 
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50, 
#'                                  treatment_order = c("E","C"), conf = 0.95, FDdistM = Bat_distM, FDcut_number = 30)
#'
#'
#' meta_inci_FD_q1 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "FD", order.q = 1, 
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50, 
#'                                  treatment_order = c("E","C"), conf = 0.95, FDdistM = Bat_distM, FDcut_number = 30)
#'
#'
#' meta_inci_FD_q2 = iNEXTmeta_beta(data = Bat_incidence_data, diversity = "FD", order.q = 2, 
#'                                  datatype = "incidence_raw", level = NULL, nboot = 50, 
#'                                  treatment_order = c("E","C"), conf = 0.95, FDdistM = Bat_distM, FDcut_number = 30)
#'
#' }
#' 
#' 
#' 
#' @examples
#' \donttest{
#'
#'
#'
#' # Plot meta forest plots for three types of diversity and four types of dissimilarity
#' # of taxonomic diversity, using incidence data with orders q = 0, 1, and 2
#' 
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "Beta")
#'
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_TD_q0, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_TD_q1, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_TD_q2, num_round = 3, range = c(-20, 15), type = "1-S")
#'
#'
#'
#' # Plot meta forest plots for three types of diversity and four types of dissimilarity
#' # of phylogenetic diversity, using incidence data with orders q = 0, 1, and 2
#' 
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_PD_q0, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_PD_q1, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_PD_q2, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' 
#' 
#' # Plot meta forest plots for three types of diversity and four types of dissimilarity
#' # of functional diversity, using incidence data with orders q = 0, 1, and 2
#' 
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_FD_q0, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' 
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_FD_q1, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' 
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "Gamma")
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "Alpha")
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "Beta")
#' 
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "1-C")
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "1-U")
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "1-V")
#' ggiNEXTmeta(meta_inci_FD_q2, num_round = 3, range = c(-20, 15), type = "1-S")
#' 
#' }