#' Implementation of the KDIGO AKI guideline SCr criteria for AKI diagnosis and staging
#'
#' @description Acute kidney injury (AKI) diagnosis and staging using the Kidney Disease: Improving Global Outcomes (KDIGO) AKI serum creatinine (SCr) guideline criteria \[1\]. Admissions are first screened for AKI; AKI is diagnosed if the SCr sample values show an absolute or relative (compared to the SCr baseline) increase within the required time-windows, as specified in the KDIGO AKI guideline (See the KDIGO AKI guideline at 2.1.1 on page 19). Next, admissions with AKI are staged; for each SCr sample value at and after AKI onset, the absolute increase (compared to the minimum value before AKI onset), relative increase (compared to the SCr baseline) and absolute SCr values are computed and associated stages assigned. See the KDIGO AKI guideline at 2.1.2, Table 2 (page 19) and the "Timeframe for diagnosis and staging" paragraph on pages 30-31.
#' @description All data presented are entirely fictional and perhaps even nonsensical. They are meant to serve an educational purpose in understanding how the code for AKI diagnosis and staging according to KIDIGO AKI guideline criteria, works. Furthermore, while very effort was made to ensure that accurate information was implemented in the code to implement AKI diagnosis and staging according to the KIDIGO AKI guideline serum creatinine criteria, we accept no liability or responsibility to any person or organization as a consequence of any reliance upon the information contained in the code or code use.
#' @description \[1\] Kellum JA, Lameire N, Aspelin P, et al. Kidney disease: Improving global outcomes (KDIGO) acute kidney injury work group. KDIGO clinical practice guideline for acute kidney injury. Kidney International Supplements. 2012;2(1):1-138.
#' @description https://kdigo.org/wp-content/uploads/2016/10/KDIGO-2012-AKI-Guideline-English.pdf
#'
#' @encoding UTF-8
#' @param data Data frame with SCr samples in rows and the following columns: Admission IDs ('admission' - Integer), admission time ('admission_time' - POSIX), SCr sample timestamps ('time' - POSIX) and SCr sample values ('value' - Numeric) in μmol/L
#' @param wrong_stage_3_crit Boolean, if TRUE, an incorrect stage 3 SCr criterion will be used to stage admissions with AKI: An increase in SCr by 4.0 mg/dl (353.6 μmol/L) instead of an increase to or beyond 4.0 mg/dl (353.6 μmol/L)
#' @return Data frame containing the original data and the added AKI diagnosis and stage status per SCr sample:
#' @return \[,1\] - 'admission': admission identification number (integer)
#' @return \[,2\] - 'admission_time': admission timestamp (POSIX)
#' @return \[,3\] - 'time': timestamp of SCr value (POSIX)
#' @return \[,4\] - 'value': SCr concentration (μmol/L, numeric)
#' @return \[,5\] - 'hours_since_admission': number of hours passed since the admission (numeric)
#' @return \[,6\] - 'creatinine_baseline': SCr baseline concentration, defined as the first SCr within the first 24h of admission (μmol/L ,numeric)
#' @return \[,7\] - 'creatinine_increase_large': any increase >= 26.5 μmol/L during the admission (logical)
#' @return \[,8\] - 'max_creatinine_increase_previous_48': max increase in previous 48h. Only calculated for non-baseline values if creatinine_increase_large == TRUE (μmol/L, numeric)
#' @return \[,9\] - 'AKI_absolute_increase': set to 1 if max_creatinine_increase_previous_48 >= 26.5 μmol/L, else set to 0. Set to NA for baseline values (numeric; \[0, 1, NA\]). See the KDIGO AKI guideline at 2.1.1 (page 19)
#' @return \[,10\] - 'AKI_relative_increase': set to 1 if SCr value / baseline >= 1.5 and increase occurred within previous 7 days, else set to 0. Set to NA for baseline values (integer; \[0, 1, NA\]). See the KDIGO AKI guideline at 2.1.1 (page 19)
#' @return \[,11\] - 'AKI': set to 1 if AKI_absolute_increase or AKI_relative_increase is set to 1. Set to NA for baseline values (integer; \[0, 1, NA\])
#' @return \[,12\] - 'stage_abs_incr': absolute increase in SCr for AKI cases. Calculated for SCr samples at or after AKI onset, by subtracting the minimum SCr before AKI onset from the SCr sample value (μmol/L, numeric)
#' @return \[,13\] - 'stage_rel_incr': relative increase in SCr for AKI cases. Calculated for SCr samples at or after AKI onset, by dividing the SCr sample value by the SCr baseline (numeric)
#' @return \[,14\] - 'AKI_stage_abs': AKI stage associated with the absolute increase in SCr and the absolute SCr value according to the KDIGO AKI SCr criteria (numeric; \[1, 2, 3, NA\]). See the KDIGO AKI guideline at 2.1.2, Table 2 (page 19) and the "Timeframe for diagnosis and staging" paragraph on pages 30-31
#' @return \[,15\] - 'AKI_stage_rel': AKI stage associated with the relative increase in SCr according to the KDIGO AKI SCr criteria (numeric; \[1, 2, 3, NA\]). See the KDIGO AKI guideline at 2.1.2, Table 2 (page 19) and the "Timeframe for diagnosis and staging" paragraph on pages 30-31
#' @return \[,16\] - 'AKI_stage': Maximum stage from the absolute increase stage and the relative increase stage (numeric; \[1, 2, 3, NA\])
#' @export
#'
#### KDIGO SCr scanner function ####
KDIGO_SCr_scanner = function(data, wrong_stage_3_crit=FALSE){

  #### Check arguments and make a cluster for parallel processing ####

  # check arguments
  if (missing(data) | !all(colnames(data) == c("admission", "admission_time", "time", "value")) | !all(class(data[["admission"]]) == "integer" & class(data[["time"]]) == c("POSIXct","POSIXt") & class(data[["admission_time"]]) == c("POSIXct","POSIXt") & class(data[["value"]]) == "numeric")){
    print.noquote("Please provide a dataframe containing four columns: Admission IDs ('admission' - Integer), admission time ('admission_time' - POSIX), SCr sample timestamps ('time' - POSIX) and SCr sample values ('value' - Numeric) in μmol/L...") & return(1)}

  # make cluster for parallel processing
  no_cores = detectCores() - 1
  cl = makeCluster(no_cores)
  clusterEvalQ(cl, library(dplyr))

  #### Sample timestamps, negative and NA SCr value filter, SCr outlier filter, sort dataframe ####

  # find the number of hours passed since the admission moment
  data$hours_since_admission = as.numeric(difftime(data$time,data$admission_time, units = "hours"))

  # remember the number of rows for a potential notification of removed rows
  nrows_before_filter = nrow(data)

  # remove NA values
  data = data[!is.na(data$value),]

  # remove negative creatinine measurements
  data = data[!data$value < 0,]

  # remove creatinine measurements > 10000 μmol/L
  data = data[!data$value > 10000,]

  # notify if any rows removed
  if (nrow(data) < nrows_before_filter){
    print.noquote(paste0("Warning: Removed ", nrows_before_filter - nrow(data), " rows due to NA SCr values, negative SCr values or extremely large SCr values (>10000 μmol/L)..."))
  }

  # sort the samples on admission ID and then on hours_since_admission
  data = arrange(data, admission, hours_since_admission)

  #### SCr baselines ####
  print.noquote("Calculating creatinine baselines: First SCr within the first 24h of admission...")

  # for every admission with creatinine measurements: get the baseline from the first SCr within the first 24h of admission
  creatinine_baselines_first = aggregate(.~admission, data = data[data$hours_since_admission<24,c("admission","value")], FUN=first)
  colnames(creatinine_baselines_first)[-1] = "creatinine_baseline"

  # add the baselines to the data dataframe
  data = join(data, creatinine_baselines_first, by="admission")

  # obtain a vector specifying which samples are SCr baseline samples
  creatinine_baselines_first_pointer = data[data$hours_since_admission<24,] %>% group_by(admission) %>% dplyr::filter(row_number()==1) %>% as.data.frame()
  creatinine_baselines_first_pointer$point = TRUE
  pointer_cols = c("admission","value","admission_time","time","hours_since_admission")
  creatinine_baselines_first_pointer = join(data[,pointer_cols], creatinine_baselines_first_pointer[,c(pointer_cols, "point")], by=pointer_cols)[["point"]]
  creatinine_baselines_first_pointer[is.na(creatinine_baselines_first_pointer)] = FALSE

  #### Find AKIs ####
  print.noquote("Finding AKI by absolute SCr increase...")

  # See the KDIGO AKI guideline at 2.1.1 (page 19)

  # find admissions with absolute increase AKI
  # to save time, first find admissions with any increase in the SCr values >= 26.5 μmol/L (without a 48h window) using the larger_increase_than function
  admissions_max_increase = data.frame("admission" = na.omit(unique(data$admission)),stringsAsFactors = FALSE)
  admissions_max_increase$larger_than = unlist(pblapply(admissions_max_increase$admission, FUN=larger_increase_than, data, 26.5, cl=cl))
  admissions_max_increase = join(data[,c("admission","value")], admissions_max_increase, by="admission")
  data$creatinine_increase_large = admissions_max_increase$larger_than

  # now add a column to find the max increase in all 48h windows if any increase >= 26.5 μmol/L using the max_creatinine_increase_48h function
  data$max_creatinine_increase_previous_48 = NA
  data[!is.na(data$creatinine_increase_large) & data$creatinine_increase_large,]$max_creatinine_increase_previous_48 = pbapply(data[!is.na(data$creatinine_increase_large) & data$creatinine_increase_large,], 1, FUN = max_creatinine_increase_48h, data[!is.na(data$creatinine_increase_large) & data$creatinine_increase_large,], cl=cl)

  # add a column to check if the max increase from each 48h window >= 26.5 μmol/L
  data$AKI_absolute_increase = as.integer(data$max_creatinine_increase_previous_48 >= 26.5)

  # set AKI_absolute increase to 0 if creatinine_increase_large == FALSE and max increase is known
  data[!is.na(data$creatinine_increase_large) & data$creatinine_increase_large==FALSE,"AKI_absolute_increase"] = 0

  # set AKI_absolute increase to NA if creatinine_increase_large == NA
  data[is.na(data$creatinine_increase_large),"AKI_absolute_increase"] = NA

  # set AKI_absolute increase to NA at baseline SCr samples
  data[creatinine_baselines_first_pointer,"AKI_absolute_increase"] = NA

  # find AKI cases by 1.5x baseline increase
  print.noquote("Finding AKI by relative SCr increase...")
  data$AKI_relative_increase = as.integer(data$value >= 1.5 * data$creatinine_baseline)

  # set AKI relative increase columns to NA at SCr baselines
  data[creatinine_baselines_first_pointer,c("AKI_relative_increase")] = NA

  # drop AKIs from >= 1.5x baseline where increase did not occur within  7 days
  data[!is.na(data$AKI_relative_increase) & data$AKI_relative_increase == 1,]$AKI_relative_increase = as.integer(apply(data[!is.na(data$AKI_relative_increase) & data$AKI_relative_increase == 1,c("admission", "hours_since_admission","value")], MARGIN = 1, FUN = relative_increase_week_check, data = data))

  # add a column to summarize the absolute and relative increase AKI columns
  data$AKI = as.integer(data$AKI_relative_increase == 1 | data$AKI_absolute_increase == 1)

  # set AKI columns to NA at SCr baselines
  data[creatinine_baselines_first_pointer,c("AKI")] = NA

  #### AKI staging ####
  print.noquote("AKI staging...")

  # See the KDIGO AKI guideline at 2.1.2, Table 2 (page 19) and the "Timeframe for diagnosis and staging" paragraph on pages 30-31

  # collect admissions with an AKI
  admissions_AKI = unique(data[data$AKI == 1 & !is.na(data$AKI),"admission"])

  # assign the AKI stages using the AKI_stager function
  data[data$admission %in% admissions_AKI,c("stage_abs_incr","stage_rel_incr","AKI_stage_abs","AKI_stage_rel","AKI_stage")] = AKI_stager(admissions_AKI,data[,c("admission","AKI","value","hours_since_admission","creatinine_baseline")], wrong_stage_3_crit = wrong_stage_3_crit)

  #### Return result ####
  stopCluster(cl)
  print.noquote("Done")
  return(data)
}

#### function to find if any SCr increase during admission >= increase argument ####
larger_increase_than = function(admission, dataframe, increase){

  # obtain all samples of the admission
  samples = dataframe[dataframe["admission"]==admission,][["value"]]
  if (length(samples)<=1){return(NA)}

  # find any increase as specified in argument
  result = max(apply(combn(samples,2),2,FUN=diff))>=increase

  return(result)

}


#### function to find max SCr increase in 48h ####
max_creatinine_increase_48h = function(row, dataframe){

  # grab the datetime
  datetime = row["time"]

  # obtain all samples within previous 48h
  samples_within_window = between(as.numeric(difftime(datetime, dataframe[dataframe["admission"]==trimws(row[["admission"]]),"time"], units = "hours")),0,48)
  samples_within_window = dataframe[dataframe["admission"]==trimws(row[["admission"]]),][samples_within_window,][["value"]]

  if (length(samples_within_window)<=1){return(NA)}

  # find max increase in window
  result = max(apply(combn(samples_within_window,2),2,FUN=diff))

  return(result)
}

#### function to check if >= 1.5x increase in SCr compared to baseline occurred within previous 7 days ####
relative_increase_week_check = function(row, data){

  # grab all previous samples within past 7 days
  samples = data[data[["admission"]] == row[["admission"]] & data[["hours_since_admission"]] < row[["hours_since_admission"]] & data[["hours_since_admission"]] > (row[["hours_since_admission"]] - 7*24),]

  # if no samples in these 7 days: return 0
  if (nrow(samples) == 0){return(0)}

  # find the minimum value within this 7 day window
  min = min(samples[["value"]])

  # if the value from the row >= 1.5 * the min value, then AKI by relative increase is within 7 days
  check = as.integer(row[["value"]] / min >= 1.5)
  return(check)

}

#### function to find AKI stage based on absolute SCr increase, relative SCr increase and absolute SCr value, optional argument to apply an incorrect stage 3 SCr criterion ####
AKI_stager = function(admissions, data, wrong_stage_3_crit=FALSE){

  # See the KDIGO AKI guideline at 2.1.2, Table 2 (page 19) and the "Timeframe for diagnosis and staging" paragraph on pages 30-31

  # create a dataframe to store the AKI stages
  aki_stages = data.frame()

  # for each admission: calculate the AKI stage
  for (admn in admissions){

    # obtain position of first AKI sample
    first_aki_index = min(which(data[data$admission==admn,"AKI"]>0))

    # get the number of SCr samples
    n_samples = nrow(data[data$admission==admn,])

    # get the min SCr before AKI onset
    min_scr = min(data[data$admission==admn,"value"][1:first_aki_index-1])

    # make a vector to store absolute increases for the AKI stages
    stage_abs_incr = rep(NA,first_aki_index-1)

    # for each SCr sample starting from the AKI onset, calculate the absolute increase for the AKI stages
    stage_abs_incr = c(stage_abs_incr,
                       data[data$admission==admn,"value"][first_aki_index:n_samples] - min_scr)

    # calculate relative SCr increase compared to baseline
    stage_rel_incr = data[data$admission==admn,"value"] / data[data$admission==admn,"creatinine_baseline"]
    stage_rel_incr[1:first_aki_index-1] = NA

    # make a dataframe for the stages
    stages_admn = data.frame("stage_abs_incr" = stage_abs_incr,
                             "stage_rel_incr" = stage_rel_incr)
    stages_admn$stage_abs = NA
    stages_admn$stage_rel = NA

    # if the wrong_stage_3_crit argument is set to TRUE, apply the incorrect stage 3 SCr criterion
    if (wrong_stage_3_crit){

      # assign stage per absolute or relative increase with the wrong stage 3 criterion
      stages_admn[stages_admn$stage_rel_incr >= 1.5 & stages_admn$stage_rel_incr < 2.0 & !is.na(stages_admn$stage_rel_incr),"stage_rel"] = 1
      stages_admn[stages_admn$stage_abs_incr >= 26.5 & stages_admn$stage_abs_incr < 353.6 & !is.na(stages_admn$stage_abs_incr),"stage_abs"] = 1
      stages_admn[stages_admn$stage_rel_incr >= 2.0 & stages_admn$stage_rel_incr < 3.0 & !is.na(stages_admn$stage_rel_incr),"stage_rel"] = 2
      stages_admn[stages_admn$stage_rel_incr >= 3.0 & !is.na(stages_admn$stage_rel_incr),"stage_rel"] = 3
      stages_admn[stages_admn$stage_abs_incr >= 353.6 & !is.na(stages_admn$stage_abs_incr),"stage_abs"] = 3

    } else{

      # correctly assign stage if wrong_stage_3_crit argument is set to FALSE
      stages_admn[stages_admn$stage_rel_incr >= 1.5 & stages_admn$stage_rel_incr < 2.0 & !is.na(stages_admn$stage_rel_incr),"stage_rel"] = 1
      stages_admn[stages_admn$stage_abs_incr >= 26.5 & !is.na(stages_admn$stage_abs_incr),"stage_abs"] = 1
      stages_admn[stages_admn$stage_rel_incr >= 2.0 & stages_admn$stage_rel_incr < 3.0 & !is.na(stages_admn$stage_rel_incr),"stage_rel"] = 2
      stages_admn[stages_admn$stage_rel_incr >= 3.0 & !is.na(stages_admn$stage_rel_incr),"stage_rel"] = 3
      stages_admn[c(rep(NA,first_aki_index-1),data[data$admission==admn,"value"][first_aki_index:n_samples]) >= 353.6 & !is.na(c(rep(NA,first_aki_index-1),data[data$admission==admn,"value"][first_aki_index:n_samples])),"stage_abs"] = 3
    }

    # get the highest stage from the absolute or relative increase stage
    stages_admn$AKI_stage = apply(stages_admn[,c("stage_abs","stage_rel")],1,FUN=maxNA)

    # add the stages to the AKI stage dataframe
    aki_stages = rbind(aki_stages,stages_admn)

  }

  # convert to matrix and return
  aki_stages = as.matrix(aki_stages)
  colnames(aki_stages) = NULL
  return(aki_stages)
}

#### function to find max of vector with na.rm=T and if all elements NA, return NA instead of -inf ####
maxNA = function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
