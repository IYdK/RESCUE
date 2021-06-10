#' Dummy dataset to test the KDIGO_SCr_scanner function
#'
#' All data presented in the KDIGO SCr scanner dummy dataset are entirely fictional and perhaps even nonsensical. They are meant to serve an educational purpose in understanding how the code for AKI diagnosis and staging according to KDIGO AKI guideline criteria, works. Furthermore, while very effort was made to ensure that accurate information was implemented in the code to implement AKI diagnosis and staging according to the KDIGO AKI guideline serum creatinine criteria, we accept no liability or responsibility to any person or organization as a consequence of any reliance upon the information contained in the code or code use.
#'
#' @format A data frame with 297 observations and 4 variables:
#'
#' @format \[,1\] - 'admission': admission identification number (integer)
#' @format \[,2\] - 'admission_time': admission timestamp (POSIX)
#' @format \[,3\] - 'time': timestamp of serum creatinine value (POSIX)
#' @format \[,4\] - 'value': serum creatinine concentration (μmol/L, numeric)

#'
#'
"KDIGO_SCr_scanner_dummy_data"
