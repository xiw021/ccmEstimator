
#' Check Arguments for CCM estimands
#'
#' @param Y
#' @param T1
#' @param T2
#' @param M
#' @param data
#'
#'
#'
#' @return A dataframe
#'
#'

checkArguments <- function(Y,T1,T2,M,data = NULL){
  if (!is.null(data)){
    if (!(is.character(Y)&is.character(T1)&is.character(T2)&is.character(M))){
      stop("Y,T1,T2,M must all be name of the outcome variable (as character/string) in data.frame if a dataframe provided.")
    }
    para.df <- subset(data, select = c(Y,T1,T2,M))
    names(para.df) <- c('Y','T1','T2','M')
  }
  else{
    para.df <- data.frame(Y,T1,T2,M)
  }
  # check if Y,T1,T2,M are numeric
  if (!(is.numeric(para.df$Y)&is.numeric(para.df$T1)&is.numeric(para.df$T2)&is.numeric(para.df$M))){
    stop("Y,T1,T2,M must all be numeric.")
  }
  # drop NA values with warning messages
  oldlength = length(para.df$Y)
  para.df = na.omit(para.df)
  if (oldlength != length(para.df$Y)){
    warning(paste0(oldlength - length(para.df$Y), " observations removed due to missing values"))
  }
  # remove T1=T2=1.
  oldlength = length(para.df$Y)
  para.df = subset(para.df,!((T1==1)&(T2==1)))
  if (oldlength != length(para.df$Y)){
    warning(paste0(oldlength - length(para.df$Y), " observations removed due to T1 = T2 = 1"))
  }
  # check if T1 and T2 are binary
  if (sum(!para.df$T1%in%c(0,1),!para.df$T2%in%c(0,1)) != 0){
    stop('T1,T2 must be binary vectors')
  }
  # T1 T2 value test.
  oneZero = subset(para.df,(T1 == 1)&(T2==0))
  ZeroOne = subset(para.df,(T1 == 0)&(T2==1))
  bothZero = subset(para.df,(T1 == 0)&(T2==0))
  if (length(oneZero$T1)*length(ZeroOne$T1)*length(bothZero$T1) == 0){
    stop('T1,T2 must include combinations of 1,0; 0,1 and 0,0.')
  }
  return(para.df)

}
