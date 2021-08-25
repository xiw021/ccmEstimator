
#' @title Subsidiary Function for Comparative Causal Mediation Analysis
#' @description Subsidiary function to check correctness of data structure and return final data for analysis
#' @usage checkData(Y,T1,T2,M,data = NULL)
#' @param Y numeric outcome variable. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param T1 binary indicator for first treatment. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param T2 binary indicator for second treatment. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param M numeric mediator variable. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param data an optional data frame containing the variables to be used in analysis
#' @return A data frame that contains the final data to be analyzed in the \code{getCCM()} function.
#' @note This function is called internally and thus should not be used directly.
#' @author Kirk Bansak and Xiaohan Wu
#' @references Bansak, K. (2020). Comparative causal mediation and relaxing the assumption of no mediator-outcome confounding: An application to international law and audience costs. Political Analysis, 28(2), 222-243.
#' @importFrom stats na.omit
#' @examples
#' data(ICAapp)
#' final.dat <- checkData(Y = "dapprp", T1 = "trt1", T2 = "trt2", M = "immorp", data = ICAapp)
#' @author Kirk Bansak and Xiaohan Wu
#' @references Bansak, K. (2020). Comparative causal mediation and relaxing the assumption of no mediator-outcome confounding: An application to international law and audience costs. Political Analysis, 28(2), 222-243.
#' @export

checkData <- function(Y,T1,T2,M,data = NULL){

  if (!is.null(data)){
    if (!is.data.frame(data)){
      stop("Data argument can only take a data.frame object.")
    }
    if (!(is.character(Y)&is.character(T1)&is.character(T2)&is.character(M))){
      stop("Y, T1, T2, and M must all be the names of the variables (as character/string) if a data frame is provided.")
    }
    para.df <- subset(data, select = c(Y,T1,T2,M))
    names(para.df) <- c('Y','T1','T2','M')
  }
  else{
    para.df <- data.frame(Y,T1,T2,M)
  }
  # check if Y,T1,T2,M are numeric
  if (!(is.numeric(para.df$Y)&is.numeric(para.df$T1)&is.numeric(para.df$T2)&is.numeric(para.df$M))){
    stop("Y, T1, T2, and M must all be numeric.")
  }
  # drop NA values with warning messages
  oldlength = length(para.df$Y)
  para.df = na.omit(para.df)
  if (oldlength != length(para.df$Y)){
    warning(paste0(oldlength - length(para.df$Y), " observations removed due to NA values."))
  }
  # remove T1=T2=1.
  oldlength = length(para.df$Y)
  para.df = subset(para.df,!((T1==1)&(T2==1)))
  if (oldlength != length(para.df$Y)){
    warning(paste0(oldlength - length(para.df$Y), " observations with T1 = T2 = 1 removed."))
  }
  # check if T1 and T2 are binary
  if (sum(!para.df$T1%in%c(0,1),!para.df$T2%in%c(0,1)) != 0){
    stop("T1 and T2 must be binary variables.")
  }
  # T1 T2 value test.
  oneZero = subset(para.df,(T1 == 1)&(T2==0))
  ZeroOne = subset(para.df,(T1 == 0)&(T2==1))
  bothZero = subset(para.df,(T1 == 0)&(T2==0))
  if (length(oneZero$T1)*length(ZeroOne$T1)*length(bothZero$T1) == 0){
    stop("Data must include units where (T1,T2) equals (1,0), (0,1), and (0,0).")
  }
  return(para.df)
}
