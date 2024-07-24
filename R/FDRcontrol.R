#' @title Sequential Hypothesis Testing with FDR controlled
#' @description
#' This function is to calculate cutoff value, and the corresponding rejection
#'     result, with p-value sequence, tuning parameter selection method, and
#'     controlled FDR input.
#' @param p_value.seq p-value sequence
#' @param tuning_type indicator to be chosen from either either 'DOS' or 'Storey'
#' @param tuning tuning parameter for DOS or Storey method
#' @param FDR the expected FDR
#' @return list of 2 elements, including the cutoff, and the rejection indicator
#' @export
#' @examples
#' data(data_list_example)
#' TStat = AtST_stat(data.list.example)$TestStat
#' p_value.seq = 1-pchisq(q=TStat,df=1)
#' SeqHyTest_FDR(p_value.seq,tuning_type='DOS',tuning=1,FDR=0.2)
SeqHyTest_FDR = function(p_value.seq,tuning_type,tuning,FDR){
  allowed_tuning_type = c('DOS','Storey')

  if (!tuning_type %in% allowed_tuning_type) {
    stop("Invalid 'tuning_type'. Please choose either 'DOS' or 'Storey'.")
  }

  d = length(p_value.seq)
  if (tuning_type=='DOS') {
    tuning_seq = c(0.5,1)
    sorted_p_value.seq = sort(p_value.seq,decreasing=F)
    d_alpha.vec = c()
    alp = tuning
    for (k in 1:floor(d/2)) {
      d_alpha.vec = c(d_alpha.vec,
                      (sorted_p_value.seq[2*k]-2*sorted_p_value.seq[k])/(k^alp))
    }
    hat_k_alpha = which.max(d_alpha.vec)
    hat_pi_1 = (hat_k_alpha/d
                - sorted_p_value.seq[hat_k_alpha])/(1-sorted_p_value.seq[hat_k_alpha])
    hat_pi_0 = 1-hat_pi_1
  }else if (tuning_type=='Storey') {
    theta = tuning
    hat_pi_0 = sum(p_value.seq>theta)/(d*(1-theta))
  }

  hat_FDR = function(delta){
    if (sum(p_value.seq<=delta)>=1) {
      res = (hat_pi_0*delta)/(sum(p_value.seq<=delta)/d)
    }else{
      res = (hat_pi_0*delta)/(1/d)
    }
    return(res)
  }
  delta_seq = sort(seq(0,0.7,0.0001)^2,decreasing=T)

  for (delta in delta_seq) {
    hat_fdr = hat_FDR(delta=delta)
    if (hat_fdr<FDR) {
      break
    }
  }

  accept_indicate = (p_value.seq>delta)
  reject_indicate = 1 - accept_indicate
  return(list(cutoff=delta,reject_indicator=reject_indicate))
}
