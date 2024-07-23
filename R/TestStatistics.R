#' @title Statistics Calculation
#' @description
#' This function is to calculate the test statistics with data input.
#' @param data.list a list with 3 elements, including X, A, Y
#' @return test statistics
#' @export
#' @importFrom ncvreg ncvreg
#' @examples
#' data(data_list_example)
#' TestStatistics(data.list.example)
TestStatistics = function(data.list){

  # data retrieving
  X = data.list$X
  A = data.list$A
  M = data.list$M
  Y = data.list$Y
  d = ncol(M)
  q = ncol(X)
  n = nrow(X)

  # tuning parameter setting
  an = 0.15*n/log(n)
  bn = 2.25*log(log(n))/(log((log(log(d)+1)+1))+1)

  # penalty factors setting
  penalty.factor_beta = c(rep(0,2),rep(1,q))
  penalty.factor_gamma = c(rep(0,1),rep(1,q))

  # do it for every M[,k] in a loop
  hat_beta_M_seq = c();hat_gamma_A_seq = c()
  hat_sigma2_beta_margin = c();hat_sigma2_gamma_margin = c()
  hat_beta_gamma_per = NULL
  Mk_per_func = function(Mk){
    # penalized regression for Y~M[,k]
    lambda_seq = sort(exp(seq(from=log(sqrt(1/n)/2),to=log(4),length.out=100)),
                      decreasing=T)
    # ncvreg: bic
    bic_beta = ncvreg::ncvreg(X=cbind(A,Mk,X),y=Y,family='gaussian',penalty='SCAD',
                      lambda=lambda_seq,
                      penalty.factor=penalty.factor_beta)
    bic_seq_beta = n*log(bic_beta$loss/n)+
      log(n)*apply(unname(bic_beta$beta),2,function(hat_beta_lambda){
        sum(hat_beta_lambda!=0)
      })
    chosen_lambda_beta = lambda_seq[which.min(bic_seq_beta)]
    hat_beta = bic_beta$beta[,which.min(bic_seq_beta)]

    # estimators assignment of beta
    hat_beta_A = hat_beta[2]
    hat_beta_Mk = hat_beta[3]
    hat_beta_X = hat_beta[4:(3+q)]
    # estimated active sets of beta_X
    hat_act_ind_beta_X = which(hat_beta_X!=0)
    hat_act_num_beta_X = length(hat_act_ind_beta_X)
    # estimated Y~M error terms and corresponding variance
    pred_Y = predict(bic_beta,cbind(A,Mk,X),lambda=chosen_lambda_beta)

    hat_eps_Yk = Y-pred_Y
    # estimated variance matrix of hat_beta_M[k]
    Bnk = cbind(rep(1,n),A,X[,hat_act_ind_beta_X])
    hat_Sigma_Bk_Mk = t(Bnk)%*%matrix(Mk,ncol=1)/n
    hat_Sigma_Bk_Bk = t(Bnk)%*%Bnk/n
    hat_Sigma_Bk_Bk_Mk = hat_Sigma_Bk_Bk-hat_Sigma_Bk_Mk%*%t(hat_Sigma_Bk_Mk)/mean(Mk^2)
    hat_sigma2_beta = mean(((1+as.numeric(t(hat_Sigma_Bk_Mk)%*%solve(hat_Sigma_Bk_Bk_Mk)%*%
                                            hat_Sigma_Bk_Mk)/mean(Mk^2))*Mk
                            -as.vector(t(hat_Sigma_Bk_Mk)%*%solve(hat_Sigma_Bk_Bk_Mk)%*%t(Bnk)))^2*
                             (hat_eps_Yk^2))/(mean(Mk^2)^2)

    # penalized regression for M[,k]~A
    # ncvreg: bic
    bic_gamma = ncvreg::ncvreg(X=cbind(A,X),y=Mk,family='gaussian',penalty='SCAD',
                       lambda=lambda_seq,
                       penalty.factor=penalty.factor_gamma)
    bic_seq_gamma = n*log(bic_gamma$loss/n)+
      log(n)*apply(unname(bic_gamma$beta),2,function(hat_gamma_lambda){
        sum(hat_gamma_lambda!=0)
      })
    chosen_lambda_gamma = lambda_seq[which.min(bic_seq_gamma)]
    hat_gamma = bic_gamma$beta[,which.min(bic_seq_gamma)]

    # estimators assignment of gamma
    hat_gamma_A = hat_gamma[2]
    hat_gamma_X = hat_gamma[3:(q+2)]
    hat_act_ind_gamma_X = which(hat_gamma_X!=0)
    hat_act_num_gamma_X = length(hat_act_ind_gamma_X)
    # estimated M[,k]~A error term[k] and corresponding variance[k]
    pred_Mk = predict(bic_gamma,cbind(A,X),lambda=chosen_lambda_gamma)

    hat_eps_Mk = Mk-pred_Mk
    # estimated variance of hat_gamma[k]
    Rnk = cbind(rep(1,n),X[,hat_act_ind_gamma_X])
    hat_Sigma_Rk_A= t(Rnk)%*%matrix(A,ncol=1)/n
    hat_Sigma_Rk_Rk = t(Rnk)%*%Rnk/n
    hat_Sigma_Rk_Rk_A = hat_Sigma_Rk_Rk-hat_Sigma_Rk_A%*%t(hat_Sigma_Rk_A)/mean(A^2)
    # estimated variance of hat_gamma_A[k]
    hat_sigma2_gamma = mean(((1+as.numeric(t(hat_Sigma_Rk_A)%*%solve(hat_Sigma_Rk_Rk_A)%*%
                                             hat_Sigma_Rk_A)/mean(A^2))*A
                             -as.vector(t(hat_Sigma_Rk_A)%*%solve(hat_Sigma_Rk_Rk_A)%*%t(Rnk)))^2*
                              (hat_eps_Mk^2))/(mean(A^2))^2
    return(c(hat_beta_Mk,hat_sigma2_beta,hat_gamma_A,hat_sigma2_gamma))
  }
  for (k in 1:d) {
    Mk = M[,k]
    hat_beta_gamma_per = cbind(hat_beta_gamma_per,Mk_per_func(Mk))
  }

  hat_beta_M_seq = hat_beta_gamma_per[1,]
  hat_sigma2_beta_margin = hat_beta_gamma_per[2,]
  hat_gamma_A_seq = hat_beta_gamma_per[3,]
  hat_sigma2_gamma_margin = hat_beta_gamma_per[4,]

  T_beta = sqrt(n)*hat_beta_M_seq/sqrt(hat_sigma2_beta_margin)
  T_gamma = sqrt(n)*hat_gamma_A_seq/sqrt(hat_sigma2_gamma_margin)

  I1 = ((an*T_gamma^2)>(T_beta^2))
  I2 = ((an*T_beta^2)>(T_gamma^2))
  T_k_2 = (n*(hat_beta_M_seq*I1+hat_gamma_A_seq*I2)^2/
             ((hat_sigma2_beta_margin*I1+hat_sigma2_gamma_margin*I2+(1/n))*(1+1/sqrt(n))) +
             I1*I2*abs(T_beta*T_gamma)/bn)

  T_choice = I1 + 2*I2
  res_vec_1 = c(T_k_2,T_beta,T_gamma,T_choice,hat_beta_M_seq,hat_gamma_A_seq,
                hat_sigma2_beta_margin,hat_sigma2_gamma_margin)
  # return(unname(res_vec_1))
  return(unname(T_k_2))
}
