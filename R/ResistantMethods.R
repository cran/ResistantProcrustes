
#'Repeated Median Method without reflection
#'Align two matrices using repeated median method that cannot deal with reflection
#'@param X1 This matrix is fixed during the transformation
#'@param X2 This matrix is transformed to fit the other matrix
#' @return
#' A list containing rotation matrix,scaling,translation of matrix X1, translation of matrix X2, objective function, and sum of square value
#' @importFrom fBasics tr
#' @importFrom stats median
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export
MM2_estimator_no_reflection <- function(X1,X2) { #X2 is transformed;
  trans_X1 = colMeans(X1)
  trans_X2 = colMeans(X2)
  k_value = dim(X2)[1]
  k_1 = rep(1,k_value)
  X1_tran = X1-k_1%*%t(trans_X1)
  X2_tran = X2-k_1%*%t(trans_X2)
  s1 = sqrt(tr(X1_tran%*%t(X1_tran)))
  s2 = sqrt(tr(X2_tran%*%t(X2_tran)))
  X1_p = X1_tran/s1
  X2_p = X2_tran/s2
  temp3 = svd(t(X1_p)%*%X2_p)
  rot = temp3$v%*%t(temp3$u)
  X2_p = X2
  X1_p = X1
  k_2 = rep(1,k_value-1)
  tau_i = rep(0,k_value)
  for(i in 1:k_value){
    X1_i = k_2%*%t(X1_p[i,])
    X1_j = X1_p[-i,]
    X1_ji = X1_j-X1_i
    X2_i = k_2%*%t(X2_p[i,])
    X2_j = X2_p[-i,]
    X2_ji = X2_j-X2_i
    tau_i[i] = median(sqrt(diag(X1_ji%*%t(X1_ji))/diag(X2_ji%*%t(X2_ji))))
  }
  tau_tild = median(tau_i)
  theta_i = rep(0,k_value)
  for(i in 1:k_value){
    X1_i = k_2%*%t(X1_p[i,])
    X1_j = X1_p[-i,]
    X1_ji = X1_j-X1_i
    X2_i = k_2%*%t(X2_p[i,])
    X2_j = X2_p[-i,]
    X2_ji = X2_j-X2_i
    ang_vec = atan2(X2_ji[,2],X2_ji[,1])-atan2(X1_ji[,2],X1_ji[,1])
    ind = ang_vec>pi
    ang_vec[ind] = ang_vec[ind]-2*pi
    ind = ang_vec<(-pi)
    ang_vec[ind] = ang_vec[ind]+2*pi
    theta_i[i] = median(ang_vec)
  }
  theta_tild = median(theta_i)
  H_tild = diag(rep(cos(theta_tild),2))
  H_tild[2,1] = sin(theta_tild)
  H_tild[1,2] = -H_tild[2,1]
  T_now = X1_p - tau_tild*X2_p%*%H_tild
  trans_1 = median(T_now[,1])
  trans_2 = median(T_now[,2])

  ##########################################
  obj_vec = (T_now[,1]-trans_1)^2 + (T_now[,2]-trans_2)^2
  obj_fun = sum(sqrt(obj_vec))
  ##########################################

  return(list(rot=H_tild,scal=tau_tild,trans_1=trans_1,trans_2=trans_2,s=obj_fun,obj_vec=obj_vec))
}


############################################################################
#'Repeated Median Method with reflection
#'Align two matrices using the repeated median method while allowing for reflection. This method is recommended for practical applications where reflected configurations may occur.
#'@param X1 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param X2 A matrix representing the configuration to be transformed to best align with X1.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1 and X2, objective function value, and the sum of squared errors.
#' @importFrom fBasics tr
#' @importFrom stats median
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export

MM2_estimator_with_reflection <- function(X1,X2) { #X2 is transformed;
  k_value = dim(X2)[1]
  MM2_no_reflect_save = MM2_estimator_no_reflection(X1,X2)
  MM2_with_reflect_save = MM2_estimator_no_reflection(X1,X2%*%diag(c(-1,1)))
  obj_trunc_no_reflect = sum(sort(MM2_no_reflect_save$obj_vec)[1:ceiling(k_value/2)])
  obj_trunc_with_reflect = sum(sort(MM2_with_reflect_save$obj_vec)[1:ceiling(k_value/2)])
  if(obj_trunc_no_reflect<obj_trunc_with_reflect){
    rot = MM2_no_reflect_save$rot
    scal = MM2_no_reflect_save$scal
    trans_1 = MM2_no_reflect_save$trans_1
    trans_2 = MM2_no_reflect_save$trans_2
    s = MM2_no_reflect_save$s
  }
  else{
    rot = diag(c(-1,1))%*%(MM2_with_reflect_save$rot)
    scal = MM2_with_reflect_save$scal
    trans_1 = MM2_with_reflect_save$trans_1
    trans_2 = MM2_with_reflect_save$trans_2
    s = MM2_with_reflect_save$s
  }

  return(list(rot=rot,scal=scal,trans_1=trans_1,trans_2=trans_2,s=s))
}



############################################################################
#'Huber and Biweight methods
#'Align two matrices using the Huber and Biweight loss functions via iterative majorization.
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param tuning_con Tuning constant beyond which the resistant loss functions become flat.
#'@param loss Specifies the loss function to use: `"Huber"` or `"Biweight"`.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1, weight matrix, objective function value, and convergence.
#' @importFrom fBasics tr
#' @importFrom stats median mad
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export

majorization <- function(X1,X2,tuning_con,loss) { #loss=1 for Huber, and 2 for biweight; X1 is transformed
  SX = mean(apply(X1,2,mad,constant=1))
  if(loss==1){
    tuning = tuning_con*SX  # SX
  }else{
    tuning = tuning_con*SX  # 4*SX
  }
  weight = diag(rep(1,dim(X1)[1]))
  rot = diag(rep(1,dim(X1)[2]))
  trans_X1 = colMeans(X1)
  trans_X2 = colMeans(X2)
  k_1 = rep(1,dim(X1)[1])
  scal = 1
  ri = rep(0,dim(X1)[1])
  obj_old = 1000000000
  obj_new = 10000000
  dead_count=0
  converge = 0
  while (abs(obj_old-obj_new)>0.000001) {
    dead_count = dead_count +1
    if(dead_count>10000){
      converge = 1
      break
    }else{
      obj_old = obj_new
    }
    for(i in 1:length(ri)){
      temp = X2[i,]-trans_X2 - scal*t(rot)%*%(X1[i,]-trans_X1)
      ri[i] = sqrt(sum(temp^2))
    }
    if(dead_count==1){
      weight = rep(1,dim(X1)[1])
      obj_new = 100000
    }else if(loss==1){
      weight = rep(1,dim(X1)[1])
      ind = (ri>=tuning)
      if(sum(ind)>0){
        weight[ind] = tuning/(ri[ind])
        obj_new = sum(tuning*ri[ind]-0.5*(tuning^2))
      }else{
        obj_new = 0
      }
      if(sum(!ind)>0){
        obj_new = obj_new + sum(0.5*(ri[!ind]^2))
      }
    }else{
      weight = rep(0,dim(X1)[1])
      ind = (ri<=tuning)
      if(sum(ind)>0){
        weight[ind] = (1-((ri[ind])/tuning)^2)^2
        obj_new = (tuning^2)*sum(1-(1-(ri[ind]^2)/(tuning^2))^3)/6
      }else{
        obj_new = 0
      }
      if(sum(!ind)>0){
        obj_new = obj_new + sum(!ind)*(tuning^2)/6
      }
    }
    weight = diag(weight)
    trans_X1 = colSums(weight%*%X1)/tr(weight)
    trans_X2 = colSums(weight%*%X2)/tr(weight)
    X1_star = X1 - k_1%*%t(trans_X1)
    X2_star = X2 - k_1%*%t(trans_X2)
    sqw = sqrt(weight)
    X1_star = sqw%*%X1_star
    X2_star = sqw%*%X2_star
    temp3 = svd(t(X2_star)%*%X1_star)
    rot = temp3$v%*%t(temp3$u)
    #scal = Re(tr(sqrtm(t(X1_star)%*%X2_star%*%t(X2_star)%*%X1_star))/tr(X1_star%*%t(X1_star)))
    scal = tr(t(X2_star)%*%X1_star%*%rot)/tr(t(X1_star)%*%X1_star)
  }
  #print(abs(obj_old-obj_new))
  trans_P = trans_X2 - scal*t(rot)%*%trans_X1    # new code
  return(list(rot=rot,scal=scal,trans_P=trans_P,weight=weight,s=obj_new,converge=converge))  #converge=0, converge
  #  return(list(rot=rot,scal=scal,trans_X1=trans_X1,trans_X2=trans_X2,weight=weight))         #converge=1, not converge
}

########################################################################################
set_tuning_constant = seq(0.5,5,0.1)
Res_initial = 100000000000000
#' Huber and Biweight Methods with Automatic Tuning Constant Selection
#' Align two matrices using Huber or biweight loss functions via iterative
#' majorization, with automatic selection of the tuning constant.
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param loss Specifies the loss function to use: `"Huber"` or `"Biweight"`.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1, weight matrix, objective function value, and convergence.
#' @importFrom fBasics tr
#' @importFrom stats median
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export

majorization_tuned <- function(X1,X2,loss) { #loss=1 for Huber, and 2 for biweight; X1 is transformed
  k_value = dim(X1)[1]
  k_4 = rep(1,k_value)
  Res_save = Res_initial
  for(tuning_con in set_tuning_constant){
    error_r=try( ( maj_result <- majorization(X1,X2,tuning_con,loss) ),
                 silent=TRUE )

    if ('try-error' %in% class(error_r)) {
      next
    }

    if (maj_result$converge==1){
      next
    }

    Res = X2 - maj_result$scal*X1%*%(maj_result$rot)- k_4%*%t(c(maj_result$trans_P[1,1],maj_result$trans_P[2,1]))
    Res = sum( sort(sqrt(rowSums(Res^2)))[1:ceiling(k_value/2)] )
    if(Res<Res_save){
      Res_save = Res
      rot=maj_result$rot
      scal=maj_result$scal
      trans_P=maj_result$trans_P
      weight=maj_result$weight
      s=maj_result$s
    }
  }
  converge=1
  if(Res_save < Res_initial){
    converge=0
  }
  return(list(rot=rot,scal=scal,trans_P=trans_P,weight=weight,s=s,s_permu=Res_save,  #converge=0, converge
              converge=converge))  #converge=1, not converge; s_permu can be used for permunova
}


########################################################################################
xi_func <- function(s,c_value,A_value,E) {
  E = E/s
  ind = (E<=c_value)
  Enow = E[ind]
  result = sum( (Enow^2)/2 - (Enow^4)/2/(c_value^2) + (Enow^6)/6/(c_value^4) )
  result = result + sum(!ind)*(c_value^2)/6
  result = (result/length(E)) - A_value
  return(result)
}

####################################################
#' @importFrom stats median dchisq
integrand_xi <- function(y,c_value,dim_m) {
  result = y/2 - (y^2)/2/(c_value^2) + (y^3)/6/(c_value^4)
  result = result*dchisq(y, dim_m)
  return(result)
}

#' @importFrom stats integrate pchisq

A_value_func <- function(c_value,dim_m) {
  A_value = integrate(integrand_xi, lower=0, upper=c_value^2, c_value=c_value, dim_m=dim_m)$value
  A_value = A_value + pchisq(c_value^2, dim_m, lower.tail=FALSE)*(c_value^2)/6
  return(A_value)
}

BDP_func <- function(c_value,dim_m,BDP_value) {
  A_value = A_value_func(c_value,dim_m)
  result = A_value - BDP_value*(c_value^2)/6
  return(result)
}

dim_m = 2
BDP_value = 0.5
c_value = uniroot(BDP_func,c(1,10),dim_m=dim_m,BDP_value=BDP_value)$root
A_value = A_value_func(c_value,dim_m)
#c_value=2.66080
#A_value=0.58999
####################################################
#'S-estimator method
#'Align two matrices using the S-estimator.
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param c_value It controls where the loss function begins to flatten.
#'@param A_value It is a normalization constant used in the S-estimator defining equation.
#'@param nR Number of randomly sampled subsets.
#'@param s_l Lower bound for the scale s.
#'@param s_u Upper bound for the scale s.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1,and scale s.
#' @importFrom fBasics tr
#' @importFrom stats median uniroot
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export

#s_l=0.000000001
#s_u=10

S_estimator <- function(X1,X2,c_value,A_value,nR,s_l=0.000000001,s_u=10) { #X1 is transformed; nR is the number of subsets;
  k_value = dim(X1)[1]                              #nR<=0 means using all subsets
  p_value = dim(X1)[2]+1                            #s_l and s_u are lower and upper bounds for s
  total_comb = combinations(k_value,p_value)        #if there is problem then increase s_u
  if(nR<=0){
    nR = dim(total_comb)[1]
    combr = total_comb
  }else{
    combr = total_comb[sample.int(dim(total_comb)[1],nR),]
  }
  k_1c = rep(1,dim(X1)[1])
  k_1 = rep(1,p_value)
  s_save = s_u
  for (i in 1:nR) {
    X1_sub = X1[combr[i,],]
    X2_sub = X2[combr[i,],]
    trans_X1 = colMeans(X1_sub)
    trans_X2 = colMeans(X2_sub)
    X1_tran = X1_sub-k_1%*%t(trans_X1)
    X2_tran = X2_sub-k_1%*%t(trans_X2)
    temp3 = svd(t(X2_tran)%*%X1_tran)
    rot = temp3$v%*%t(temp3$u)
    #scal = Re(tr(sqrtm(t(X1_tran)%*%X2_tran%*%t(X2_tran)%*%X1_tran))/tr(X1_tran%*%t(X1_tran)))
    scal = tr(t(X2_tran)%*%X1_tran%*%rot)/tr(X1_tran%*%t(X1_tran))
    ####################################
    X1_tran = X1-k_1c%*%t(trans_X1)
    X2_tran = X2-k_1c%*%t(trans_X2)
    E_i = X2_tran - scal*X1_tran%*%rot
    E_i = sqrt(rowSums(E_i*E_i))
    test = xi_func(s_save,c_value,A_value,E_i)
    if(test<0){
      s_now = uniroot(xi_func,c(s_l,s_u),c_value=c_value,A_value=A_value,E=E_i)$root
      if(s_now<s_save){
        s_save = s_now
        trans_X1_save = trans_X1
        trans_X2_save = trans_X2
        rot_save = rot
        scal_save = scal
      }
    }
    ####################################
  }
  trans_P = trans_X2_save - scal_save*t(rot_save)%*%trans_X1_save    # new code
  return(list(rot=rot_save,scal=scal_save,trans_P=trans_P,s=s_save))  #new code
  #  return(list(rot=rot_save,scal=scal_save,trans_X1=trans_X1_save,trans_X2=trans_X2_save,s=s_save))
}
########################################################################################

xi_func_LMS <- function(s,A_value,E) {
  ind = ((E/s)>=1)
  result = sum(ind)/length(E) - A_value
  return(result)
}

#'Least median of squares (LMS) method
#'Align two matrices using the LMS.
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param nR Number of randomly sampled subsets.
#'@param s_l_LMS Lower bound for the scale s.
#'@param s_u_LMS Upper bound for the scale s.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1,and scale s.
#' @importFrom fBasics tr
#' @importFrom stats median uniroot
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export

#s_l_LMS=0.000000001
#s_u_LMS=10
LMS_estimator <- function(X1,X2,nR,s_l_LMS=0.000000001,s_u_LMS=10) { #X1 is transformed; nR is the number of subsets;
  k_value = dim(X1)[1]                              #nR<=0 means using all subsets
  p_value = dim(X1)[2]+1                            #s_l_LMS and s_u_LMS are lower and upper bounds for s
  total_comb = combinations(k_value,p_value)        #if there is problem then increase s_u_LMS
  if(nR<=0){
    nR = dim(total_comb)[1]
    combr = total_comb
  }else{
    combr = total_comb[sample.int(dim(total_comb)[1],nR),]
  }
  k_1c = rep(1,dim(X1)[1])
  k_1 = rep(1,p_value)
  A_value = floor((k_value+1)/2)/k_value
  s_save = s_u_LMS
  for (i in 1:nR) {
    X1_sub = X1[combr[i,],]
    X2_sub = X2[combr[i,],]
    trans_X1 = colMeans(X1_sub)
    trans_X2 = colMeans(X2_sub)
    X1_tran = X1_sub-k_1%*%t(trans_X1)
    X2_tran = X2_sub-k_1%*%t(trans_X2)
    temp3 = svd(t(X2_tran)%*%X1_tran)
    rot = temp3$v%*%t(temp3$u)
    #scal = Re(tr(sqrtm(t(X1_tran)%*%X2_tran%*%t(X2_tran)%*%X1_tran))/tr(X1_tran%*%t(X1_tran)))
    scal = tr(t(X2_tran)%*%X1_tran%*%rot)/tr(X1_tran%*%t(X1_tran))
    ####################################
    X1_tran = X1-k_1c%*%t(trans_X1)
    X2_tran = X2-k_1c%*%t(trans_X2)
    E_i = X2_tran - scal*X1_tran%*%rot
    E_i = sqrt(rowSums(E_i*E_i))
    test = xi_func_LMS(s_save,A_value,E_i)
    if(test<0){
      s_now = uniroot(xi_func_LMS,c(s_l_LMS,s_u_LMS),A_value=A_value,E=E_i)$root
      if(s_now<s_save){
        s_save = s_now
        trans_X1_save = trans_X1
        trans_X2_save = trans_X2
        rot_save = rot
        scal_save = scal
      }
    }
    ####################################
  }
  trans_P = trans_X2_save - scal_save*t(rot_save)%*%trans_X1_save    # new code
  return(list(rot=rot_save,scal=scal_save,trans_P=trans_P,s=s_save))  #new code
  #  return(list(rot=rot_save,scal=scal_save,trans_X1=trans_X1_save,trans_X2=trans_X2_save,s=s_save))
}
########################################################################################
#'Improved Least median of squares (ILMS) method
#'Align two matrices using the ILMS.
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param nR Number of randomly sampled subsets.
#'@param s_l_LMS Lower bound for the scale s.
#'@param s_u_LMS Upper bound for the scale s.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1,and scale s.
#' @importFrom fBasics tr
#' @importFrom stats median uniroot
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @examples
#' X1 <- matrix(c(
#'   0.3068822, -1.24752606,
#'   -2.1150143, 4.24520347,
#'   0.5271354, 1.64468795,
#'   -0.1264126, -0.08113392,
#'   -1.4749896, -2.26016700,
#'   -0.2347269, -1.27047483,
#'   -0.4833630, 1.43158341,
#'   -2.3499156, 0.85892382,
#'   -2.1021129, 0.40267502,
#'   0.5812934, -0.57210742,
#'   0.7397990, 1.54898604,
#'   -1.2173015, -0.65486036,
#'   -0.6093984, -1.16937177,
#'   0.3833907, -0.72611343,
#'   -0.8423056, 1.22427834
#' ), ncol = 2, byrow = TRUE)
#'
#' X2 <- matrix(c(
#'   1.96990726, -1.67860936,
#'   1.16952104, -3.19996198,
#'   -1.64680983, 3.40009837,
#'   1.51340633, 0.90762220,
#'   6.02678651, -1.29463718,
#'   2.56525182, -0.70441947,
#'   0.08388971, 3.87254957,
#'   4.11361309, 4.77780158,
#'   -0.91329668, 0.01870814,
#'   0.73714253, -0.50391759,
#'   -1.60842305, 3.07780435,
#'   3.90039741, 1.29926458,
#'   3.31127451, -0.28114418,
#'   0.85909201, -0.86365691,
#'   1.52112556, 4.02668249
#' ), ncol = 2, byrow = TRUE)
#' fit <- LMS_LTS_estimator(X1,X2,nR=0,s_l_LMS=0.000000001,s_u_LMS=10)
#' @export

#s_l_LMS=0.000000001
#s_u_LMS=10
LMS_LTS_estimator <- function(X1,X2,nR=0,s_l_LMS=0.000000001,s_u_LMS=10) { #X1 is transformed; nR is the number of subsets;
  k_value = dim(X1)[1]                              #nR<=0 means using all subsets
  p_value = dim(X1)[2]+1                            #s_l_LMS and s_u_LMS are lower and upper bounds for s
  total_comb = combinations(k_value,p_value)        #if there is problem then increase s_u_LMS
  if(nR<=0){
    nR = dim(total_comb)[1]
    combr = total_comb
  }else{
    combr = total_comb[sample.int(dim(total_comb)[1],nR),]
  }
  k_1c = rep(1,dim(X1)[1])
  k_1 = rep(1,p_value)
  A_value = floor((k_value+1)/2)/k_value
  s_save = s_u_LMS
  for (i in 1:nR) {
    X1_sub = X1[combr[i,],]
    X2_sub = X2[combr[i,],]
    trans_X1 = colMeans(X1_sub)
    trans_X2 = colMeans(X2_sub)
    X1_tran = X1_sub-k_1%*%t(trans_X1)
    X2_tran = X2_sub-k_1%*%t(trans_X2)
    temp3 = svd(t(X2_tran)%*%X1_tran)
    rot = temp3$v%*%t(temp3$u)
    #scal = Re(tr(sqrtm(t(X1_tran)%*%X2_tran%*%t(X2_tran)%*%X1_tran))/tr(X1_tran%*%t(X1_tran)))
    scal = tr(t(X2_tran)%*%X1_tran%*%rot)/tr(X1_tran%*%t(X1_tran))
    ####################################
    X1_tran = X1-k_1c%*%t(trans_X1)
    X2_tran = X2-k_1c%*%t(trans_X2)
    E_i = X2_tran - scal*X1_tran%*%rot
    E_i = sqrt(rowSums(E_i*E_i))
    test = xi_func_LMS(s_save,A_value,E_i)
    if(is.na(test)) {next}
    if(test<0){
      s_now = uniroot(xi_func_LMS,c(s_l_LMS,s_u_LMS),A_value=A_value,E=E_i)$root
      if(s_now<s_save){
        s_save = s_now
        trans_X1_save = trans_X1
        trans_X2_save = trans_X2
        rot_save = rot
        scal_save = scal
      }
    }
    ####################################
  }
  #################################################################################
  n_trim = ceiling(k_value/2)
  X1_tran = X1-k_1c%*%t(trans_X1_save)
  X2_tran = X2-k_1c%*%t(trans_X2_save)
  E_i = X2_tran - scal_save*X1_tran%*%rot_save
  E_i = rowSums(E_i*E_i)
  ind = order(E_i)[1:n_trim]
  ##
  X1_sub = X1[ind,]
  X2_sub = X2[ind,]
  trans_X1_save = colMeans(X1_sub)
  trans_X2_save = colMeans(X2_sub)
  k_1_trim = rep(1,n_trim)
  X1_tran = X1_sub-k_1_trim%*%t(trans_X1_save)
  X2_tran = X2_sub-k_1_trim%*%t(trans_X2_save)
  temp3 = svd(t(X2_tran)%*%X1_tran)
  rot_save = temp3$v%*%t(temp3$u)
  #scal = Re(tr(sqrtm(t(X1_tran)%*%X2_tran%*%t(X2_tran)%*%X1_tran))/tr(X1_tran%*%t(X1_tran)))
  scal_save = tr(t(X2_tran)%*%X1_tran%*%rot_save)/tr(X1_tran%*%t(X1_tran))
  ##
  E_i = X2_tran - scal_save*X1_tran%*%rot_save
  s_save = sum(E_i*E_i)
  #################################################################################
  trans_P = trans_X2_save - scal_save*t(rot_save)%*%trans_X1_save    # new code
  return(list(rot=rot_save,scal=scal_save,trans_P=trans_P,s=s_save))  #new code
  #  return(list(rot=rot_save,scal=scal_save,trans_X1=trans_X1_save,trans_X2=trans_X2_save,s=s_save))
}


############################################################################
#############################################################################'Least trimmed squares (LTS) method
#'Align two matrices using the ILMS.
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#'@param nR Number of randomly sampled subsets.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1, and least trimmed squares.
#' @importFrom fBasics tr
#' @importFrom stats median
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @export

LTS_estimator <- function(X1,X2,nR) { #X1 is transformed; nR is the number of subsets;
  k_value = dim(X1)[1]                #nR<=0 means using all subsets
  p_value = dim(X1)[2]+1
  total_comb = combinations(k_value,p_value)
  if(nR<=0){
    nR = dim(total_comb)[1]
    combr = total_comb
  }else{
    combr = total_comb[sample.int(dim(total_comb)[1],nR),]
  }
  k_1c = rep(1,dim(X1)[1])
  k_1 = rep(1,p_value)
  n_trim = ceiling(k_value/2)
  s_save = 10^10
  for (i in 1:nR) {
    X1_sub = X1[combr[i,],]
    X2_sub = X2[combr[i,],]
    trans_X1 = colMeans(X1_sub)
    trans_X2 = colMeans(X2_sub)
    X1_tran = X1_sub-k_1%*%t(trans_X1)
    X2_tran = X2_sub-k_1%*%t(trans_X2)
    temp3 = svd(t(X2_tran)%*%X1_tran)
    rot = temp3$v%*%t(temp3$u)
    #scal = Re(tr(sqrtm(t(X1_tran)%*%X2_tran%*%t(X2_tran)%*%X1_tran))/tr(X1_tran%*%t(X1_tran)))
    scal = tr(t(X2_tran)%*%X1_tran%*%rot)/tr(X1_tran%*%t(X1_tran))
    ####################################
    X1_tran = X1-k_1c%*%t(trans_X1)
    X2_tran = X2-k_1c%*%t(trans_X2)
    E_i = X2_tran - scal*X1_tran%*%rot
    E_i = sort(rowSums(E_i*E_i))
    s_now = sum(E_i[1:n_trim])
    if(s_now<s_save){
      s_save = s_now
      trans_X1_save = trans_X1
      trans_X2_save = trans_X2
      rot_save = rot
      scal_save = scal
    }
    ####################################
  }
  trans_P = trans_X2_save - scal_save*t(rot_save)%*%trans_X1_save    # new code
  return(list(rot=rot_save,scal=scal_save,trans_P=trans_P,s=s_save))  #new code
  #  return(list(rot=rot_save,scal=scal_save,trans_X1=trans_X1_save,trans_X2=trans_X2_save,s=s_save))
}


########################################################################################
#'The standard Procrustes analysis
#'Align two matrices using Procrustes analysis
#'@param X1 A matrix representing the configuration to be transformed to best align with X2.
#'@param X2 A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @return
#' A list containing the estimated rotation matrix, scaling factor,translation vectors for X1, and the sum of square values.
#' @importFrom fBasics tr
#' @importFrom stats median
#' @importFrom gllvm gllvm
#' @importFrom gtools combinations
#' @examples
#' X1 <- matrix(c(
#'   0.3068822, -1.24752606,
#'   -2.1150143, 4.24520347,
#'   0.5271354, 1.64468795,
#'   -0.1264126, -0.08113392,
#'   -1.4749896, -2.26016700,
#'   -0.2347269, -1.27047483,
#'   -0.4833630, 1.43158341,
#'   -2.3499156, 0.85892382,
#'   -2.1021129, 0.40267502,
#'   0.5812934, -0.57210742,
#'   0.7397990, 1.54898604,
#'   -1.2173015, -0.65486036,
#'   -0.6093984, -1.16937177,
#'   0.3833907, -0.72611343,
#'   -0.8423056, 1.22427834
#' ), ncol = 2, byrow = TRUE)
#'
#' X2 <- matrix(c(
#'   1.96990726, -1.67860936,
#'   1.16952104, -3.19996198,
#'   -1.64680983, 3.40009837,
#'   1.51340633, 0.90762220,
#'   6.02678651, -1.29463718,
#'   2.56525182, -0.70441947,
#'   0.08388971, 3.87254957,
#'   4.11361309, 4.77780158,
#'   -0.91329668, 0.01870814,
#'   0.73714253, -0.50391759,
#'   -1.60842305, 3.07780435,
#'   3.90039741, 1.29926458,
#'   3.31127451, -0.28114418,
#'   0.85909201, -0.86365691,
#'   1.52112556, 4.02668249
#' ), ncol = 2, byrow = TRUE)
#' fit <- Procrustes_fun(X1,X2)
#' @export

Procrustes_fun <- function(X1,X2) { #X1 is transformed;
  k_1 = rep(1,dim(X1)[1])
  trans_X1 = colMeans(X1)
  trans_X2 = colMeans(X2)
  X1_tran = X1-k_1%*%t(trans_X1)
  X2_tran = X2-k_1%*%t(trans_X2)
  temp3 = svd(t(X2_tran)%*%X1_tran)
  rot = temp3$v%*%t(temp3$u)
  scal = tr(t(X2_tran)%*%X1_tran%*%rot)/tr(X1_tran%*%t(X1_tran))
  E_i = X2_tran - scal*X1_tran%*%rot
  s_save = sum(E_i*E_i)
  trans_P = trans_X2 - scal*t(rot)%*%trans_X1
  return(list(rot=rot,scal=scal,trans_P=trans_P,s=s_save))
}





#######################################################
#######################################################
####permutation test for the eight methods
#' Repeated Median Method Permutation Test
#' Performs a permutation test for the repeated median method.
#' @param P_mat A matrix representing the fixed configuration that remains
#'unchanged during the transformation.
#' @param Q_mat A matrix representing the configuration to be transformed to
#'best align with `P_mat`.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the repeated median method.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @export

MM2_permute <- function(P_mat,Q_mat,P_N,S0,verbose = FALSE){


  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= MM2_estimator_with_reflection(Q_mat,P_mat_new)
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S))
}
##----------------------------------------------

#' Huber Method Permutation Test
#' Performs a permutation test for the Huber method.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the Huber method.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @export

Huber_permute <- function(P_mat,Q_mat,P_N,S0,verbose = FALSE){
  converge_vec = rep(0,P_N)

  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= majorization_tuned(P_mat_new,Q_mat,1)
    converge_vec[i] = result_1$converge
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S,converge_vec=converge_vec))
}

##------------------------------------------------------
#' Biweight Method Permutation Test
#' Performs a permutation test for the Biweight method.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the Biweight method.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @export

Biweight_permute <- function(P_mat,Q_mat,P_N,S0,verbose = FALSE){
  converge_vec = rep(0,P_N)

  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<=P_N) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= majorization_tuned(P_mat_new,Q_mat,2)
    converge_vec[i] = result_1$converge
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S,converge_vec=converge_vec))
}


##----------------------------------------------------------

#' Standard Procrustes Analysis Permutation Test
#' Performs a permutation test for the standard Procrustes analysis.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the Procrustes analysis.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @examples
#' P_mat <- matrix(c(
#'   0.3068822, -1.24752606,
#'   -2.1150143, 4.24520347,
#'   0.5271354, 1.64468795,
#'   -0.1264126, -0.08113392,
#'   -1.4749896, -2.26016700,
#'   -0.2347269, -1.27047483,
#'   -0.4833630, 1.43158341,
#'   -2.3499156, 0.85892382,
#'   -2.1021129, 0.40267502,
#'   0.5812934, -0.57210742,
#'   0.7397990, 1.54898604,
#'   -1.2173015, -0.65486036,
#'   -0.6093984, -1.16937177,
#'   0.3833907, -0.72611343,
#'   -0.8423056, 1.22427834
#' ), ncol = 2, byrow = TRUE)
#'
#' Q_mat <- matrix(c(
#'   1.96990726, -1.67860936,
#'   1.16952104, -3.19996198,
#'   -1.64680983, 3.40009837,
#'   1.51340633, 0.90762220,
#'   6.02678651, -1.29463718,
#'   2.56525182, -0.70441947,
#'   0.08388971, 3.87254957,
#'   4.11361309, 4.77780158,
#'   -0.91329668, 0.01870814,
#'   0.73714253, -0.50391759,
#'   -1.60842305, 3.07780435,
#'   3.90039741, 1.29926458,
#'   3.31127451, -0.28114418,
#'   0.85909201, -0.86365691,
#'   1.52112556, 4.02668249
#' ), ncol = 2, byrow = TRUE)
#' result_Procrustes <- Procrustes_fun(P_mat,Q_mat)
#' S0 <- result_Procrustes$s
#' P_N <- 10
#' permute_fit <- Procrustes_permute(P_mat,Q_mat,P_N,S0)
#' @export

Procrustes_permute <- function(P_mat,Q_mat,P_N,S0,verbose = FALSE){


  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result = Procrustes_fun(P_mat_new,Q_mat)
    S[i] = result$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S))
}
##-----------------------------------------------------------

#' S-estimator Permutation Test
#' Performs a permutation test for the S-estimator method.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the S-estimator method.
#' @param nR Number of randomly sampled subsets. If `nR <= 0`, all possible subsets are used.
#' @param s_l Lower bound for the scale s.
#' @param s_u Upper bound for the scale s.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @export

S_permute <- function(P_mat,Q_mat,P_N,S0,nR=0,s_l=0.000000001,s_u=10,verbose = FALSE){


  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= S_estimator(P_mat_new,Q_mat,c_value,A_value,nR,s_l,s_u)
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S))
}

##------------------------------------------------------------

#' The Least Median of Squares (LMS) Permutation Test
#' Performs a permutation test for the LMS method.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the LMS method.
#' @param nR Number of randomly sampled subsets. If `nR <= 0`, all possible subsets are used.
#' @param s_l_LMS Lower bound for the scale s.
#' @param s_u_LMS Upper bound for the scale s.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @export

LMS_permute <- function(P_mat,Q_mat,P_N,S0,nR=0,s_l_LMS=0.000000001,s_u_LMS=10,verbose = FALSE){


  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= LMS_estimator(P_mat_new,Q_mat,nR,s_l_LMS,s_u_LMS)
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S))
}

##---------------------------------------------------------------------

#' The Least Trimmed Squares (LTS) Permutation Test
#' Performs a permutation test for the LTS method.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the LTS method.
#' @param nR Number of randomly sampled subsets. If `nR <= 0`, all possible subsets are used.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @export

LTS_permute <- function(P_mat,Q_mat,P_N,S0,nR=0,verbose = FALSE){


  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= LTS_estimator(P_mat_new,Q_mat,nR)
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S))
}



##------------------------------------------------------------

#' The Improved Least Median of Squares (ILMS) Permutation Test
#' Performs a permutation test for the ILMS method.
#'@param P_mat A matrix representing the configuration to be transformed to best align with Q_mat.
#'@param Q_mat A matrix representing the fixed configuration that remains unchanged during the transformation.
#' @param P_N Number of permutations.
#' @param S0 Observed objective function value from the ILMS method.
#' @param nR Number of randomly sampled subsets. If `nR <= 0`, all possible subsets are used.
#' @param s_l_LMS Lower bound for the scale s.
#' @param s_u_LMS Upper bound for the scale s.
#' @param verbose Logical; if \code{TRUE}, displays progress messages during iterations. Defaults to \code{FALSE}.
#' @return
#' A list containing:
#' \item{count_n}{The number of permuted objective function values less than or
#'   equal to the observed value.}
#' \item{p_value}{The permutation test p-value.}
#' \item{S}{A vector of objective function values obtained from the
#'   permutations.}
#' @examples
#' P_mat <- matrix(c(
#'   0.3068822, -1.24752606,
#'   -2.1150143, 4.24520347,
#'   0.5271354, 1.64468795,
#'   -0.1264126, -0.08113392,
#'   -1.4749896, -2.26016700,
#'   -0.2347269, -1.27047483,
#'   -0.4833630, 1.43158341,
#'   -2.3499156, 0.85892382,
#'   -2.1021129, 0.40267502,
#'   0.5812934, -0.57210742,
#'   0.7397990, 1.54898604,
#'   -1.2173015, -0.65486036,
#'   -0.6093984, -1.16937177,
#'   0.3833907, -0.72611343,
#'   -0.8423056, 1.22427834
#' ), ncol = 2, byrow = TRUE)
#'
#' Q_mat <- matrix(c(
#'   1.96990726, -1.67860936,
#'   1.16952104, -3.19996198,
#'   -1.64680983, 3.40009837,
#'   1.51340633, 0.90762220,
#'   6.02678651, -1.29463718,
#'   2.56525182, -0.70441947,
#'   0.08388971, 3.87254957,
#'   4.11361309, 4.77780158,
#'   -0.91329668, 0.01870814,
#'   0.73714253, -0.50391759,
#'   -1.60842305, 3.07780435,
#'   3.90039741, 1.29926458,
#'   3.31127451, -0.28114418,
#'   0.85909201, -0.86365691,
#'   1.52112556, 4.02668249
#' ), ncol = 2, byrow = TRUE)
#' result_ILMS <- LMS_LTS_estimator(P_mat,Q_mat,nR=0,s_l_LMS=0.000000001,s_u_LMS=10)
#' S0 <- result_ILMS$s
#' P_N <- 10
#' permute_fit <- ILMS_permute(P_mat,Q_mat,P_N,S0,nR=0,s_l_LMS=0.000000001,s_u_LMS=10)
#' @export

ILMS_permute <- function(P_mat,Q_mat,P_N,S0,nR=0,s_l_LMS=0.000000001,s_u_LMS=10,verbose = FALSE){

  count_n = 0
  i=1
  S= rep(0,P_N)
  while (i<P_N+1) {

    P_mat_new = P_mat[sample.int(n=dim(P_mat)[1],size=dim(P_mat)[1],replace=FALSE),]
    result_1= LMS_LTS_estimator(P_mat_new,Q_mat,nR,s_l_LMS,s_u_LMS)
    S[i] = result_1$s
    if(S[i]<=S0){
      count_n = count_n+1
    }
    if (verbose) {
      message(paste0("Permutation ", i, " of ", P_N))
    }
    i =i+1
  }
  p_value = (count_n+1)/(P_N+1)

  return(list(count_n=count_n,p_value=p_value,S=S))
}



