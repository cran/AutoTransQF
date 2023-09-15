#' @importFrom VGAM qgumbel
#' @importFrom stats IQR
#' @importFrom stats median
#' @importFrom stats sd
#' @export

autotransfuncQF <-
function(vari, istat, paraindex) {
  #The function is used for transformation of each individual variable
  #The variable paraindex is used to calculate transformation parameter beta
  i = paraindex
  row_num = dim(matrix(vari, nrow = 1))[1]
  col_num = dim(matrix(vari, nrow = 1))[2]
  beta = sign(i)*(exp(abs(i))-1)
  alpha = abs(1/beta)
  
  #get the interquartile range of the values in the data vector
  iqrange_vari = IQR(vari, type = 2)
  
  #use the interquartile range as the data range if not zero 
  #otherwise use the full data range
  if(iqrange_vari == 0) {
    range_vari = range(vari)
  }else{
    range_vari = iqrange_vari
  }
  
  #transform based on the sign of the beta
  if (beta  == 0) {
    vark1 = vari
  }else if (beta > 0) {
    vark1 = log(vari-min(vari)+alpha*range_vari)
  }else {
    vark1 = -log(max(vari)-vari+alpha*range_vari)
  }
  
  #standardize transformed data vector
  mabd = mean(abs(vark1 - median(vark1))) 
  mabd = mabd* sqrt(pi / 2)
  
  if (mabd == 0) {
    final_vari = matrix(0, row_num, col_num)
  }else {
    vark1 = (vark1 - median(vark1)) / mabd
    # winsorise transformed data vector
    p95 = -qgumbel(0.05)
    ak = (2*log(row_num*col_num))^(-0.5)
    bk = (2*log(row_num*col_num) - log(log(row_num*col_num)) - log(4*pi))^0.5
    TH = p95 * ak + bk

    
    
    vark1[vark1 > TH] = rep(TH, sum(vark1 > TH))
    vark1[vark1 < -TH] = rep(-TH, sum(vark1 < -TH))
    
    #re-standardize by substract the mean and divide by standard deviation
    final_vari  = (vark1 - mean(vark1)) / sd(vark1)
  }
  
  #calculate evaluation stat
  #both AD statistics and skewness test how far away data is from normality
  if(istat == 1) {
    eval = ADStatQF(final_vari) 
  }else {
    eval = skewness(final_vari) 
  }
  return(final_vari)
}
