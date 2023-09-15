#' @importFrom stats pnorm
#' @importFrom stats sd
#' @export
ADStatQF <-
function(x) {
    #The function is used to calculate the Anderson Darling Test Statistics of Standard Normal Distribution
    #vectors tested need lengths greater than 7
    n = length(x)
    if (n < 7) {
        warning('Sample Size must be greater than 7.')
    }
    else {
        test1 = sort(x)
        f = pnorm(test1, mean = mean(test1), sd = sd(test1))
        i = 1:n
        S = sum(((2*i)-1) / n * (log(f)+log(1-f[n+1-i])))
        ADval = -n-S
        return(ADval)
    }
}
