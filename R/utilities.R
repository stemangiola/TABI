sigmoid_4_param = function(x, A, y_cross, slope, inflection){
  
  if(y_cross<0) stop("y_cross must be > 0")
  
  x = matrix(x, ncol = 1)
  
  A + 
    (
      (y_cross) * 
        exp(   
          log1p_exp(inflection * slope[1]) -
            log1p_exp(   -( x %*% slope ) + inflection * slope[1]  ) 
          
        )
    ) 
  
}

log1p_exp = function(x){ log(1+exp(x)) }


