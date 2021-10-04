#' Bandwidth Golden Selection
#'
#' @description This is golden selection function to select the optimal bandwidth. Created by MC, edited by BL, further improved by CL.
#'
#' @author M.C.
#'
#' @param fun       The function to calcalate CV score or AIC score of a specific bandwidth.
#' @param xL        The lower boundary of the potential bandwidth range.
#' @param xU        The upper boundary of the potential bandwidth range.
#' @param adapt.bw  If TRUE, apadtive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param ...       The input parameters of the "fun"
#'
#' @return Optimal bandwidth for GWPR
#' @noRd
gold <- function(fun,xL,xU,adapt.bw=F,...)
{
  lower.limitation <- xL
  upper.limitation <- xU
  if(adapt.bw)
  {
    eps = 1e-4
  }
  else
  {
    eps = (upper.limitation - lower.limitation) * 1e-6 * 2
  } #even global scale the eps is roughly 1km
  R <- (sqrt(5)-1)/2  # R <- 0.61803398....
  iter <- 1
  d <- R*(xU-xL)
  if (adapt.bw)
  {
    x1 <- floor(xL+d)
    x2 <- round(xU-d)
  }
  else
  {
    x1 <- xL+d
    x2 <- xU-d
  }
  f1 <- eval(fun(x1,...))
  f2 <- eval(fun(x2,...))
  d1<-f2-f1
  # Establish initial value of xopt:
  if (f1 < f2)
    xopt  <-  x1
  else xopt  <-  x2
  # start main loop
  ea <- 100
  while ((abs(d) > eps) & (abs(d1) > eps) )
  {
    d <- R*d
    if   (f1 < f2)
    {
      xL  <-  x2
      x2  <-  x1
      if (adapt.bw)
      {
        x1 <- round(xL+d)
      }
      else
      {
        x1 <- xL+d
      }
      #x1  <-  xL + d
      if(x1 < upper.limitation)
      {
        f2  <-  f1
        f1 <- eval(fun(x1,...))
      }
      else
      {
        f2  <-  f1
        x1 <- upper.limitation
        f1 <- eval(fun(x1,...))
      } # avoid the xopt spikes the upper boundary
    }
    else
    {
      xU  <-  x1
      x1  <-  x2
      if (adapt.bw)
        x2 <- floor(xU - d)
      else
        x2  <-  xU - d
      if(x2 > lower.limitation)
      {
        f1  <-  f2
        f2 <- eval(fun(x2,...))
      }
      else
      {
        f1  <-  f2
        x2  <- lower.limitation
        f2 <- eval(fun(x2,...))# avoid the xopt spikes the lower boundary
      }
    }
    iter  <-  iter + 1 # iter is the test parameter in the first version (CL)
    # Establish value of xopt after iteration:
    if    (f1 < f2)
    {
      xopt  <-  x1
    }
    else
    {
      xopt  <-  x2
    }
    d1<-f2-f1
  }
  xopt
}
