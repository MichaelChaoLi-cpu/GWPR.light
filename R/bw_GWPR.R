#' Bandwidth selection for basic GWPR
#'
#' @description A function for automatic bandwidth selection to calibrate a GWPR model
#'
#' @usage bw.GWPR(formula, data, index, SDF, adaptive = FALSE, p = 2, bigdata = FALSE,
#'                upperratio = 0.25, effect = "individual",
#'                model = c("pooling", "within", "random"),
#'                random.method = "swar", approach = c("CV","AIC"), kernel = "bisquare",
#'                longlat = FALSE, doParallel = FALSE, cluster.number = 2,
#'                human.set.range = FALSE, h.upper = NULL, h.lower = NULL,
#'                gradientIncrement = FALSE, GI.step = NULL, GI.upper = NULL,
#'                GI.lower = NULL)
#'
#' @param formula            The regression formula: : Y ~ X1 + ... + Xk
#' @param data               data.frame for the Panel data
#' @param index              A vector of the two indexes: (c("ID", "Time"))
#' @param SDF                Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param adaptive           If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param p                  The power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param bigdata            TRUE or FALSE, if the dataset exceeds 40,000, we strongly recommend set it TRUE
#' @param upperratio         Set the ratio between upper boundary of potential bandwidth range and the forthest distance of SDF, if bigdata = T. (default value: 0.25)
#' @param effect             The effects introduced in the model, one of "individual" (default) , "time", "twoways", or "nested"
#' @param model              Panel model transformation: (c("within", "random", "pooling"))
#' @param random.method      Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param approach           Score used to optimize the bandwidth, c("CV", "AIC")
#' @param kernel             bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                            gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                            exponential: wgt = exp(-vdist/bw);
#'                            tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                            boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param longlat            If TRUE, great circle distances will be calculated
#' @param doParallel         If TRUE, "cluster": multi-process technique with the parallel package would be used.
#' @param cluster.number     The number of the clusters that user wants to use
#' @param human.set.range    If TRUE, the range of bandwidth selection for golden selection could be set by the user
#' @param h.upper            The upper boundary of the potential bandwidth range for golden selection.
#' @param h.lower            The lower boundary of the potential bandwidth range for golden selection.
#' @param gradientIncrement  The bandwidth selection become gradient increment, if TRUE
#' @param GI.step            The step length of the increment.
#' @param GI.upper           The upper boundary of the gradient increment selection.
#' @param GI.lower           The lower boundary of the gradient increment selection.
#'
#' @return The optimal bandwidth
#'
#' @export
#'
#' @import dplyr
#' @importFrom data.table setDT
#' @importFrom sp coordinates bbox
#' @import parallel
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#'
#' @examples
#' \donttest{
#' data(TransAirPolCalif)
#' data(California)
#' formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
#'    Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
#'    Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
#'    Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
#'    Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
#'    pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax
#'
#' bw.CV.Fix <- bw.GWPR(formula = formula.GWPR, data = TransAirPolCalif,
#'                      index = c("GEOID", "year"),
#'                      SDF = California, adaptive = FALSE, p = 2, bigdata = FALSE,
#'                      effect = "individual", model = "within", approach = "CV",
#'                      kernel = "bisquare", longlat = FALSE,
#'                      gradientIncrement = TRUE, GI.step = 0.5, GI.upper = 5,
#'                      GI.lower = 1.5)
#'
#' bw.CV.Fix
#'
#' bw.AIC.Fix <- bw.GWPR(formula = formula.GWPR, data = TransAirPolCalif,
#'                       index = c("GEOID", "year"),
#'                       SDF = California, adaptive = FALSE, p = 2, bigdata = FALSE,
#'                       effect = "individual", model = "within", approach = "AIC",
#'                       kernel = "bisquare", longlat = FALSE, doParallel = FALSE)
#' bw.AIC.Fix
#' }
bw.GWPR <- function(formula, data, index, SDF, adaptive = FALSE, p = 2, bigdata = FALSE, upperratio = 0.25,
                    effect = "individual", model = c("pooling", "within", "random"), random.method = "swar",
                    approach = c("CV","AIC"), kernel = "bisquare", longlat = FALSE, doParallel = FALSE,
                    cluster.number = 2, human.set.range = FALSE, h.upper = NULL, h.lower = NULL,
                    gradientIncrement = FALSE, GI.step = NULL, GI.upper = NULL, GI.lower = NULL)
{
  if(length(index) != 2)
  {
    stop("The \"index\" have included \"ID\" or/and \"time\" index.")
  }
  if(!((index[1] %in% colnames(data)) & (index[2] %in% colnames(data))))
  {
    stop("The data.frame(data) does not have the index columns.")
  }
  if(!(index[1] %in% colnames(SDF@data)))
  {
    stop("The SDF does not have the \"ID\" columns.")
  }
  if(!(model %in% c("pooling", "within", "random")))
  {
    stop("This version GWPR only accept \"pooling\", \"within\", or \"random\"")
  }
  if(!(approach %in% c("CV", "AIC")))
  {
    stop("This version GWPR only accept \"CV\"and \"AIC\" approach ")
  }

  # Data preparation
  varibale_name_in_equation <- all.vars(formula)
  #data <- dplyr::select(data, index, varibale_name_in_equation)
  data <- dplyr::select(data, dplyr::all_of(index), dplyr::all_of(varibale_name_in_equation))
  ##### 22.6.17 we change this, it is the problem on linux
  data$raw_order_data <- 1:nrow(data)
  raw_id <- index[1]
  colnames(data)[1] <- "id"
  index[1] <- "id"

  # Assuming unbalanced panel, get individuals' ID and max record number of individuals
  ID <- dplyr::select(data, "id")
  .N <- 0
  ID_num <- data.table::setDT(ID)[,list(Count = .N),names(ID)]
  if(model == "within" )
  {
    data <- drop_ID_with_single_observation(data, ID_num)
    ID <- dplyr::select(data, "id")
    ID_num <- data.table::setDT(ID)[,list(Count = .N),names(ID)]
  }

  # Judge the datasize of calculation
  if ((nrow(ID_num) > 1000) & !doParallel)
  {
    message("Dear my friend, more than 1,000 individuals in your dataset.\n",
        "It would be time-consuming. We use \"-\" and \"|\" to inform you where\n",
        "we are. A \"-\" is equal to 2.5% in once score calculation. A \"|\" is\n",
        "25%. Do not feel nervous or boring! Your research is on the way!\n",
        ".............................................................\n")
    huge_data_size <- TRUE
  }
  else
  {
    huge_data_size <- FALSE
  }

  # Panel SDF preparation
  SDF@data <- dplyr::select(SDF@data, dplyr::all_of(raw_id))
  colnames(SDF@data)[1] <- "id"
  dp.locat <- sp::coordinates(SDF)
  coord <- cbind(as.data.frame(dp.locat), SDF@data$id)
  SDF <- 0 # drop the SDF
  colnames(coord) <- c("X", "Y", "id")
  data <- dplyr::left_join(data, coord, by = "id")
  lvl1_data <- data # data put into calculation
  #### data have been prepared

  #### set range of selection
  if (gradientIncrement)
  {
    if (is.null(GI.upper) | is.null(GI.lower) | is.null(GI.step))
    {
      stop("Please input upper, lower boundaries (GI.upper and GI.lower) and step length (GI.step) of GI")
    }
    message("Since GI method is used, so the GI.upper is the real upper boundary: ",
            GI.upper, " lower boundary: ", GI.lower," step length: ", GI.step)
    lower <- GI.lower
    upper <- GI.upper
  } ### gradientIncrement selection
  else
  {
    if(human.set.range)
    {
      message("...............................................................................................\n",
              "Now, the range of bandwidth selection for golden selection is set by the user\n",
              "We assume that the user is familiar with bandwidth selection, and uses this setting to reduce calculation time.\n",
              "If not, please stop the current calculation, and set \"human.set.range\" as FALSE\n",
              "...............................................................................................\n")
      if(is.null(h.lower))
      {
        stop("Need to set the lower boundary (h.lower)!")
      }
      if(is.null(h.upper))
      {
        stop("Need to set the upper boundary (h.upper)!")
      }
      if(h.lower > h.upper)
      {
        stop("h.lower should be smaller than h.upper")
      }
      lower <- h.lower
      upper <- h.upper
    } ### set human range for golden selection
    else
    {
      # decide upper and lower boundary
      lower.freedom <- protect_model_with_enough_freedom(formula = formula, data = lvl1_data, ID_list = ID_num,
                                                         index = index, p = p, longlat = longlat)
      message("To make sure every subsample have enough freedom, the minimum number of individuals is ",lower.freedom, "\n")
      if(adaptive)
      {
        upper <- nrow(ID_num)
        lower <- lower.freedom + 1
        if((model == "random")&(random.method == "swar"))
        {
          lower <- length(varibale_name_in_equation) + 1
        }
      }
      else
      {
        b.box <- sp::bbox(dp.locat)
        upper <- sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower <- upper/5000
        if ((model == "random")&(random.method == "swar"))
        {
          lower <- protect_model_with_least_individuals(lvl1_data, ID_num, index, kernel, p, longlat,
                                                        bw_panel = (length(varibale_name_in_equation) + 1))
        }
        else
        {
          lower <- protect_model_with_least_individuals(lvl1_data, ID_num, index,
                                                        kernel, p, longlat, bw_panel = lower.freedom)
        }
      }

      if(bigdata)
      {
        upper <- upper * upperratio
        lower <- lower
      }

      if(bigdata)
      {
        message("You set the \"bigdata\" is: ", bigdata, ". The ratio is: ", upperratio, "\n",
                "Now the lower boundary of the bandwidth selection is ", lower, ", and upper boundary is ", upper,".\n",
                "Note: if the optimal bandwidth is close to the upper boundary, you need to increase the ratio.\n",
                "However, you should also know that the larger ratio requires more memory. Please, balance them and enjoy your research.\n")
      }
      if(huge_data_size)
      {
        message("Data Prepared! Go!............................................\n")
      }
    }
    message("The upper boundary is ", upper,", and the lower boundary is ", lower,"\n")
  } ### automatically golden select bandwidth

  #### begin to select
  if(doParallel)
  {
    message("..................................................................................\n")
    message("You use parallel process, so be careful about your memory usage. Cluster number: ", cluster.number,"\n")
    if(cluster.number > parallel::detectCores())
    {
      stop("The cluster number exceeds the cores you have")
    }
    else
    {
      if(cluster.number > (parallel::detectCores()-2))
      {
        warning("You might use too many cluster, only one left for other task.")
      }
    }
    #0.2.0
    if (gradientIncrement)
    {
      if (is.null(GI.upper) | is.null(GI.lower) | is.null(GI.step))
      {
        stop("Please input upper, lower boundaries (GI.upper and GI.lower) and step length (GI.step) of GI")
      }
      BandwidthVector <- c()
      ScoreVector <- c()
      if (adaptive)
      {
        if(approach == "CV")
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- CV_A_para(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                               formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                               model = model, index = index, kernel = kernel, effect = effect,
                               random.method = random.method,  cluster.number = cluster.number)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
        else
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- AIC_A_para(bw = bw.now, data_input = lvl1_data, ID_list = ID_num,
                               formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                               model = model, index = index, kernel = kernel, effect = effect,
                               random.method = random.method,  cluster.number = cluster.number)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
      }
      else
      {
        if(approach == "CV")
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- CV_F_para(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                               formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                               model = model, index = index, kernel = kernel, effect = effect,
                               random.method = random.method,  cluster.number = cluster.number)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
        else
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- AIC_F_para(bw = bw.now, data_input = lvl1_data, ID_list = ID_num,
                               formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                               model = model, index = index, kernel = kernel, effect = effect,
                               random.method = random.method,  cluster.number = cluster.number)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
      }
      bw <- BandwidthSocreTable
    }
    else
    {
      if (adaptive)
      {
        if(approach == "CV")
        {
          bw <- gold(CV_A_para, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method,  cluster.number = cluster.number)
        }
        else
        {
          bw <- gold(AIC_A_para, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, cluster.number = cluster.number)
        }
      }
      else
      {
        if(approach == "CV")
        {
          bw <- gold(CV_F_para, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method,  cluster.number = cluster.number)
        }
        else
        {
          bw <- gold(AIC_F_para, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method,  cluster.number = cluster.number)
        }
      }
      return(bw)
    }
  } #### do parallel
  else
  { #### not do parallel
    #0.2.0
    if (gradientIncrement)
    {
      if (is.null(GI.upper) | is.null(GI.lower) | is.null(GI.step))
      {
        stop("Please input upper, lower boundaries (GI.upper and GI.lower) and step length (GI.step) of GI")
      }
      BandwidthVector <- c()
      ScoreVector <- c()
      if (adaptive)
      {
        if(approach == "CV")
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- CV_A(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                          formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                          model = model, index = index, kernel = kernel, effect = effect,
                          random.method = random.method, huge_data_size = huge_data_size)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
        else
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- AIC_A(bw = bw.now, data_input = lvl1_data, ID_list = ID_num,
                          formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                          model = model, index = index, kernel = kernel, effect = effect,
                          random.method = random.method, huge_data_size = huge_data_size)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
      }
      else
      {
        if(approach == "CV")
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- CV_F(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                          formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                          model = model, index = index, kernel = kernel, effect = effect,
                          random.method = random.method, huge_data_size = huge_data_size)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
        else
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- AIC_F(bw = bw.now, data_input = lvl1_data, ID_list = ID_num,
                          formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                          model = model, index = index, kernel = kernel, effect = effect,
                          random.method = random.method, huge_data_size = huge_data_size)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        }
      }
      bw <- BandwidthSocreTable
    } ### gradient increment

    #0.1.2 /|\
    else
    {
      if (adaptive)
      {
        if(approach == "CV")
        {
          bw <- gold(CV_A, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
        else
        {
          bw <- gold(AIC_A, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
      }
      else
      {
        if(approach == "CV")
        {
          bw <- gold(CV_F, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
        else
        {
          bw <- gold(AIC_F, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
      }
    }
  }
  #v0.1.1
  return(bw)
}
