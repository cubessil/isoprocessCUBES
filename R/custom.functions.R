#' Multiple plot function
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @export
multiplot <- function(..., plotlist=NULL, file=NULL, cols=1, layout=NULL, width=8.5, height=11) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    if(!is.null(file)) pdf(file=file,width,height)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
    if(!is.null(file)) dev.off()
  }
}

#############################################
# calculate means of subsets of a dataframe #
#############################################

#' Calculate means using group_by type
#'
#' @export
calc_means <- function(df, x_col) {
  df$x <- df[[x_col]]
  df %>%
    group_by(type) %>%
    summarise(
      n = n(),
      mean = mean(x),
      `mean + sigma` = mean + sd(x),
      `mean - sigma` = mean - sd(x),
      `mean + 2 sigma` = mean + 2 * sd(x),
      `mean - 2 sigma` = mean - 2 * sd(x)
    ) %>%
    gather(linetype, yintercept, -type, -n) %>%
    filter(!is.na(yintercept))
}

#############################################
# calculate means of subsets of a dataframe #
#############################################

#' Calculate means using groub_by(`Identifier 2`)
#'
#' @export
calc_means_dxf <- function(df, x_col) {
  df$x <- df[[x_col]]
  df %>%
    group_by(`Identifier 2`) %>%
    summarise(
      n = n(),
      mean = mean(x),
      `mean + sigma` = mean + sd(x),
      `mean - sigma` = mean - sd(x),
      `mean + 2 sigma` = mean + 2 * sd(x),
      `mean - 2 sigma` = mean - 2 * sd(x)
    ) %>%
    gather(linetype, yintercept, -type, -n) %>%
    filter(!is.na(yintercept))
}
#########################################
# functions for clumped data processing #
#########################################

############
### Calculating external errors after initial york regressions
###########

#' Calc ext err
#'
#' Description
#'
#' @param df this will usually be a datafram with heated gas data
#' @param mass either 47 or 48
#' @param slope york slope from the initial regression
#' @param intercept york intercept from the initial regression
#' @param N number of heated gases
#' @export
calc_exterr<-function(df, mass, slope, intercept, N){
  Chi.value<-(df[[paste0("D",mass)]]-(slope * df[[paste0("d",mass)]] + intercept))/df[[paste0("D",mass,".sterr")]]
  df[[paste0("d",mass,"exterr")]]<-sqrt(df[[paste0("d",mass,".sterr")]]*abs(sum(Chi.value))/(N-2))
  df[[paste0("D",mass,"exterr")]]<-sqrt(df[[paste0("D",mass,".sterr")]]*abs(sum(Chi.value))/(N-2))
  return(df)
}

###############
## calculating confidence and predictive intervals for york regressions
###############

#' York conf predictio interval
#'
#' Description
#'
#' @param reg.df dataframe that the york regression is based on - might be HG only or same as all data
#' @param all.df all the data of interest for constructing confidence intervals; might be same as reg.df
#' @param mass typically either 47 or 48
#' @param slope york slope
#' @param intercept york intercept
#' @param conf.level level of confidence interval, e.g. use 0.95 for 95\% confidence
#' @export
york.conf.pred.interval<-function(reg.df, all.df, mass, slope, intercept, conf.level) {
  line.predicted.Y<-slope * reg.df[[paste0("d",mass)]] + intercept #calculates predicted D48 value for just the regression data
  W <- sqrt( 2 * qf(conf.level, 2, length(line.predicted.Y)-2) )   #the Working-Hotelling multiplier
  SE.est.pred <- sqrt(sum((reg.df[[paste0("D",mass)]]-line.predicted.Y)^2)/(length(line.predicted.Y)-2))*sqrt(1+1/length(line.predicted.Y)+(all.df[[paste0("d",mass)]]-mean(reg.df[[paste0("d",mass)]]))^2/sum((reg.df[[paste0("d",mass)]]-mean(reg.df[[paste0("d",mass)]]))^2)) #calculates predictive interval parameter
  SE.est.conf <- sqrt(sum((reg.df[[paste0("D",mass)]] - line.predicted.Y)^2)/(length(line.predicted.Y) - 2)) * sqrt(1/length(line.predicted.Y)+(all.df[[paste0("d",mass)]]-mean(reg.df[[paste0("d",mass)]]))^2/sum((reg.df[[paste0("d",mass)]]-mean(reg.df[[paste0("d",mass)]]))^2)) #calculates confidence interval parameter
  conf.int.df<-all.df
  conf.int.df[[paste0("D",mass,".conf.int")]] <- W*SE.est.conf #calculates the amount of D48 excess that falls within the confidence interval
  conf.int.df[[paste0("D",mass,".pred.int")]] <- W*SE.est.pred
  conf.int.df[["all.predicted.Y"]]<-slope * all.df[[paste0("d",mass)]] + intercept #predicted D48 value from line for all data points, needed for plotting data
  return(conf.int.df)
}


#####################
## calculate temps and se from D47 values
######################

#' Original Ghosh calibration, using CIT reference frame values
#' @export
convert_CIT.D47_to_CIT.Ghosh.temp <- function(D47) {
  round((59200/(D47+0.02))^0.5-273.15, 1)
}

#' Ghosh calibration in ARF ref frame, using ARF D47 values
#' @export
convert_ARF.D47_to_ARF.Ghosh.temp <- function(D47) {
  round((63600/(D47+0.0047))^0.5-273.15, 1)
}

#' Dennis Calibration in ARF ref frame, using ARF D47 values
#' @export
convert_ARF.D47_to_ARF.Dennis.temp <- function(D47) {
  round((36200/(D47-0.292))^0.5-273.15, 1)
}


#' Calc D47 temp se
#' @export
calc_D47.Temp_se <- function(Temp, D47, D47se) {
  round(sqrt(abs((((Temp + 273.15)^6)/4)*(6.35164845416249E-13 + 2 * D47 * - 1.03937857088367E-12 +
                                                          ((D47)^2) * 1.7284993474114E-12 + 0.0000169461014527562^2 * D47se^2))), 1)
}

#' Calc d18
#' @export
calc_d18Ow <- function(d18Om, Temp) {
  round((((1.03092 * d18Om + 30.92) + 1000)/(exp((((18.03 * 10^3)/(Temp + 273.15)) - 32.42) / 10^3))) - 1000, 1)
}

#' Calc d18O
#' @export
calc_d18Owse <- function (d18Om, d18Omse, Temp, Tempse) {
  round(sqrt((((18.03 * ((1.03092 * d18Om + 30.92) + 1000) * exp(0.03242)) / (exp(18.03 / (Temp + 273.15)) * (Temp + 273.15)^2) * Tempse)^2) + ((exp(0.03242) / exp(18.03/(Temp + 273.15)) * d18Omse)^2)), 1)
}

calc_d18Ow_daeron <- function (d18Om, Temp) {
  A = 17.57
  B = 29.13
  T_c = Temp + 273.15
  d18O_vsmow = (1.03092 * d18Om + 30.92)
  
  numerator = d18O_vsmow + 1000
  alpha  = exp((A * 1000 / T_c - B)/1000)
  d18Ow = round(((numerator/alpha) - 1000),1)
  
  return (d18Ow)
}


calc_d18Owse_daeron <- function (d18Om, d18Omse, Temp, Tempse) {
  A = 17.57
  B = 29.13
  T_c = Temp + 273.15
  d18O_vsmow = (1.03092 * d18Om + 30.92)
  
  first_part = (((A *((d18O_vsmow) + 1000) * exp(B/1000)) / (exp(A /T_c) * T_c^2) * Tempse)^2)  
  second_part = ((exp(B/1000) / exp(A/(T_c)) * d18Omse)^2)
  d18Owse = round(sqrt(first_part +  second_part),1)
  
  return(d18Owse)
}
#####################################
## York Regression in general form ##
#####################################

#' York reg general form
#'
#' script for generating a york regression of a dataset with errors in both x and y
#' @export
york.regression<-function(X,x.err,Y,y.err,error.corr)  #for clumps, X is d47, Y is D47, and the term "error.corr" is for how correlated the y.err are with the x.err - for my purposes, I usually use a value of 0
{  					#opens the function, everything with one "<"in here is internal, two "<" means I can call the value by naming it
  weightsX<-1/(x.err^2)	#set up to input error rather than variance; could use either stdev or sterr in this depending on how your statistics are being used
  weightsY<-1/(y.err^2)
  alpha<-sqrt(weightsX*weightsY)
  initial.slope<-(coef(lm(Y~X))[[2]]) #uses typical least squares linear model to calculate an initial slope to start the iterations
  initial.intercept<-(coef(lm(Y~X))[[1]])	#uses typical least squares linear model to calculate an initial intercept to start the interations
  Wi<-(weightsX*weightsY/(weightsX+initial.slope^2*weightsY-2*initial.slope*error.corr*alpha))
  Ui<-X-(sum(Wi*X)/sum(Wi))
  Vi<-Y-(sum(Wi*Y)/sum(Wi))
  beta.i<-Wi*(Ui/weightsY+initial.slope*Vi/weightsX-(initial.slope*Ui+Vi)*error.corr/alpha)
  york.slope<-sum(Wi*beta.i*Vi)/sum(Wi*beta.i*Ui)
  york.intercept<-sum(Wi*Y)/sum(Wi)-york.slope*sum(Wi*X)/sum(Wi)

  york.slope.temp<-numeric(1000)    #defines the intermediate slope values, for each step of iteration. Neccesary if I want R to store these values so i can monitor the changes
  york.int.temp<-numeric(1000)  #defines the intermediate intercept values, for each step of iteration. Neccesary if I want to R to store these values so i can monitor the changes

  for (i in 1:1000) {     #defines loop for subsequent iterations, brackets open the part that contains the looped equations
    Wi<-(weightsX*weightsY/(weightsX+york.slope^2*weightsY-2*york.slope*error.corr*alpha))
    Ui<-X-(sum(Wi*X)/sum(Wi))
    Vi<-Y-(sum(Wi*Y)/sum(Wi))
    beta.i<-Wi*(Ui/weightsY+york.slope*Vi/weightsX-(york.slope*Ui+Vi)*error.corr/alpha)
    york.slope<-sum(Wi*beta.i*Vi)/sum(Wi*beta.i*Ui)   #this gives your final slope value
    york.intercept<-sum(Wi*Y)/sum(Wi)-york.slope*sum(Wi*X)/sum(Wi)  #final york intercept value
    york.slope.temp[i]<-york.slope
    york.int.temp[i]<-york.intercept
  }                        #closes loop
  xi<-sum(Wi*X)/sum(Wi)+beta.i
  ui<-xi-sum(Wi*xi)/sum(Wi)
  york.slopevar<-1/sum(Wi*ui^2)	  #variance of the slope
  york.intvar<-1/sum(Wi)+(sum(Wi*xi)/sum(Wi))^2*york.slopevar   #variance of the intercept
  N<-length(X)
  york.fit<-sum(Wi*(Y-york.slope*X-york.intercept)^2)/(N-2)      #goodness of fit.
  S<-sum(1/(x.err^2))
  Sx<-sum(X/(x.err^2))
  Sxx<-sum(X^2/(x.err^2))
  Cab<- -Sx/(S*Sxx-(Sx^2)) #measure of covariance of slope and intecept, maybe lifted from Least squares approach?
  return(lapply(setNames(ls(), ls()), function(i) get(i, -1)))
  }      #closes the function.  Everything after this runs the function or calls the results


 ########################################################
 # functions for adding regression equations to ggplots #
 # use with the geom_text(x=, y= label=lm_eqn(...))     #
 ########################################################

#' lm poly eq on plot
#' for just a second-order polynomial
#' @export
lm_poly2eqn = function(x,y){
  m = lm(y ~ poly(x,2, raw=TRUE)) #need raw=TRUE to get same coefficients as excel and same fit done by GGPLOT2
  eq <- substitute(italic(y) == c %.% italic(x)^2 + b %.% italic(x) + a * "," ~ italic(r)^2 == r2,
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),
                        c = format(coef(m)[3], digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq))
}

# lm eqn. on plots
# for just a basic linear regression
# @export
# lm_eqn = function(x,y){
#   m = lm(y ~ x)
#   eq <- substitute(italic(y) == b %.% italic(x) + a * "," ~ italic(r)^2 == r2,
#                    list(a = format(coef(m)[1], digits = 2),
#                         b = format(coef(m)[2], digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 2)))
#   as.character(as.expression(eq))
# }

#' lm eqn. on plots
#' form that deals with negatives better, normal linear regression
#' for just a basic linear regression
#' @export
lm_eqn = function(x,y) {
  m = lm(y ~ x)
  l <- list(a = format(abs(coef(m)[1]), digits = 2),
            b = format(coef(m)[2], digits = 2),
            r2 = format(summary(m)$r.squared, digits = 2));

  if (coef(m)[1] >= 0)  {
    eq <- substitute(italic(y) == b %.% italic(x) + a *","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == b %.% italic(x) - a *","~~italic(r)^2~"="~r2,l)
  }

  as.character(as.expression(eq));
}

#' plot eqns
#'
#' form that does this for york regression results
#' @export
york_eqn = function(slope, intercept){

  l<-list(m=round(slope, 4),
          b=round(abs(intercept), 4)  );

  if (intercept >= 0)  {
    eq <- substitute(italic(y) == m %.% italic(x) + b, l)
  } else {
    eq <- substitute(italic(y) == m %.% italic(x) - b, l)
  }

  as.character(as.expression(eq));
}


###########################
#' function for ascribing a color from the wheel
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


  ######################
  ## coding examples: expressions
  #how to code strings of text, useful for using "parse" eg in geom_text
  #ggplot(data.frame(x=1:3, y=1:3), aes(x,y)) + geom_point() + labs(x=expression(Delta^2 + a == x^2 + 5 ~ "hello" ~ sqrt(x) ~ "" ~ Delta), y=expression(Delta[47]^"wurst" ~ "C"))

  ######################
  #math expressions that can be used in plots to create d18O min and fluid contours, if x is temp, and y is 18O
  #e.g. fluid contour is y=0.97006*(((z+1000)*exp((((18.03*10^3)/(x+273.15))-32.42)/10^3))-1000)-29.94 for a given z (d18O fluid)
  #e.g. 1.03091*(((z+1000)*1/(exp((((18.03*10^3)/(x+273.15))-32.42)/10^3)))-1000) + 30.91 for a given z (d18O rock)

#' contours
#' @export
contours<-function(x, k, fun) {
  xN<-length(x)
  kN<-length(k)
  allX<-rep(x, times=kN)
  allK<-rep(k, each=xN)
  y<-do.call(fun, args=list(x=allX, k=allK))
  data.frame(x=allX, y=y, k=allK)
}

#example code for how to call the contour functions
#  modelDF.fluid<-contours(
#   x=seq(from=0, to=40, by=0.1),
#  k=seq(from=-9, to=3, by=1),
#    fun=function(x, k) 0.97006*(((k+1000)*exp((((18.03*10^3)/(x+273.15))-32.42)/10^3))-1000)-29.94)

#  modelDF.rock<-contours(
#    x=seq(from=20, to=60, by=0.1),
#    k=seq(from=-9, to=3, by=1),
#    fun=function(x, k) 1.03091*(((k+1000)*1/(exp((((18.03*10^3)/(x+273.15))-32.42)/10^3)))-1000) + 30.91)

#example code to plot the contours
#ggplot(modelDF, aes(x, y, group=k, label=k)) + geom_line(colour="red") + geom_text(data=subset(modelDF, x==-5), aes(y=y+0.3), colour="red")

# excel export ---------

#' add worksheet with data
#' @export
add_ws_with_data <- function(wb, sheet, data) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet=sheet, data)
  return(wb)
}

##########################################
# code for color-blind friendly palettes
#########################

.onLoad <- function(libname, pkgname) {
  # make global variables for palettes
  # The palette with grey:
  cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # The palette with black:
  cbbPalette <<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}


# To use for fills, add
#scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)





#######################################
#Function for defining legend size#
# Create a function with three arguments
# p : pre-defined plot
# legend_size : legend size, expressed as percentage of window width
# legend_text_size : font size in points

    ggplot_size_legend <- function(p, legend_size=0.1, legend_text_size=10)
    {
      Layout <- grid.layout(nrow = 1,
                            ncol = 2,
                            widths = unit(c(1-legend_size, legend_size),
                                          c("null", "null")),
                            heights = unit(1, "null"))

      vplayout <- function(...) {
        grid.newpage()
        pushViewport(viewport(layout = Layout))
      }

      subplot <- function(x, y) viewport(layout.pos.row = x,
                                         layout.pos.col = y)

      #create two plots, one with legend, one without

      pl <- p + opts(legend.position = "none")
      pn <- p + theme_grey(legend_text_size) + opts(keep = "legend_box")

      # print the plot

      vplayout()
      print(pl, vp = subplot(1, 1))
      print(pn, vp = subplot(1, 2))
    }

######################
# expand data frames #
######################

# (same as expand grid but paramters can also be whole data frames)

expand.df <- function(...) {
  # convert all params to data.frames
  l <- list(...)
  dfs <- lapply(1:length(l), function(i) { if(is.data.frame(l[[i]])) l[[i]] else as.data.frame(l[i])})

  # get indices grid
  indices <- lapply(rev(dfs), function(df) seq_len(nrow(df)))
  ind.grid <- do.call(expand.grid, indices)

  #use subsetting and cbind
  exp.dfs <- lapply(1:length(dfs), function(i) dfs[[i]][ind.grid[,length(dfs)+1-i], , drop = F])
  do.call(cbind, exp.dfs)
}

############################
# multiple legends display #
############################

#code from Seb to get ggplot legends to display correctly when using separate variables for multiple aesthetics, e.g. fill and shape are determines by different variable; note, this can often be better displayed with facet

guides_fill_shape <- function(...){
  guides(fill = guide_legend(override.aes = list(colour = "white", size=8, shape = 22), ...),
         shape = guide_legend(override.aes = list(fill = "gray"), ...))
}

                    
                    
#' STEYX function
#'
#' Translation of STEYX function Excel > R
#' Returns the standard error of the predicted y-value for each x in the regression.
#' http://office.microsoft.com/en-au/excel-help/steyx-function-HP010062545.aspx
#' @param x known x
#' @param y known y
#' @keywords steyx
#' @export
#' @examples 
#' steyx(10,150)

steyx <- function(x,y) {
# check parameters
if (missing(x))
        stop("Specify x!")
if (missing(y))
        stop("Specify y!")
if (length(x)!=length(y))
	stop("X and Y not of the same length...")

# calculations
n <- length(x)
a <- sqrt((1/(n-2))*(sum((y-mean(y))^2)-((sum((x-mean(x))*(y-mean(y))))^2)/(sum((x-mean(x))^2))))

# result
return(a)
}


