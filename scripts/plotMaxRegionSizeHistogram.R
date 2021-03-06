#!/usr/bin/Rscript

# Plot histogram of max region size
# from Axelrod experiments, for critical value of q at 
# set of fixed radius values
# this verision usies the output from collect_results.sh collecting results
# from the model e.g. /lattice-jointactivity-simcoop-social-noise-cpp-end/model
# run by
# lattice-python-mpi/src/axelrod/geo/expphysicstimeline/multiruninitmain.py
# but with 'end' on the command line so writes stats at start and end only.
#
# Fixed value of F and L and zero noise only. Limited set of
# radius values coded here.
#
# The name of the CSV file with the data is given on stdin
#
# Usage:
#
# Rscript plotMaxRegionSizeHistogram.R data.csv outputprefix
#
#
# E.g. Rscript plotMaxRegionSizeHistogram.R results.csv end_radius
#
# Output is output.eps files
# given on command line, with -max_region_size and -num_regions suffix  etc.
# e.g. outputprefix-max_region_size and outputprefix-num_regions.eps
#
# ADS April 2016
#
# $Id: plotMaxRegionSizeHistogram.R 1043 2016-10-13 04:07:51Z stivalaa $


library(ggplot2)
library(doBy)
library(reshape)
library(scales)
library(grid)
library(gridExtra)


bla <- function(variable, value) {

    # note use of bquote(), like a LISP backquote, evaluates only inside
    # .()
# http://stackoverflow.com/questions/4302367/concatenate-strings-and-expressions-in-a-plots-title
# but actually it didn't work, get wrong value for beta_s (2, even though no
   # such beta_s exists ?!? some problem with efaluation frame [couldn't get
# it to work by chaning where= arugment either), so using
# substitute () instead , which also doesn't work, gave up.
# --- turns out problem was I'd forgotten that all these vsariables have
# been converted to factors, so have to do .(levels(x)[x]) not just .(x)
    sapply      (value, FUN=function(xval ) 
        if (variable == "beta_s") {
          print(xval)
          bquote(beta[s]==.(levels(xval)[xval]))
        }
        else if (variable == "n") {
          bquote(N/L^2==.(levels(xval)[xval]))
        }
        else if (variable == "n_immutable") {
          bquote(F[I] == .(levels(xval)[xval]))
        }
        else if (variable == "m") {
          bquote(L == .(levels(xval)[xval]))
        }
        else if (variable == "pool_multiplier") {
          bquote(plain("pool multiplier") ~~ m == .(levels(xval)[xval]))
        }
        else if (variable == "mpcr") {
          bquote(plain("MPCR") == .(xval))
        }
        else {
          bquote(.(variable) == .(levels(xval)[xval]))
        }
      )
}

# http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

my_scientific_10 <- function(x) {
  # also remove "1 x", just have the exponent

  parse(text=gsub("1e", " 10^", scientific(as.numeric(x))))
}



Fvalue <- 5   # uses only one value of F
Lvalue <- 100 # uses only one value of L

radii <- c(1, 2, 3, 4)

responses <- c('max_region_size', 'num_regions')
response_names <- c(expression(symbol("\074")~S[max]~symbol("\076")/L^2), 'Number of regions' )


if (length(commandArgs(trailingOnly=TRUE)) != 2) {
    cat("Usage: Rscript plotMaxRegionSizeHistogram.R results.csv outputprefix\n")
    quit(save='no')
}
results_filename <- commandArgs(trailingOnly=TRUE)[1]
output_prefix <- commandArgs(trailingOnly=TRUE)[2]

experiment <- read.table(results_filename, sep=',',header=T,stringsAsFactors=F)

experiment <- experiment[which(experiment$time > 0),]
experiment <- experiment[which(experiment$noise == 0),]
experiment <- experiment[which(experiment$F == Fvalue),]
experiment <- experiment[which(experiment$m == Lvalue),]
experiment <- experiment[which(experiment$radius  %in% radii),]

experiment$m <- as.factor(experiment$m)

  
D <- melt(experiment, id=c('n','m','F','q','radius','pool_multiplier','noise_rate','run') )
E0 <- D
D<-summaryBy(value ~ n + m + F + q + radius + pool_multiplier + noise_rate +  variable, data=D, FUN=c(mean, sd))
  
for (i in 1:length(responses)) {

  plotlist <- list()

  response <- responses[i]
  if (!(response %in% colnames(experiment))) {
      print(paste('skipping response ', response, ' not in data'))
      next
  }

  response_name <- response_names[i]
  Dst <- D[which(D$variable == response), ]
  E <- E0[which(E0$variable == response), ]
  if (all(is.nan(Dst$value.mean))) {
      print(paste('skipping response ', response, ' all values are NaN'))
      next
  }

  # get critical value of q where variance of order parmeter
  # (max region size or number of regions) is largest
  if (response == 'max_region_size' || response == 'num_regions') {
      for (m in unique(Dst$m)) {
          for (radius in unique(Dst$radius)) {
              Dst0 <- Dst[which(Dst$noise_rate == 0 & Dst$m == m & Dst$radius == radius),]
              q_c <- Dst0[which(Dst0$value.sd >= max(Dst0$value.sd)), 'q']

              if (length(q_c) < 1) {
                  q_c <- NA
              }
              cat(response, ' m = ', m, ' radius = ', radius, ' q_c = ', q_c, '\n')              
              Dst[which(Dst$m == m & Dst$radius == radius),'q_c'] <- q_c
          }
      }
  }

  E$q_c <- NA
  for (radius in unique(E$radius)) {
    q_c <- Dst[which(Dst$radius == radius), 'q_c']
    stopifnot(all(q_c[1] == q_c)) 
    q_c <- q_c[1]
    E[which(E$radius == radius),'q_c'] <- q_c
  }
  Dst <- E


  Dst0 <- Dst
  for (radius in radii) {
      print(radius)#XXX
      Dst <- Dst0[which(Dst0$radius == radius),]
      p <- ggplot(Dst, aes(x = value))
      p <- p + theme_bw()
      p <- p + theme(plot.margin = 	    unit(c(0,0,0,0), "lines"),
                    axis.text.x =       element_text(size = 12, colour = "black", lineheight = 0.2),
                    axis.text.y =       element_text(size = 12, colour = "black", lineheight = 0.2),
                    axis.title.x =       element_text(size = 12, colour = "black", lineheight = 0.2),
                    axis.title.y =       element_text(angle = 90, size = 12, colour = "black", lineheight = 0.2),

                  strip.text.x =	element_text(size = 12, colour = "black"),
                  strip.text.y =	element_text(angle = -90, size = 10, colour = "black"),
            legend.text = 	element_text(size = 12, colour = "black"),
                  legend.key =  	element_rect(fill = "white", colour = "white"),

                    axis.ticks =        element_line(colour = "black"),
                    axis.ticks.length = unit(0.1, "cm"),
                    strip.background =  element_rect(fill = "white", colour = "white"),
                    panel.grid.minor =  element_blank(),
                    panel.grid.major =  element_blank(),
                    panel.border =      element_rect(colour = "black"),
                    axis.ticks.margin = unit(0.1, "cm")
      )
      # http://stackoverflow.com/questions/11766856/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion
      p <- p + geom_histogram(binwidth=diff(range(Dst$value))/10,
                              aes(y=..count../sum(..count..)),
                              colour="black", fill="white")
      p <- p + xlim(c(0, 1))
      p <- p + xlab(response_name)
      p <- p + ylab("Relative frequency")

    #  p <- p + ggtitle(bquote(list(F == .(Fvalue), L == .(Lvalue) )))
      p <- p + ggtitle(bquote(list(R == .(radius), q[c] == .(Dst$q_c) )))

      plotlist <- c(plotlist, list(p))
    }

  # EPS suitable for inserting into LaTeX
  postscript(paste(paste(output_prefix, response, sep='-'), 'eps',sep='.'),
             onefile=FALSE,paper="special",horizontal=FALSE, 
             width = 9, height = 6)
  do.call(grid.arrange, plotlist)
  dev.off()
  
}



