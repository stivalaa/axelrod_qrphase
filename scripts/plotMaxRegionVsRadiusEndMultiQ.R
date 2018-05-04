#!/usr/bin/Rscript

# Plot radius on x axis and max cluster size on y axis for results
# from Axelrod / joint activity experiments. 
# this verision usies the output from collect_results.sh collecting results
# from the model e.g. /lattice-jointactivity-simcoop-social-noise-cpp-end/model
# run by
# lattice-python-mpi/src/axelrod/geo/expphysicstimeline/multiruninitmain.py
# but with 'end' on the command line so writes stats at start and end only.
# Also create spearate file with num clusters on y axis.
# Also plots value of q where variance of max region size and
# number regions is maximum.
#
# Fxied value of F and zero noise only. Limited set of q values coded here.
#
# The name of the CSV file with the data is given on stdin
#
# Usage:
#
# Rscript plotMaxRegionVsRadiusEndMultiQ.R data.csv outputprefix
#
#
# E.g. Rscript plotMaxRegionVsRadiusEndMultiQ.R results.csv end_radius
#
# Output is output.eps files
# given on command line, with -max_region_size and -num_regions suffix  etc.
# e.g. outputprefix-max_region_size and outputprefix-num_regions.eps
#
# ADS April 2016
#
# $Id: plotMaxRegionVsRadiusEnd.R 903 2016-08-12 00:17:52Z stivalaa $


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



zSigma <- 1.96 # number of sd for 95% CI

Fvalue <- 5   # uses only one value of F

qs <- c(15, 30, 45, 60, 75, 90)
max_radius <- 7 # spread graph out more, all monoculture after this anyway

responses <- c('max_region_size', 'num_regions',  'num_cultures')
response_names <- c(expression(symbol("\074")~S[max]~symbol("\076")/L^2), 'Number of regions', 'Number of cultures')


if (length(commandArgs(trailingOnly=TRUE)) != 2) {
    cat("Usage: Rscript plotMaxRegionVsRadiusEndMultiQ.R results.csv outputprefix\n")
    quit(save='no')
}
results_filename <- commandArgs(trailingOnly=TRUE)[1]
output_prefix <- commandArgs(trailingOnly=TRUE)[2]

experiment <- read.table(results_filename, sep=',',header=T,stringsAsFactors=F)

experiment <- experiment[which(experiment$time > 0),]
experiment <- experiment[which(experiment$noise == 0),]
experiment <- experiment[which(experiment$F == Fvalue),]

experiment <- experiment[which(experiment$q %in% qs),]
experiment <- experiment[which(experiment$radius  <= max_radius),]

  experiment$m <- as.factor(experiment$m)
  experiment$q <- as.factor(experiment$q)
  
  D <- melt(experiment, id=c('n','m','F','q','radius','pool_multiplier','noise_rate','run') )
  D<-summaryBy(value ~ n + m + F + q + radius + pool_multiplier + noise_rate + variable, data=D, FUN=c(mean, sd))

  
for (i in 1:length(responses)) {
  response <- responses[i]
  if (!(response %in% colnames(experiment))) {
      print(paste('skipping response ', response, ' not in data'))
      next
  }

  response_name <- response_names[i]
  Dst <- D[which(D$variable == response), ]
  if (all(is.nan(Dst$value.mean))) {
      print(paste('skipping response ', response, ' all values are NaN'))
      next
  }

  # get critical value of radius where variance of order parmeter
  # (max region size or number of regions) is largest
  if (response == 'max_region_size' || response == 'num_regions') {
      stopifnot(D$noise == 0)
      for (m in unique(Dst$m)) {
          for (q in qs) {
              Dst0 <- Dst[which(Dst$m == m & Dst$q == q),]
              critical_radius <- Dst0[which(Dst0$value.sd >= max(Dst0$value.sd)), 'radius']

              if (length(critical_radius) < 1) {
                  critical_radius <- NA
              }
              cat(response, ' m = ', m, ' q = ', q, ' critical_radius = ', critical_radius, '\n')              
              Dst[which(Dst$m == m & Dst$q == q),'critical_radius'] <- critical_radius
          }
      }
  }

#print(Dst)#XXX
      
#  print('before ggplot'); print(lsos()) #XXX

  p <- ggplot(Dst, aes(x = radius, y = value.mean,
                       colour = as.factor(q), linetype=as.factor(q), shape = as.factor(q))   )
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
  if (response == "num_cooperators") {
    p <- p + theme(legend.position = c(0, 1), legend.justification = c(0, 1)) # top left (for num_cooperators plot)
  }
  # else leave legnd outside plot region

  p <- p + geom_point()
  p <- p + geom_errorbar(aes(ymin=value.mean - zSigma*value.sd, ymax=value.mean + zSigma*value.sd) )
  p <- p + geom_line()
  # plot vertical line at critical value of radius
  #p <- p + geom_vline(aes(xintercept = as.numeric(critical_radius)), colour = "gray55", linetype = "dashed")
#  p <- p + geom_vline(aes(xintercept = as.numeric(critical_radius), colour = as.factor(q), linetype=as.factor(q)) )

  #p <- p + xlab('von Neumann radius R')
  p <- p + xlab(expression(R))
  p <- p + ylab(response_name)
  p <- p + facet_grid(m ~ ., labeller=bla)
  p <- p + ggtitle (bquote(list(F == .(Fvalue))))

  p <- p + scale_colour_brewer('q', palette = "Dark2")
  p <- p + scale_linetype('q')
  p <- p + scale_shape_manual('q', values=1:nlevels(as.factor(Dst$q))) # needed for more than 6 shapes


# EPS suitable for inserting into LaTeX
  postscript(paste(paste(output_prefix, response, sep='-'), 'eps',sep='.'),
             onefile=FALSE,paper="special",horizontal=FALSE, 
             width = 9, height = 6)

  print( p)
  dev.off()
}

