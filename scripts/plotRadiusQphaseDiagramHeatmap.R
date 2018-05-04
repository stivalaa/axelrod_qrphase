#!/usr/bin/Rscript

# Plot heat map of q and radius with heat color of max region size as 
# q-radius phase diagram using results from  from Axelrod experiments.
# this verision usies the output from collect_results.sh collecting results
# from the model e.g. /lattice-jointactivity-simcoop-social-noise-cpp-end/model
# run by
# lattice-python-mpi/src/axelrod/geo/expphysicstimeline/multiruninitmain.py
# but with 'end' on the command line so writes stats at start and end only.
#
# Fxied value of F and zero noise only. Facet on  m (lattice size).
#
# The name of the CSV file with the data is given on stdin
#
# Usage:
#
# Rscript plotRadiusQphaseDiagramHeatmap.R data.csv outputfilenameprefix
#
#
# E.g. Rscript plotRadiusQphaseDiagramHeatmap.R results.csv phase_diagram
#
# Output is output.eps files with -m<latticesize> appended e.g.
# phase_daigram-m50.eps
#
# ADS August 2016
#
# $Id: plotCriticalRadiusVsQend.R 892 2016-08-11 00:02:05Z stivalaa $


library(doBy)
library(reshape2)
#library(gplots)
library(ggplot2)
library(grid)


Fvalue <- 5   # uses only one value of F


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


if (length(commandArgs(trailingOnly=TRUE)) != 2) {
    cat("Usage: Rscript plotRadiusQphaseDiagramHeatmap.R results.csv outputprefix\n")
    quit(save='no')
}
results_filename <- commandArgs(trailingOnly=TRUE)[1]
base_output_prefix <- commandArgs(trailingOnly=TRUE)[2]

orig_experiment <- read.table(results_filename, sep=',',header=T,stringsAsFactors=F)

orig_experiment <- orig_experiment[which(orig_experiment$time > 0),]
orig_experiment <- orig_experiment[which(orig_experiment$noise == 0),]
orig_experiment <- orig_experiment[which(orig_experiment$F == Fvalue),]


experiment <- orig_experiment

stopifnot(length(unique(experiment$pool_multiplier)) <= 1)
stopifnot(length(unique(experiment$num_joint_activities)) <= 1)

output_prefix <- base_output_prefix
 
experiment$m <- factor(experiment$m)
 
D <- melt(experiment, id=c('n','m','F','q','radius','pool_multiplier','noise_rate','run') )
D<-summaryBy(value ~ n + m + F + q + radius + pool_multiplier + noise_rate + variable, data=D, FUN=c(mean, sd))
  
  
Dst <- D[which(D$variable == 'max_region_size'), ]

#  X <- acast(Dst, q ~ radius, value.var="value.mean"  ,fun.aggregate=mean)
#  heatmap.2(X , dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none',
#            xlab='radius', ylab='q', col=heat.colors, revC=TRUE)

p <- ggplot(Dst, aes(radius, q))
p <- p + geom_tile(aes(fill = value.mean))
p <- p + xlab(expression(R))
p <- p + ylab(expression(q))
#p <- p + theme(aspect.raio=1) # doen't work ('not a valid theme element name')
p <- p + coord_fixed(ratio = max(Dst$radius)/max(Dst$q))
#p <- p + scale_fill_gradient(expression(group("<", S[max], ">")/L^2)) # does not work
#p <- p + scale_fill_gradient(expression(group("[", S[max], "]")/L^2))  # works but not what we need
p <- p + scale_fill_distiller(expression(symbol("\074")~S[max]~symbol("\076")/L^2))
p <- p + facet_grid(. ~ m, labeller=bla)

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

# EPS suitable for inserting into LaTeX
postscript(paste(output_prefix, 'eps', sep='.'),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
print( p )
dev.off()

