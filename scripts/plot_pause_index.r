library(PEPPROr)
args<-commandArgs(trailingOnly = T)
labelCuts = function(breakPoints, digits=1, collapse="-", infBins=FALSE) {
  labels <- 
    apply(round(cbind(breakPoints[-length(breakPoints)],   
                      breakPoints[-1]),digits), 1, paste0,
          collapse=collapse) 
  
  if (infBins) {
    labels[1] <- paste0("<", breakPoints[2])
    labels[length(labels)] <- paste0(">", breakPoints[length(breakPoints)-1])
  }
  return(labels)
}
cutDists = function(vec, divisions = c(-Inf, -1e6, -1e4, -1000, -100, 0,
                                       100, 1000, 10000, 1e6, Inf)) {
  if (is.list(vec)) {
    x = lapply(vec, cutDists)
    
    # To accommodate multiple lists, we'll need to introduce a new 'name'
    # column to distinguish them.
    nameList = names(vec)
    if(is.null(nameList)) {
      nameList = 1:length(query) # Fallback to sequential numbers
    }
    
    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, sapply(x, nrow))
    
    return(xb)
  }
  divisions <- unique(divisions)
  labels <- labelCuts(signif(divisions, 3), collapse=" to ", infBins=TRUE)
  #message(paste0("breaks: ", paste0(divisions, collapse=" ")))
  cuts   <- cut(vec, divisions, labels)
  return(as.data.frame(table(cuts)))
}
theme_PEPPRO <- function(base_family = "sans", ...){
  theme_classic(base_family = base_family, base_size = 14, ...) +
    theme(
      axis.line = element_line(size = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      aspect.ratio = 1,
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
}
plotPI <- function(pi, name='pause indicies',
                   raw=FALSE,
                   type=c("histogram", "boxplot", "violin"),
                   annotate=TRUE) {
    # TODO: make summary plot of these that IS boxplots
    if (exists(pi)) {
      PI <- data.table(get(pi))
    } else if (file.exists(pi)) {
      PI <- fread(pi)
    } else {
      stop(paste0("FileExistsError: ", pi, " could not be found."))
      quit(save = "no", status = 1, runLast = FALSE)
    }
    # colnames(PI) <- c("chr", "start", "end", "name", "pi", "strand")
    colnames(PI) <- c('name','Gene','ppc','ppm','ppd','pps','gbc','gbm','gbd','pi')


    div <- c(-Inf, seq(from=-1, to=3, by=0.25), Inf)

    # ensure breaks are not duplicated
    div <- unique(div)
    lowerLabel <- paste0(round(
        (nrow(PI[log10(PI$pi) < -1, ]) / nrow(PI)) * 100, 2), '%')
    upperLabel <- paste0(round(
        (nrow(PI[log10(PI$pi) > 3, ]) / nrow(PI)) * 100, 2), '%')


    
    if (length(div) <= 3) {
        base_plot <- ggplot(data = PI, aes(x=log10(pi)))
    } else {
        # calculate a frequency table with the specified divisions
        pi_table  <- cutDists(log10(PI$pi), divisions = div)
        base_plot <- ggplot(data = pi_table,  aes(x=cuts, y=Freq))
    }
    if (length(div) <= 3) {
        q = base_plot +
            geom_histogram(col="black", fill=I("transparent")) +
            geom_vline(aes(xintercept=median(log10(PI$pi))),
                       color="gray", linetype="dashed", size=1) +
            geom_vline(aes(xintercept=mean(log10(PI$pi))),
                       color="light gray", linetype="dotted", size=1) +
            labs(x="Pause indicies", y="Frequency") +
            scale_x_log10(limits = c(0.001, 50),
                          expand = expand_scale(mult = c(0, 0)),
                          labels=fancyNumbers,
                          breaks=prettyLogs) +
            annotation_logticks(sides = c("rl"))
    } else {
        q <- base_plot +
            geom_bar(stat="identity",
                     fill = c("gray80",
                              rep("gray40", (length(div)-3)),
                              "gray80")) + 
            labs(x="Pause indicies", y="Frequency") +
            geom_text(aes(label=c(lowerLabel, upperLabel)),
                      size=theme_get()$text[["size"]]/4,
                      data=pi_table[c(1,length(pi_table$Freq)),],
                      vjust=0.5, hjust=-0.1, angle=90)
    }
    max_x <- length(layer_scales(q)$x$range$range)
    max_y <- suppressMessages(
      suppressWarnings(layer_scales(q)$y$range$range[2]))

    if (is.na(max_x)) {max_x <- Inf}


    label1 <- paste("'median'[raw]", ":", round(median(PI$pi), 2))
    label2 <- paste("'mean'[raw]", ":", round(mean(PI$pi), 2))
    label3 <- paste("'median'[log[10]]", ":",
                    round(median(log10((PI$pi))), 2))
    label4 <- paste("'mean'[log[10]]", ":", round(mean(log10(PI$pi)), 2))


    q <- q + annotate("text", x = floor(max_x), y = floor(max_y),
                      size=theme_get()$text[["size"]]/4,
                      hjust="right", vjust=1.05,
                      label = label1, parse=TRUE) +
        annotate("text", x = floor(max_x), y = floor(max_y),
                 size=theme_get()$text[["size"]]/4,
                 hjust="right", vjust=2.25,
                 label = label2, parse=TRUE) +
        annotate("text", x = floor(max_x), y = floor(max_y),
                size=theme_get()$text[["size"]]/4,
                 hjust="right", vjust=3.25,
                 label = label3, parse=TRUE) +
        annotate("text", x = floor(max_x), y = floor(max_y),
                 size=theme_get()$text[["size"]]/4,
                 hjust="right", vjust=4.30,
                 label = label4, parse=TRUE) +
        theme_PEPPRO()

    return(q)
}
sample_name        <- sampleName(args[1])
name               <- basename(sample_name)
suppressWarnings(p <- plotPI(pi=args[1], name=name, raw=FALSE,
                             type=tolower('histogram'),
                             annotate=TRUE))
# Save plot to pdf file
pdf(file=paste0(sample_name, "_pause_index.pdf"),
    width = 4, height = 4, useDingbats=F)
suppressWarnings(print(p))
invisible(dev.off())

# Save plot to png file
png(filename = paste0(sample_name, "_pause_index.png"),
    width = 275, height = 275)
suppressWarnings(print(p))
invisible(dev.off())

if (exists("p")) {
  write("Pause index plot completed!\n", stdout())
} else {
  write("Unable to produce pause index plot!\n", stdout())
}