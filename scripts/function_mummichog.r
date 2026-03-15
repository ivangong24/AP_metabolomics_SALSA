library(MetaboAnalystR)
library(ggrepel)

mummichog <- function(wd_mum, input, sub_dir1, sub_dir2) {
  wd <- getwd()
  setwd(wd_mum)
  
  main_dir <- sub("(.*metaboAnalyst).*", "\\1", normalizePath(getwd()))
  input_path <- file.path(list.dirs(here::here(main_dir, "Input", sub_dir1, sub_dir2)), 
                          input)
  # input_path <- paste("../Input/", input, sep = "")
  
  mSet<-InitDataObjects("mass_all", "mummichog", FALSE, default.dpi = 300)
  mSet<-SetPeakFormat(mSet, "rmp")
  mSet<-UpdateInstrumentParameters(mSet, 10.0, "mixed", "yes", 0.02);
  mSet<-Read.PeakListData(mSet, input_path);
  mSet<-SanityCheckMummichogData(mSet)
  add.vec <- c("M-H [1-]","M-2H [2-]","M-H2O-H [1-]","M [1+]","M+H [1+]","M+Na [1+]")
  mSet<-Setup.AdductData(mSet, add.vec);
  mSet<-PerformAdductMapping(mSet, "mixed")
  mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
  # pval <- sort(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])[ceiling(length(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])*0.1)]
  mSet<-SetMummichogPval(mSet, 0.1)
  mSet<-PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)
  mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 300, width = 10)
  
  setwd(wd)

  # change plot
  library(ggrepel)
  mat <- mSet$mummi.resmat
  mat <- mat %>%
    as.data.frame()
  #  filter(Hits.sig >= 3)

  y <- -log10(mat[, 5])
  x <- mat[, 8]/mat[, 4]
  pathnames <- rownames(mat)
  Hits.sig <- mat$Hits.sig

  inx <- order(y, decreasing = TRUE)
  y <- y[inx]
  x <- x[inx]
  pathnames <- pathnames[inx]
  Hits.sig <- Hits.sig[inx]
  df_plot <- data.frame(y, x, pathnames, Hits.sig)
  radi.vec <- sqrt(abs(x))

  p <- df_plot %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(size = radi.vec,
                   color = y),
               stroke = 0.5) +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_gradient(low = "yellow", high = "red",
                         name = "-log10(p)") +
    xlab("Enrichment Factor") +
    ylab("-log10(p)") +
    theme_minimal() +
    theme(legend.position = "none")

  top_indices <- (df_plot$y > -log10(0.05)) & (df_plot$Hits.sig >= 3)

  p <- p +
    geom_text_repel(aes(label = pathnames),
                    data = df_plot[top_indices, ], size = 3)
  return(p)
}

