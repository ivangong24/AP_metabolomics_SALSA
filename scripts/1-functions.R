## ---------------------------
##
## Script name: 1-functions.R
## Purpose of script: To create functions for analysis in the project
##
## Author: Yufan Gong
##
## Date Created: 2025-11-12
##
## Date Modified: 2025-11-12
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------

# 1. Check the directory --------------------------------------------------

#please open all the R scripts in AP_Metabolomics_SALSA.Rproj

{
  getwd()
  options(mc.cores = 9)
}


# 2. Loading packages -----------------------------------------------------

{
  pacman::p_load(
    #For creating tables
    "kableExtra",     #create amazing tables: kbl()
    "skimr",          #summary statistics: skim()
    #"dataxray",      #explore data: report_xray()
    "arsenal",        #create tables: tableby()
    "expss",          #create contingency tables: calc_cro_cpct()  
    "huxtable",       #as_hux_table()
    "flextable",      #as_flex_table()
    "gtsummary",      #create amazing tables: tbl_summary()
    
    #For loading and exporting data
    "readxl",         #read in excel data: read_xlsx
    "writexl",        #write in excel data: write_xlsx()
    "haven",          #read in sas data: read_sas()
    "here",           #setting the directory in the project: here()
    
    #For manipulating data
    "rlang",          #for Non-standard evaluation: eval(), expr(), ensym(), caller_env(), exec(), !!
    "magrittr",       #for the pipe operator: %>% and %<>%
    "lubridate",      #for manipulating dates: intervals(), durations()
    "labelled",       #labelleling the data: set_variable_labels(), set_value_labels()
    
    # Enhancing plots
    "scales",         #makes easy to format percent, dollars, comas: percent()
    "ggalt",          #makes easy splines: geom_xsplines()
    "ggeasy",         #applies labels among other things: easy_labs()
    "gridExtra",      #combining plots and tables on plots: grid.arrange(), tableGrob()
    "ggpubr",         #combines plots: ggarrange()
    "ggthemes",       #blind colors
    "ggVennDiagram",  #venn diagram
    "Amelia",         #check missing pattern: missmap()
    "pheatmap",       #create pretty heatmap: pheatmap(),
    "ggnewscale",     #create multiple scales: new_scale()
    "ggrepel",        #avoid overlapping labels: geom_text_repel()
    "DataExplorer",   #create a report of the data: create_report()
    "patchwork",      #combines plots: wrap_plots(), plot_layout()
    "plotly",         #create interactive plots: ggplotly()
    
    # Other great packages
    "glue",           #replaces paste: glue()
    "Hmisc",          #explore the data: describe()
    "mice",           #imput missing data: mice()
    "gmodels",        #create contigency table: CrossTable()
    "meta",           #meta models: metagen()
    "codebook",       #amazing package to set labels: dict_to_list()
    "foreach",        #executing R code repeatedly: foreach(), %do%, %dopar%
    "doParallel",     #Provides a parallel backend for the %dopar% function: registerDoParallel(),
    "doFuture",       #Provides a parallel backend for the %dofuture% function: registerDoFuture()
    "future",         #Provides a future backend for the %dofuture% function: plan()
    "furrr",          #Provides a purrr like syntax for future: future_map(), future_pmap()
    
    # For analysis
    "lme4",           #linear mixed-effects model: lmer()
    "limma",          #linear models for microarray data: lmFit() 
    "nlme",           #linear and nonlinear mixed effects:
    "survival",       #survival analysis: Surv(), coxph()
    "WGCNA",          #weighted correlation network analysis: GOenrichmentAnalysis()
    "mixOmics",       #block.plsda(),
    "MetaboAnalystR", # pathway analysis: PerformPSEA()
    
    #For data cleaning
    "tidyverse"       #data manipulation and visualization:select(), mutate()
  )
}

# 3. Create sub folders ---------------------------------------------------

# c("scripts", "qmd", "tables", "figures", "data") %>%
#   map(dir.create)

# 4. clean env ------------------------------------------------------------

rm(list = ls())

# 5. create functions -----------------------------------------------------

## negate function for `%in%`
`%notin%` <- Negate(`%in%`)

## convert plain text in a list to text within quotes
quote_all <- function(...) {
  args <- rlang::enquos(...)  # capture all inputs as quosures
  # applies transformation and returns character vector
  purrr::map_chr(args, ~{
    expr <- rlang::get_expr(.x) # extracts the expression from each quosure
    if (rlang::is_symbol(expr)) {
      rlang::as_string(expr) 
    } else {
      as.character(expr)
    }
  }) |> 
    stringr::str_c(sep = "")
}

# create table1
table1 <- function(table) {
  
  table |> 
    gtsummary::as_hux_table() -> hux
  
  table |> 
    gtsummary::as_flex_table() -> flex
  
  return(list(hux=hux, flex=flex))
  
}

# read in files in different formats
read_file <- function(file) {

  # get lowercase file extension
  ext <- tolower(tools::file_ext(file))

  switch(ext,
         
    # CSV
    csv   = readr::read_csv(file),
    
    # Excel
    xlsx  = readxl::read_xlsx(file),
    
    # TXT (tab or general delimited)
    txt   = readr::read_delim(file),
    
    # SAS
    sas7bdat = haven::read_sas(file),

    # Stata
    dta = haven::read_dta(file),

    # R serialized objects
    rds = readr::read_rds(file),
    
    # RData (assume it contains exactly one data.frame/tibble)
    rdata = {
      env  <- new.env()
      load(file, envir = env)
      objs <- as.list(env)

      objs_df <- objs |> 
        map(as.data.frame)

      if (length(objs_df) == 1 && inherits(objs_df[[1]], c("data.frame", "tbl_df", "matrix"))) {
        objs_df[[1]]
      } else {
        stop(
          "RData file `", file, "` does not contain exactly one data.frame/tibble/matrix.\n",
          "It contains: ", paste(names(objs_df), collapse = ", "), ".\n",
          "Please load and handle this file manually."
        )
      }
    },

    # default if no match
    stop("Unsupported file format: ", ext)
  )
}

# Merge all files within a single subdirectory
sdir_merge <- function(sdir) {
  files <- list.files(
    sdir,
    pattern = "\\.(csv|xlsx?|txt|sas7bdat)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    message("No data files found in: ", sdir)
    return(tibble())
  }

  # read, then bind
  files |> 
    purrr::map(read_file) |> 
    purrr::map(~rename_all(.x, str_to_lower)) |>    # standardize column names to lower case
    purrr::list_rbind() |>                          # fast row-bind (dplyr >= 1.1.0)
    dplyr::mutate(.source_dir = basename(sdir))      # optional: keep source info
}

# winsorize function
extreme_remove_percentile_win <- function(x) {
  Q = quantile(x, c(0.01, 0.99), na.rm = TRUE)
  x = dplyr::case_when(
    # x < Q[1] ~ Q[1],
    x > Q[2] ~ Q[2],
    TRUE ~ x)
}

create_dir <- function(dir){
  if(!dir.exists(dir)){
    dir.create(dir)
  }
}

#--------------------------------End of the code--------------------------------