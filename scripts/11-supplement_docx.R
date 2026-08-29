## ---------------------------
##
## Script name: 11-supplement_docx.R
## Purpose of script: Build Supplement/Supporting Information.docx by
##                    appending the supplement tables from
##                    Supplement/Supplement tables_composite_meta.xlsx to a
##                    copy of Supplement/Supporting Information_old.docx.
##
##                    Sheet titles become paragraphs in Word's "caption"
##                    style so a List of Tables can be inserted from
##                    References > Insert Table of Figures (Caption label =
##                    Table) without further setup.
##
##                    All sheet titles, ordering, and table content come
##                    from the xlsx itself; nothing is hard-coded here.
##
## Author: Yufan Gong
##
## Date Created: 2026-06-10
##
## ---------------------------

pacman::p_load(officer, flextable, openxlsx, here, dplyr, purrr, stringr)

OLD_DOCX <- here::here("Supplement", "Supporting Information_old.docx")
XLSX     <- here::here("Supplement", "Supplement tables_composite_meta.xlsx")
OUT_DOCX <- here::here("Supplement", "Supporting Information.docx")

# Maximum rows rendered per table in the docx. Word docs with tens of
# thousands of rows trigger an xml2 segfault during officer's print step
# and become impractical to open. Any sheet longer than this cap is
# truncated and gets a small italic footer pointing to the full Excel
# workbook. Set to Inf to disable; the xlsx workbook always holds the
# unabridged data.
TABLE_ROW_LIMIT <- 200

# 1. Inject a "caption" paragraph style into the source docx ---------------
# The original docx has no Caption style, so officer's body_add_par(style =
# "caption") would silently drop the style assignment and Word's TOC tools
# wouldn't be able to find the captions. We patch styles.xml inside the zip
# to add a Caption style based on Normal, formatted as Times New Roman 11
# bold, and then load that patched copy with officer.

inject_caption_style <- function(src_docx, dest_docx) {
  tmp <- tempfile("docx_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  utils::unzip(src_docx, exdir = tmp)
  styles_path <- file.path(tmp, "word", "styles.xml")
  txt <- paste(readLines(styles_path, warn = FALSE, encoding = "UTF-8"),
               collapse = "")

  if (!grepl('w:styleId="Caption"', txt, fixed = TRUE)) {
    caption_xml <- paste0(
      '<w:style w:type="paragraph" w:styleId="Caption">',
        '<w:name w:val="caption"/>',
        '<w:basedOn w:val="Normal"/>',
        '<w:next w:val="Normal"/>',
        '<w:qFormat/>',
        '<w:pPr><w:spacing w:before="120" w:after="120"/></w:pPr>',
        '<w:rPr>',
          '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/>',
          '<w:b/><w:sz w:val="22"/>',
        '</w:rPr>',
      '</w:style>'
    )
    txt <- sub("</w:styles>",
               paste0(caption_xml, "</w:styles>"),
               txt, fixed = TRUE)
    writeLines(txt, styles_path, useBytes = TRUE)
  }

  # Re-zip preserving the docx ZIP layout.
  prev_wd <- setwd(tmp); on.exit(setwd(prev_wd), add = TRUE)
  files <- list.files(".", recursive = TRUE, all.files = TRUE)
  if (file.exists(dest_docx)) file.remove(dest_docx)
  utils::zip(dest_docx, files, flags = "-Xq")
}

base_docx <- tempfile(fileext = ".docx")
inject_caption_style(OLD_DOCX, base_docx)

doc <- officer::read_docx(base_docx)

# 2. Resolve the actual heading + caption style names ---------------------
style_info <- officer::styles_info(doc, type = "paragraph")
pick_style <- function(...) {
  candidates <- unlist(list(...))
  hits <- candidates[candidates %in% style_info$style_name]
  if (length(hits) == 0)
    stop("No matching style found: ", paste(candidates, collapse = ", "))
  hits[1]
}
heading_style <- pick_style("heading 1", "Heading 1")
caption_style <- pick_style("caption", "Caption", "Table Caption")

# 3. Helpers --------------------------------------------------------------
# Mirror the numeric / scientific formatting from 10-supplement_tables.R so
# the docx tables render the same precision as the xlsx.
number_cols <- c("logFC", "VIP",
                 "FE estimate", "I^2 (%)",
                 "pathway size", "# of peaks", "# of sig peaks",
                 "empirical")
sci_cols    <- c("p", "FDR",
                 "FE p", "FE FDR", "Cochran Q p",
                 "adj.p")

fmt_col <- function(v, cn) {
  if (cn %in% number_cols && is.numeric(v))
    return(ifelse(is.na(v), "", formatC(v, format = "f", digits = 2)))
  if (cn %in% sci_cols && is.numeric(v))
    return(ifelse(is.na(v), "", formatC(v, format = "e", digits = 2)))
  v <- as.character(v)
  v[is.na(v)] <- ""
  v
}

read_sheet <- function(s) {
  title <- as.character(
    openxlsx::readWorkbook(wb, sheet = s, rows = 1, colNames = FALSE)[1, 1])
  data  <- openxlsx::readWorkbook(wb, sheet = s, startRow = 2,
                                  colNames = TRUE, check.names = FALSE)
  list(title = title, data = data)
}

# Build a flextable that approximates the xlsx visual: TNR 10 (one size
# down so wide tables fit on the page), bold header, horizontal rules on
# the header (top + bottom) and the last body row.
build_flextable <- function(d) {
  if (nrow(d) == 0)
    d <- tibble::as_tibble(setNames(
           as.list(rep("", length(colnames(d)))), colnames(d))) |>
         dplyr::slice(0)

  for (cn in colnames(d)) d[[cn]] <- fmt_col(d[[cn]], cn)

  border <- officer::fp_border(width = 1, color = "black")

  ft <- flextable::flextable(d) |>
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "left", part = "all") |>
    flextable::border_remove() |>
    flextable::hline_top(part = "header", border = border) |>
    flextable::hline_bottom(part = "header", border = border) |>
    flextable::hline_bottom(part = "body", border = border)

  vcenter_cols <- which(colnames(d) %in%
                          c("met", "confidence level", "reference database"))
  if (length(vcenter_cols))
    ft <- flextable::valign(ft, j = vcenter_cols, valign = "center",
                            part = "body")

  flextable::set_table_properties(ft, layout = "autofit", width = 1)
}

# 4. Append the tables section --------------------------------------------
wb <- openxlsx::loadWorkbook(XLSX)
sheet_names <- openxlsx::sheets(wb)

doc <- doc |>
  officer::body_add_break() |>
  officer::body_add_par("Supplemental Tables", style = heading_style)

trunc_note_style <- officer::fp_text(font.family = "Times New Roman",
                                     font.size = 9, italic = TRUE)

purrr::walk(sheet_names, function(s) {
  payload <- read_sheet(s)
  n_total <- nrow(payload$data)
  truncated <- n_total > TABLE_ROW_LIMIT
  if (truncated) payload$data <- payload$data[seq_len(TABLE_ROW_LIMIT), ]

  message("Adding ", s, " (", nrow(payload$data),
          if (truncated) sprintf(" of %d", n_total) else "",
          " rows) ...")

  doc <<- doc |>
    officer::body_add_par(payload$title, style = caption_style) |>
    flextable::body_add_flextable(build_flextable(payload$data))

  if (truncated) {
    note <- sprintf(
      "Showing top %d rows of %d total; full table available in the accompanying Excel workbook 'Supplement tables_composite_meta.xlsx'.",
      TABLE_ROW_LIMIT, n_total)
    doc <<- officer::body_add_fpar(doc,
              officer::fpar(officer::ftext(note, prop = trunc_note_style)))
  }

  doc <<- officer::body_add_par(doc, "")
})

print(doc, target = OUT_DOCX)
message("Wrote ", OUT_DOCX)

#--------------------------------End of the code--------------------------------
