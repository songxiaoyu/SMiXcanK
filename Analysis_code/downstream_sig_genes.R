genes <- c("PSMB1","SLC4A7","CASP8","FAR2","EPN2","EP300","CELF4","CRHR1","GTF3A",
           "PDLIM4","USP37","POC1B","RMC1","WDPCP","CD200R1","TRMT61A","STXBP4","CTSW",
           "SNX32","ZNF169","DPP7","RMI1","RPRML","C8orf33","FES","PLEKHM1","AC008267.5",
           "DPP3","AC073912.1","IQCH-AS1","AC125494.2","AC016355.1")

# Install if needed:
# install.packages(c("httr2","jsonlite","dplyr","tibble","stringr"))

library(httr2)
library(jsonlite)
library(dplyr)
library(tibble)
library(stringr)

endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

# 1) First resolve gene symbols -> Ensembl IDs using OT "search"
resolve_ensembl <- function(symbol){
  q <- '
  query($q: String!) {
    search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
      hits { id entity }
    }
  }'
  body <- list(query = q, variables = list(q = symbol))
  res <- request(endpoint) |>
    req_method("POST") |>
    req_headers(`Content-Type`="application/json") |>
    req_body_json(body) |>
    req_perform()
  dat <- resp_body_json(res)
  if (length(dat$data$search$hits) == 0) return(NA_character_)
  dat$data$search$hits[[1]]$id
}

ensembl <- setNames(vapply(genes, resolve_ensembl, character(1)), genes)

# 2) For each Ensembl ID, ask OT for disease associations and keep breast cancer–related ones
fetch_bc_assoc <- function(ensg){
  if (is.na(ensg)) return(tibble())
  q <- '
  query($ensg: String!) {
    target(ensemblId: $ensg) {
      associatedDiseases(page: {index: 0, size: 200}) {
        rows {
          disease { id name }
          score
          datasourceScores { id score }
        }
      }
    }
  }'
  body <- list(query = q, variables = list(ensg = ensg))
  res <- request(endpoint) |>
    req_method("POST") |>
    req_headers(`Content-Type`="application/json") |>
    req_body_json(body) |>
    req_perform()
  dat <- resp_body_json(res)
  rows <- dat$data$target$associatedDiseases$rows
  if (is.null(rows) || length(rows) == 0) return(tibble())
  tibble(
    disease = vapply(rows, \(x) x$disease$name, character(1)),
    score   = vapply(rows, \(x) x$score, numeric(1))
  ) |>
    filter(str_detect(tolower(disease), "breast"))
}

ot_hits <- lapply(ensembl, fetch_bc_assoc)
ot_summary <- tibble(
  gene_name = names(ot_hits),
  ensembl = unname(ensembl),
  has_breast_assoc = vapply(ot_hits, nrow, integer(1)) > 0,
  top_breast_assoc = vapply(ot_hits, \(df) if(nrow(df)) df$disease[which.max(df$score)] else NA_character_, character(1)),
  top_score = vapply(ot_hits, \(df) if(nrow(df)) max(df$score) else NA_real_, numeric(1))
)
merged_ot <-merge(ot_summary, out_ct2_final, on='gene_name')

ot_summary |> arrange(desc(has_breast_assoc), desc(top_score))
write.csv(merged_ot, '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/sig_genes_literature_ct2.csv', row.names = FALSE)

# 3 cell type
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")
merged$fdr_p_join_ct3 <- p.adjust(merged$p_join_ct3, method = "BH")
genes <- merged[merged$fdr_p_join_ct3  < fdr_cutoff,'gene_name']
ensembl <- setNames(vapply(genes, resolve_ensembl, character(1)), genes)


ot_hits <- lapply(ensembl, fetch_bc_assoc)
ot_summary <- tibble(
  gene_name = names(ot_hits),
  ensembl = unname(ensembl),
  has_breast_assoc = vapply(ot_hits, nrow, integer(1)) > 0,
  top_breast_assoc = vapply(ot_hits, \(df) if(nrow(df)) df$disease[which.max(df$score)] else NA_character_, character(1)),
  top_score = vapply(ot_hits, \(df) if(nrow(df)) max(df$score) else NA_real_, numeric(1))
)
merged_ot <-merge(ot_summary, out_ct3_final, on='gene_name')

write.csv(merged_ot, '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/sig_genes_literature_ct3.csv', row.names = FALSE)

