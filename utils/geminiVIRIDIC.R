## Libraries----------------------------------------------------------------
suppressPackageStartupMessages(
    {
library(magrittr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(seqinr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(IRanges, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(pheatmap, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(fastcluster, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(parallelDist, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(furrr, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(future, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
    })
## Initialize project-------------------------------------------------------
input0 <- as.character(commandArgs(trailingOnly = TRUE))
in1 <- str_which(input0, "blastres=")
in2 <- str_which(input0, "out=")
input <- ""
input[1] <- str_remove(input0[in1], "blastres=")
input[2] <- str_remove(input0[in2], "out=")
blastn_out_path <- input[1]
outfmt<- ' -outfmt "6 qseqid sseqid evalue bitscore qlen slen qstart qend sstart send qseq sseq nident gaps"'
## functions  -------------------------------------------------------
nident_fun <- function(DF, gen_len)
{
  if(nrow(DF)== 1)
  {
    dfs <- DF%>%
      mutate(nident_recalc = nident)%>%
      mutate(alig_q = abs(qstart-qend)+1)
  }else
  {
    dfs <- DF %>%
      filter(qstart == 1, qend == gen_len)
    if(nrow(dfs) == 1)
    {
      dfs <- dfs%>%
        mutate(nident_recalc = nident)%>%
        mutate(alig_q = abs(qstart-qend)+1)
    } else
    {
      if(nrow(dfs) == 0)
      {
        rm(dfs)
        ir <- IRanges(start = as.numeric(DF$qstart), end = as.numeric(DF$qend), names = seq(from = 1, to = nrow(DF)))
        cov <- coverage(ir) %>%
          as.vector()
        cov2 <- cov[cov >1]
        if(length(cov2) == 0)
        {
          dfs <- DF %>%
            mutate(nident_recalc = nident)%>%
            mutate(alig_q = abs(qstart-qend)+1)
        }else
        {
          rm(ir, cov, cov2)
          dfs <- DF %>%
            mutate(alig_q = abs(qstart-qend)+1) %>%
            arrange(desc(alig_q))
          ##remove all hits incuded in another hit (can I maybe use the code below with the cov? because then I get rid of 2 for loops)
          for(a in nrow(dfs):2)
          {
            for(b in 1:(a-1))
            {
              if(dfs$qstart[b] <= dfs$qstart[a] & dfs$qend[a] <= dfs$qend[b])
              {
                dfs <- dfs[-a, ]
                break()
              }
            }
            rm(b)
          }
          rm(a)
          ir <- IRanges(start = dfs$qstart, end = dfs$qend, names = seq(from = 1, to = nrow(dfs)))
          cov <- coverage(ir) %>%
            as.vector()
          cov2 <- cov[cov >1]
          if(length(cov2) > 0)
          {
            #remove hits that are completely overlapping two other ranges
            for(i in nrow(dfs):1)
            {
              cov_r <- cov[dfs$qstart[i]:dfs$qend[i]] %>%
                min()
              if(cov_r > 1)
              {
                dfs <- dfs[-i, ]
                ir <- IRanges(start = dfs$qstart, end = dfs$qend, names = seq(from = 1, to = nrow(dfs)))
                cov <- coverage(ir) %>%
                  as.vector()
              }
              rm(cov_r)
            }
            rm(i)
            cov2 <- cov[cov >1]
            if(length(cov2) == 0)
            {
              dfs <- dfs %>%
                mutate(nident_recalc = nident)
            }else
            {
              dfs <- dfs %>%
                arrange(dplyr::desc(qstart)) %>%
                mutate(nident_recalc = nident) %>%
                mutate(qstart_recalc = qstart) %>%
                mutate(qend_recalc = qend)
              for(a in nrow(dfs):2)
              {
                for(b in (a-1):1)
                {
                  if((dfs$qend_recalc[a] >= dfs$qstart_recalc[b]) & (dfs$qstart_recalc[a] < dfs$qstart_recalc[b]))
                  {
                    overlap <- dfs$qend_recalc[a] - dfs$qstart_recalc[b] + 1
                    q_over_a <- dfs$qseq[a] %>%
                      str_replace_all("-", "") %>%
                      str_sub(start = -overlap, end = -1) %>%
                      s2c %>%
                      paste0("-*") %>%
                      c2s()
                    q_over_a_ext <- dfs$qseq[a] %>%
                      str_extract(q_over_a) %>%
                      s2c()
                    s_over_a <- dfs$sseq[a] %>%
                      str_sub(start = -length(q_over_a_ext), end = -1) %>%
                      s2c()
                    diffe_a <- str_match(q_over_a_ext, s_over_a) %>%
                      is.na() %>%
                      sum()
                    dfs[a,"nident_recalc"] <- dfs$nident[a] - (length(q_over_a_ext) - diffe_a)
                    dfs[a, "qend_recalc"] <- dfs$qstart_recalc[b]-1
                    rm(overlap, q_over_a, q_over_a_ext, s_over_a, diffe_a)
                  }
                }
              }
            }
          }else
          {
            dfs <- dfs %>%
              mutate(nident_recalc = nident)
          }
        }
      }else
      {
        dfs <- dfs%>%
          mutate(nident_recalc = 100000)%>%
          mutate(alig_q = abs(qstart-qend)+1)
      }
    }
  }
  nident_o <- sum(dfs$nident_recalc)
  align_q_o <- sum(dfs$alig_q)
  outp <- c(nident_o, align_q_o)
  return(outp)
}
s_nident_fun <- function(q, s, DF)
{
  DF2 <- DF %>%
    filter(sseqid == q, qseqid == s)
  if(nrow(DF2) == 1)
  {
    ret <- c(DF2$q_nident, DF2$q_fract_aligned)
  }else
  {
    if(nrow(DF2) == 0)
    {
      ret <- c(0, 0)
    }else
    {ret <- c(-100, -100)}
  }
  return(ret)
}
det_dist_fun <- function(simDF, th, sim_dist) ##th = 95 or 70
{
  simDF_ma <- as.matrix(simDF)
  th <- as.numeric(th)
  if(sim_dist == "distance")
  {
    th <- 100 - th
    elems <- base::which(simDF <= th)
    intergenome <- simDF_ma[elems] %>%
      max()
  }else
  {
    elems <- base::which(simDF >= th)
    intergenome <- simDF_ma[elems] %>%
      min()
  }
  sel_elem <- which(simDF == intergenome)
  return(sel_elem)
}
## -CALCULATE Intergenomic similarities/distances-------------------------------------------------
print("Read blast results")
blastn_DF <- data.table::fread(blastn_out_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, drop = c(3, 4))
colnames(blastn_DF) <- c("qseqid", "sseqid", "qlen", "slen", "qstart", "qend", "sstart", "send", "qseq", "sseq", "nident", "gaps")
print("Make group_by")
blastn_DF_g <- blastn_DF %>%
  group_by(qseqid, sseqid, qlen, slen) %>%
  nest()
rm(blastn_DF)
plan(multisession)
print("Analysis")
options(future.globals.maxSize= 256 * 1024 * 1024^2)
outp_ls <- future_map2(.x = blastn_DF_g$data, .y = blastn_DF_g$qlen, .f = nident_fun)
blastn_DF_g[, "q_nident"] <- lapply(outp_ls, "[[", 1) %>%
  unlist()
blastn_DF_g[, "q_aligned"] <- lapply(outp_ls, "[[", 2) %>%
  unlist()
blastn_DF_g <- blastn_DF_g %>%
  mutate(q_fract_aligned = round(q_aligned/qlen, digits = 2))
rm(outp_ls)
options(future.globals.maxSize= 256 * 1024 * 1024^2)
ret_ls <- future_map2(.x = blastn_DF_g$qseqid, .y = blastn_DF_g$sseqid, .f = s_nident_fun, DF = blastn_DF_g)
blastn_DF_g[, "s_nident"] <- lapply(ret_ls, "[", 1) %>%
  unlist
blastn_DF_g[, "s_fract_aligned"] <- lapply(ret_ls, "[", 2) %>%
  unlist
rm(ret_ls)
blastn_DF_gbkp <- blastn_DF_g
blastn_DF_g <- blastn_DF_gbkp %>%
  select(-data) %>%
  ungroup() %>%
  mutate(qs_nident = (as.numeric(q_nident) + as.numeric(s_nident))) %>%
  mutate(qs_len = (qlen + slen)) %>%
  mutate(interg_sim = (qs_nident/qs_len)*100) %>%
  mutate(interg_sim = round(interg_sim, digits =3)) %>%
  mutate(fract_qslen = (pmin(qlen, slen)/pmax(qlen,slen))) %>%
  mutate(fract_qslen = round(fract_qslen, digits = 1)) %>%
  select(qseqid, sseqid, qlen, slen, q_nident, q_aligned, s_nident, qs_nident, qs_len, interg_sim, fract_qslen, q_fract_aligned, s_fract_aligned)
blastn_DF_g <- blastn_DF_g %>%
  as.data.frame(stringsAsFactors = FALSE)
## SAVE outputs ---------------------------------------------------
#save calculated BlastDF as csv
write.table(x = blastn_DF_g, file = paste0(input[2]), sep = "\t", row.names = FALSE, col.names = TRUE)