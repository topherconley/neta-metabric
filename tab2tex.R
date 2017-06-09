
xhubs2textab <- function(xhubs, cap, tfile) { 
  library(xtable)
  colnames(xhubs)[1] <- "CNA hub"
  colnames(xhubs)[3] <- "pot. # cis"
  txhubs <- xtable(
    x = xhubs, 
    caption = cap,
    label = "tab:MatrixData",
    align = c("l|", rep("c", ncol(xhubs) - 1), "p{1.5in}"),
    scalebox = 0.75,
    digits = c(rep(3, ncol(xhubs)+1))
  )
  print(
    txhubs,
    table.placement = "H",
    caption.placement = "top",
    include.rownames = FALSE,
    include.colnames = TRUE,
    size = "small",
    tabular.environment = 'longtable',
    floating = FALSE,
    add.to.row = list(pos = list(0),command = "\\hline \\endhead "),
    file = tfile
  )  
}

mod2textab <- function(etab, cap, tfile) { 
  te <- xtable(
    x = etab,
    caption = cap,
    label = "tab:MatrixData",
    align = c("l|", "l", "p{2in}", "l", "p{1in}"),scalebox = 0.75,
    digits = c(rep(3, ncol(xhubs)+1))
  )
  print(
    te,
    table.placement = "H",
    caption.placement = "top",
    include.rownames = FALSE,
    include.colnames = TRUE,
    size = "small",
    tabular.environment = 'longtable',
    floating = FALSE,
    add.to.row = list(pos = list(0),command = "\\hline \\endhead "),
    file = tfile
  )
}