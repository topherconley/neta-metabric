
probe2gene <- function()  { 
  
  
  #GET DIMENSIONS OF Y-Y NETWORK
  #probes without entrez id
  noeg_idx <- which(is.na(yinfo$EntrezReannotated))
  probe_noeg <- yinfo$id[noeg_idx]
  #entrez id with at least one probe
  egtab <- table(yinfo$EntrezReannotated)
  eg2probeid <- lapply(names(egtab), function(eg) {
    yinfo$id[yinfo$EntrezReannotated %in% eg]  
  })
  names(eg2probeid) <- names(egtab)
  #FORM NEW NETWORK
  qnew <- length(noeg_idx) + length(eg2probeid)
  cnet <- list(yy = matrix(data = 0, nrow = qnew, ncol = qnew),
               xy =  matrix(data = 0, nrow = nrow(net$xy), ncol = qnew))
               
  #NEW ID ORDER
  #(1) entrez with single probe
  id1 <- names(eg2probeid)[sapply(eg2probeid, length) == 1]
  #(2) entrez with multiple probes
  id2 <- names(eg2probeid)[sapply(eg2probeid, length) > 1]
  #(3) probe with no entrez
  id3 <- probe_noeg
  #new id vector 
  idnew <- c(id1, id2, id3)
  
  #entrez genes with 
  meg <- eg2probeid[sapply(eg2probeid, length) > 1]
  #take union of edges between probes with same Entrez gene id
  
  for(m in meg) { 
    uyedge <- apply(X = net$yy[m,], MARGIN = 2, FUN = max)
    uxedge <- apply(X = net$xy[,m], MARGIN = 1, FUN = max)
    net$yy[m[1],] <- uyedge
    net$xy[,m[1]] <- uxedge
  }
  #drop secondary listing of probes
  secprobes <- unlist(sapply(meg, function(m) m[2:length(m)]))
  keepers <- setdiff(rownames(net$yy), secprobes)
  net$yy <- net$yy[keepers,keepers]
  net$xy <- net$xy[,keepers]
  #adjust order of yinfo
  yinfo <- yinfo[match(keepers,yinfo$id),]
  bdeg$yy <- bdeg$yy[,keepers]
  
  
  
  
  probe2eg <- yinfo[,c("id", "EntrezReannotated")]
  #all duplicates should be missing values now. 
  eg <- probe2eg$EntrezReannotated
  stopifnot(all(is.na(eg[which(duplicated(eg))])))
  #replace ID with Entrez gene
  yinfo$probeid <- yinfo$id
  yinfo$id <- ifelse(is.na(eg), yinfo$id, eg)
  stopifnot(!anyDuplicated(yinfo$id))
  
  
  rownames(net$xy) <- xinfo$id; colnames(net$xy) <- yinfo$id;
  rownames(net$yy) <- yinfo$id; colnames(net$yy) <- yinfo$id;
  colnames(bdeg$xy) <- xinfo$id; colnames(bdeg$yy) <- yinfo$id; 
  
}