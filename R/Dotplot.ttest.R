DotPlot.ssc <- function (input, gene, group.by, split.by = NULL, point.by = "Patient", 
                         num.cutoff = 10, title = "", number_labels = T, text_sizes = c(20, 
                                                                                        10, 5, 10, 5, 5, 3), spread = NULL, theme = "classic", 
                         jitter.alpha = 0.5, jitter.width = 0.2, jitter.size = 2, 
                         col.trans = "", bar = T, link = FALSE, ttest = F) 
{
  group.by.levels <- levels(input@meta.data[, group.by])
  if (!is.null(split.by)) split.by.levels <- levels(input@meta.data[, split.by])
  meta <- input@meta.data[, c(group.by, split.by, point.by)]
  tmp <- tapply(input@assays$RNA@data[gene, ], apply(meta, 
                                                     1, function(x) paste0(x, collapse = ":")), function(x) c(mean(x), 
                                                                                                              length(x), sum(x > 0)))
  df <- data.frame(gene = sapply(tmp, function(x) x[1]), num = sapply(tmp, 
                                                                      function(x) x[2]), num_above0 = sapply(tmp, function(x) x[3]))
  df$num_cutoff <- df$num
  df$num_cutoff[df$num < num.cutoff] <- NA
  df.meta <- as.data.frame(t(as.data.frame(sapply(names(tmp), 
                                                  function(x) strsplit(x, ":")))))
  colnames(df.meta) <- c(group.by, split.by, point.by)
  df <- cbind(df, df.meta)
  tmp <- tapply(input@assays$RNA@data[gene, ], 
                apply(meta[, c(group.by, split.by), drop=F], 1, function(x) paste0(x, collapse = ":")), 
                function(x) c(mean(x), round(sum(x > 0)/length(x), digits = 3)))
  df.stat <- data.frame(bar_mean = sapply(tmp, function(x) x[1]), 
                        frac_sum = sapply(tmp, function(x) x[2]))
  df.meta <- as.data.frame(t(as.data.frame(sapply(names(tmp), 
                                                  function(x) strsplit(x, ":")))))
  colnames(df.meta) <- c(group.by, split.by)
  df.stat <- cbind(df.stat, df.meta)
  if (any(!is.null(spread))) {
    others <- setdiff(unique(df[, spread[1]]), spread[2])
    ind <- which(df[, spread[1]] == spread[2])
    rmdf <- df[ind, ]
    df <- df[-ind, ]
    ind.stat <- which(df.stat[, spread[1]] == spread[2])
    rmdf.stat <- df.stat[ind.stat, ]
    df.stat <- df.stat[-ind.stat, ]
    for (i in 1:length(others)) {
      rmdf[, spread[1]] <- others[i]
      df <- rbind(df, rmdf)
      rmdf.stat[, spread[1]] <- others[i]
      df.stat <- rbind(df.stat, rmdf.stat)
    }
  }
  df[, group.by] <- factor(df[, group.by], levels = group.by.levels)
  df.stat[, group.by] <- factor(df.stat[, group.by], levels = group.by.levels)
  if (!is.null(split.by)){
    df[, split.by] <- factor(df[, split.by], levels = split.by.levels)
    df.stat[, split.by] <- factor(df.stat[, split.by], levels = split.by.levels)
  }
  
  g <- ggplot(df)
  if (number_labels == T) {
    g <- g + stat_summary(aes_string(x = group.by, y = "num"), 
                          fun.data = function(x) {
                            return(c(y = -max(df$gene)/25, label = sum(x)))
                          }, colour = "black", geom = "text", size = text_sizes[7])
    g <- g + geom_text(data = df.stat, aes_string(x = group.by, 
                                                  y = -max(df$gene)/10, label = "frac_sum"), size = text_sizes[7])
  }
  if (theme == "bw") {
    g <- g + theme_bw()
  }else {
    g <- g + theme_classic()
  }
  if (title == "") {
    title <- gene
    g <- g + labs(title = title, y = gene)
  }else {
    g <- g + labs(title = title, y = gene)
  }
  g <- g + theme(plot.title = element_text(size = text_sizes[1]), 
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), 
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + geom_jitter(aes_string(x = group.by, y = "gene", 
                                  col = "num_cutoff"), width = jitter.width, size = jitter.size, 
                       alpha = jitter.alpha)
  g <- g + scale_color_gradientn(colors = c("blue", "red", 
                                            "orange"))
  if (bar) {
    g <- g + geom_col(data = df.stat, aes_string(x = group.by, 
                                                 y = "bar_mean", fill = group.by), alpha = 0.3)
  }
  if (link) {
    g <- g + geom_line(aes_string(x = group.by, y = "gene", 
                                  group = point.by), alpha = 0.6)
  }
  if (length(split.by) == 1) {
    g <- g + facet_grid(facets = reformulate(split.by), 
                        scales = "free_x", space = "free_x")
  }else if (length(split.by) == 2) {
    g <- g + facet_grid(facets = reformulate(split.by[1], 
                                             split.by[2]), scales = "free_x", space = "free_x")
  }else if (length(split.by) > 2) {
    stop("Parameter split.by needs to be a string with equal or less than two variables.")
  }
  g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
  
  if (ttest == F){
    return(g)
  }else{
    # do t-test on dots pass filter (with more than a certain number of cells)
    df.ttest <- df[!is.na(df$num_cutoff),]
    list.ttest <- list()
    if (is.null(split.by)){
      #  t-test all
      list.ttest[[1]] <- df.ttest
    }else{
      # plot split in multiple panels. Do t-test across groups with in each panel
      list.ttest <- split(df.ttest, df.ttest[,split.by])
    }
    
    
    ttestresult <- lapply(list.ttest, function(tb){
      groups <- unique(as.character(tb[,group.by]))
      comparisons <- as.data.frame(t(combn(groups,2)))
      colnames(comparisons) <- c("group1","group2")
      t.res <- data.frame()
      for (i in 1:nrow(comparisons)){
        tb2 <- split(tb, as.character(tb[,group.by]))
        if (nrow(tb2[[1]]) == 1 | nrow(tb2[[2]]) == 1) next
        if (link){
          # if link is true, then do paired t-test with connected points.
          pointpaired <- intersect(tb2[[1]][,point.by], tb2[[2]][,point.by])
          if (length(pointpaired) < 2) message("Less than 3 points paired")
          tb2[[1]] <- tb2[[1]][match(pointpaired, tb2[[1]][,point.by]),]
          tb2[[2]] <- tb2[[2]][match(pointpaired, tb2[[2]][,point.by]),]
          tmp <- t.test(x=tb2[[comparisons[i,1]]][,"gene"], y=tb2[[comparisons[i,2]]][,"gene"], paired = T)
          
        }else{
          # do with all samples with regular t-test
          tmp <- t.test(x=tb2[[comparisons[i,1]]][,"gene"], y=tb2[[comparisons[i,2]]][,"gene"], paired = F)
        }
        t.res <- rbind(t.res,
                     data.frame(t = tmp$statistic,
                                df = tmp$parameter,
                                p = tmp$p.value,
                                confint.low = tmp$conf.int[1],
                                confint.high = tmp$conf.int[2],
                                mean.1 = tmp$estimate[1],
                                mean.2 = tmp$estimate[2],
                                stderr = tmp$stderr))
      }
      if (nrow(t.res) > 0){
        res <- cbind(comparisons, t.res)
        return(res)
      }else{
        return(NULL)
      }

    })
    output <- list(g=g, ttest=ttestresult)
    return(output)
    
  }
  
  
  
}
