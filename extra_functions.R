
#### HELPER FUNCTIONS

# transform for count data
# using sqrt with minimum value of 1, as per advice at
# http://www.webcitation.org/query?url=http%3A%2F%2Fpareonline.net%2Fgetvn.asp%3Fv%3D8%26n%3D6&date=2010-02-11
tr = function(x) return(sqrt(1 + x))

log.tr = function(x) return(log(1 + x))

undo.tr = function(y) return(y^2 - 1)

undo.log.tr = function(y) return(exp(y) - 1)

#### Some of my variables are extracted by looking to see if MEDLINE records have a MeSH term.
# This is fine, but sometimes MEDLINE records are incomplete, they are in process or aren't on the list
# of journals indexed by MeSH indexers.  This means a lack of relevant MeSH term does not mean the MeSH
# term does not apply.  Thus, must replace these 0s with NAs to indicate the data is indeed missing
medline.values = function(raw_values, medline_status) {
    values = raw_values
    values[medline_status!="indexed for MEDLINE"] = NA
    return(values)
}

print.thresh = function(x, thresh, decreasing=TRUE) {
    response = ""
    if (decreasing) {
        x.thresh = x[which(x>=thresh)]
    } else {
        x.thresh = x[which(x<=thresh)]
    }
    x.sorted = sort(x.thresh, decreasing)   
    for (ii in 1:length(x.sorted)) {
        response = cat(response, sprintf("%7.2f %s\n", x.sorted[ii], names(x.sorted)[ii]))
    }
    if (length(response>0)) {
        print(response)
    }
}

####Get scores:
# based on middle of factanal()
factor.scores.bartlett = function(x, fa.fit, na.action=NULL) {
    Lambda <- fa.fit$loadings
    z <- as.matrix(x)
    if (!is.numeric(z)) 
#        stop("factor analysis applies only to numerical variables")
        z = as.matrix(colwise(as.numeric)(x))
    zz <- scale(z, TRUE, TRUE)
    d <- 1/fa.fit$uniquenesses
    tmp <- t(Lambda * d)
    scores <- t(solve(tmp %*% Lambda, tmp %*% t(zz)))
    rownames(scores) <- rownames(z)
    colnames(scores) <- colnames(Lambda)
    if (!is.null(na.action)) 
        scores <- napredict(na.act, scores)
    scores
}

# Extracted then modified from hetcor.data.frame from polycor
adjust.to.positive.definite = function(inputcor) {
    min.eigen = min(eigen(inputcor, only.values=TRUE)$values)
    if (min.eigen < 0){
        print("will try to make correlation matrix positive-definite")
        cor.corrected <- nearcor(inputcor)  # also could try nearPD
        if (!cor.corrected$converged) {
            stop("attempt to make correlation matrix positive-definite failed")
        }
        print("the correlation matrix has been adjusted to make it positive-definite")
        outputcor = cor.corrected$cor
        rownames(outputcor) <- rownames(inputcor)
        colnames(outputcor) <- colnames(inputcor)
    } else {
        print("Correlation matrix already positive-definite, so no adjustment needed")
        outputcor = inputcor
    }
    outputcor
}


# Based on hetcore.data.frame from polycor
# Modified slightly to call different correlation function for continuous correlations
# and print out updates
"hetcor.modified" <-
function(data, ML=FALSE, std.err=TRUE, use=c("complete.obs", "pairwise.complete.obs"),
  bins=4, pd=TRUE, ...){
  se.r <- function(r, n){
    rho <- r*(1 + (1 - r^2)/(2*(n - 3))) # approx. unbiased estimator
    v <- (((1 - rho^2)^2)/(n + 6))*(1 + (14 + 11*rho^2)/(2*(n + 6)))
    sqrt(v)
    }
  use <- match.arg(use)
  if (class(data) != "data.frame") stop("argument must be a data frame.")
  if (use == "complete.obs") data <- na.omit(data)
  p <- length(data)
  if (p < 2) stop("fewer than 2 variables.")
  R <- matrix(1, p, p)
  Type <- matrix("", p, p)
  SE <- matrix(0, p, p)
  N <- matrix(0, p, p)
  Test <- matrix(0, p, p)
  diag(N) <- if (use == "complete.obs") nrow(data)
             else sapply(data, function(x) sum(!is.na(x)))
  for (i in 2:p) {
    print(i)
    for (j in 1:(i-1)){
      x <- data[[i]]
      y <- data[[j]]
      if (inherits(x, c("numeric", "integer")) && inherits(y, c("numeric", "integer"))) {
#         r <- cor(x, y, use="complete.obs")
         r <- rcorr(x, y, type="pearson")$r[1,2]
#         Type[i, j] <- Type[j, i] <- "Pearson"
         Type[i, j] <- Type[j, i] <- "rcorr Pearson"
         R[i, j] <- R[j, i] <- r
         if (std.err) {
           n <- sum(complete.cases(x, y))
           SE[i, j] <- SE[j, i] <- se.r(r, n)
           N[i, j] <- N[j, i] <- n
           Test[i, j] <- pchisq(chisq(x, y, r, bins=bins), bins^2 - 2, lower.tail=FALSE)
           }
         }
      else if (inherits(x, "factor") && inherits(y, "factor")) {
         Type[i, j] <- Type[j, i] <- "Polychoric"
         result <- polychor(x, y, ML=ML, std.err=std.err)
         if (std.err){
           n <- sum(complete.cases(x, y))
           R[i, j] <- R[j, i] <- result$rho
           SE[i, j] <- SE[j, i] <- sqrt(result$var[1,1])
           N[i, j] <- N[j, i] <- n
           Test[i, j] <- if (result$df > 0)
                pchisq(result$chisq, result$df, lower.tail=FALSE)
                else NA
           }
         else R[i, j] <- R[j, i] <- result
         }
       else {
         if (inherits(x, "factor") && inherits(y, c("numeric", "integer")))
           result <- polyserial(y, x, ML=ML, std.err=std.err, bins=bins)
         else if (inherits(x, c("numeric", "integer")) && inherits(y, "factor"))
           result <- polyserial(x, y, ML=ML, std.err=std.err, bins=bins)
         else {
             stop("columns must be numeric or factors.")
             }
         Type[i, j] <- Type[j, i] <- "Polyserial"
         if (std.err){
           n <- sum(complete.cases(x, y))
           R[i, j] <- R[j, i] <- result$rho
           SE[i, j] <- SE[j, i] <- sqrt(result$var[1,1])
           N[i, j] <- N[j, i] <- n
           Test[i, j] <- pchisq(result$chisq, result$df, lower.tail=FALSE)
           }
         else R[i, j] <- R[j, i] <- result
         }
       }
     }
   if (pd) {
       if (min(eigen(R, only.values=TRUE)$values) < 0){
            cor <- nearcor(R)
            if (!cor$converged) stop("attempt to make correlation matrix positive-definite failed")
            warning("the correlation matrix has been adjusted to make it positive-definite")
            R <- cor$cor
        }
    }
   rownames(R) <- colnames(R) <- names(data)
   result <- list(correlations=R, type=Type, NA.method=use, ML=ML)
   if (std.err) {
     rownames(SE) <- colnames(SE) <- names(data)
     rownames(N) <- colnames(N) <- names(N)
     rownames(Test) <- colnames(Test) <- names(data)
     result$std.errors <- SE
     result$n <- if (use == "complete.obs") n else N
     result$tests <- Test
     }
   class(result) <- "hetcor"
   result
   }


# based on heatmap.2
# Changed so that labels are on a different axis, close to clustering
heatmap.3 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    # Also put them on the other side    
    axis(3, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(xlab)) 
        mtext(xlab, side = 3, line = margins[1] - 1.25)
    # Also put them on the top axis
    axis(2, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 2, line = margins[2] - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x))
            tmpbreaks[length(tmpbreaks)] <- max(abs(x))
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


bargraph.CI.ordered = function (x.factor, response, group = NULL, split = FALSE, col = NULL, 
    angle = NULL, density = NULL, lc = TRUE, uc = TRUE, legend = FALSE, 
    ncol = 1, leg.lab = NULL, x.leg = NULL, y.leg = NULL, cex.leg = 1, 
    bty = "n", bg = "white", space = if (split) c(-1, 1), err.width = if (length(levels(as.factor(x.factor))) > 
        10) 0 else 0.1, err.col = "black", err.lty = 1, fun = function(x) mean(x, 
        na.rm = TRUE), ci.fun = function(x) c(fun(x) - se(x), 
        fun(x) + se(x)), ylim = NULL, xpd = FALSE, data = NULL, 
    subset = NULL, ...) 
{
    subset <- eval(substitute(subset), envir = data)
    if (!is.null(data)) {
        if (!is.null(subset)) 
            data <- subset(data, subset)
        x.factor <- eval(substitute(x.factor), envir = data)
        response <- eval(substitute(response), envir = data)
        group <- eval(substitute(group), envir = data)
    }
    subset = NULL
    if (split) {
        if (length(group[[1]]) > 1) {
            print("Error: Can't split for > 2 group levels")
            stop()
        }
        split.fn <- function(x, y) {
            val <- levels(group)[[1]]
            ifelse(x == val, y * -1, y)
        }
        response <- split.fn(group, response)
    }
    if (is.null(group)) 
        groups = x.factor
    else {
        if (length(group[[1]]) > 1) {
            group <- do.call("paste", c(group, sep = "."))
        }
        groups <- list(group, x.factor)
    }
    mn.data.unordered <- tapply(response, groups, fun)
    mn.data <- sort(mn.data.unordered, decreasing=TRUE)
    CI.dat.unordered <- tapply(response, groups, ci.fun)
    CI.dat <- CI.dat.unordered[order(mn.data.unordered, decreasing=TRUE)]
    null.fn <- function(x) if (is.null(x[[1]])) 
        rep(NaN, 2)
    else x
    if (!is.null(group)) 
        CI.dat <- apply(CI.dat, c(1, 2), null.fn)
    CI.data <- array(unlist(CI.dat), c(2, if (is.null(group)) 1 else nrow(mn.data), 
        length(levels(as.factor(x.factor)))))
    CI.L <- CI.data[1, , ]
    CI.H <- CI.data[2, , ]
    replace.NA <- function(x) if (is.na(x)) 
        0
    else x
    if (!is.null(group)) {
        CI.L <- apply(CI.L, c(1, 2), replace.NA)
        CI.H <- apply(CI.H, c(1, 2), replace.NA)
    }
    else {
        CI.L <- as.vector(unlist(lapply(CI.L, replace.NA)))
        CI.H <- as.vector(unlist(lapply(CI.H, replace.NA)))
    }
    if (is.null(ylim)) 
        ylim <- c(min(0, CI.L), max(CI.H))
    xvals <- barplot(mn.data, ylim = ylim, beside = TRUE, col = col, 
        density = density, angle = angle, space = space, xpd = xpd, 
        ...)
    nlevels.x <- dim(mn.data)[1]
    if (is.null(group)) 
        arrows(xvals, if (lc) 
            CI.L
        else mn.data, xvals, if (uc) 
            CI.H
        else mn.data, angle = 90, length = err.width, code = 3)
    else {
        nlevels.y <- dim(mn.data)[2]
        for (i in 1:nlevels.y) arrows(xvals[, i], if (lc) 
            CI.L[, i]
        else mn.data[, i], xvals[, i], if (uc) 
            CI.H[, i]
        else mn.data[, i], angle = 90, code = 3, col = err.col, 
            lty = err.lty, length = err.width)
    }
    if (!is.null(group) & legend) {
        legend(x = if (is.null(x.leg)) 
            0.8 * max(xvals)
        else x.leg, y = if (is.null(y.leg)) 
            max(ylim)
        else y.leg, legend = if (is.null(leg.lab)) 
            levels(as.factor(group))
        else leg.lab, bty = bty, bg = bg, ncol = ncol, fill = if (is.null(col)) 
            gray.colors(nlevels.x)
        else col, density = density, angle = angle, cex = cex.leg)
    }
    invisible(list(xvals = xvals, vals = mn.data, CI = CI.data))
}


plot.summary.rms.norangelabels = function (x, at, log = FALSE, q = c(0.7, 0.8, 0.9, 0.95, 0.99), 
    xlim, nbar, cex = 1, nint = 10, cex.c = 0.5, cex.t = 1, clip = c(-1e+30, 
        1e+30), main, ...) 
{
    scale <- attr(x, "scale")
    adjust <- attr(x, "adjust")
    Type <- x[, "Type"]
    x <- x[Type == 1, , drop = FALSE]
    lab <- dimnames(x)[[1]]
    effect <- x[, "Effect"]
    se <- x[, "S.E."]
    if (!log && any(Type == 2)) {
        fun <- exp
        tlab <- scale[2]
    }
    else {
        fun <- function(x) x
        if (log) {
            if (length(scale) == 2) 
                tlab <- scale[2]
            else tlab <- paste("exp(", scale[1], ")", sep = "")
        }
        else tlab <- scale[1]
    }
    if (!length(scale)) 
        tlab <- ""
    if (!missing(main)) 
        tlab <- main
    augment <- if (log | any(Type == 2)) 
        c(0.1, 0.5, 0.75, 1)
    else 0
    n <- length(effect)
    out <- qnorm((max(q) + 1)/2)
    if (missing(xlim) && !missing(at)) 
        xlim <- range(if (log) logb(at) else at)
    else if (missing(xlim)) {
        xlim <- fun(range(c(effect - out * se, effect + out * 
            se)))
        xlim[1] <- max(xlim[1], clip[1])
        xlim[2] <- min(xlim[2], clip[2])
    }
    else augment <- c(augment, if (log) exp(xlim) else xlim)
    fmt <- function(k) {
        m <- length(k)
        f <- character(m)
        for (i in 1:m) f[i] <- format(k[i])
        f
    }
    #lb <- ifelse(is.na(x[, "Diff."]), lab, paste(lab, " - ", 
    #    fmt(x[, "High"]), ":", fmt(x[, "Low"]), sep = ""))
    lb <- lab
    plot.new()
    par(new = TRUE)
    mxlb <- 0.1 + max(strwidth(lb, units = "inches", cex = cex))
    tmai <- par("mai")
    on.exit(par(mai = tmai))
    par(mai = c(tmai[1], mxlb, 1.5 * tmai[3], tmai[4]))
    outer.widths <- fun(effect + out * se) - fun(effect - out * 
        se)
    if (missing(nbar)) 
        nbar <- n
    npage <- ceiling(n/nbar)
    is <- 1
    for (p in 1:npage) {
        ie <- min(is + nbar - 1, n)
        plot(1:nbar, rep(0, nbar), xlim = xlim, ylim = c(1, nbar), 
            type = "n", axes = FALSE, xlab = "", ylab = "")
        if (cex.t > 0) 
            title(tlab, cex = cex.t)
        lines(fun(c(0, 0)), c(nbar - (ie - is), nbar), lty = 2)
        if (log) {
            pxlim <- pretty(exp(xlim), n = nint)
            pxlim <- sort(unique(c(pxlim, augment)))
            pxlim <- pxlim[pxlim >= exp(xlim[1])]
            if (!missing(at)) 
                pxlim <- at
            axis(3, logb(pxlim), lab = format(pxlim))
        }
        else {
            pxlim <- pretty(xlim, n = nint)
            pxlim <- sort(unique(c(pxlim, augment)))
            pxlim <- pxlim[pxlim >= xlim[1]]
            if (!missing(at)) 
                pxlim <- at
            axis(3, pxlim)
        }
        imax <- (is:ie)[outer.widths[is:ie] == max(outer.widths[is:ie])][1]
        for (i in is:ie) {
            confbar(nbar - (i - is + 1) + 1, effect[i], se[i], 
                q = q, type = "h", fun = fun, cex = cex.c, labels = i == 
                  imax, clip = clip, ...)
            mtext(lb[i], 2, 0, at = nbar - (i - is + 1) + 1, 
                cex = cex, adj = 1, las = 1)
        }
        if (adjust != "") {
            adjto <- paste("Adjusted to:", adjust, sep = "")
            xx <- par("usr")[2]
            if (nbar > ie) 
                text(xx, nbar - (ie - is + 1), adjto, adj = 1, cex = cex)
#            else title(sub = adjto, adj = 1, cex = cex)
        }
        is <- ie + 1
    }
    invisible()
}
