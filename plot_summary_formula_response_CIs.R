dotchart2.CIs = 
function (data, labels, groups = NULL, gdata = NA, horizontal = TRUE, 
    pch = 16, xlab = "", ylab = "", auxdata, auxgdata = NULL, 
    auxtitle, lty = if (.R.) 1 else 2, lines = TRUE, dotsize = 0.8, 
    cex = par("cex"), cex.labels = cex, cex.group.labels = cex.labels * 
        1.25, sort. = TRUE, add = FALSE, dotfont = par("font"), 
    groupfont = 2, reset.par = add, xaxis = TRUE, width.factor = 1.1, 
    lcolor = if (.R.) "gray" else par("col"), ...) 
{
    if (.R. && !add) {
        plot.new()
        par(new = TRUE)
    }
    ieaux <- if (missing(auxdata)) 
        FALSE
    else is.expression(auxdata)
    mtextsrt <- function(..., srt = 0) if (.R.) 
        mtext(..., las = 1)
    else mtext(..., srt = srt)
    ndata <- length(data)
    if (missing(labels)) {
        if (!is.null(names(data))) 
            labels <- names(data)
        else labels <- paste("#", seq(along = ndata))
    }
    else labels <- rep(as.character(labels), length = ndata)
    if (missing(groups)) {
        glabels <- NULL
        gdata <- NULL
    }
    else {
        if (!sort.) {
            ug <- unique(as.character(groups))
            groups <- factor(as.character(groups), levels = ug)
        }
        groups <- oldUnclass(groups)
        glabels <- levels(groups)
        gdata <- rep(gdata, length = length(glabels))
        ord <- order(groups, seq(along = groups))
        groups <- groups[ord]
        data <- data[ord]
        labels <- labels[ord]
        if (!missing(auxdata)) 
            auxdata <- auxdata[ord]
    }
    alldat <- c(data, gdata)
    if (!missing(auxdata)) {
        auxdata <- c(auxdata, auxgdata)
        if (!ieaux) 
            auxdata <- format(auxdata)
    }
    alllab <- paste(c(labels, glabels), "")
    tcex <- par("cex")
    tmai <- par("mai")
    oldplt <- par("plt")
    if (reset.par) 
        on.exit(par(mai = tmai, cex = tcex, usr = tusr))
    par(cex = cex)
    mxlab <- 0.1 + max(strwidth(labels, units = "inches", cex = cex.labels), 
        if (length(glabels)) strwidth(glabels, units = "inches", 
            cex = cex.group.labels)) * width.factor
    if (horizontal) {
        tmai2 <- tmai[3:4]
        if (!missing(auxdata)) 
            tmai2[2] <- 0.2 + width.factor * max(strwidth(if (ieaux) auxdata else format(auxdata), 
                units = "inches", cex = cex.labels))
        par(mai = c(tmai[1], mxlab, tmai2))
        if (!add) 
            plot(alldat, seq(along = alldat), type = "n", ylab = "", 
                axes = FALSE, xlab = "", ...)
            plot(alldat, seq(along = alldat), type = "n", ylab = "", 
                axes = FALSE, xlab = "", ...)
        logax <- par("xaxt") == "l"
    }
    else {
        par(mai = c(mxlab, tmai[2:4]))
        if (!add) 
            plot(seq(along = alldat), alldat, type = "n", xlab = "", 
                axes = FALSE, ylab = "", ...)
        logax <- par("yaxt") == "l"
    }
    tusr <- par("usr")
    if (!add && logax) {
        if (horizontal) 
            abline(v = 10^tusr[1:2], h = tusr[3:4])
        else abline(v = tusr[1:2], h = 10^tusr[3:4])
    }
    else if (!add) 
        abline(v = tusr[1:2], h = tusr[3:4])
    den <- ndata + 2 * length(glabels) + 1
    if (horizontal) {
        if (!add && xaxis) 
            mgp.axis(1, axistitle = xlab)
        delt <- (-(tusr[4] - tusr[3]))/den
        ypos <- seq(tusr[4], by = delt, length = ndata)
    }
    else {
        if (!add) 
            mgp.axis(2, axistitle = xlab)
        delt <- (tusr[2] - tusr[1])/den
        ypos <- seq(tusr[1], by = delt, length = ndata)
    }
    if (!missing(groups)) {
        ypos1 <- ypos + 2 * delt * (if (length(groups) > 1) 
            cumsum(c(1, diff(groups) > 0))
        else 1)
        diff2 <- c(3 * delt, diff(ypos1))
        ypos2 <- ypos1[abs(diff2 - 3 * delt) < abs(0.001 * delt)] - 
            delt
        ypos <- c(ypos1, ypos2) - delt
    }
    ypos <- ypos + delt
    nongrp <- 1:ndata
    if (horizontal) {
        xmin <- par("usr")[1]
        if (!add && lines) 
            abline(h = ypos[nongrp], lty = lty, lwd = 1, col = lcolor)
        points(alldat, ypos, pch = pch, cex = dotsize * cex, 
            font = dotfont)
        for (i in nongrp) {  ##HAP (note, not in vertical plots yet)
            N = as.numeric(auxdata[i])
            CIs = binconf(alldat[i]*N, N, method="wilson")
            lines(c(CIs[2], CIs[3]), c(ypos[i],ypos[i]), col="blue", lwd=1)
          } 
  
        if (!add && !missing(auxdata)) {
            faux <- if (ieaux) 
                auxdata
            else paste(" ", format(auxdata), sep = "")
            upedge <- par("usr")[4]
            outerText(faux, ypos[nongrp], adj = 1, cex = cex.labels)
            if (!missing(auxtitle)) {
                auxtitle <- paste(" ", auxtitle, sep = "")
                outerText(auxtitle, upedge + strheight(auxtitle, 
                  cex = cex.labels)/2, adj = 1, cex = cex.labels, 
                  setAside = faux[1])
            }
        }
        if (!add) {
            labng <- alllab[nongrp]
            bracket <- substring(labng, 1, 1) == "[" | substring(labng, 
                nchar(labng), nchar(labng)) == "]"
            yposng <- ypos[nongrp]
            s <- !bracket
            if (!is.na(any(s)) && any(s)) 
                mtextsrt(paste(labng[s], ""), 2, 0, at = yposng[s], 
                  srt = 0, adj = 1, cex = cex.labels)
            s <- bracket
            if (!is.na(any(s)) && any(s)) {
                if (.R.) 
                  text(rep(par("usr")[1], sum(s)), yposng[s], 
                    labng[s], adj = 1, cex = cex.labels, srt = 0, 
                    xpd = NA)
                else if (.SV4. && under.unix) 
                  text(rep(par("usr")[1], sum(s)), yposng[s], 
                    labng[s], adj = 1, cex = cex.labels, srt = 0)
                else {
                  xmin <- par("usr")[1] - max(nchar(labng[s])) * 
                    0.5 * cex.labels * par("1em")[1]
                  text(rep(xmin, sum(s)), yposng[s], labng[s], 
                    adj = 0, cex = cex.labels, srt = 0)
                }
            }
            if (!missing(groups)) 
                mtextsrt(paste(alllab[-nongrp], ""), 2, 0, at = ypos[-nongrp], 
                  srt = 0, adj = 1, cex = cex.group.labels, font = groupfont)
        }
    }
    else {
        if (!add && lines) 
            abline(v = ypos[nongrp], lty = lty, lwd = 1, col = lcolor)
        points(ypos, alldat, pch = pch, cex = dotsize * cex, 
            font = dotfont)

        if (!add) 
            mtextsrt(alllab[nongrp], 1, 0, at = ypos[nongrp], 
                srt = 90, adj = 1, cex = cex.labels)
        if (!add && !missing(groups)) 
            mtextsrt(alllab[-nongrp], 1, 0, at = ypos[-nongrp], 
                srt = 90, adj = 1, cex = cex.group.labels, font = groupfont)
    }
    plt <- par("plt")
    if (horizontal) {
        frac <- (oldplt[2] - oldplt[1])/(oldplt[2] - plt[1])
        umin <- tusr[2] - (tusr[2] - tusr[1]) * frac
        tusr <- c(umin, tusr[2:4])
    }
    else {
        frac <- (oldplt[4] - oldplt[3])/(oldplt[4] - plt[3])
        umin <- tusr[4] - (tusr[4] - tusr[3]) * frac
        tusr <- c(tusr[1:2], umin, tusr[4])
    }
    invisible()
}


plot.summary.formula.response.CIs = 
function (x, which = 1, vnames = c("labels", "names"), xlim, 
    xlab, pch = c(16, 1, 2, 17, 15, 3, 4, 5, 0), superposeStrata = TRUE, 
    dotfont = 1, add = FALSE, reset.par = TRUE, main, subtitles = TRUE, 
    ...) 
{
    stats <- x
    stats <- oldUnclass(stats)
    vnames <- match.arg(vnames)
    ul <- vnames == "labels"
    at <- attributes(stats)
    ns <- length(at$strat.levels)
    if (ns > 1 && length(which) > 1) 
        stop("cannot have a vector for which if > 1 strata present")
    if (ns < 2) 
        superposeStrata <- FALSE
    vn <- if (ul) 
        at$vlabel
    else at$vname
    Units <- at$vunits
    vn <- ifelse(Units == "", vn, paste(vn, " [", Units, "]", 
        sep = ""))
    vn <- vn[vn != ""]
    d <- dim(stats)
    n <- d[1]
    nstat <- d[2]/ns
    vnd <- factor(rep(vn, at$nlevels))
    dn <- dimnames(stats)
    if (missing(xlim)) 
        xlim <- range(stats[, nstat * ((1:ns) - 1) + 1 + which], 
            na.rm = TRUE)
    if (missing(main)) 
        main <- at$funlab
    nw <- length(which)
    pch <- rep(pch, length = if (superposeStrata) ns else nw)
    dotfont <- rep(dotfont, length = nw)
    opar <- if (.R.) 
        par(no.readonly = TRUE)
    else par()
    if (reset.par) 
        on.exit(par(opar))
    if (superposeStrata) 
        Ns <- apply(stats[, nstat * ((1:ns) - 1) + 1], 1, sum)
    for (is in 1:ns) {
        for (w in 1:nw) {
            js <- nstat * (is - 1) + 1 + which[w]
            z <- stats[, js]
            if (missing(xlab)) 
                xlab <- if (nw > 1) 
                  dn[[2]][js]
                else at$ylabel
            dotchart2.CIs(z, groups = vnd, xlab = xlab, xlim = xlim, 
                auxdata = if (superposeStrata) 
                  Ns
                else stats[, js - which[w]], auxtitle = "N", 
                sort = FALSE, pch = pch[if (superposeStrata) 
                  is
                else w], dotfont = dotfont[w], add = add | w > 
                  1 | (is > 1 && superposeStrata), reset.par = FALSE, 
                ...)
            lines(0.3, 0.7)
            if (ns > 1 && !superposeStrata) 
                title(paste(paste(main, if (main != "") 
                  "   "), at$strat.levels[is]))
            else if (main != "") 
                title(main)
            if (ns == 1 && subtitles) {
                title(sub = paste("N=", at$n, sep = ""), adj = 0, 
                  cex = 0.6)
                if (at$nmiss > 0) 
                  title(sub = paste("N missing=", at$nmiss, sep = ""), 
                    cex = 0.6, adj = 1)
            }
        }
    }
    if (superposeStrata) {
        Key <- if (.R.) 
            function(x = NULL, y = NULL, lev, pch) {
                oldpar <- par(usr = c(0, 1, 0, 1), xpd = NA)
                on.exit(par(oldpar))
                if (is.list(x)) {
                  y <- x$y
                  x <- x$x
                }
                if (!length(x)) 
                  x <- 0
                if (!length(y)) 
                  y <- 1
                rlegend(x, y, legend = lev, pch = pch, ...)
                invisible()
            }
        else function(x = NULL, y = NULL, lev, pch, ...) {
            if (length(x)) {
                if (is.list(x)) {
                  y <- x$y
                  x <- x$x
                }
                key(x = x, y = y, text = list(lev), points = list(pch = pch), 
                  transparent = TRUE, ...)
            }
            else key(text = list(lev), points = list(pch = pch), 
                transparent = TRUE, ...)
            invisible()
        }
        formals(Key) <- list(x = NULL, y = NULL, lev = at$strat.levels, 
            pch = pch)
        storeTemp(Key)
    }
    invisible()
}