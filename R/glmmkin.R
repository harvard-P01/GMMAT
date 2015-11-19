glmmkin <- function(fixed, data = parent.frame(), kins, family = binomial(link = "logit"), method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE, ...) {
	call <- match.call()
	if(!class(kins) %in% c("matrix", "list"))
		stop("Error: \"kins\" must be a matrix or a list.")
	if(class(kins) == "list" && length(table(sapply(kins, dim))) != 1)
		stop("Error: when \"kins\" is a list, all its elements must be square matrices of the same size.")
	if(!method %in% c("REML", "ML"))
		stop("Error: \"method\" must be \"REML\" or \"ML\".")
	method.optim <- try(match.arg(method.optim, c("AI", "Brent", "Nelder-Mead")))
	if(class(method.optim) == "try-error")
		stop("Error: \"method.optim\" must be \"AI\", \"Brent\" or \"Nelder-Mead\".")
	if(method.optim == "AI" && method == "ML")
		stop("Error: method \"ML\" not available for method.optim \"AI\", use method \"REML\" instead.")
	if(method.optim == "Brent" && class(kins) == "list")
		stop("Error: method.optim \"Brent\" can only be applied in one-dimensional optimization, use a matrix for \"kins\".")
	if(class(family) != "family")
		stop("Error: \"family\" must be an object of class \"family\".")
	if(!family$family %in% c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"))
		stop("Error: family must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
	if(method.optim != "Brent" && class(kins) == "matrix") kins <- list(kins1 = kins)
	fit0 <- glm(formula = fixed, data = data, family = family, ...)
	idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
	if(class(kins) == "matrix") kins <- kins[idx, idx]
	else {
	        for(i in 1:length(kins)) kins[[i]] <- kins[[i]][idx, idx]
	}
	fit <- glmmkin.fit(fit0, kins, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
	fit$call <- call
	class(fit) <- "glmmkin"
	return(fit)
}

glmmkin.fit <- function(fit0, kins, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE) {
	if(method.optim == "Brent") {
		fit <- glmmkin.brent(fit0, kins, method = method, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
		if(fit$theta[2]/fit$theta[1] < 1.01 * tol) {
			warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
			fit <- glmmkin.brent(fit0, kins, method = method, tau = 0, fixtau = 1, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
		}
	} else {
		names(kins) <- paste("kins", 1:length(kins), sep="")
		if(method.optim == "AI") {
			fixtau.old <- rep(0, length(kins)+1)
			fit <- glmmkin.ai(fit0, kins, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau.new <- 1*(fit$theta < 1.01 * tol)
			while(any(fixtau.new != fixtau.old)) {
				warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
				fixtau.old <- fixtau.new
				fit <- glmmkin.ai(fit0, kins, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
				fixtau.new <- 1*(fit$theta < 1.01 * tol)
			}
		} else {
			fixtau.old <- rep(0, length(kins))
			fit <- glmmkin.nm(fit0, kins, method = method, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
			while(any(fixtau.new != fixtau.old)) {
				warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
				fixtau.old <- fixtau.new
				tau <- rep(1, length(kins))
				tau[which(fixtau.old == 1)] <- 0
				fit <- glmmkin.nm(fit0, kins, method = method, tau = tau, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
				fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
			}
		}
	}
	return(fit)
}

glmmkin.ai <- function(fit0, kins, tau = rep(0, length(kins)+1), fixtau = rep(0, length(kins)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {
	y <- fit0$y
	n <- length(y)
	offset <- fit0$offset
	if(is.null(offset)) offset <- rep(0, n)
	family <- fit0$family
	eta <- fit0$linear.predictors
	mu <- fit0$fitted.values
	mu.eta <- family$mu.eta(eta)
	Y <- eta - offset + (y - mu)/mu.eta
	sqrtW <- mu.eta/sqrt(fit0$family$variance(mu))
	X <- model.matrix(fit0)
	alpha <- fit0$coef
	if(verbose) {
		cat("Fixed-effect coefficients:\n")
		print(alpha)
	}
	if(family$family %in% c("poisson", "binomial")) {
		tau[1] <- 1
		fixtau[1] <- 1
	}
	q <- length(kins)
	idxtau <- which(fixtau == 0)
	q2 <- sum(fixtau == 0)
	if(q2 > 0) {
	        tau[fixtau == 0] <- rep(var(Y)/(q+1), q2)
		Sigma <- tau[1]*diag(1/sqrtW^2)
		for(i in 1:q) Sigma <- Sigma + tau[i+1]*kins[[i]]
		Sigma_i <- chol2inv(chol(Sigma))
		Sigma_iX <- crossprod(Sigma_i, X)
		P <- Sigma_i - tcrossprod(tcrossprod(Sigma_iX, chol2inv(chol(crossprod(X, Sigma_iX)))), Sigma_iX)
		PY <- crossprod(P, Y)
		tau0 <- tau
		for(i in 1:q2) {
		        if(i == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/sqrtW)^2) - sum(diag(P)/sqrtW^2))/n)
			else {
	        	        PAPY <- crossprod(P, crossprod(kins[[idxtau[i]-1]], PY))
				tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (crossprod(Y, PAPY) - sum(P*kins[[idxtau[i]-1]]))/n)
			}
		}
	}
	for (i in seq_len(maxiter)) {
		if(verbose) cat("\nIteration ", i, ":\n")
		alpha0 <- alpha
		tau0 <- tau
		fit <- .Call("fitglmm_ai", Y, X, length(kins), kins, sqrtW^2, tau, fixtau, maxiter, tol)
		tau <- as.numeric(fit$tau)
		cov <- as.matrix(fit$cov)
		alpha <- as.numeric(fit$alpha)
		eta <- as.numeric(fit$eta) + offset
		if(verbose) {
			cat("Variance component estimates:\n")
			print(tau)
			cat("Fixed-effect coefficients:\n")
			print(alpha)
		}
		mu <- family$linkinv(eta)
		mu.eta <- family$mu.eta(eta)
		Y <- eta - offset + (y - mu)/mu.eta
		sqrtW <- mu.eta/sqrt(family$variance(mu))
		if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
	}
	converged <- ifelse(i < maxiter, TRUE, FALSE)
	res <- y - mu
	P <- fit$P
	return(list(theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, P=P, residuals=res, cov=cov, converged=converged))
}

glmmkin.brent <- function(fit0, kins, method = "REML", tau = 1, fixtau = 0, maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE) {
	y <- fit0$y
	n <- length(y)
	offset <- fit0$offset
	if(is.null(offset)) offset <- rep(0, n)
	family <- fit0$family
	eta <- fit0$linear.predictors
	mu <- fit0$fitted.values
	mu.eta <- family$mu.eta(eta)
	Y <- eta - offset + (y - mu)/mu.eta
	sqrtW <- mu.eta/sqrt(fit0$family$variance(mu))
	X <- model.matrix(fit0)
	alpha <- fit0$coef
	if(verbose) {
		cat("Fixed-effect coefficients:\n")
		print(alpha)
	}
	dispersion <- ifelse(family$family %in% c("poisson", "binomial"), "N", "Y")
	method <- ifelse(method=="REML", "R", "L")
	for (i in seq_len(maxiter)) {
		if(verbose) cat("\nIteration ", i, ":\n")
		alpha0 <- alpha
		tau0 <- tau
		fit <- .Call("fitglmm_brent", Y, X, kins, sqrtW, method, dispersion, tau, fixtau, maxiter, tol, taumin, taumax, tauregion)
		tau <- as.numeric(fit$tau)
		cov <- as.matrix(fit$cov)
		alpha <- as.numeric(fit$alpha)
		eta <- as.numeric(fit$eta) + offset
		if(verbose) {
			cat("Variance component estimates:\n")
			print(tau)
			cat("Fixed-effect coefficients:\n")
			print(alpha)
		}
		mu <- family$linkinv(eta)
		mu.eta <- family$mu.eta(eta)
		Y <- eta - offset + (y - mu)/mu.eta
		sqrtW <- mu.eta/sqrt(family$variance(mu))
		if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
		if(family$family == "gaussian") break
	}
	converged <- ifelse(i < maxiter, TRUE, FALSE)
	res <- y - mu
	fit$eval <- as.numeric(fit$eval)
	P <- fit$U %*% (t(fit$U) * fit$eval) - fit$U %*% (fit$UtX * fit$eval) %*% (cov %*% t(fit$UtX) %*% (t(fit$U) * fit$eval))
	if(dispersion=="N") {
		theta <- c(1, tau)
	} else {
		phi <- ifelse(method == "R", as.numeric(t(Y) %*% P %*% Y)/(n - ncol(X)), as.numeric(t(Y) %*% P %*% Y)/n)
		theta <- phi * c(1, tau)
		P <- P/phi
		cov <- phi*cov
	}
	return(list(theta=theta, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, P=P, residuals=res, cov=cov, converged=converged))
}

glmmkin.nm <- function(fit0, kins, method = "REML", tau = rep(1, length(kins)), fixtau = rep(0, length(kins)), maxiter = 500, tol = 1e-5, verbose = FALSE) {
	y <- fit0$y
	n <- length(y)
	offset <- fit0$offset
	if(is.null(offset)) offset <- rep(0, n)
	family <- fit0$family
	eta <- fit0$linear.predictors
	mu <- fit0$fitted.values
	mu.eta <- family$mu.eta(eta)
	Y <- eta - offset + (y - mu)/mu.eta
	sqrtW <- mu.eta/sqrt(fit0$family$variance(mu))
	X <- model.matrix(fit0)
	alpha <- fit0$coef
	if(verbose) {
		cat("Fixed-effect coefficients:\n")
		print(alpha)
	}
	dispersion <- ifelse(family$family %in% c("poisson", "binomial"), "N", "Y")
	method <- ifelse(method=="REML", "R", "L")
	for (i in seq_len(maxiter)) {
		if(verbose) cat("\nIteration ", i, ":\n")
		alpha0 <- alpha
		tau0 <- tau
		fit <- .Call("fitglmm_nm", Y, X, length(kins), kins, sqrtW^2, method, dispersion, tau, fixtau, maxiter, tol)
		tau <- as.numeric(fit$tau)
		cov <- as.matrix(fit$cov)
		alpha <- as.numeric(fit$alpha)
		eta <- as.numeric(fit$eta) + offset
		if(verbose) {
			cat("Variance component estimates:\n")
			print(tau)
			cat("Fixed-effect coefficients:\n")
			print(alpha)
		}
		mu <- family$linkinv(eta)
		mu.eta <- family$mu.eta(eta)
		Y <- eta - offset + (y - mu)/mu.eta
		sqrtW <- mu.eta/sqrt(family$variance(mu))
		if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
		if(family$family == "gaussian") break
	}
	converged <- ifelse(i < maxiter, TRUE, FALSE)
	res <- y - mu
	fit$eval <- as.numeric(fit$eval)
	P <- fit$U %*% (t(fit$U) * fit$eval) - fit$U %*% (fit$UtX * fit$eval) %*% (cov %*% t(fit$UtX) %*% (t(fit$U) * fit$eval))
	if(dispersion=="N") {
		theta <- c(1, tau)
	} else {
		phi <- ifelse(method == "R", as.numeric(t(Y) %*% P %*% Y)/(n - ncol(X)), as.numeric(t(Y) %*% P %*% Y)/n)
		theta <- phi * c(1, tau)
		P <- P/phi
		cov <- phi*cov
	}
	return(list(theta=theta, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, P=P, residuals=res, cov=cov, converged=converged))
}

# res = res/phi
glmm.score.text <- function(res, P, infile, outfile, tol = 1e-5, center = T, select = NULL, missing.method = "impute2mean", infile.nrow = NULL, infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", infile.ncol.skip = 1, infile.ncol.print = 1, infile.header.print = "SNP", nperbatch = 100) {
	if(is.null(infile.nrow)) {
		if(grepl("\\.gz$", infile)) infile.nrow <- as.integer(system(paste("zcat", infile, "| wc -l | awk '{print $1}'"), intern = T))
		else if(grepl("\\.bz2$", infile)) infile.nrow <- as.integer(system(paste("bzcat", infile, "| wc -l | awk '{print $1}'"), intern = T))
		else infile.nrow <- as.integer(system(paste("wc -l", infile, "| awk '{print $1}'"), intern = T))
	}
	if(!is.numeric(infile.nrow) | infile.nrow < 0)
		stop("Error: number of rows of the input file is incorrect!")
	if(!is.numeric(infile.nrow.skip) | infile.nrow.skip < 0)
	        stop("Error: number of skipped rows of the input file is incorrect!")
	if(!is.numeric(infile.ncol.skip) | infile.ncol.skip < 0)
		stop("Error: number of skipped cols of the input file is incorrect!")
	if(length(infile.ncol.print) != length(infile.header.print))
		stop("Error: number of cols selected to print does not match number of header names!")
	if(is.null(infile.ncol.print))
		infile.ncol.print <- 0
	if(is.null(infile.header.print))
		infile.header.print <- infile.na
	if(any(!is.numeric(infile.ncol.print)) | any(infile.ncol.print < 0) | any(infile.ncol.print > infile.ncol.skip))
		stop("Error: cols selected to print have incorrect indices!")
	if(any(infile.ncol.print != sort(infile.ncol.print)))
		stop("Error: col indices must be sorted increasingly in infile.ncol.print!")
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)
	if(is.null(select)) select <- 1:length(res)
	select2 <- select[select > 0]
	if(length(select2) != length(res) | any(sort(select2) != 1:length(res))) stop("Error: select is a vector of orders, individuals not in res should be coded 0!")
	if(center) {
		time <- .Call("glmm_score_text", res, P, infile, outfile, tol, 'c', miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, infile.header.print, nperbatch, select)
	} else {
		time <- .Call("glmm_score_text", res, P, infile, outfile, tol, 'n', miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, infile.header.print, nperbatch, select)
	}
	print(sprintf("Computational time: %.2f seconds", time))
	invisible(time)
}

glmm.wald.text <- function(fixed, data = parent.frame(), kins, family = binomial(link = "logit"), infile, snps, snp.col = 1, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, center = T, select = NULL, missing.method = "impute2mean", infile.nrow = NULL, infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", infile.ncol.skip = 1, infile.ncol.print = 1, infile.header.print = "SNP", verbose = FALSE, ...) {
	if(!class(kins) %in% c("matrix", "list"))
		stop("Error: \"kins\" must be a matrix or a list.")
	if(class(kins) == "list" && length(table(sapply(kins, dim))) != 1)
		stop("Error: when \"kins\" is a list, all its elements must be square matrices of the same size.")
	if(!method %in% c("REML", "ML"))
		stop("Error: \"method\" must be \"REML\" or \"ML\".")
	method.optim <- try(match.arg(method.optim, c("AI", "Brent", "Nelder-Mead")))
	if(class(method.optim) == "try-error")
		stop("Error: \"method.optim\" must be \"AI\", \"Brent\" or \"Nelder-Mead\".")
	if(method.optim == "AI" && method == "ML")
		stop("Error: method \"ML\" not available for method.optim \"AI\", use method \"REML\" instead.")
	if(method.optim == "Brent" && class(kins) == "list")
		stop("Error: method.optim \"Brent\" can only be applied in one-dimensional optimization, use a matrix for \"kins\".")
	if(class(family) != "family")
		stop("Error: \"family\" must be an object of class \"family\".")
	if(!family$family %in% c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"))
		stop("Error: family must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
        if(method.optim != "Brent" && class(kins) == "matrix") kins <- list(kins1 = kins)
	if(is.null(infile.nrow)) {
		if(grepl("\\.gz$", infile)) infile.nrow <- as.integer(system(paste("zcat", infile, "| wc -l | awk '{print $1}'"), intern = T))
		else if(grepl("\\.bz2$", infile)) infile.nrow <- as.integer(system(paste("bzcat", infile, "| wc -l | awk '{print $1}'"), intern = T))
		else infile.nrow <- as.integer(system(paste("wc -l", infile, "| awk '{print $1}'"), intern = T))
	}
	if(!is.numeric(infile.nrow) | infile.nrow < 0)
		stop("Error: number of rows of the input file is incorrect!")
	if(!is.numeric(infile.nrow.skip) | infile.nrow.skip < 0)
	        stop("Error: number of skipped rows of the input file is incorrect!")
	if(!is.numeric(infile.ncol.skip) | infile.ncol.skip <= 0)
		stop("Error: number of skipped cols of the input file is incorrect!")
	if(length(infile.ncol.print) != length(infile.header.print))
		stop("Error: number of cols selected to print does not match number of header names!")
	if(!is.numeric(infile.ncol.print) | (!snp.col %in% infile.ncol.print))
		stop("Error: snp.col is not in infile.ncol.print!")
	if(any(!is.numeric(infile.ncol.print)) | any(infile.ncol.print < 0) | any(infile.ncol.print > infile.ncol.skip))
		stop("Error: cols selected to print have incorrect indices!")
	if(any(infile.ncol.print != sort(infile.ncol.print)))
		stop("Error: col indices must be sorted increasingly in infile.ncol.print!")
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)
	nn <- ifelse(class(kins) == "matrix", nrow(kins), nrow(kins[[1]]))
	if(is.null(select)) select <- 1:nn
	select2 <- select[select > 0]
	if(length(select2) != nn | any(sort(select2) != 1:nn)) stop("Error: select is a vector of orders, individuals not in kins should be coded 0!")
	idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
	nn <- length(idx)
	data <- data[idx, ]
	select2 <- match(idx, select)
	select[-select2] <- 0
	select[select2] <- 1:nn
	if(class(kins) == "matrix") kins <- kins[idx, idx]
	else {
	        for(i in 1:length(kins)) kins[[i]] <- kins[[i]][idx, idx]
	}
	snpinfo <- matrix(NA, length(snps), length(infile.header.print))
	N <- AF <- BETA <- SE <- PVAL <- converged <- rep(NA, length(snps))
	for(ii in 1:length(snps)) {
	        snp <- snps[ii]
		if(verbose) cat("\nAnalyze SNP ", ii, ": ", snp, "\n")
		if(center) {
			readfile <- .Call("glmm_wald_text", nn, snp, infile, tol, 'c', miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, snp.col, select)
		} else {
			readfile <- .Call("glmm_wald_text", nn, snp, infile, tol, 'n', miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, snp.col, select)
		}
		if(readfile$skip == 2) { # snp not found in the file
			snpinfo[ii, which(infile.ncol.print == snp.col)] <- snp
		} else {
			snpinfo[ii, ] <- readfile$snpinfo
			N[ii] <- readfile$N
			AF[ii] <- readfile$AF
			if(readfile$skip != 1) { # snp
				data$SNP__ <- as.numeric(readfile$G)
				data$SNP__[data$SNP__ < (-999)] <- NA
				fit0 <- glm(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, family = family, ...)
				idx <- match(rownames(model.frame(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, na.action = na.omit)), rownames(model.frame(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, na.action = na.pass)))
				tmpkins <- kins
				if(class(tmpkins) == "matrix") tmpkins <- tmpkins[idx, idx]
				else {
				        for(i in 1:length(tmpkins)) tmpkins[[i]] <- tmpkins[[i]][idx, idx]
				}
				fit <- glmmkin.fit(fit0, tmpkins, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
				BETA[ii] <- fit$coefficients[length(fit$coefficients)]
				SE[ii] <- sqrt(diag(fit$cov)[length(fit$coefficients)])
				PVAL[ii] <- pchisq((BETA[ii]/SE[ii])^2, 1, lower.tail=F)
				converged[ii] <- fit$converged
			}
		}
	}
	res <- data.frame(snpinfo, N, AF, BETA, SE, PVAL, converged)
	names(res)[1:length(infile.header.print)] <- infile.header.print
	return(res)
}

glmm.score.bimbam <- function(res, P, infile, outfile, tol = 1e-5, center = T, select = NULL, missing.method = "impute2mean", infile.nrow = NULL, infile.na = "NA", nperbatch = 100) {
	glmm.score.text(res, P, infile, outfile, tol, center, select, missing.method, infile.nrow, infile.nrow.skip=0, infile.sep=",", infile.na, infile.ncol.skip=3, infile.ncol.print=1:3, infile.header.print=c("SNP", "ALLELE1", "ALLELE2"), nperbatch=nperbatch)
}

glmm.wald.bimbam <- function(fixed, data = parent.frame(), kins, family = binomial(link = "logit"), infile, snps, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, center = T, select = NULL, missing.method = "impute2mean", infile.nrow = NULL, infile.na = "NA", verbose = FALSE, ...) {
	glmm.wald.text(fixed, data, kins, family, infile, snps, snp.col = 1, method, method.optim, maxiter, tol, taumin, taumax, tauregion, center, select, missing.method, infile.nrow, infile.nrow.skip = 0, infile.sep = ",", infile.na, infile.ncol.skip = 3, infile.ncol.print = 1:3, infile.header.print = c("SNP", "ALLELE1", "ALLELE2"), verbose = verbose, ...)
}

glmm.score.bed <- function(res, P, plinkfiles, outfile, center = T, select = NULL, missing.method = "impute2mean", nperbatch = 100) {
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)
	if(is.null(select)) select <- 1:length(res)
	select2 <- select[select > 0]
	if(length(select2) != length(res) | any(sort(select2) != 1:length(res))) stop("Error: select is a vector of orders, individuals not in res should be coded 0!")
	bimfile <- paste(plinkfiles, "bim", sep=".")
	bedfile <- paste(plinkfiles, "bed", sep=".")
	famfile <- paste(plinkfiles, "fam", sep=".")
	if(length(select) != as.integer(system(paste("wc -l", famfile, "| awk '{print $1}'"), intern = T))) stop("Error: number of individuals in plink fam file incorrect!")
	if(center) {
		time <- .Call("glmm_score_bed", res, P, bimfile, bedfile, outfile, 'c', miss.method, nperbatch, select)
	} else {
		time <- .Call("glmm_score_bed", res, P, bimfile, bedfile, outfile, 'n', miss.method, nperbatch, select)
	}
	print(sprintf("Computational time: %.2f seconds", time))
	invisible(time)
}

glmm.wald.bed <- function(fixed, data = parent.frame(), kins, family = binomial(link = "logit"), plinkfiles, snps, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, center = T, select = NULL, missing.method = "impute2mean", verbose = FALSE, ...) {
	if(!class(kins) %in% c("matrix", "list"))
		stop("Error: \"kins\" must be a matrix or a list.")
	if(class(kins) == "list" && length(table(sapply(kins, dim))) != 1)
		stop("Error: when \"kins\" is a list, all its elements must be square matrices of the same size.")
	if(!method %in% c("REML", "ML"))
		stop("Error: \"method\" must be \"REML\" or \"ML\".")
	method.optim <- try(match.arg(method.optim, c("AI", "Brent", "Nelder-Mead")))
	if(class(method.optim) == "try-error")
		stop("Error: \"method.optim\" must be \"AI\", \"Brent\" or \"Nelder-Mead\".")
	if(method.optim == "AI" && method == "ML")
		stop("Error: method \"ML\" not available for method.optim \"AI\", use method \"REML\" instead.")
	if(method.optim == "Brent" && class(kins) == "list")
		stop("Error: method.optim \"Brent\" can only be applied in one-dimensional optimization, use a matrix for \"kins\".")
	if(class(family) != "family")
		stop("Error: \"family\" must be an object of class \"family\".")
	if(!family$family %in% c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"))
		stop("Error: family must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
        if(method.optim != "Brent" && class(kins) == "matrix") kins <- list(kins1 = kins)
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)
	nn <- ifelse(class(kins) == "matrix", nrow(kins), nrow(kins[[1]]))
	if(is.null(select)) select <- 1:nn
	select2 <- select[select > 0]
	if(length(select2) != nn | any(sort(select2) != 1:nn)) stop("Error: select is a vector of orders, individuals not in kins should be coded 0!")
	idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
	nn <- length(idx)
	data <- data[idx, ]
	select2 <- match(idx, select)
	select[-select2] <- 0
	select[select2] <- 1:nn
	if(class(kins) == "matrix") kins <- kins[idx, idx]
	else {
	        for(i in 1:length(kins)) kins[[i]] <- kins[[i]][idx, idx]
	}
	bimfile <- paste(plinkfiles, "bim", sep=".")
	bedfile <- paste(plinkfiles, "bed", sep=".")
	famfile <- paste(plinkfiles, "fam", sep=".")
	if(length(select) != as.integer(system(paste("wc -l", famfile, "| awk '{print $1}'"), intern = T))) stop("Error: number of individuals in plink fam file incorrect!")
	snpinfo <- matrix(NA, length(snps), 6)
	N <- AF <- BETA <- SE <- PVAL <- converged <- rep(NA, length(snps))
	for(ii in 1:length(snps)) {
	        snp <- snps[ii]
		if(verbose) cat("\nAnalyze SNP ", ii, ": ", snp, "\n")
		if(center) {
			readfile <- .Call("glmm_wald_bed", nn, snp, bimfile, bedfile, 'c', miss.method, select)
		} else {
			readfile <- .Call("glmm_wald_bed", nn, snp, bimfile, bedfile, 'n', miss.method, select)
		}
		if(readfile$skip == 2) { # snp not found in the file
			snpinfo[ii, 2] <- snp
		} else {
			snpinfo[ii, ] <- readfile$snpinfo
			N[ii] <- readfile$N
			AF[ii] <- readfile$AF
			if(readfile$skip != 1) { # snp
				data$SNP__ <- as.numeric(readfile$G)
				data$SNP__[data$SNP__ < (-999)] <- NA
				fit0 <- glm(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, family = family, ...)
				idx <- match(rownames(model.frame(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, na.action = na.omit)), rownames(model.frame(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, na.action = na.pass)))
				tmpkins <- kins
				if(class(tmpkins) == "matrix") tmpkins <- tmpkins[idx, idx]
				else {
				        for(i in 1:length(tmpkins)) tmpkins[[i]] <- tmpkins[[i]][idx, idx]
				}
				fit <- glmmkin.fit(fit0, tmpkins, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
				BETA[ii] <- fit$coefficients[length(fit$coefficients)]
				SE[ii] <- sqrt(diag(fit$cov)[length(fit$coefficients)])
				PVAL[ii] <- pchisq((BETA[ii]/SE[ii])^2, 1, lower.tail=F)
				converged[ii] <- fit$converged
			}
		}
	}
	res <- data.frame(snpinfo, N, AF, BETA, SE, PVAL, converged)
	names(res)[1:6] <- c("CHR", "SNP", "cM", "POS", "A1", "A2")
	return(res)
}
