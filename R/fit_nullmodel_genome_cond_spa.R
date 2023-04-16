#' Fitting conditional generalized linear mixed models with known relationship matrices
#' for conditional analysis in imbalanced case-control setting.
#'
#' The \code{fit_nullmodel_genome_cond_spa} function fit regression models for conditional analysis in imbalanced 
#' case-control setting, which provides the preliminary step for subsequent conditional variant-set tests in
#' conditional analysis. Each chromosome has a separate null model for conditional analysis. See \code{fit_nullmodel} for more details.
#' @param fixed an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the fixed effects model to be fitted. For multiple phenotype analysis,
#' \code{\link{formula}} recognized by \code{\link{lm}}, such as \code{cbind(y1,y2,y3) ~ x1 + x2},
#' can be used in \code{fixed} as fixed effects.
#' @param data a data frame or list (or object coercible by \code{as.data.frame} to a data frame)
#' containing the variables in the model.
#' @param kins a known positive semi-definite relationship matrix
#' (e.g. kinship matrix in genetic association studies) or a list of known
#' positive semi-definite relationship matrices. The rownames and colnames of
#' these matrices must at least include all samples as specified in the \code{id} column
#' of the data frame \code{data}. If \code{kins} is NULL, \code{fit_nullmodel}
#' will switch to the generalized linear model with no random effects.
#' @param use_sparse a logical switch of whether the provided dense \code{kins} matrix should be
#' transformed to a sparse matrix (default = TRUE).
#' @param use_SPA a logical switch determines if the null model fitting occurs in an imbalanced case-control setting (default = TRUE).
#' @param agds_dir file directory of annotated GDS (aGDS) files for all chromosomes (1-22)
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis and should
#' contain 4 columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT)  (default = NULL).
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param MAC_cutoff the cutoff of the minimum minor allele count of known variants adjusted in
#' conditional analysis (default = 20).
#' @param output_path the directory for the output files.
#' @param cond_null_model_name the file name of conditional null models (default = NULL).
#' @param phenotype_id id of samples.
#' @param phenotype outcome in regression.
#' @param kins_cutoff the cutoff value for clustering samples to make the output matrix sparse block-diagonal
#' (default = 0.022).
#' @param id a column in the data frame \code{data}, indicating the id of samples.
#' When there are duplicates in \code{id}, the data is assumed to be longitudinal with repeated measures.
#' @param random.slope an optional column indicating the random slope for time effect used
#' in a mixed effects model for longitudinal data. It must be included in the names of \code{data}.
#' There must be duplicates in \code{id} and \code{method.optim} must be "AI" (default = NULL).
#' @param groups an optional categorical variable indicating the groups used in a
#' heteroscedastic linear mixed model (allowing residual variances in different groups to be different).
#' This variable must be included in the names of \code{data}, and \code{family} must be "gaussian"
#' and \code{method.optim} must be "AI" (default = NULL).
#' @param family a description of the error distribution and link function to be used
#' in the model. This can be a character string naming a family function, a family
#' function or the result of a call to a family function. (See \code{\link{family}} for details of family functions).
#' @param method method of fitting the generalized linear mixed model. Either "REML" or "ML" (default = "REML").
#' @param method.optim optimization method of fitting the generalized linear mixed model.
#' Either "AI", "Brent" or "Nelder-Mead" (default = "AI").
#' @param maxiter a positive integer specifying the maximum number of iterations when
#' fitting the generalized linear mixed model (default = 500).
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates below which iterations should be stopped (default = 1e-5).
#' @param taumin the lower bound of search space for the variance component parameter \eqn{\tau} (default = 1e-5),
#' used when \code{method.optim} = "Brent". See Details.
#' @param taumax the upper bound of search space for the variance component parameter \eqn{\tau} (default = 1e5),
#' used when \code{method.optim} = "Brent". See Details.
#' @param tauregion the number of search intervals for the REML or ML estimate of the variance component
#' parameter \eqn{\tau} (default = 10), used when \code{method.optim} = "Brent". See Details.
#' @param verbose a logical switch for printing detailed information (parameter estimates in each iteration)
#' for testing and debugging purpose (default = FALSE).
#' @param ... additional arguments that could be passed to \code{\link{glm}}.
#' @return The function returns objects of the null models fit from \code{\link{fit_nullmodel}} 
#' and whether the \code{kins} matrix is sparse when fitting the null model, each chromosome has one output. 
#' See \code{\link{fit_nullmodel}} for more details.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Chen, H., et al. (2019). Efficient variant set mixed model association tests for continuous and
#' binary traits in large-scale whole-genome sequencing studies. \emph{The American Journal of Human Genetics}, \emph{104}(2), 260-274.
#' (\href{https://doi.org/10.1016/j.ajhg.2018.12.012}{pub})
#' @references Chen, H. (2021). GMMAT: Generalized linear Mixed Model Association Tests Version 1.3.2.
#' (\href{https://cloud.r-project.org/web/packages/GMMAT/vignettes/GMMAT.pdf}{web})
#' @export

fit_nullmodel_genome_cond_spa <- function(fixed, data = parent.frame(), kins, use_sparse = TRUE, use_SPA = TRUE,
								   agds_dir,known_loci,geno_missing_imputation=c("mean","minor"),MAC_cutoff=20,
								   output_path,cond_null_model_name = NULL, phenotype_id, phenotype,
								   kins_cutoff = 0.022, id, random.slope = NULL, groups = NULL,
								   family = binomial(link = "logit"), method = "REML",
								   method.optim = "AI", maxiter = 500, tol = 1e-5,
								   taumin = 1e-5, taumax = 1e5, tauregion = 10,
								   verbose = FALSE, ...){

	phenotype.id <- phenotype_id

	for(chr in 1:22)
	{
		if(sum(known_loci[,1]==chr)>0)
		{
			## get known loci genotype
			known_loci_chr <- known_loci[known_loci[,1]==chr,,drop=FALSE]

			gds.path <- agds_dir[chr]
			genofile <- seqOpen(gds.path)

			## Genotype of Adjusted Variants
			position <- as.numeric(seqGetData(genofile, "position"))
			REF <- as.character(seqGetData(genofile, "$ref"))
			ALT <- as.character(seqGetData(genofile, "$alt"))
			variant.id <- seqGetData(genofile, "variant.id")

			rs_num_in <- c()
			for(i in 1:dim(known_loci_chr)[1])
			{
				rs_num_in <- c(rs_num_in,which((position==known_loci_chr[i,2])&(REF==known_loci_chr[i,3])&(ALT==known_loci_chr[i,4])))
			}

			variant.id.in <- variant.id[rs_num_in]
			seqSetFilter(genofile,variant.id=variant.id.in,sample.id=phenotype.id)

			## genotype id
			id.genotype <- seqGetData(genofile,"sample.id")
			id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
			phenotype.id.merge <- data.frame(phenotype.id)
			phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
			id.genotype.match <- phenotype.id.merge$index

			Geno_adjusted <- seqGetData(genofile, "$dosage")
			Geno_adjusted <- Geno_adjusted[id.genotype.match,,drop=FALSE]

			## impute missing
			if(!is.null(dim(Geno_adjusted)))
			{
				if(dim(Geno_adjusted)[2]>0)
				{
					if(geno_missing_imputation=="mean")
					{
						Geno_adjusted <- matrix_flip_mean(Geno_adjusted)$Geno
					}
					if(geno_missing_imputation=="minor")
					{
						Geno_adjusted <- matrix_flip_minor(Geno_adjusted)$Geno
					}
				}
			}

			if(class(Geno_adjusted)[1]=="numeric")
			{
				Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
			}

			AF <- apply(Geno_adjusted,2,mean)/2
			MAF <- AF*(AF<0.5) + (1-AF)*(AF>=0.5)

			MAF_cutoff <- MAC_cutoff/dim(Geno_adjusted)[1]

			case_number_carrier <- crossprod(Geno_adjusted,phenotype)

			if(sum((MAF>MAF_cutoff)&(case_number_carrier>0))>0)
			{
				print(chr)

				Geno_adjusted <- Geno_adjusted[,(MAF>MAF_cutoff)&(case_number_carrier>0),drop=FALSE]

				if(class(Geno_adjusted)[1]=="numeric")
				{
					Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
				}

				colnames(Geno_adjusted)=c(paste("Geno_adjusted",1:(ncol(Geno_adjusted)),sep=""))
				data_cond <- cbind(data,Geno_adjusted)

				nullmodel_formula_cond <- fixed
				for(kk in 1:ncol(Geno_adjusted))
				{
					nullmodel_formula_cond <- paste0(nullmodel_formula_cond,"+Geno_adjusted",kk)
				}

				### fit null model
				obj.nullmodel.cond <- fit_nullmodel(nullmodel_formula_cond, data = data_cond, kins = kins, id = id, use_sparse = use_sparse,  use_SPA = use_SPA, family = binomial(link = "logit"), verbose=T)

				save(obj.nullmodel.cond, file = paste0(output_path,cond_null_model_name,".chr",chr,".Rdata"))

				seqClose(genofile)
			}
		}
	}
}

