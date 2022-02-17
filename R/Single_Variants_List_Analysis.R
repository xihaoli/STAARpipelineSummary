#' Calculate individual-variant p-values of a list of variants
#'
#' The \code{Single_Variants_List_Analysis} function takes in a list of variants to calculate the p-values and effect sizes of the input variants
#' @param agds_dir file directory of annotated GDS (aGDS) files for all chromosomes (1-22).
#' @param single_variants_list name a data frame containing the information of variants to be functionally annotated. The data frame must include 4 columns with
#' the following names: "CHR" (chromosome number), "POS" (position), "REF" (reference allele), and "ALT" (alternative allele).
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{fit_nullmodel} function in the \code{STAARpipeline} package,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{genesis2staar_nullmodel} function in the \code{STAARpipeline} package.
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @return a data frame containing the basic information (chromosome, position, reference allele and alternative allele)
#' the score test p-values, and the effect sizes for the input variants.
#' @export

Single_Variants_List_Analysis <- function(agds_dir,single_variants_list,obj_nullmodel,
                                          QC_label="annotation/filter",geno_missing_imputation=c("mean","minor")){

	## evaluate choices
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)

	Sigma_i <- obj_nullmodel$Sigma_i
	Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
	cov <- obj_nullmodel$cov

	residuals.phenotype <- obj_nullmodel$scaled.residuals

	single_variants_list_info <- single_variants_list[,c("CHR","POS","REF","ALT")]

	single_variants_list_annotation <- c()
	for(chr in 1:22)
	{
		print(chr)
		if(sum(single_variants_list_info$CHR==chr)>0)
		{
			single_variants_list_info_chr <- single_variants_list_info[single_variants_list_info$CHR==chr,,drop=FALSE]

			gds.path <- agds_dir[chr]
			genofile <- seqOpen(gds.path)

			position <- as.numeric(seqGetData(genofile, "position"))
			REF <- as.character(seqGetData(genofile, "$ref"))
			ALT <- as.character(seqGetData(genofile, "$alt"))
			variant_id <- seqGetData(genofile, "variant.id")

			chr_info <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT,variant_id=variant_id)

			single_variants_list_info_chr <- dplyr::left_join(single_variants_list_info_chr,chr_info,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
			variant.id.in <- single_variants_list_info_chr$variant_id[!is.na(single_variants_list_info_chr$variant_id)]

			seqSetFilter(genofile,variant.id=variant.id.in)

			## genotype id
			id.genotype <- seqGetData(genofile,"sample.id")

			id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
			phenotype.id.merge <- data.frame(phenotype.id)
			phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
			id.genotype.match <- phenotype.id.merge$index


			Geno <- seqGetData(genofile, "$dosage")
			Geno <- Geno[id.genotype.match,,drop=FALSE]

			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)
			}

			MAF <- Geno$MAF
			ALT_AF <- 1 - Geno$AF

			CHR <- as.numeric(seqGetData(genofile, "chromosome"))
			position <- as.numeric(seqGetData(genofile, "position"))
			REF <- as.character(seqGetData(genofile, "$ref"))
			ALT <- as.character(seqGetData(genofile, "$alt"))
			filter <- seqGetData(genofile, QC_label)

			N <- rep(samplesize,length(CHR))

			Score_test <- Individual_Score_Test(Geno$Geno, Sigma_i, Sigma_iX, cov, residuals.phenotype)

			single_variants_list_annotation_chr <-data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,ALT_AF=ALT_AF,QC_label=filter,MAF=MAF,N=N,
			                                                 pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
			                                                 Score=Score_test$Score,Score_se=Score_test$Score_se,
			                                                 Est=Score_test$Est,Est_se=Score_test$Est_se)

			single_variants_list_annotation <- rbind(single_variants_list_annotation,single_variants_list_annotation_chr)
			seqClose(genofile)
		}

	}
	single_variants_list_info_annotation <- dplyr::left_join(single_variants_list,single_variants_list_annotation,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))

	return(single_variants_list_info_annotation)
}

