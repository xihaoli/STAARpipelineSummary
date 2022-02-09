#' Functionally annotate rare variants in a coding mask
#'
#' The \code{Gene_Centric_Coding_Info} function takes in a coding mask of a gene to functionally annotate the rare variants in the mask.
#' @param category the coding functional category of rare variants to be functionally annotated. Choices include
#' \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, \code{synonymous} (default = \code{plof}).
#' @param chr chromosome.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{fit_nullmodel} function in the \code{STAARpipeline} package,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{genesis2staar_nullmodel} function in the \code{STAARpipeline} package.
#' @param gene_name name of the gene to be annotated.
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis and should
#' contain 4 columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT) (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{known_loci}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{known_loci}
#' and taking the residuals (default = \code{optimal}).
#' @param QC_label channel name of the QC label in the GDS/aGDS file.
#' @param variant_type type of variant included in the conditional analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @return a data frame containing the basic information (chromosome, position, reference allele and alternative allele),
#' unconditional and conditional the score test p-values,
#' and annotation scores for the rare variants of the input coding mask.
#' @export

Gene_Centric_Coding_Info <- function(category=c("plof","plof_ds","missense","disruptive_missense","synonymous"),
                                     chr,genofile,obj_nullmodel,gene_name,known_loci=NULL,rare_maf_cutoff=0.01,
                                     method_cond=c("optimal","naive"),
                                     QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                     Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,Annotation_name=NULL){

	## evaluate choices
	category <- match.arg(category)
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	if(is.null(known_loci))
	{
		known_loci <- data.frame(chr=logical(0),pos=logical(0),ref=character(0),alt=character(0))
	}

	if(category=="plof")
	{
		results <- info_plof(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
		                     method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)
	}

	if(category=="plof_ds")
	{
		results <- info_plof_ds(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
		                        method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)
	}

	if(category=="missense")
	{
		results <- info_missense(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
		                         method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)
	}

	if(category=="disruptive_missense")
	{
		results <- info_disruptive_missense(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
		                                    method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)
	}

	if(category=="synonymous")
	{
		results <- info_synonymous(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
		                           method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)
	}

	return(results)
}

