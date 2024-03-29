#' Functionally annotate rare variants in a genetic region
#'
#' The \code{Sliding_Window_Info} function takes in the location of a genetic region to functionally annotate the rare variants in the region.
#' @param chr chromosome.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{fit_nullmodel} function in the \code{STAARpipeline} package,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{genesis2staar_nullmodel} function in the \code{STAARpipeline} package.
#' @param start_loc starting location (position) of the genetic region to be annotated.
#' @param end_loc ending location (position) of the genetic region to be annotated.
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
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type variants include in the conditional analysis. Choices include "variant", "SNV", or "Indel" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Annotation_name a vector of qualitative/quantitative annotation names user wants to extract.
#' @return A data frame containing the basic information (chromosome, position, reference allele and alternative allele),
#' unconditional and conditional the score test p-values,
#' and annotation scores for the input variants.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

Sliding_Window_Info <- function(chr,genofile,obj_nullmodel,start_loc,end_loc,known_loci=NULL,rare_maf_cutoff=0.01,
                                method_cond=c("optimal","naive"),
                                QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,Annotation_name){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	geno_missing_imputation <- match.arg(geno_missing_imputation)
	variant_type <- match.arg(variant_type)

	## SPA status
	if(!is.null(obj_nullmodel$use_SPA))
	{
		use_SPA <- obj_nullmodel$use_SPA
	}else
	{
		use_SPA <- FALSE
	}
	
	phenotype.id <- as.character(obj_nullmodel$id_include)

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")

	## Position
	position <- as.numeric(seqGetData(genofile, "position"))

	is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
	seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	if(!use_SPA)
	{
		position <- as.numeric(seqGetData(genofile, "position"))
		Indiv_Uncond <- Indiv_Score_Test_Region(Geno,obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff)

		position <- position[!is.na(Indiv_Uncond$Score)]
		Individual_p <- Indiv_Uncond$pvalue[!is.na(Indiv_Uncond$pvalue)]
	}


	## MAF
	AF <- colMeans(Geno)/2
	MAF <- pmin(AF,1-AF)

	########################################################
	#           Annotation
	########################################################

	CHR <- as.numeric(seqGetData(genofile, "chromosome"))
	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	filter <- seqGetData(genofile, QC_label)

	Anno.Int.PHRED.sub <- NULL
	Anno.Int.PHRED.sub.name <- NULL

	for(k in 1:length(Annotation_name))
	{
		if(Annotation_name[k]%in%Annotation_name_catalog$name)
		{
			Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
			Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

			Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
		}
	}
	Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
	colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name

	if(!use_SPA)
	{
		Info_Basic <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,QC_label=filter,MAF=MAF,Score_Stat=Indiv_Uncond$Score,SE_Score=Indiv_Uncond$SE,pvalue=Indiv_Uncond$pvalue)
		Info_Basic <- Info_Basic[!is.na(Indiv_Uncond$pvalue),]

		Anno.Int.PHRED <- Anno.Int.PHRED.sub[!is.na(Indiv_Uncond$pvalue),]
	}else
	{
		Info_Basic <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,QC_label=filter,MAF=MAF)	
		Info_Basic <- Info_Basic[MAF>0,]
		
		Anno.Int.PHRED <- Anno.Int.PHRED.sub[MAF>0,]
	}

	seqResetFilter(genofile)

	##########################################################
	#                 Conditional Analysis
	##########################################################
	
	if(!use_SPA)
	{
		### known SNV Info
		if(is.null(known_loci))
		{
			known_loci <- data.frame(CHR=logical(0),POS=logical(0),REF=character(0),ALT=character(0))
		}
		known_loci_chr <- known_loci[known_loci[,1]==chr,,drop=FALSE]
		known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),,drop=FALSE]

		POS <- Info_Basic$POS

		position <- as.numeric(seqGetData(genofile, "position"))
		REF <- as.character(seqGetData(genofile, "$ref"))
		ALT <- as.character(seqGetData(genofile, "$alt"))
		variant.id <- seqGetData(genofile, "variant.id")

		sub_start_loc <- min(POS)
		sub_end_loc <- max(POS)
		### known variants needed to be adjusted
		known_loci_chr_region <- known_loci_chr[(known_loci_chr[,2]>=sub_start_loc-1E6)&(known_loci_chr[,2]<=sub_end_loc+1E6),]

		## Genotype of Adjusted Variants
		rs_num_in <- c()
		for(i in 1:dim(known_loci_chr_region)[1])
		{
			rs_num_in <- c(rs_num_in,which((position==known_loci_chr_region[i,2])&(REF==known_loci_chr_region[i,3])&(ALT==known_loci_chr_region[i,4])))
		}

		variant.id.in <- variant.id[rs_num_in]

		if(length(rs_num_in)>=1)
		{
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

			Geno_adjusted <- Geno_adjusted[,MAF>0]
			if(class(Geno_adjusted)[1]=="numeric")
			{
				Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
			}

			seqResetFilter(genofile)

			Indiv_Cond <- Indiv_Score_Test_Region_cond(Geno,Geno_adjusted,obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,method_cond=method_cond)
			pvalue_cond <- Indiv_Cond$pvalue_cond[!is.na(Indiv_Cond$pvalue_cond)]
			Score_Stat_cond <- Indiv_Cond$Score_cond[!is.na(Indiv_Cond$pvalue_cond)]
			SE_Score_cond <- Indiv_Cond$SE_cond[!is.na(Indiv_Cond$pvalue_cond)]
		}else
		{
			pvalue_cond <- Info_Basic$pvalue
			Score_Stat_cond <- Info_Basic$Score_Stat
			SE_Score_cond <- Info_Basic$SE_Score
			seqResetFilter(genofile)
		}

		Info_Basic <- cbind(Info_Basic,Score_Stat_cond,SE_Score_cond,pvalue_cond)

		### Variants_in_conditional_analysis
		if(dim(known_loci_chr)[1]>=1)
		{
			Variants_in_Cond <- rep(0,dim(known_loci_chr)[1])
			known_loci_chr <- cbind(known_loci_chr,Variants_in_Cond)
			known_loci_chr <- known_loci_chr[,-1]
			Info_Basic <- dplyr::left_join(Info_Basic,known_loci_chr,by=c("POS"="POS","REF"="REF","ALT"="ALT"))
			Info_Basic$Variants_in_Cond[is.na(Info_Basic$Variants_in_Cond)] <- 1
		}else
		{
			Variants_in_Cond <- rep(1,dim(Info_Basic)[1])
			Info_Basic <- cbind(Info_Basic,Variants_in_Cond)
		}
	}
	
	Info_Basic_Anno <- cbind(Info_Basic,Anno.Int.PHRED)

	seqResetFilter(genofile)
	return(Info_Basic_Anno)
}

