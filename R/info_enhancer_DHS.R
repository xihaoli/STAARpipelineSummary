info_enhancer_DHS <- function(chr,genofile,obj_nullmodel,gene_name,known_loci,rare_maf_cutoff=0.01,
                              method_cond=c("optimal","naive"),
                              QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                              Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,Annotation_name){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)
	
	## SPA status
	if(!is.null(obj_nullmodel$use_SPA))
	{
		use_SPA <- obj_nullmodel$use_SPA
	}else
	{
		use_SPA <- FALSE
	}

	n_pheno <- obj_nullmodel$n.pheno	
	phenotype.id <- as.character(obj_nullmodel$id_include)

	## Enhancer
	varid <- seqGetData(genofile, "variant.id")

	# Now extract the GeneHancer with rOCRs Signal Overlay
	genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
	genehancer <- genehancerAnno!=""

	rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
	rOCRs <- rOCRsAnno!=""
	rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
	rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
	seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])

	# Variants that covered by whole GeneHancer without rOCRs overlap
	genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
	enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
	enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
	enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
	enhancervpos <- as.numeric(seqGetData(genofile,"position"))
	enhancervref <- as.character(seqGetData(genofile,"$ref"))
	enhancervalt <- as.character(seqGetData(genofile,"$alt"))
	dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

	rm(varid)
	gc()

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

	rm(filter)
	gc()

	variant.id <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- variant.id[SNVlist]

	dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
	dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
	dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
	dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

	seqResetFilter(genofile)

	rm(dfHancerrOCRsVarGene)
	gc()

	### Gene
	is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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
		
		if(n_pheno==1)
		{
			Indiv_Uncond <- Indiv_Score_Test_Region(Geno,obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff)

			position <- position[!is.na(Indiv_Uncond$Score)]
			Individual_p <- Indiv_Uncond$pvalue[!is.na(Indiv_Uncond$pvalue)]
		}else
		{
			Indiv_Uncond <- Indiv_Score_Test_Region_multi(Geno,obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff)	
			
			position <- position[!is.na(Indiv_Uncond$pvalue)]
			Individual_p <- Indiv_Uncond$pvalue[!is.na(Indiv_Uncond$pvalue)]
		}
	}

	## MAF
	AF <- colMeans(Geno)/2
	MAF <- pmin(AF,1-AF)

	########################################################
	#                 Annotation
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
		if(n_pheno==1)
		{
			Info_Basic <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,QC_label=filter,MAF=MAF,Score_Stat=Indiv_Uncond$Score,SE_Score=Indiv_Uncond$SE,pvalue=Indiv_Uncond$pvalue)
			Info_Basic <- Info_Basic[!is.na(Indiv_Uncond$pvalue),]

			Anno.Int.PHRED <- Anno.Int.PHRED.sub[!is.na(Indiv_Uncond$pvalue),]
		}else
		{
			Info_Basic <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,QC_label=filter,MAF=MAF,pvalue=Indiv_Uncond$pvalue)
			Info_Basic <- Info_Basic[!is.na(Indiv_Uncond$pvalue),]

			Anno.Int.PHRED <- Anno.Int.PHRED.sub[!is.na(Indiv_Uncond$pvalue),]
		
		}
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

			if(n_pheno==1)
			{
				Indiv_Cond <- Indiv_Score_Test_Region_cond(Geno,Geno_adjusted,obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,method_cond=method_cond)
				pvalue_cond <- Indiv_Cond$pvalue_cond[!is.na(Indiv_Cond$pvalue_cond)]
				Score_Stat_cond <- Indiv_Cond$Score_cond[!is.na(Indiv_Cond$pvalue_cond)]
				SE_Score_cond <- Indiv_Cond$SE_cond[!is.na(Indiv_Cond$pvalue_cond)]
			}else
			{
				Indiv_Cond <- Indiv_Score_Test_Region_multi_cond(Geno,Geno_adjusted,obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,method_cond=method_cond)
				pvalue_cond <- Indiv_Cond$pvalue_cond[!is.na(Indiv_Cond$pvalue_cond)]
			}
		}else
		{
			if(n_pheno==1)
			{
				pvalue_cond <- Info_Basic$pvalue
				Score_Stat_cond <- Info_Basic$Score_Stat
				SE_Score_cond <- Info_Basic$SE_Score
			}else
			{
				pvalue_cond <- Info_Basic$pvalue
			}
			seqResetFilter(genofile)
		}
	}

	Gene <- rep(gene_name,dim(Info_Basic)[1])
	
	if(!use_SPA)
	{
		if(n_pheno==1)
		{
			Info_Basic <- cbind(Gene,Info_Basic,Score_Stat_cond,SE_Score_cond,pvalue_cond)
		}else
		{
			Info_Basic <- cbind(Gene,Info_Basic,pvalue_cond)
		}
		
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
	}else
	{
		Info_Basic <- cbind(Gene,Info_Basic)
	}

	Info_Basic_Anno <- cbind(Info_Basic,Anno.Int.PHRED)

	seqResetFilter(genofile)
	return(Info_Basic_Anno)
}

