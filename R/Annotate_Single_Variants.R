#' Functionally annotate a list of variants
#'
#' The \code{Annotate_Single_Variants} function takes in a list of variants to functionally annotate the input variants
#' @param agds_dir file directory of annotated GDS (aGDS) files for all chromosomes (1-22).
#' @param single_variants_list a data frame containing the information of variants to be functionally annotated. The data frame must include 4 columns with
#' the following names: "CHR" (chromosome number), "POS" (position), "REF" (reference allele), and "ALT" (alternative allele).
#' @param QC_label channel name of the QC label in the GDS/aGDS file  (default = "annotation/filter").
#' @param Annotation_dir channel name of the annotations in the aGDS file (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the annotation names and the corresponding channel names in the aGDS file.
#' @param Annotation_name a vector of qualitative/quantitative annotation names user wants to extract.
#' @return a data frame containing the basic information (chromosome, position, reference allele and alternative allele)
#' and annotation scores for the input variants.
#' @export

Annotate_Single_Variants <- function(agds_dir,single_variants_list,
                                     QC_label="annotation/filter",Annotation_dir="annotation/info/FunctionalAnnotation",
                                     Annotation_name_catalog,Annotation_name){

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

			single_variants_list_annotation_chr <- data.frame(CHR=CHR,POS=position,REF=REF,ALT=ALT,QC_label=filter)
			single_variants_list_annotation_chr <- cbind(single_variants_list_annotation_chr,Anno.Int.PHRED.sub)

			single_variants_list_annotation <- rbind(single_variants_list_annotation,single_variants_list_annotation_chr)
			seqClose(genofile)
		}

	}
	single_variants_list_info_annotation <- dplyr::left_join(single_variants_list,single_variants_list_annotation,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))

	return(single_variants_list_info_annotation)
}

