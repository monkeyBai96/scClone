# functions for scClone
# bsh2024


#' @title Convert SNV to signature
#' @description Extract trinucleotide context and signature from SNV table.
#' @param dat SNV data frame.
#' @param ref Genome build hg19/hg38/mm10.
#' @return Data frame with context, alteration and signature columns.
#' @export
#' @importFrom SomaticSignatures mutationContext
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom IRanges IRanges
#' @importFrom VariantAnnotation VRanges
snv_to_sig <- function(dat = dat, ref = ref){
  # column need
  colnames(dat)[1:4] <- c("Chr","Pos","Ref","Alt")
  # delete INDEL and chrUn
  dat <- dat %>% filter(Ref %in% c("A","T","C","G") & Alt %in% c("A","T","C","G"))
  dat <- dat %>% filter(Chr %in% paste0("chr", c(1:22,"X","Y")))
  dat <- dat[!duplicated(dat),]
  sca_vr = VariantAnnotation::VRanges(
    seqnames =  dat$Chr,
    ranges = IRanges::IRanges(start = as.numeric(dat$Pos),end = as.numeric(dat$Pos)),
    ref = dat$Ref,
    alt = dat$Alt)
  if(ref=="mm10"){
    sca_motifs = SomaticSignatures::mutationContext(sca_vr, BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  }else if(ref=="hg19"){
    sca_motifs = SomaticSignatures::mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  }else if(ref=="hg38"){
    sca_motifs = SomaticSignatures::mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  }else{
    cat("Error ref provided.\n")
  }
  dat1 <- as.data.frame(sca_motifs)
  rownames(dat1) <- rownames(dat)
  dat1$signature <- paste0(substring(dat1$context,1,1),"[",substring(dat1$alteration,1,1),">",substring(dat1$alteration,2,2),"]",substring(dat1$context,3,3))
  dat1 <- dat1[,c(1:3,6,7,12:14)]
  return(dat1)
}


#' @title Plot 96-trinucleotide signature
#' @description Generate COSMIC-style barplot for 96 signatures and return table.
#' @param array Named numeric vector of signature counts.
#' @param name Sample name for plot title.
#' @return List of signature table and ggplot object.
#' @export
show_signature <- function(array = array, name = name){
  context_96 <- as.array(rep(0,96))
  names(context_96) <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
                         "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
                         "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
                         "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
                         "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
                         "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  array <- array[which(names(array) %in% names(context_96))]
  context_96[match(names(array), names(context_96))] <- as.numeric(array)
  p <- sigminer::show_catalogue(context_96, mode = "SBS", style = "cosmic", samples_name = name)
  list <- list(table = context_96, plot = p)
  return(list)
}
