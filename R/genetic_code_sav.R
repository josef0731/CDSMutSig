#' CDSMutSig: Mapping Coding Sequences between DNA Mutational signatures and Amino Acid Changes
#'
#' A set of functions which use Ensembl and Bioconductor genome objects
#' to map between given protein sequences to corresponding coding DNA
#' sequences, extract single-amino acid variants (SNVs) and corresponding
#' nucleotide substitutions.
#'
#' @docType package
#' @name CDSMutSig
NULL
#> NULL


################################
#' Get possible missense substitutions arising from single-nucleotide variants
#'
#' This function removes missense substitutions that are impossible to achieve
#' without accumulating multiple single-nucleotide variants, based on the standard
#' genetic code.
#'
#' @return A data.frame with four columns:
#'   \code{MUT}: codon corresponding to mutant amino acid
#'   \code{WT}: codon corresponding to wild-type amino acid
#'   \code{WT_AA}: wild-type amino acid (one-letter code)
#'   \code{MUT_AA}: wild-type amino acid (one-letter code)
#'
#' @importFrom Biostrings GENETIC_CODE stringDist reverseComplement DNAString
#' @importFrom reshape2 melt
#' @export
getPossibleMissenseSNVs <- function()
{
  genetic_code <- data.frame( Biostrings::GENETIC_CODE )
  genetic_code$code <- rownames( genetic_code )
  colnames( genetic_code ) <- c( "AA", "code" )
  # edit distance
  subs_dist <- as.matrix( Biostrings::stringDist( genetic_code$code, diag = FALSE ) )
  rownames( subs_dist ) <- genetic_code$code; colnames( subs_dist ) <- genetic_code$code
  subs_dist <- reshape2::melt( subs_dist )
  colnames( subs_dist ) <- c("WT", "MUT", "dist")
  snvs <- subs_dist[ which(subs_dist$dist == 1), ]
  snvs <- merge( snvs, genetic_code, by.x = "WT", by.y = "code",
                 all.x = TRUE, all.y = FALSE, sort = FALSE)
  colnames( snvs ) <- c("WT", "MUT", "dist", "WT_AA")
  snvs <- merge( snvs, genetic_code, by.x = "MUT", by.y = "code",
                 all.x = TRUE, all.y = FALSE, sort = FALSE)
  colnames( snvs ) <- c("MUT", "WT", "dist", "WT_AA", "MUT_AA")
  snvs_missense <- snvs[ which(snvs$WT_AA != "*"), ]
  snvs_missense <- snvs[ which(snvs$MUT_AA != "*"), ]
  snvs_missense <- snvs[ which(snvs$WT_AA != snvs$MUT_AA), ]
  out <- snvs_missense[, -3] # remove the distance column
  out$WT <- as.character( out$WT )
  out$MUT <- as.character( out$MUT )
  out$WT_AA <- as.character( out$WT_AA )
  out$MUT_AA <- as.character( out$MUT_AA )
  # reverse complement the WT and MUT columns
  out2 <- out
  out2$WT <- sapply(out2$WT, function(x) as.character( Biostrings::reverseComplement( Biostrings::DNAString(x) )))
  out2$MUT <- sapply(out2$MUT, function(x) as.character( Biostrings::reverseComplement( Biostrings::DNAString(x) )))
  out <- do.call("rbind", list(out, out2))
  out
}

#' Get amino acid sequence ranges residing in each exon
#'
#' Given a GRanges object of exon boundaries, this function determines
#' the range of amino acids residing in each exon. Split AA will appear in both exons.
#'
#' @param GRangesObj A \code{GRangeObj} containing boundaries of each exon in the coding sequence. Generated internally in \code{mapCodon}.
#'
#' @return A data.frame with six columns:
#'   \code{exon}: exon number
#'   \code{lengths}: length (in nucleotides) of the exon
#'   \code{prot_start}: amino acid position at the beginning of the exon
#'   \code{prot_end}: amino acid position at the end of the exon
#'   \code{carried_forfward}: the number of nucleotides from the previous exon which does not yield a complete codon.
#'   \code{accu_length}: running sum of length (in nucleotides) of the coding DNA sequence.
#'
#' @importFrom GenomicRanges width
getExonProtBoundary <- function(GRangesObj){
  out <- data.frame(exon = 1:length(GRangesObj), lengths = GenomicRanges::width(GRangesObj))
  out$prot_start <- NA; out$prot_end <- NA; out$carried_forward <- NA
  # accumulated length in nt
  out$accu_length <- NA
  for(i in 1:length(GRangesObj)){
    out[i, "accu_length"] <- sum(out$lengths[1:i])
  }
  for(i in 1:length(GRangesObj)){
    if(i == 1){
      real_length <- out[1, "lengths"]
      out[1, "prot_start"] <- 1
      out[1, "carried_forward"] <- 0
    } else {
      # figure out effective length after accommodating the carried-forward AA from previous exon
      real_length <- out[i, "lengths"] - carry_forward
      out[i, "prot_start"] <- out[i-1, "prot_end"]
      out[i, "carried_forward"] <- carry_forward
      # adjust the start AA position based on the number of nt carried forward from previous exon
      if( carry_forward %in% c(0, 3) ) out[i, "prot_start"] <- out[i, "prot_start"] + 1
    }
    real_codons <- floor( real_length / 3)
    # figure out the number of nt remains after fitting in maximum possible number of AA
    if(i == 1) remainder <- real_length - real_codons * 3
    else if(i > 1) remainder <- out[i, "accu_length"] - out[i, "prot_start"] * 3 - real_codons * 3
    # get end AA position after fitting in maximum possible numbers of AA
    out[i, "prot_end"] <- out[i, "prot_start"] + real_codons - 1
    if( remainder > 0 ) out[i, "prot_end"] <- out[i, "prot_end"] + 1
    if( i > 1 ){
      # adjust end AA position based on whether there is a split AA  carried from the previous exon
      if( carry_forward > 0) out[i, "prot_end"] <- out[i, "prot_end"] + 1
    }
    # set the number of nt to carry forward to next exon
    carry_forward <- 3 - remainder
  }
  out#[, c("exon", "prot_start", "prot_end")]
}

#' Map codons for each amino acid for a given protein.
#'
#' @param prot_id character, Ensembl protein ID
#' @param ensembldbObj Ensembl DB Object from the \code{EnsDb} package (see example)
#' @param genomeObj Bioconductor genome object from which coding sequences are extracted (see example)
#' @param seqLevelStyle Either "UCSC" or "NCBI". Denote conventions used in genome object to name chromosomes. Need changing to match the genome object. (Default: "NCBI")
#'
#' @return A data.frame with four columns:
#'   \code{AA_pos}: integer, amino acid position of the given Ensembl protein
#'   \code{AA}: amino acid (one-letter code)
#'   \code{codon}: codon corresponding to \code{AA}
#'   \code{exon}: the exon(s) in which the given \code{AA} can be found. Amino acids which span across exons are denoted as a semicolon-separated triplets of exon IDs, with each entry corresponding to each position of the codon.
#'
#' @examples
#' require(EnsDb.Hsapiens.v75)
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' id <- "ENSP00000233146"
#' codons <- mapCodon(id, EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
#'                    BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
#'                    seqLevelStyle = "UCSC")
#' @importFrom ensembldb proteins proteinToGenome seqlevelsStyle
#' @importFrom AnnotationFilter ProteinIdFilter
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom GenomeInfoDb genome
#' @importFrom BSgenome getSeq
#' @export
mapCodon <- function(prot_id, ensembldbObj, genomeObj, seqLevelStyle = "NCBI")
{
  if( class( ensembldbObj ) != "EnsDb" )
    stop("ensembldbObj should be a EnsDb object from EnsDb packages. Check documentation of mapCodon function.")
  if( class( genomeObj ) != "BSgenome" )
    stop("genomeObj should be a BSgenome object. Check documentation of mapCodon function.")
  if( !seqLevelStyle %in% c("NCBI", "UCSC") )
    stop("seqLevelStyle must be either 'NCBI' or 'UCSC'. Check documentation of mapCodon function.")
  sequence <- ensembldb::proteins( ensembldbObj, filter = AnnotationFilter::ProteinIdFilter(prot_id),
                                   return.type = "AAString")
  sequence <- as.character( sequence )
  seqrange <- IRanges::IRanges(start = 1, end = nchar(sequence), names = prot_id)
  g_range <- ensembldb::proteinToGenome(seqrange, ensembldbObj)[[1]]
  if( length(g_range) == 0 ){
    return( NULL )
  }
  # split AA sequence into a vector of 1-letter codes
  wt <- unlist(strsplit(sequence, split = ''))
  pos <- 1:length(wt)
  exons <- getExonProtBoundary( g_range ) # mapping of AA pos to exons
  # for split AA positions, parse how many nt is mapped to either exon
  exons$split_aa_ntUsed <- c(3 - exons$carried_forward[-1], 0)
  split_aa <- intersect(exons$prot_end, exons$prot_start) # list of split AA positions
  exon_map <- IRanges::IRanges(start = exons$prot_start, end = exons$prot_end)
  # get CDS
  if( seqLevelStyle == "UCSC" ){
    ensembldb::seqlevelsStyle(g_range) <- "UCSC"; GenomeInfoDb::genome(g_range) <- "hg19"
  }
  if( any( ! as.character(GenomicRanges::seqnames(g_range)) %in%
           GenomeInfoDb::seqnames(genomeObj) ) )
    stop("Failed to fetch CDS. Are you sure you have set the correct seqLevelStyle? (Check package vignette for details.)")
  g_seq <- BSgenome::getSeq(genomeObj, g_range)
  cds_seq <- paste(g_seq, collapse = "")
  # split CDS in codons
  out <- data.frame( AA_pos = pos, AA = wt, stringsAsFactors = FALSE )#, codon = unlist(cds_codons) )
  out$codon <- sapply(out$AA_pos, function(x) substr(cds_seq, 3*x-2, 3*x))
  out$exon <- sapply(out$AA_pos, function(x){
    if(x %in% split_aa){
      n <- exons[exons$prot_end == x, "split_aa_ntUsed"][1]
      paste(c(rep(exons[exons$prot_end == x, "exon"], n),
              rep(exons[exons$prot_start == x, "exon"], 3 - n))[1:3], collapse = ";")
    } else {
      as.data.frame(IRanges::findOverlaps(IRanges::IRanges(x), exon_map))$subjectHits
    }
  })
  rownames(out) <- NULL
  out
}

#' Map mutations given a table of codons for each amino acid for a given protein.
#'
#' @param codonMap data.frame of codons for each amino acid obtained from \code{mapCodon} (see example)
#' @param possibleMuts data.frame of possible mutations arising froms single-nucleotide susbtitutions obtained from \code{getPossibleMissenseSNVs} (see example)
#'
#' @return A data.frame with six columns:
#'   \code{AA}: wild-type amino acid (one-letter code)
#'   \code{WT_codon}: codon corresponding to wild-type amino acid
#'   \code{AA_pos}: integer, amino acid position of the given Ensembl protein
#'   \code{exon}: the exon(s) in which the given \code{AA} can be found. Amino acids which span across exons are denoted as a semicolon-separated triplets of exon IDs, with each entry corresponding to each position of the codon.
#'   \code{MUT_codon}: codon corresponding to mutant amino acid
#'   \code{MUT_AA}: mutant amino acid (one-letter code)
#'
#' @examples
#' require(EnsDb.Hsapiens.v75)
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' id <- "ENSP00000233146"
#' codons <- mapCodon(id, EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
#'                    BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
#'                    seqLevelStyle = "UCSC")
#' muts <- mapMut(codons, getPossibleMissenseSNVs())
#' @export
mapMut <- function(codonMap, possibleMuts)
{
  out <- merge(codonMap, possibleMuts, by.x = c("AA", "codon"),
               by.y = c("WT_AA", "WT"), all.x = TRUE, sort = FALSE)
  colnames(out) <- c("AA", "WT_codon", "AA_pos", "exon", "MUT_codon", "MUT_AA")
  out[order(out$AA_pos), ]
}

#' Get the codon position mutated to yield a given missense substitution
#'
#' @param wt_c wild-type codon
#' @param mut_c mutant codon
#'
#' @return A integer with value of 1, 2 or 3, corresponding to the position within the codon which is mutated from \code{wt_c} to \code{mut_c}.
#'
#' @examples
#' getMutCodonPos( "ATG", "CTG" ) # return 1
#' getMutCodonPos( "ATG", "ACG" ) # return 2
#' getMutCodonPos( "ATG", "ATA" ) # return 3
#'
#' @export
getMutCodonPos <- function(wt_c, mut_c){
  # get which position of the codon (i.e. 1, 2 or 3) which is substituted
  # comparing the wt codon and mut codon
  out <- sapply(1:3, function(x) substr(wt_c, x, x) == substr(mut_c, x, x))
  which(!out)
}

#' Get 96 mutational contexts used in defining mutational signatures
#'
#' Used as a general function to obtain all 96 possibilities for counting, graphing etc.
#'
#' @return A character vector containing the 96 mutational contexts defined with the single-base substitution (C>A, C>T, C>G, T>A, T>C, T>G) together with 1 nt 5' and 1 nt 3'.
#'
#' @export
get96Contexts <- function()
{
  # 96 contexts (3bp, with 6 substitutions)
  bases <- c("A", "C", "G", "T")
  subs <- c(sapply(c("C", "T"),
                   function(y) paste(y, bases[ -(which(bases == y)) ], sep =">") ))
  c(sapply(subs, function(y){
    c(sapply(bases, function(x) paste0(x, "[", y, "]", bases)))
  }))
}

#' Get 32 DNA motifs used in defining mutational signatures
#'
#' Used as a general function to obtain all 32 possibilities for counting, graphing etc.
#'
#' @return A character vector containing the 32 DNA motifs defined with the wild-type nucleotide (either C or T) at which single-base substitution occurs, together with 1 nt 5' and 1 nt 3'
#'
#' @export
get32Contexts <- function()
{
  # 32 contexts (3bp with C or T in the middle, without substitutions)
  bases <- c("A", "C", "G", "T")
  subs <- c("C", "T")
  c(sapply(subs, function(y){
    c(sapply(bases, function(x) paste0(x,  y, bases)))
  }))
}

#' Map mutations given a table of codons for each amino acid for a given protein.
#'
#' @param prot_id character, Ensembl protein ID
#' @param prot_length integer, length (in amino acids) of the given protein
#' @param mutMap data.frame list of mutations obtained from \code{mapMut}
#' @param ensembldbObj Ensembl DB Object from the \code{EnsDb} package (see example)
#' @param genomeObj Bioconductor genome object from which coding sequences are extracted (see example)
#' @param flank integer, the number of flanking nucleotides with which to define mutational signatures. Only the value of 1 has been catered for at the moment. (Default: 1)
#' @param contexts character vector, all possible DNA motifs on which mutational signatures are defined. Obtained from \code{get32Contexts()} (see example)
#' @param seqLevelStyle Either "UCSC" or "NCBI". Denote conventions used in genome object to name chromosomes. Need changing to match the genome object. (Default: "NCBI")
#'
#' @return A data.frame with eight columns:
#'   \code{chr}: chromosome
#'   \code{g_pos}: genomic position
#'   \code{AA_pos}: integer, amino acid position of the given Ensembl protein
#'   \code{WT_AA}: wild-type amino acid (one-letter code)
#'   \code{MUT_AA}: mutant amino acid (one-letter code)
#'   \code{WT_codon}: codon corresponding to wild-type amino acid
#'   \code{MUT_codon}: codon corresponding to mutant amino acid
#'   \code{MutSig}: the corresponding mutational signature (96-contexts) corresponding to the given AA substitution
#'
#' @examples
#' require(EnsDb.Hsapiens.v75)
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' id <- "ENSP00000233146"
#' codons <- mapCodon(id, EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
#'                    BSgenome.Hsapiens.UCSC.hg19::Hsapiens, seqLevelStyle = "UCSC")
#' codons_muts <- mapMut(codons, getPossibleMissenseSNVs())
#' codons_mutsigs <- mapMutSig(id, nrow(codons), codons_muts,
#'                             EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
#'                             BSgenome.Hsapiens.UCSC.hg19::Hsapiens, seqLevelStyle = "UCSC")
#' @importFrom ensembldb proteinToGenome
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb genome
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings DNAString reverseComplement
#' @export
mapMutSig <- function(prot_id, prot_length, mutMap, ensembldbObj, genomeObj,
                      flank = 1, contexts = get32Contexts(), seqLevelStyle = "NCBI")
{
  if( class( ensembldbObj ) != "EnsDb" )
    stop("ensembldbObj should be a EnsDb object from EnsDb packages. Check documentation of mapCodon function.")
  if( class( genomeObj ) != "BSgenome" )
    stop("genomeObj should be a BSgenome object. Check documentation of mapCodon function.")
  if( !seqLevelStyle %in% c("NCBI", "UCSC") )
    stop("seqLevelStyle must be either 'NCBI' or 'UCSC'. Check documentation of mapCodon function.")
  if( class( all.equal.character( colnames( mutMap ),
                                  c("AA", "WT_codon", "AA_pos",  "exon", "MUT_codon", "MUT_AA")) ) == "character" )
    stop("mutMap appears to be different from output of the mapMut function. Please ensure you directly pass the output from mapMut for this parameter.")
  seqrange <- IRanges::IRanges(start = 1, end = prot_length, names = prot_id)
  g_range <- ensembldb::proteinToGenome(seqrange, ensembldbObj)[[1]]
  mutMap <- mutMap[ which(!is.na(mutMap$MUT_codon)), ]
  # which position of the codon (1, 2 or 3) is mutated?
  mutMap$codon_pos <- apply(mutMap[, c("WT_codon", "MUT_codon")], MARGIN = 1,
                            function(x) getMutCodonPos(x[1], x[2]))
  mutMap$exon <- apply(mutMap[, c("exon", "codon_pos")], MARGIN = 1, function(x){
    if( grepl(";", x[1])){
      # if split exon, identify which exon the mutated nt is in, and substitute
      exon_list <- unlist(strsplit(x[1], split = ";"))
      as.numeric( exon_list[ as.numeric(x[2]) ] )
    } else as.numeric( x[1] )
  })
  # get genomic position of the substitution
  # = start(r) + 3n + m - 4 - sum(l_x) where x in [1, ..., r - 1]
    # n = AA position
    # m = mutated position in codon (i.e. either 1, 2 or 3), and
    # this mutated position is in the rth exon
  g_pos <- data.frame(g_range[, 1:4])
  g_pos$seqnames <- as.character(g_pos$seqnames)
  mutMap_gpos <- apply(mutMap, MARGIN = 1, function(x){
    n <- as.numeric(x[3]); m <- as.numeric(x[7]); r <- as.numeric(x[4])
    pos <- as.numeric(g_pos[r, "start"]) + 3 * n + m - 4
    if( r > 1 ) pos <- pos - sum( g_pos[1:(r-1), "width"])
    return(c(g_pos[1, "seqnames"], pos))
  })
  mutMap_gpos <- data.frame(t(mutMap_gpos), stringsAsFactors = FALSE)
  colnames(mutMap_gpos) <- c("chr", "pos")
  mutMap_gpos$pos <- as.numeric(mutMap_gpos$pos)
  # get the CDS with flanking lengths added to starts and ends of each exon;
  # subsetting this sequence at correct positions will give the sequence context of subs
  GenomicRanges::start(g_range) <- GenomicRanges::start(g_range) - flank
  GenomicRanges::end(g_range) <- GenomicRanges::end(g_range) + flank
  if( seqLevelStyle == "UCSC" ){
    ensembldb::seqlevelsStyle(g_range) <- "UCSC"; GenomeInfoDb::genome(g_range) <- "hg19"
  }
  if( any( ! as.character(GenomicRanges::seqnames(g_range)) %in%
           GenomeInfoDb::seqnames(genomeObj) ) )
    stop("Failed to fetch CDS. Are you sure you have set the correct seqLevelStyle? (Check package vignette for details.)")
  context_seq <- BSgenome::getSeq(genomeObj, g_range)
  context_seq <- paste(context_seq, collapse = "")
  # extract WT and MUT nucleotides
  mutMap$WT_nt <- substr(mutMap$WT_codon, mutMap$codon_pos, mutMap$codon_pos)
  mutMap$MUT_nt <- substr(mutMap$MUT_codon, mutMap$codon_pos, mutMap$codon_pos)
  # transform the codon_pos to actual pos in the assembled context_seq
  mutMap$codon_pos <- apply(mutMap[, c("AA_pos", "exon", "codon_pos")], 1, function(x){
    # the transformed position =
    # 3n + k + m - 3 where:
    # n = AA position
    # m = mutated position in codon (i.e. either 1, 2 or 3), and
    # k = f + 2f(r - 1) if this mutated position is in the rth exon and f = flank
    x <- as.numeric(x)
    3*x[1] + flank + 2 * flank * (x[2] - 1) + x[3] - 3
  })
  mutMap$context <- sapply(mutMap$codon_pos, function(x){
    substr(context_seq, x - flank, x + flank)
  })
  for(i in 1:nrow(mutMap)){
    if( !mutMap[i, "context"] %in% contexts ){
      context <- Biostrings::reverseComplement( Biostrings::DNAString(mutMap[i, "context"]) )
      mut <- Biostrings::reverseComplement( Biostrings::DNAString(mutMap[i, "MUT_nt"]) )
      mutMap[i, "context"] <- as.character(context)
      mutMap[i, "MUT_nt"] <- as.character(mut)
    }
  }
  mutMap$MutSig <- apply(mutMap[, c("context", "MUT_nt")], 1, function(x){
    paste0(substr(x[1], 1, flank), "[", substr(x[1], 1 + flank, 1 + flank),
           ">", x[2], "]", substr(x[1], 2 + flank, nchar(x[1])))
  })
  out <- mutMap[, c("AA_pos", "AA", "MUT_AA", "WT_codon", "MUT_codon", "MutSig")]
  out <- data.frame( mutMap_gpos, out, stringsAsFactors = FALSE )
  colnames(out) <- c("chr", "g_pos", "AA_pos", "WT_AA", "MUT_AA", "WT_codon", "MUT_codon", "MutSig")
  out
}
