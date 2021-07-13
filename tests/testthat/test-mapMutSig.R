test_that("input check works", {
  id <- "ENSP00000233146" # canonical MSH2
  codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                     genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                     seqLevelStyle = "UCSC")
  codons_muts <- mapMut(codons, getPossibleMissenseSNVs())
  codons_muts2 <- codons_muts
  colnames( codons_muts2 )[1] <- "random"
  expect_error(
    codons_mutsigs <- mapMutSig(id, prot_length = nrow(codons), mutMap = codons_muts,
                                ensembldbObj = "EnsDb.Hsapiens.v75",
                                genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                seqLevelStyle = "UCSC"),
    "ensembldbObj should be a EnsDb object from EnsDb packages. Check documentation of mapCodon function."
  )
  expect_error(
    codons_mutsigs <- mapMutSig(id, prot_length = nrow(codons), mutMap = codons_muts,
                                ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                genomeObj = "BSgenome.Hsapiens.UCSC.hg19",
                                seqLevelStyle = "UCSC"),
    "genomeObj should be a BSgenome object. Check documentation of mapCodon function."
  )
  expect_error(
    codons_mutsigs <- mapMutSig(id, prot_length = nrow(codons), mutMap = codons_muts,
                                ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                seqLevelStyle = "abc"),
    "seqLevelStyle must be either 'NCBI' or 'UCSC'. Check documentation of mapCodon function."
  )
  expect_error(
    codons_mutsigs <- mapMutSig(id, prot_length = nrow(codons), mutMap = codons_muts2,
                                ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                seqLevelStyle = "UCSC"),
    "mutMap appears to be different from output of the mapMut function. Please ensure you directly pass the output from mapMut for this parameter."
  )
})

test_that("complain if sequence level style is wrongly set", {
  id <- "ENSP00000233146" # canonical MSH2
  codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                     genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                     seqLevelStyle = "UCSC")
  codons_muts <- mapMut(codons, getPossibleMissenseSNVs())
  expect_error(
    codons_mutsigs <- mapMutSig(id, prot_length = nrow(codons), mutMap = codons_muts,
                                ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                seqLevelStyle = "NCBI"),
    "Failed to fetch CDS"
  )
})

test_that("check correct output", {
  id <- "ENSP00000233146"
  id <- "ENSP00000233146" # canonical MSH2
  codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                     genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                     seqLevelStyle = "UCSC")
  codons_muts <- mapMut(codons, getPossibleMissenseSNVs())
  codons_mutsigs <- mapMutSig(id, prot_length = nrow(codons), mutMap = codons_muts,
                              ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                              genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                              seqLevelStyle = "UCSC")
  expect_is(codons_mutsigs, "data.frame")
  expect_equal(nrow(codons_mutsigs), 6590)
  expect_equal(ncol(codons_mutsigs), 8)
  expect_equal( codons_mutsigs[3, "g_pos"], 47630332 )
  expect_equal( codons_mutsigs[273, "AA_pos"], 40 )
  expect_equal( codons_mutsigs[273, "WT_AA"], "G" )
  expect_equal( codons_mutsigs[273, "MUT_AA"], "C" )
  expect_equal( codons_mutsigs[504, "MutSig"], "A[T>A]C" )
})
