test_that("input check works", {
  id <- "ENSP00000233146" # canonical MSH2
  expect_error(
    codons <- mapCodon(id, ensembldbObj = "EnsDb.Hsapiens.v75",
                       genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                       seqLevelStyle = "UCSC"),
    "ensembldbObj should be a EnsDb object from EnsDb packages. Check documentation of mapCodon function."
  )
  expect_error(
    codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                       genomeObj = "BSgenome.Hsapiens.UCSC.hg19",
                       seqLevelStyle = "UCSC"),
    "genomeObj should be a BSgenome object. Check documentation of mapCodon function."
  )
  expect_error(
    codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                       genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                       seqLevelStyle = "abc"),
    "seqLevelStyle must be either 'NCBI' or 'UCSC'. Check documentation of mapCodon function."
  )
})

test_that("complain if sequence level style is wrongly set", {
  id <- "ENSP00000233146"
  expect_error(
    codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                       genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                       seqLevelStyle = "NCBI"),
    "Failed to fetch CDS"
  )
})

test_that("check correct output", {
  id <- "ENSP00000233146"
  codons <- mapCodon(id, ensembldbObj = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                     genomeObj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                     seqLevelStyle = "UCSC")
  expect_is(codons, "data.frame")
  expect_equal(nrow(codons), 934)
  expect_equal(ncol(codons), 4)
  expect_equal( codons[3, "codon"], "GTG" )
  expect_equal( codons[273, "AA_pos"], 273 )
  expect_equal( codons[273, "AA_pos"], 273 )
  expect_equal( codons[273, "AA"], "V" )
  expect_equal( codons[504, "exon"], "9;10;10" )
})


