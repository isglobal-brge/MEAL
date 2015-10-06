context("Create Methylation Set")

test_that("Conversion from beta matrix", {
  beta_matrix <- matrix(runif(4), nrow = 2)
  colnames(beta_matrix) <- c("H", "M")
  rownames(beta_matrix) <- c("cg00050873", "cg00212031")
  phenotypes <- data.frame(age = c(12, 23))
  rownames(phenotypes) <- c("H", "M")
  expect_match(class(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes)),
               "MethylationSet")
  beta_frame <- data.frame(beta_matrix)
  expect_match(class(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes)),
               "MethylationSet")
  annot <- data.frame(chr = c("chr1", "chr2"), pos = c(12490, 124109), stringsAsFactors = FALSE) 
  rownames(annot) <- c("cg00050873", "cg00212031")
  set <- prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes,
                                      annotation = annot)
  expect_match(class(set),"MethylationSet")
  expect_match(annotation(set), "custom")
  beta_FRAME <- DataFrame(beta_matrix)
  expect_error(prepareMethylationSet(matrix = beta_FRAME, phenotypes = phenotypes),
               "matrix is not a minfi class nor a data.frame nor a matrix")
  beta_vector <- beta_matrix[ , 1]
  expect_error(prepareMethylationSet(matrix = beta_vector, phenotypes = phenotypes), 
               "matrix is not a minfi class nor a data.frame nor a matrix")
})

test_that("Check phenotype", {
  beta_matrix <- matrix(runif(10), nrow = 5)
  colnames(beta_matrix) <- c("VAL0156", "VAL0372")
  rownames(beta_matrix) <- c("cg00050873", "cg00212031", "cg00213748", "cg00214611", "cg00455876")
  phenotypes <- c(12, 23)
  names(phenotypes) <- c("VAL0156", "VAL0372")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes), 
                     "phenotypes must be a data.frame or an AnnotatedDataFrame.")  
  phenotypes <- DataFrame(a = c(12, 23))
  rownames(phenotypes) <- c("VAL0156", "VAL0372")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes), "phenotypes must be a data.frame or an AnnotatedDataFrame.")  
})

test_that("Empty arguments", {
  beta_matrix <- matrix(nrow = 0, ncol = 0)
  phenotypes <- data.frame()
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
               "Matrix is empty.")
  phenotypes <- data.frame(age = c(12, 23))
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
                                          "Matrix is empty.")
  beta_matrix <- matrix(runif(10), nrow = 5)
  phenotypes <- data.frame()
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
                                            "Rownames of matrix must contain probe names.")
  rownames(beta_matrix) <- c("cg00050873", "cg00212031", "cg00213748", "cg00214611", "cg00455876")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
               "Colnames of matrix must contain sample names.")
  colnames(beta_matrix) <- c("VAL0156", "VAL0372")
    expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
               "Phenotypes is empty.")
  })

test_that("No names", {
  beta_matrix <- matrix(runif(4), nrow = 2)
  phenotypes <- data.frame(age = c(12, 23))
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
               "Rownames of matrix must contain probe names.")
  rownames(beta_matrix) <- c("cg00050873", "cg00212031")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
               "Colnames of matrix must contain sample names.")
  colnames(beta_matrix) <- c("H", "M")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes),
               "There are no common samples between the beta matrix and the phenotypes table. Please, check sample names.")
  rownames(phenotypes) <- c("H", "M")
  annot <- data.frame(chr = c("chr1", "chr2"), pos = c(12490, 124109)) 
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes,
                                            annotation = annot), "There are no probes with annotation.")
  rownames(annot) <- c("cg00050873", "cg00212031")
  colnames(annot) <- c("c", "pos")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes,
                                            annotation = annot), 
               "Annotation data.frame must contain a column name equal to \"chromosome\" argument.")
  colnames(annot) <- c("chr", "p")
  expect_error(prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes,
                                            annotation = annot), 
               "Annotation data.frame must contain a column name equal to \"position\" argument.")
  
  })