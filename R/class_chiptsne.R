.datatable.aware=TRUE
setOldClass(c('data.frame'))
setOldClass(c('data.table', 'data.frame'))

setClass(Class = "ChIPtsne",

         slots = c(
             config_dt = "data.table",
             col_attrib = "character",
             row_attrib = "character",
             id_attrib = "character",
             profile_dt = "data.table",
             query_gr = "GRanges",
             tsne_res = "data.table",
             hic_1d = "data.table"

         ),

         validity = function(object){
             errors <- character()
             # mat_cnames = c("i", "j", "val")
             # if (length(intersect(colnames(object@hic_2d), mat_cnames)) != length(mat_cnames)){
             #     msg <- "colnames of hic_2d must be c(i, j, val)"
             #     errors <- c(errors, msg)
             # }
             # reg_cnames = c("seqnames", "start", "end", "index")
             # if (length(intersect(colnames(object@hic_1d), reg_cnames)) != length(reg_cnames)){
             #     msg <- "colnames of hic_1d must be c(seqnames, start, end, index)"
             #     errors <- c(errors, msg)
             # }
             # if (length(errors) == 0) TRUE else errors
             TRUE
         }
)

setMethod("initialize", "ChIPtsne", function(.Object, matrix_file, regions_file, parameters) {
    if(missing(matrix_file) & missing(regions_file) & missing(parameters)){
        return(.Object)
    }

    ### set provided
    .Object@matrix_file = matrix_file
    .Object@regions_file = regions_file
    .Object@parameters = parameters

    ### do stuff



    validObject(.Object)
    .Object

})

ChIPtsne = function(config_dt){
    #process args
    new("ChIPtsne",
        matrix_file = matrix_file,
        regions_file = regions_file,
        parameters = parameters)
}

setMethod("show", "ChIPtsne",
          function(object) {
              print(object@config_dt)
          }
)
