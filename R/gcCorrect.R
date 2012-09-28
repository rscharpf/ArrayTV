setGeneric("gcCorrect", function(object, ...) {standardGeneric("gcCorrect")})

setMethod("gcCorrect", signature(object="matrix"), function(object, ...){
    gcCorrectMain(object, ...)
})
