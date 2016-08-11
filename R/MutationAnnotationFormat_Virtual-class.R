setClass(
    Class="MutationAnnotationFormat_Virtual",
    representation=representation(position="data.table",
                                  mutation="data.table",
                                  sample="data.table",
                                  meta="data.table", "VIRTUAL")
)

setMethod(
    f="getPosition",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@position)
    }
)

setMethod(
    f="getMutation",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@mutation)
    }
)

setMethod(
    f="getSample",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@Sample)
    }
)

setMethod(
    f="getMeta",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@meta)
    }
)
