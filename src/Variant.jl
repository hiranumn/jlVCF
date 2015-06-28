export Variant, getChrom, getPos, getId, getRef, getAlt, getQual, getFilter, getInfo, getFormat

#Stores data for each variant
type Variant
    CHROM::ASCIIString
    POS::Int64
    ID::Array{ASCIIString,1}
    REF::ASCIIString
    ALT::Array{ASCIIString,1}
    #QUAL::Array{Float64, 1}
    QUAL
    FILTER::Array{ASCIIString,1} 
    INFO::Dict{ASCIIString, Any}
    FORMAT::Dict{ASCIIString, Dict{ASCIIString, Any}}
end

#Get methods for variants
getChrom(v::Variant) = v.CHROM
getPos(v::Variant) = v.POS
getId(v::Variant) = v.ID
getRef(v::Variant) = v.REF
getAlt(v::Variant) = v.ALT
getQual(v::Variant) = v.QUAL
getFilter(v::Variant) = v.FILTER
getInfo(v::Variant) = v.INFO
getFormat(v::Variant) = v.FORMAT

