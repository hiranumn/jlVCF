#jlVCF
VCF file parser module for Julia language.

##Sample usage:
```
vc = VCFIterator("jlVCFtest/somaticsniper.vcf")
while !eof(vc)
    #v is a variable for variant class.
    v = next(vc)
	print(v.POS, v.CHROM
end
```

##Classes and their fields
See vcf file format specification for the detailed explanation of what each field means. 
* Variant
  * CHROM::ASCIIString
  * POS::Int64
  * ID::Array{ASCIIString,1}
  * REF::ASCIIString
  * ALT::Array{ASCIIString,1}
  * QUAL
  * FILTER::Array{ASCIIString,1}
  * INFO::Dict{ASCIIString, Any}; keys of this dict are infokey and the values are infovalue.
  * FORMAT::Dict{ASCIIString, Dict{ASCIIString, Any}}. The first layer of keys are sample names, and the second layer keys are fo matkeys.
* VCFIterator
  * version::ASCIIString; version of this VCF file (not functional yet)
  * filename::ASCIIString; name of the source vcf file
  * samples::Array{ASCIIString, 1}; list of sample names
  * contigs::Dict{ASCIIString, Int64}; dict of contigs with contig names as keys and positions as values 
  * filters::Dict{ASCIIString, ASCIIString}; dict of filters with filter names as keys and their descriptions as values 
  * reference::ASCIIString; reference genome name

##Available functions
*VCFIterator related
  * VCFIterator(filename::ASCIIString)
  * getVersion(vc::VCFIterator)
  * getFilename(vc::VCFIterator)
  * getSampleNames(vc::VCFIterator)
  * getINFOProperties(vc::VCFIterator)
  * getFORMATProperties(vc::VCFIterator)
  * getContigs(vc::VCFIterator)
  * getFilters(vc::VCFIterator)
  * getReference(vc::VCFIterator)
  * next(vc::VCFIterator)
  * eof(vc::VCFIterator)
  * close(vc::VCFIterator)
  * reset(vc::VCFIterator)
*Variant related
  * etChrom(v::Variant)
  * getPos(v::Variant)
  * getId(v::Variant)
  * getRef(v::Variant)
  * getAlt(v::Variant)
  * getQual(v::Variant)
  * getFilter(v::Variant)
  * getInfo(v::Variant)
  * getFormat(v::Variant)

Version 0.01

Naozumi Hiranuma (hiranumn at cs dot washington dot edu)

University of Washington Computer Science and Engineering
