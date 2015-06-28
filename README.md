# jlVCF
VCF file parser module for Julia language.

Sample usage:
vc = VCFIterator("jlVCFtest/somaticsniper.vcf")
while !eof(vc)
    #v is a variable for variant class.
    v = next(vc)
    println(v)
end

















Version 0.01
Naozumi Hiranuma (hiranumn at cs dot washington dot edu)
University of Washington Computer Science and Engineering
