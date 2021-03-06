using jlVCF

#test loading
vc = VCFIterator("data/somaticsniper.vcf")
println(vc.samples)
count = 0
v = "dmy"
while !eof(vc)
    if count%10000 == 0
        println(count)
    elseif count == 1
        println(v.POS)
    end
    count += 1 
    v = next(vc)
    println(v)
end

#reseting first time
reset(vc)
v = next(vc)
println(v)

#reseting twice 
reset(vc)
v = next(vc)
println(v)

close(vc)