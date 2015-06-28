export VCFIterator, getVersion, getFilename, getSampleNames, getINFOProperties, getFORMATProperties, getContigs, getFilters, getReference, next, eof, close, reset
import Base.eof, Base.close, Base.reset

#Main class for VCFiterator
type VCFIterator
    version::ASCIIString
    filename::ASCIIString
    samples::Array{ASCIIString, 1}
    infoTypes::Dict{ASCIIString, Any}
    infoFlags::Array{ASCIIString, 1}
    formatTypes::Dict{ASCIIString, Any}
    contigs::Dict{ASCIIString, Int64}
    filters::Dict{ASCIIString, ASCIIString}
    reference::ASCIIString
    file 
end

#Stores INFO and FORMAT section properties
type InfoField
    Number::ASCIIString 
    Type::ASCIIString
    Description::ASCIIString
end

#Standard get methods
getVersion(vc::VCFIterator) = vc.version
getFilename(vc::VCFIterator) = vc.filename
getSampleNames(vc::VCFIterator) = vc.samples
getINFOProperties(vc::VCFIterator) = vc.infoTypes
getFORMATProperties(vc::VCFIterator) = vc.formatTypes
getContigs(vc::VCFIterator) = vc.contigs
getFilters(vc::VCFIterator) = vc.filters
getReference(vc::VCFIterator) = vc.reference

#Gets next variant.
function next(vc::VCFIterator)
    if !eof(vc.file)
        line = readline(vc.file)
        return generateVariant(vc, line)
    else
        return false
    end
end

#Checks if it is eof.
function eof(vc::VCFIterator)
    eof(vc.file)
end

#Closes the parser
function close(vc::VCFIterator)
    close(vc.file)
end

#resets the parser to the first variant of the VCF file.
function reset(vc::VCFIterator)
    close(vc.file)
    vc.file = open(vc.filename)
    
    while !eof(vc.file)
        line = chomp(readline(vc.file))
        if line[1:2] != "##"
            break
        end
    end
end

#Parses info/format line
function readInfoLine!(typesdict::Dict{ASCIIString, Any}, line::String, flaglist)
    #parse the line
    entries = split(line[2:end-1], ',')
    Id = split(entries[1], '=')[2]
    Number = split(entries[2], '=')[2]
    Type = split(entries[3], '=')[2]
    Description = split(entries[4], '=')[2]
    temp = InfoField(Number, Type, Description[2:end-1])
    
    #add to dict
    typesdict[Id] = temp
    
    if Type == "Flag"
        push!(flaglist, Id)
    end
end

#Parses contig line
function readContigLine!(contigs::Dict{ASCIIString, Int64}, line::String)
    entries = split(line[2:end-1], ',')
    Id = split(entries[1], '=')[2]
    Length = int(split(entries[2], '=')[2])
    contigs[Id] = Length
end

#Parses filter line
function readFilterLine!(filters::Dict{ASCIIString, ASCIIString}, line::String)
    entries = split(line[2:end-1], ',')
    Id = split(entries[1], '=')[2]
    Description = split(entries[2], '=')[2]
    filters[Id] = Description[2:end-1]
end

#Constructor for VCFIterator
function VCFIterator(filename::ASCIIString)
    f=open(filename)
    vc = VCFIterator("VCFv4.1", filename, [], Dict{ASCIIString, Any}(), [], Dict{ASCIIString, Any}(), Dict{ASCIIString, Int64}(), Dict{ASCIIString, ASCIIString}(), "", f)
    parseHeader(vc)
    
    #return at the end 
    vc 
end

#Parses through header lines
function parseHeader(vc::VCFIterator)
    while !eof(vc.file)
        line = chomp(readline(vc.file))
        if line[1:2] != "##"
            samples = split(line, "\t")[10:end]
            vc.samples = samples
            break
        end
        
        #read the info line
        if length(line) >= 6 && line[1:6]=="##INFO"
            line = line[8:end]
            readInfoLine!(vc.infoTypes, line, vc.infoFlags)
        elseif length(line) >= 8 && line[1:8]=="##FORMAT"
            line = line[10:end]
            readInfoLine!(vc.formatTypes, line, ASCIIString[])
        elseif length(line) >= 8 && line[1:8]=="##contig"
            line = line[10:end]
            readContigLine!(vc.contigs, line)
        elseif length(line) >= 8 && line[1:8]=="##FILTER"
            line = line[10:end]
            readFilterLine!(vc.filters, line)
        elseif length(line) >= 11 && line[1:11]=="##reference"
            line = line[13:end]
            vc.reference = line
        end    
            
    end 
end

#Helper function that convertss x to float except when x is a period.
function floatField(x)
    if x == "."
        return "."
    else
        return float(x)
    end
end

#Helper function that convertss x to int except when x is a period.
function intField(x)
    if x == "."
        return "."
    else
        return int(x)
    end
end

#Constructor for variant given a line of vcf file
function generateVariant(vc::VCFIterator, line::ASCIIString)
    #remove newline and split by tabs. 
    line = split(chomp(line), "\t")
    
    #Assigining variables
    chr = convert(ASCIIString, line[1])
    pos = intField(line[2])
    id_ = split(line[3])
    ref = convert(ASCIIString, line[4])
    alt = split(line[5],",")
    if line[6] == "."
        qual = "."
    else
        qual = float(line[6])
    end
    filter = split(line[7],";")
    
    #creating INFO field
    information = Dict{ASCIIString, Any}()
    tempINFO = split(line[8],";")
    
    #instantiate all flags to be false
    if tempINFO[1] != "."
        for i in 1:length(vc.infoFlags)
            information[vc.infoFlags[i]] = false
        end
        for i in 1:length(tempINFO)
            entry = tempINFO[i]
            id = split(entry,'=')[1]

            #check if it is flag or not
            contentType = (vc.infoTypes[id]).Type

            if contentType != "Flag"
                contents = split(entry,'=')[2]
                #parse through contents and convert types
                contents = split(contents, ",")

                if contentType == "Integer"
                    contents = map(intField, contents)
                elseif contentType == "Float"
                    contents = map(floatField, contents)
                elseif contentType != "String" && contentType != "Char"
                    println("INVALID TYPE")
                end

                information[id] = contents
            else
                #set variable to true if the flag exists
                information[id] = true 
            end
        end
    end
    
    #creating FORMAT fields if they exist
    if length(line)>8
        format = Dict{ASCIIString, Dict{ASCIIString, Any}}()
        formatEntries = split(line[9],":")
        for i in 10:length(line)
            tempFORMAT = split(line[i],":")
            sample = vc.samples[i-9]
            tempFORMATdict = Dict{ASCIIString, Any}()
            
            for j in 1:length(tempFORMAT)
                id = formatEntries[j]
                contents = split(tempFORMAT[j], ",")
                contentType = vc.formatTypes[id].Type

                if contentType == "Integer"
                    contents = map(intField, contents)
                elseif contentType == "Float"
                    contents = map(floatField, contents)
                elseif contentType != "String" && contentType != "Char"
                    println("INVALID TYPE")
                end
            
                tempFORMATdict[id] = contents
            end
            format[sample] = tempFORMATdict
        end
    else
        format = Dict{ASCIIString,Dict{ASCIIString,Any}}()
    end

    #construct a variant
    Variant(chr, pos, id_, ref, alt, qual, filter, information, format)
end