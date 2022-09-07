# epiread Record
# ==========

mutable struct Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # NOME-seq or not
    nome::Int
    # indexes
    chrom::UnitRange{Int}
    chromstart::UnitRange{Int}
    chromend::UnitRange{Int}
    name::UnitRange{Int}
    readnum::Int
    strand::Int
    cg_rle::UnitRange{Int}
    gc_rle::UnitRange{Int}
end

"""
    EPIREAD.Record()

Create an unfilled EPIREAD record.
"""
function Record()
    return Record(
		# data, filled, nome
        UInt8[], 1:0, 0,
        # chrom, start, end, name, readnum
        1:0, 1:0, 1:0, 1:0, 0,
        # strand, cg_rle, gc_rle
        0, 1:0, 1:0
	)
end

"""
    EPIREAD.Record(data::Vector{UInt8})

Create a EPIREAD record object from `data`.

This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(
		# data, filled, nome
        data, 1:0, 0,
        # chrom, start, end, name, readnum
        1:0, 1:0, 1:0, 1:0, 0,
        # strand, cg_rle, gc_rle
        0, 1:0, 1:0
    )
    index!(record)
    return record
end

"""
    EPIREAD.Record(str::AbstractString)

Create a EPIREAD record object from `str`.

This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return convert(Record, Vector{UInt8}(str))
end

function Base.empty!(record::Record)
    record.filled = 1:0
    record.nome = 0
    record.chrom = 1:0
    record.chromstart = 1:0
    record.chromend = 1:0
    record.name = 1:0
    record.readnum = 0
    record.strand = 0
    record.cg_rle = 1:0
    record.gc_rle = 1:0
    return record
end

function GenomicFeatures.Interval(record::Record)
    name = BioGenerics.seqname(record)
    lpos = BioGenerics.leftposition(record)
    rpos = BioGenerics.rightposition(record)
    strd = hasstrand(record) ? GenomicFeatures.strand(record) : GenomicFeatures.STRAND_BOTH
    return GenomicFeatures.Interval(name, lpos, rpos, strd, record)
end

function Base.convert(::Type{GenomicFeatures.Interval}, record::Record)
    return GenomicFeatures.Interval(record)
end

function Base.convert(::Type{GenomicFeatures.Interval{Record}}, record::Record)
    return convert(GenomicFeatures.Interval, record)
end

function isfilled(record::Record)
    return !isempty(record.filled)
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = record1.filled
        r2 = record2.filled
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    end

    return isfilled(record1) == isfilled(record2) == false
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.nome,
        record.chrom,
        record.chromstart,
        record.chromend,
        record.name,
        record.readnum,
        record.strand,
		record.cg_rle,
		record.gc_rle
	)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "    chromosome: ", chrom(record))
        println(io, "         start: ", chromstart(record))
        println(io, "           end: ", chromend(record))
        println(io, "     name/read: ", "$(name(record))/$(readnum(record))")
        println(io, "        strand: ", strand(record))
        println(io, "CpG RLE string: ", cg_rle(record))
        if is_nome(record)
            println(io)
			print(io, "GpC RLE string: ", gc_rle(record))
        end
    else
        print(io, " <not filled>")
    end
end


# Accessor functions
# ------------------

"""
    chrom(record::Record)::String

Get the chromosome name of `record`.
"""
function chrom(record::Record)::String
    checkfilled(record)
    return String(record.data[record.chrom])
end

function BioGenerics.seqname(record::Record)
    return chrom(record)
end

function BioGenerics.hasseqname(record::Record)
    return isfilled(record)
end

"""
    chromstart(record::Record)::Int

Get the starting position of `record`.

Note that the first base is numbered 1.
"""
function chromstart(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.chromstart) + 1
end

function haschromstart(record::Record)
    return isfilled(record)
end

function BioGenerics.leftposition(record::Record)
    return chromstart(record)
end

function BioGenerics.hasleftposition(record::Record)
    return haschromstart(record)
end

"""
    chromend(record::Record)::Int

Get the end position of `record`.
"""
function chromend(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.chromend)
end

function haschromend(record::Record)
    return isfilled(record)
end

function BioGenerics.rightposition(record::Record)
    return chromend(record)
end

function BioGenerics.hasrightposition(record::Record)
    return haschromend(record)
end

"""
    name(record::Record)::String

Get the name of `record`.
"""
function name(record::Record)::String
    checkfilled(record)
    return String(record.data[record.name])
end

"""
    readnum(record::Record)::Int

Get the read number.
"""
function readnum(record::Record)::Char
    checkfilled(record)
    return Char(record.data[record.readnum])
end

"""
    strand(record::Record)::GenomicFeatures.Strand

Get the strand of `record`.
"""
function strand(record::Record)::GenomicFeatures.Strand
    checkfilled(record)
    return convert(GenomicFeatures.Strand, Char(record.data[record.strand]))
end

function GenomicFeatures.strand(record::Record)
    return strand(record)
end

"""
	cg_rle(record::Record)::String

Get the CpG run-length encoded string from `record`.
"""
function cg_rle(record::Record)::String
	checkfilled(record)
	return String(record.data[record.cg_rle])
end

"""
	gc_rle(record::Record)::String

Get the GpC run-length encoded string from `record`.
"""
function gc_rle(record::Record)::String
	checkfilled(record)
	if !hasgc_rle(record)
		missingerror(:gc_rle)
	end
	return String(record.data[record.gc_rle])
end

function is_nome(record::Record)
	return Bool(record.nome)
end


function unsafe_parse_byte(data::Vector{UInt8}, range::UnitRange{Int})
    val::UInt8 = 0x00
    for i in range
        val = val * 0x0a + (data[i] - UInt8('0'))
    end
    return val
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled EPIREAD record"))
    end
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int}) where T <: Signed
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end