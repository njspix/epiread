# EPIREAD Reader
# ==========
import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

mutable struct Reader <: BioGenerics.IO.AbstractReader
    state::BioGenerics.Automa.State
    index::Union{Indexes.Tabix,Nothing}

    function Reader(stream::TranscodingStream, index=nothing)
        return new(BioGenerics.Automa.State(stream, 1, 1, false), index)
    end
end

"""
    EPIREAD.Reader(input::IO; index=nothing)
    EPIREAD.Reader(input::AbstractString; index=:auto)

Create a data reader of the EPIREAD file format.

The first argument specifies the data source.
When it is a filepath that ends with *.bgz*, it is considered to be block 
compression file format (BGZF) and the function will try to find a tabix 
index file (<filename>.tbi) and read it if any.
See <http://www.htslib.org/doc/tabix.html> for bgzip and tabix tools.

Arguments
---------
- `input`: data source
- `index`: path to a tabix file
"""
function Reader(input::IO; index=nothing)
    if isa(index, AbstractString)
        index = Indexes.Tabix(index)
    end

    stream = TranscodingStreams.NoopStream(input)

    return Reader(stream, index)

end

function Reader(filepath::AbstractString; index=:auto)
    if isa(index, Symbol) && index != :auto
        throw(ArgumentError("invalid index argument: ':$(index)'"))
    end
    if endswith(filepath, ".bgz")
        input = BGZFStreams.BGZFStream(filepath)
        if index == :auto
            index = Indexes.findtabix(filepath)
        end
    else
        input = open(filepath)
    end
    return Reader(input, index = index)
end

function Base.eltype(::Type{Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.iterate(reader::Reader, nextone = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return copy(nextone), empty!(nextone) # Empty record for inplace reading and reuse of array allocations.
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::GenomicFeatures.Interval)
    if reader.index === nothing
        throw(ArgumentError("index is null"))
    end
    return Indexes.TabixOverlapIterator(reader, interval)
end

const record_machine, file_machine = (function ()
    alt = Automa.RegExp.alt
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    opt = Automa.RegExp.opt

    record = let
        chrom = re"[ -~]*"
        chrom.actions[:enter] = [:pos]
        chrom.actions[:exit] = [:record_chrom]

        chromstart = re"[0-9]+"
        chromstart.actions[:enter] = [:pos]
        chromstart.actions[:exit] = [:record_chromstart]

        chromend = re"[0-9]+"
        chromend.actions[:enter] = [:pos]
        chromend.actions[:exit] = [:record_chromend]

        name = re"[ -~]*"
        name.actions[:enter] = [:pos]
        name.actions[:exit] = [:record_name]

        readnum = re"[12]+"
        readnum.actions[:enter] = [:record_readnum] # single byte

        strand = re"[+\-.?]"
        strand.actions[:enter] = [:record_strand] #Note: single byte.

        cg_rle = re"[FDPMUATCGatcgx0-9]+"
        cg_rle.actions[:enter] = [:pos]
        cg_rle.actions[:exit] = [:record_cg_rle]

		gc_rle = re"[FDPOSATCGatcgx0-9]+"
        gc_rle.actions[:enter] = [:pos]
        gc_rle.actions[:exit] = [:record_gc_rle]

        cat(
            chrom, '\t',
            chromstart, '\t',
            chromend, '\t',
            name, '\t',
            readnum, '\t',
            strand, '\t',
            cg_rle,
            opt(cat('\t', gc_rle))
		)
    end
    record.actions[:enter] = [:mark]
    record.actions[:exit] = [:record]

    hspace = re"[ \t\v]"

    blankline = rep(hspace)

    comment = re"#.*"

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    file = rep(alt(
        cat(record, newline),
        cat(blankline, newline),
        cat(comment, newline),
    ))

    return map(Automa.compile, (record, file))
end)()

write("EPIREAD.dot", Automa.machine2dot(file_machine))
#=
run(`dot -Tsvg -o EPIREAD.svg EPIREAD.dot`)
=#

const record_actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(),
    :record_chrom => :(record.chrom = (pos:@relpos(p-1))),
    :record_chromstart => :(record.chromstart = (pos:@relpos(p-1))),
    :record_chromend => :(record.chromend = (pos:@relpos(p-1))),
    :record_name => :(record.name = (pos:@relpos(p-1))),
    :record_readnum => :(record.readnum = @relpos(p)),
    :record_strand => :(record.strand = @relpos(p)),
	:record_cg_rle => :(record.cg_rle = (pos:@relpos(p-1))),
	:record_gc_rle => :(record.gc_rle = (pos:@relpos(p-1)); record.nome==1),
    :record => :(record.filled = 1:@relpos(p-1))
)

Automa.Stream.generate_reader(
    :index!,
    record_machine,
    arguments = (:(record::Record),),
    actions = record_actions,
    # context = :(),
    initcode = :(pos = 0),
    # loopcode = :()
    # returncode = :()
) |> eval


const initcode = quote
    pos = 0
    linenum = 0
    found_record=false
    # empty!(record)
    cs, linenum = state
end

const loopcode = quote
    if found_record
        @goto __return__
    end
end

Automa.Stream.generate_reader(
    :readrecord!,
    file_machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = merge(record_actions, Dict(
        :record => quote
            appendfrom!(record.data, 1, data, @markpos, p-@markpos)
            record.filled = 1:(p-@markpos)
            found_record = true
            @escape
        end,
        :countline => :(linenum += 1),
    )),
    initcode = initcode,
    loopcode = loopcode,
    returncode = :(return cs, linenum, found_record)
) |> eval


function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs = index!(stream, record)
    if cs != 0
        throw(ArgumentError("Invalid EPIREAD record. Machine failed to transition from state $(cs)."))
    end
    return record
end

"""
    read!(rdr::Reader, rec::Record)
Read a `Record` into `rec`; overwriting or adding to existing field values.
It is assumed that `rec` is already initialized or empty.
"""
function Base.read!(rdr::Reader, record::Record)

    cs, ln, found = readrecord!(rdr.state.stream, record, (rdr.state.state, rdr.state.linenum))

    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = found

    if found
        return record
    end

    if cs == 0 || eof(rdr.state.stream)
        throw(EOFError())
    end

    throw(ArgumentError("Malformed EPIREAD file record at line $(ln). Machine failed to transition from state $(cs)."))
end