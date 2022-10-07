include("./EPIREAD.jl")
using SparseArrays
allsame(x) = all(y->y==first(x),x)
struct IndelError <: Exception
    msg::String
end

const lookup_str = "FPxMUOSACGTNacgtnD_"
const lookup_matrix = [
    'F' 'F' 'x' 'M' 'U' 'O' 'S' 'A' 'C' 'G' 'T' 'N' 'a' 'c' 'g' 't' 'n' 'D' 'F';
    'F' 'P' 'x' 'M' 'U' 'O' 'S' 'A' 'C' 'G' 'T' 'N' 'a' 'c' 'g' 't' 'n' 'D' 'P';
    'x' 'x' 'x' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'X' 'X' 'X' 'X' 'X' 'X' 'x';
    'M' 'M' 'Q' 'M' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'M' 'X' 'X' 'X' 'X' 'X' 'X' 'M';
    'U' 'U' 'Q' 'Q' 'U' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'U' 'X' 'X' 'X' 'X' 'X' 'X' 'U';
    'O' 'O' 'Q' 'Q' 'Q' 'O' 'Q' 'Q' 'Q' 'Q' 'Q' 'O' 'X' 'X' 'X' 'X' 'X' 'X' 'O';
    'S' 'S' 'Q' 'Q' 'Q' 'Q' 'S' 'Q' 'Q' 'Q' 'Q' 'S' 'X' 'X' 'X' 'X' 'X' 'X' 'S';
    'A' 'A' 'Q' 'Q' 'Q' 'Q' 'Q' 'A' 'Q' 'Q' 'Q' 'A' 'X' 'X' 'X' 'X' 'X' 'X' 'A';
    'C' 'C' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'C' 'Q' 'Q' 'C' 'X' 'X' 'X' 'X' 'X' 'X' 'C';
    'G' 'G' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'G' 'Q' 'G' 'X' 'X' 'X' 'X' 'X' 'X' 'G';
    'T' 'T' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'Q' 'T' 'T' 'X' 'X' 'X' 'X' 'X' 'X' 'T';
    'N' 'N' 'Q' 'M' 'U' 'O' 'S' 'A' 'C' 'G' 'T' 'N' 'X' 'X' 'X' 'X' 'X' 'X' 'N';
    'a' 'a' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'a' 'Q' 'Q' 'Q' 'Q' 'X' 'a';
    'c' 'c' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'Q' 'c' 'Q' 'Q' 'Q' 'X' 'c';
    'g' 'g' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'Q' 'Q' 'g' 'Q' 'Q' 'X' 'g';
    't' 't' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'Q' 'Q' 'Q' 't' 'Q' 'X' 't';
    'n' 'n' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'Q' 'Q' 'Q' 'Q' 'n' 'X' 'n';
    'D' 'D' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'X' 'D' 'D';
    'F' 'P' 'x' 'M' 'U' 'O' 'S' 'A' 'C' 'G' 'T' 'N' 'a' 'c' 'g' 't' 'n' 'D' '_'
]

function get_leftmost_chromend(cache::Vector{EPIREAD.Record})::Int64
    if length(cache) > 0
        return minimum(EPIREAD.chromend.(cache))
    else
        return Int64(0)
    end
end

# Find any reads in the cache that are too far to the left, or on a different chr than provided
function clean_record_cache!(cache::Vector{EPIREAD.Record}, max_end::Int64, new_chr::String)::Vector{EPIREAD.Record}
    reads_to_analyze = Vector{EPIREAD.Record}()
    for (i, record) in enumerate(cache)
        if (EPIREAD.chromend(record) < max_end) || (EPIREAD.chrom(record) ≠ new_chr)
            push!(reads_to_analyze, popat!(cache, i))
        end
    end
    return reads_to_analyze
end

function parse_rle(rle::String)
	reps = parse.(Int16, map(x-> x == "" ? "1" : x, split(rle, r"[[:alpha:]]")[2:end]))
	bases = only.(split(replace(rle, r"[0-9]" => ""), ""))
	return zip(bases, reps)
end 

#Note that if we call EPIREAD.chromstart(record) to get the start position of the read, it is 1-based, end-inclusive (*on-base*). Try to match this here. 
function get_feature_pos_in_read(rle::String; mode = 'M')
	J = Vector{Int64}()
	V = Vector{Int}()
	counter = 0

    if mode == 'M'# methylation + SNP
        ignore_base = "FPx_OSQ" # sequencing errors and open/shut don't indicate different epialleles
        feature_base = "MUACGTN"
    end

    try
        for (base, rep) in parse_rle(rle)
            if occursin(base, ignore_base) 
                counter += rep
            elseif occursin(base, feature_base) 
                for _ in 1:rep
                    counter += 1
                    push!(J, copy(counter))
                    push!(V, Int(base))
                end
            elseif occursin(base, "acgtn") # insertion base
                # handle this better in the future, but for right now keep everything aligned with reference
            elseif base == 'D'
                # handle this better in the future, but for right now keep everything aligned with reference
                counter += rep
            else
                @assert base == 'X'
                throw(IndelError("insertion/deletion mismatch in overlapping reads"))
            end
        end
    catch error
        if isa(error, IndelError)
            return(J, V, IndelError)
        else
            throw(e)
        end
    end

	return (J, V, nothing)
end

function get_feature_pos_absolute(rle::String, read_start::Int64)
	J, V, error = get_feature_pos_in_read(rle)
	return (J.+read_start.-1, V, error) # we have to subtract one, as this is 1-start, fully closed
end

function update_vecs(master_I::Vector{Int64}, master_J::Vector{Int64}, master_V::Vector{Int}, J::Vector{Int64}, V::Vector{Int})
	read_num = length(master_I) == 0 ? Int64(1) : master_I[end]+1
	append!(master_I, fill(read_num, length(J)))
	append!(master_J, J)
	append!(master_V, V)
	return nothing
end

function analyze_read(read::EPIREAD.Record, master_I::Vector{Int64}, master_J::Vector{Int64}, master_V::Vector{Int})
	# println("analyzing read $(EPIREAD.name(read))/$(EPIREAD.readnum(read)) ($(EPIREAD.chrom(read))): $(EPIREAD.cg_rle(read))")
	J, V, error = get_feature_pos_absolute(EPIREAD.cg_rle(read), EPIREAD.chromstart(read))
    # This method is only used for analyzing single reads, so we should never get an IndelError. If we do get an error, we need to throw it:
    if isa(error, Exception)
        throw(error)
    end
	update_vecs(master_I, master_J, master_V, J, V)
	return nothing
end

function analyze_read(rle::String, chromstart::Int64, read_name::String, master_I::Vector{Int64}, master_J::Vector{Int64}, master_V::Vector{Int}, indel_error_reads::Vector{String})
	# println("analyzing read $(EPIREAD.name(read))/$(EPIREAD.readnum(read)) ($(EPIREAD.chrom(read))): $(EPIREAD.cg_rle(read))")
	J, V, error = get_feature_pos_absolute(rle, chromstart)
    # this method is used for analyzing merged paired reads, so we might get an IndelError. If we do, add the read name to the IndelError Vector
    if isa(error, Exception)
        push!(indel_error_reads, read_name)
    end
	update_vecs(master_I, master_J, master_V, J, V)
	return nothing
end

function are_overlapped(read1::EPIREAD.Record, read2::EPIREAD.Record)
	rightmost_start = maximum([EPIREAD.chromstart(read1), EPIREAD.chromstart(read2)])
	leftmost_end = minimum([EPIREAD.chromend(read1), EPIREAD.chromend(read2)])
	return leftmost_end ≥ rightmost_start # since were' in 1-based, fully-closed coords (aka on-base coords), reads are overlapped if the start coord of one equals the end coord of the other
end

function expand_rle(rle::String)
	temp = Vector{Char}()
	for (base, rep) in parse_rle(rle)
		for i in 1:rep
			push!(temp, base)
		end
	end
	return(String(temp))
end

function compare_bases(base1::Char, base2::Char)
	if base1==base2
		return base1
	else
		return lookup_matrix[findfirst(base1, lookup_str), findfirst(base2, lookup_str)]
	end
end

function merge_rle(rle1::String, rle2::String)
	merged = Vector{Char}()
	for (base1, base2) in zip(rle1, rle2)
		push!(merged, compare_bases(base1, base2))
	end
	return(String(merged))
end

function merge_rle(read1::EPIREAD.Record, read2::EPIREAD.Record)
	# We are only guaranteed that the reads overlap by ≥ 1 base. 
	left_read, right_read = (EPIREAD.chromstart(read1) ≤ EPIREAD.chromstart(read2)) ? (read1, read2) : (read2, read1)
	left_rle, right_rle = expand_rle.(EPIREAD.cg_rle.([left_read, right_read]))

	# We don't need to pad the left side of the left_read, but we might need to pad all 3 other sides:
	# We use the character '_' (underscore) as a pad.
	# First, line up the left sides:
	left_pad = EPIREAD.chromstart(right_read) - EPIREAD.chromstart(left_read)
	if left_pad > 0
		right_rle = "_"^left_pad*right_rle
	end
	# Now line up the right sides:
	long_rle, short_rle = length(right_rle) ≥ length(left_rle) ? (right_rle, left_rle) : (left_rle, right_rle)
	diff = length(long_rle) - length(short_rle)
	if diff > 0
		# pad right end of the shorter one
		short_rle = short_rle*"_"^abs(diff)
	end
	return (merge_rle(long_rle, short_rle), EPIREAD.chromstart(left_read))
end

function analyze_record_pair(read1::EPIREAD.Record, read2::EPIREAD.Record, master_I::Vector{Int64}, master_J::Vector{Int64}, master_V::Vector{Int}, indel_error_reads::Vector{String})
	if !allsame(EPIREAD.chrom.((read1, read2)))
		for each in (read1, read2)
			analyze_read(each, master_I, master_J, master_V)
		end
	else 
		if !are_overlapped(read1, read2)
			for each in (read1, read2)
				analyze_read(each, master_I, master_J, master_V)
			end
		else
			# handle overlapping reads
			rle, chromstart = merge_rle(read1, read2)
			analyze_read(rle, chromstart, EPIREAD.name(read1), master_I, master_J, master_V, indel_error_reads)
		end 
	end
end

function analyze_chr(reader::EPIREAD.Reader, max_isize::Int, leftover_record::EPIREAD.Record)

    i = 0 # counter for number of reads analyzed
    I,J,V = Vector{Int64}(), Vector{Int64}(), Vector{Int}() # vectors for sparse matrix
    record_cache = Vector{EPIREAD.Record}() # Reads that could possibly be paired with an upcoming read
    indel_error_reads = Vector{String}() # Read names that have an indel mismatch in their overlapping portions
    last_chr = String("") # Current chromosome being analyzed

    # Put the 'leftover' record from the previous chromosome into the record_cache for analysis
    if EPIREAD.isfilled(leftover_record)
        push!(record_cache, leftover_record)
        last_chr = EPIREAD.chrom(leftover_record)
        i = 1
    end

    # Pre-allocate record
    record = EPIREAD.Record()
    while !eof(reader)
        empty!(record)
        read!(reader, record) 
        name = EPIREAD.name(record)
        chromstart = EPIREAD.chromstart(record)
        current_chr = EPIREAD.chrom(record)

        # Are we on a new chromosome?
        if (last_chr ≠ "" && last_chr ≠ current_chr)
            # Clean out the cache
            for each in record_cache
                analyze_read(each, I, J, V)
            end
            println(Base.stderr, "Completed analyis of $i reads in $(last_chr)")
            return (sparse(I,J,V), indel_error_reads, last_chr, record)
        else
            i += 1
            # Find the furthest-left *end* coordinate of all reads in the cache (if it is too far away from the beginning of this read, we'll never be able to pair to it)
            leftmost_chromend = get_leftmost_chromend(record_cache)
            # Test as above, and also check if we're on a new chr.
            if (leftmost_chromend < (chromstart - max_isize)) && (leftmost_chromend > 0) 
                # If either of the above is true, get all the un-pair-able reads out of the cache and analyze them:
                # println("cleaning cache because min_end $leftmost_chromend is more than $max_isize away from $chromstart or because last-seen $(last_chr) is different from current $(EPIREAD.chrom(record)) (read $(name))")
                for each in clean_record_cache!(record_cache, chromstart - max_isize, EPIREAD.chrom(record))
                    analyze_read(each, I, J, V)
                end
            end
            # Ok, any read still in the cache is fair game to pair with our current read. Let's see if we can find a match:
            cached_names = EPIREAD.name.(record_cache)
            if name in cached_names 
                pair_record = popat!(record_cache, findfirst(isequal(name), cached_names))
                analyze_record_pair(record, pair_record, I, J, V, indel_error_reads)
            # Ok, we couldn't find a match. Save this read for later (maybe its match just hasn't been read yet!)
            else
                push!(record_cache, copy(record))
            end
            last_chr = current_chr
        end
    end

    # Clean out the cache
    for each in record_cache
        analyze_read(each, I, J, V)
    end
    println(Base.stderr, "Completed analyis of $i reads in $(last_chr)")
    return (sparse(I,J,V), indel_error_reads, last_chr, EPIREAD.Record())
end


function parse_infile(filename::String; max_isize = 1000)

    results = Dict{String, SparseMatrixCSC{Int, Int64}}()
    indel_error_reads = Vector{String}()
    leftover_record = EPIREAD.Record()

    reader = EPIREAD.Reader(filename)
    while !eof(reader)
        (sparse_results, new_indel_errors, last_chr, leftover_record) = analyze_chr(reader, max_isize, leftover_record)
        results[last_chr] = sparse_results
        append!(new_indel_errors, indel_error_reads)
    end

    if length(indel_error_reads) > 0
        println(Base.stderr, "WARNING: $(length(indel_error_reads)) reads had mismatching indels and were not analyzed:")
        println(Base.stderr, indel_error_reads)
    end
    close(reader)

    return results
end

#==========================================
Epiallele parsing and counting functions
==========================================#

#= We need to set up this struct to allow for quick comparisons between epireads. 
	If they only differ at 'N' bases, they should be treated as identical =#

    struct Epiallele <: AbstractVector{Int64}
        x::Vector{Int64}
    end
    
    function Base.size(x::Epiallele)
        return size(x.x)
    end
    
    function Base.getindex(x::Epiallele, i::Int64)
        return x.x[i]
    end
    
    function Base.:(==)(x::Epiallele, y::Epiallele) 
        # if the epialleles are different lengths, return false
        if length(x.x) != length(y.x)
            return false
        end
        # compare elementwise; if only difference is at N's (78), return true
        for i in 1:length(x.x)
            if x.x[i] != y.x[i] && x.x[i] != 78 && y.x[i] != 78
                return false
            end
        end
        return true
    end
    
    function Base.isequal(x::Epiallele, y::Epiallele)
        return x == y
    end
    
    # Zero initializer for Epiallele
    function Base.zeros(::Type{Epiallele}, n::Int64)
        return Epiallele(zeros(Int64, n))
    end

function get_unsparse_nonzero_rows(matrix::SparseMatrixCSC{Int64, Int64})
	nzrows = unique(matrix.rowval)
	if !isempty(nzrows)
		return Matrix(matrix)[nzrows,:]
	else 
		return zeros(Int64, (0,0))
	end
end

function matrix_to_epiallele_vector(mat::Matrix{Int64})
	vec = Vector{Epiallele}(undef, size(mat, 1))
	for i in 1:size(mat, 1)
		vec[i] = Epiallele(mat[i,:])
	end
	return vec
end

function diff_from_all_epialleles(needle::Epiallele, haystack::Vector{Epiallele})
	for i in 1:length(haystack)
		if needle == haystack[i]
			return false
		end
	end
	return true
end

function count_epialleles(epiallele_vec::Vector{Epiallele})
	seen = fill(zeros(Epiallele, length(epiallele_vec[1])), length(epiallele_vec))
	unique = 0
	for i in 1:length(epiallele_vec)
		diff_from_all_epialleles(epiallele_vec[i], seen) && (unique += 1)
		seen[i] = epiallele_vec[i]
	end
	return unique
end

function tally_epialles(matrix::SparseMatrixCSC{Int64, Int64}, chr::String; window_size = 4::Int)
    n_bases = size(matrix, 2)
    println("fixedStep\tchrom=$(chr)\tstart=1\tstep=$(window_size)")
    current_bed_line = (0::Int, 0::Int) # (start, n_alleles)
    for i in 1:window_size:n_bases
        if i % 10000 == 1
            println(Base.stderr, "Tallying epialleles on $(chr) at position $(i) of $(n_bases)")
            flush(stdout)
        end
    
        upper_limit = minimum([i+window_size-1, n_bases])
        a = get_unsparse_nonzero_rows(matrix[:,i:upper_limit])
        if !isempty(a)
            n = count_epialleles(matrix_to_epiallele_vector(a))
            if n != current_bed_line[2]
                println("$(chr)\t$(current_bed_line[1])\t$(upper_limit)\t$(n)")
                current_bed_line = (i-1, n)
            end
        end
    end
end

function output_stats(results::Dict{String, SparseMatrixCSC{Int64, Int64}})
    println("track type=bedGraph")
    chrs = collect(keys(results))
    for chr in sort(chrs)
        println(Base.stderr, "Processing $chr")
        tally_epialles(results[chr], chr)
    end
end

function process_file(filename::String; max_isize = 1000)
    results = parse_infile(filename, max_isize = max_isize)
    output_stats(results)
end

process_file("./testing/test_epiread.bed")
