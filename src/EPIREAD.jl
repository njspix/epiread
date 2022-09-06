module EPIREAD

import Automa
import Automa.RegExp: @re_str
import BGZFStreams
import BioGenerics

using GenomicFeatures
using Indexes
using TranscodingStreams

include("record.jl")
include("reader.jl")

end
