using ArgParse
include("./parse_file.jl")

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_tbl s begin
		"--window", "-w"
			help="window size for epialele calculation (default 10 bp)", default=10
		"--insert-size", "-i"
			help="consider reads unpaired if greater than this distance (default: 1000 bp)"
			default=1000
		"--mode", "-m"
			help="mode of epialele calculation (default: M, methylation+SNP, \n\toptions: M, S (SNP only), Q (sequencing errors only))"
			default='M'
			arg_type = Char
		"<file>"
			help="input file"
	end

	return parse_args(s)
end

function main()
	args = parse_commandline()
	process_file(
		args["<file>"],
		max_isize = args["--insert-size"],
		window_size = args["--window"],
		mode = args["--mode"]
	)
end

main()