function fastphase()
    return nothing # todo
end

"""
    fastphase_estim_param(xdata; ...)

Runs fastPHASE to estimate `r`, `θ`, `α` parameters. Will run different initial
EM starts in parallel if Julia is started with multiple threads. `θ` (emission)
probabilities) values will automatically be flipped so that `θ[i, j]` is the
probability of observing allele A2 (usually major) in the PLINK bam file at
SNP `i` of haplotype `j`.

# Inputs
`xdata`: A `String` for binary PLINK file (without `.bed/bim/fam` extensions) 
    or a `SnpData` (see SnpArrays.jl). 

# Optional inputs
`n`: Number of samples used to fit HMM in fastPHASE. Defaults to sample size in `xdata`.
`T`: Number of different initial conditions for EM. Different initial conditions
    will be run in parallel in `Threads.nthreads()` number of threads.
`K`: Number of haplotype clusters. Defaults to 12
`C`: Number of EM iterations before convergence. Defaults to 10.
`outfile`: Extension of output alpha, theta, and r file names. Defaults to `fastphase_out`
`outdir`: Output directory. By default all output will be stored in new folder in 
    `knockoffs` in the current directory.
`fastphase_infile`: Filename of fastPHASE's input, which is the decompressed PLINK
    genotypes readable by fastPHASE. Defaults to `fastphase.inp`. 
"""
function fastphase_estim_param(
    xdata::SnpData;
    n::Int = size(xdata.snparray, 1),
    T::Int = Threads.nthreads(),
    K::Int = 12,
    C::Int = 10,
    outfile::AbstractString = "fastphase_out",
    outdir::AbstractString = "knockoffs",
    fastphase_infile::AbstractString = joinpath(outdir, "fastphase.inp"),
    )
    isdir(outdir) || mkdir(outdir)
    x = xdata.snparray
    n ≤ size(x, 1) || error("n must be smaller than the number of samples!")
    sampleid = xdata.person_info[!, :iid]
    # create input format for fastPHASE software
    p = size(x, 2)
    open(fastphase_infile, "w") do io
        println(io, n)
        println(io, p)
        for i in 1:n
            println(io, "ID ", sampleid[i])
            # print genotypes for each sample on 2 lines. The "1" for heterozygous
            # genotypes will always go on the 1st line.
            for j in 1:p
                if x[i, j] == 0x00
                    print(io, 0)
                elseif x[i, j] == 0x02 || x[i, j] == 0x03
                    print(io, 1)
                else
                    print(io, '?')
                end
            end
            print(io, "\n")
            for j in 1:p
                if x[i, j] == 0x00 || x[i, j] == 0x02
                    print(io, 0)
                elseif x[i, j] == 0x03
                    print(io, 1)
                else
                    print(io, '?')
                end
            end
            print(io, "\n")
        end
    end
    # run fastPHASE, each thread runs a different start
    Threads.@threads for i in 1:T
        tmp_out = joinpath(outdir, "tmp$(i)")
        seed = i
        s = 0 # no simulation
        U = 0 # no simulation
        H = 0
        run(`$fastPHASE_EXE -T1 -S$seed -K$K -C$C -o$tmp_out -Pp -H$H -s$s -U$U $fastphase_infile`)
    end
    # aggregate results into final output
    r, θ, α = zeros(p), zeros(p, K), zeros(p, K)
    for i in 1:T
        rtmp, θtmp, αtmp = process_fastphase_output(joinpath(outdir, "tmp$(i)"); T=1)
        r .+= rtmp
        θ .+= θtmp
        α .+= αtmp
    end
    r ./= T
    θ ./= T
    α ./= T
    α ./= sum(α, dims = 2) # normalize rows to sum to 1
    # save averaged results
    writedlm(joinpath(outdir, outfile * "_rhat.txt"), r, ' ')
    writedlm(joinpath(outdir, outfile * "_thetahat.txt"), θ, ' ')
    writedlm(joinpath(outdir, outfile * "_alphahat.txt"), α, ' ')
    # save orichar file, which determines which allele was "allele 1" in θ file
    open(joinpath(outdir, outfile * "_origchars"), "w") do io
        println(io, p)
        for i in 1:p
            println(io, "2\t01") # since all theta are flipped, all origchar will be 01
        end
    end
    # clean up
    # for i in 1:T
    #     rm(joinpath(outdir, "tmp$(i)_rhat.txt"), force=true)
    #     rm(joinpath(outdir, "tmp$(i)_thetahat.txt"), force=true)
    #     rm(joinpath(outdir, "tmp$(i)_alphahat.txt"), force=true)
    #     rm(joinpath(outdir, "tmp$(i)_hapsfrommodel.out"), force=true)
    #     rm(joinpath(outdir, "tmp$(i)_hapguess_switch.out"), force=true)
    #     rm(joinpath(outdir, "tmp$(i)_origchars"), force=true)
    #     rm(joinpath(outdir, "tmp$(i)_sampledHgivG.txt"), force=true)
    # end
    return r, θ, α
end

fastphase_estim_param(xdata::AbstractString; args...) = 
    fastphase_estim_param(SnpData(plinkfile); args...)

"""
    process_fastphase_output(filename; [T])

Reads r, θ, α into memory, averaging over `T` simulations. θ (emission probabilities)
must be flipped sometimes depending on which allele was defined as "allele 1".
fastPHASE simply uses whichever allele was observed first in sample 1 haplotype 1.
Thus, in sample 1, genotypes that start with "10" or "1?" must be flipped, and
this info is provided in the "_origchars" file. 

# Inputs
`filename`: Path to the `rhat.txt`, `thetahat.txt`, `alphahat.txt` and
    `_origchars` files are stored. E.g. Use `out`
    if your output files are `out_rhat.txt`, `out_thetahat.txt`, `out_alphahat.txt`
`T`: Number of different runs excuted by fastPHASE. This is the number of different 
    initial conditions used for EM algorithm. All files `rhat.txt`, `thetahat.txt`,
    `alphahat.txt` would therefore have `T × p` rows
"""
function process_fastphase_output(
    filename::AbstractString;
    T::Int=1,
    )
    # read full data 
    rfile = "$(filename)_rhat.txt" # T*p × 1
    θfile = "$(filename)_thetahat.txt" # T*p × K
    αfile = "$(filename)_alphahat.txt" # T*p × K
    charfile = "$(filename)_origchars" # p × 2 file
    isfile(rfile) && isfile(θfile) && isfile(αfile) && isfile(charfile) ||
        error("Files not found!")
    r_full = readdlm(rfile, comments=true, comment_char = '>', header=false)
    θ_full = readdlm(θfile, comments=true, comment_char = '>', header=false)
    α_full = readdlm(αfile, comments=true, comment_char = '>', header=false)
    flip_idx = flip_θ_index(charfile)
    # compute averages across T simulations as suggested by Scheet et al 2006
    p = Int(size(r_full, 1) / T)
    K = size(θ_full, 2)
    r, θ, α = zeros(p), zeros(p, K), zeros(p, K)
    for i in 1:T
        rows = (i - 1) * p + 1:p*i
        r .+= @view(r_full[rows])
        α .+= @view(α_full[rows, :])
        # make sure to flip θ values if the allele 1 and 2 was flipped by fastPHASE
        θtmp = @view(θ_full[rows, :])
        for i in eachindex(flip_idx)
            if flip_idx[i]
                for j in 1:K
                    θtmp[i, j] = 1 - θtmp[i, j]
                end
            end
        end
        θ .+= θtmp
    end
    r ./= T
    θ ./= T
    α ./= T
    α ./= sum(α, dims = 2) # normalize rows to sum to 1
    return r, θ, α
end

function flip_θ_index(charfile::AbstractString)
    cfile = CSV.read(charfile, DataFrame, delim='\t', header=false, skipto=2)
    flip_idx = falses(size(cfile, 1))
    for i in 1:length(flip_idx)
        if cfile[i, 2] == "10" || cfile[i, 2] == "1?" || cfile[i, 2] == 10
            flip_idx[i] = true
        end
    end
    return flip_idx
end

# function merge_knockoffs_with_original(
#     xdata::SnpData,
#     x̃data::SnpData,
#     des::AbstractString
#     )
#     d = Dict{AbstractString, SnpData}("original"=>xdata, "knockoff"=>x̃data)
#     write_plink(des, merge_plink(d))
# end

function merge_knockoffs_with_original(
    xdata::SnpData,
    x̃data::SnpData;
    des::AbstractString = "knockoff"
    )
    n, p = size(xdata)
    x, x̃ = xdata.snparray, x̃data.snparray
    xfull = SnpArray(des * ".bed", n, 2p)
    original, knockoff = sizehint!(Int[], p), sizehint!(Int[], p)
    for i in 1:p
        # decide which of original or knockoff SNP comes first
        orig, knoc = rand() < 0.5 ? (2i - 1, 2i) : (2i, 2i - 1)
        copyto!(@view(xfull[:, knoc]), @view(x̃[:, i]))
        copyto!(@view(xfull[:, orig]), @view(x[:, i]))
        push!(original, orig)
        push!(knockoff, knoc)
    end
    # copy fam files
    cp(xdata.srcfam, des * ".fam", force=true)
    # copy bim file, knockoff SNPs end in ".k"
    new_bim = copy(xdata.snp_info)
    empty!(new_bim)
    for i in 1:p
        if original[i] < knockoff[i]
            push!(new_bim, xdata.snp_info[i, :])
            push!(new_bim, x̃data.snp_info[i, :])
        else
            push!(new_bim, x̃data.snp_info[i, :])
            push!(new_bim, xdata.snp_info[i, :])
        end
    end
    CSV.write(des * ".bim", new_bim, delim='\t', header=false)
end
