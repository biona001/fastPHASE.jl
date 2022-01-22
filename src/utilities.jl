function fastphase_estim_param(
    xdata::SnpData;
    n::Int = size(xdata.snparray, 1), # number of samples used to fit HMM
    T::Int = Threads.nthreads(), # number of different initial conditions for EM
    K::Int = 12, # number of clusters
    C::Int = 10, # number of EM iterations
    out::AbstractString = "fastphase_out",
    fastphase_infile::AbstractString = "fastphase.inp",
    outdir::AbstractString = pwd()
    )
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
        rtmp, θtmp, αtmp = process_fastphase_output(outdir; T=1, extension="tmp$(i)")
        r .+= rtmp
        θ .+= θtmp
        α .+= αtmp
    end
    r ./= T
    θ ./= T
    α ./= T
    α ./= sum(α, dims = 2) # normalize rows to sum to 1
    # save averaged results
    writedlm(out * "_rhat.txt", r, ' ')
    writedlm(out * "_thetahat.txt", θ, ' ')
    writedlm(out * "_alphahat.txt", α, ' ')
    # clean up
    for i in 1:T
        rm(joinpath(outdir, "tmp$(i)_rhat.txt"), force=true)
        rm(joinpath(outdir, "tmp$(i)_thetahat.txt"), force=true)
        rm(joinpath(outdir, "tmp$(i)_alphahat.txt"), force=true)
        rm(joinpath(outdir, "tmp$(i)_hapsfrommodel.out"), force=true)
        rm(joinpath(outdir, "tmp$(i)_hapguess_switch.out"), force=true)
        rm(joinpath(outdir, "tmp$(i)_origchars"), force=true)
        rm(joinpath(outdir, "tmp$(i)_sampledHgivG.txt"), force=true)
    end
    return r, θ, α
end

"""
    process_fastphase_output(datadir, T; [extension])

Reads r, θ, α into memory, averaging over `T` simulations. 

# Inputs
`datadir`: Directory that the `rhat.txt`, `thetahat.txt`, `alphahat.txt` files are stored
`T`: Number of different runs excuted by fastPHASE. This is the number of different 
    initial conditions used for EM algorithm. All files `rhat.txt`, `thetahat.txt`,
    `alphahat.txt` would therefore have `T × p` rows
`extension`: Name in front of the output fastPHASE files. E.g. Use `extension=out_`
    if your output files are `out_rhat.txt`, `out_thetahat.txt`, `out_alphahat.txt`
"""
function process_fastphase_output(
    datadir::AbstractString;
    T::Int=1,
    extension="out"
    )
    # read full data 
    rfile = joinpath(datadir, "$(extension)_rhat.txt") # T*p × 1
    θfile = joinpath(datadir, "$(extension)_thetahat.txt") # T*p × K
    αfile = joinpath(datadir, "$(extension)_alphahat.txt") # T*p × K
    isfile(rfile) && isfile(θfile) && isfile(αfile) || error("Files not found!")
    r_full = readdlm(rfile, comments=true, comment_char = '>', header=false)
    θ_full = readdlm(θfile, comments=true, comment_char = '>', header=false)
    α_full = readdlm(θfile, comments=true, comment_char = '>', header=false)

    # compute averages across T simulations as suggested by Scheet et al 2006
    p = Int(size(r_full, 1) / T)
    K = size(θ_full, 2)
    r, θ, α = zeros(p), zeros(p, K), zeros(p, K)
    for i in 1:T
        rows = (i - 1) * p + 1:p*i
        r .+= @view(r_full[rows])
        θ .+= @view(θ_full[rows, :])
        α .+= @view(α_full[rows, :])
    end
    r ./= T
    θ ./= T
    α ./= T
    α ./= sum(α, dims = 2) # normalize rows to sum to 1
    return r, θ, α
end
