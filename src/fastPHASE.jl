module fastPHASE

using Pkg.Artifacts

export fastPHASE_EXE

if Sys.isapple()
    fastPHASE_EXE = joinpath(artifact"fastPHASE", "fastPHASE")
    run(`chmod +x $fastPHASE_EXE`) # fastPHASE executable by default is read-only
elseif Sys.islinux()
    fastPHASE_EXE = joinpath(artifact"fastPHASE", "fastPHASE")
    run(`chmod +x $fastPHASE_EXE`) # fastPHASE executable by default is read-only
else
    error("fastPHASE only supports Linux and MacOS")
end

end # module
