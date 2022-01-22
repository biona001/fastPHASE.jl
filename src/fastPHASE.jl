module fastPHASE

using Pkg.Artifacts

export fastPHASE_EXE

if Sys.isapple()
    fastPHASE_EXE = joinpath(artifact"fastPHASE", "fastPHASE")
    chmod(fastPHASE_EXE, 755)
elseif Sys.islinux()
    fastPHASE_EXE = joinpath(artifact"fastPHASE", "fastPHASE")
    chmod(fastPHASE_EXE, 755)
else
    error("fastPHASE only supports Linux and MacOS")
end

end # module
