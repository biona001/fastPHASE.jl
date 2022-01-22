module fastPHASE

using Pkg.Artifacts
using SnpArrays
using DelimitedFiles

export fastPHASE_EXE

if Sys.isapple()
    fastPHASE_EXE = joinpath(artifact"fastPHASE", "fastPHASE")
elseif Sys.islinux()
    fastPHASE_EXE = joinpath(artifact"fastPHASE", "fastPHASE")
else
    error("fastPHASE only supports Linux and MacOS")
end
run(`chmod +x $fastPHASE_EXE`) # fastPHASE executable by default is read-only

include("utilities.jl")

end # module
