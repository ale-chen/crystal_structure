using Plots
using LinearAlgebra
using SpecialFunctions
using LaTeXStrings
using SphericalHarmonics # Probably should switch to using FastSphericalHarmonics at some point
using Statistics

readSliceXYZ(filepath) = (
    # function to read slice file and return array of particle positions
        read(filepath,String)
        |> (y -> replace(y,"\n" => " ; ")) 
        |> (y -> string("[",y,"]"))
        |> Meta.parse
        |> eval
    )

function readNeighbors_gpt(filepath)
    # chatgpt readslicexyz that doesn't throw a stackoverflow error
    lines = readlines(filepath)
    return [parse.(Int, split(strip(line))) 
            for line in lines if !isempty(strip(line))]
end

directory = """/Users/alec/Desktop/S26/Spaepen+Weitz/\
Project/data/example_data/"""
positions_filepath = directory*"example_slice.xyz"
neighbor_lists = directory*"neighbor_list.txt"

σ = 1 # particle diameter

σ_c = 5σ # cell size: hard spheres would be >1σ

d = 3


poss = readSliceXYZ(positions_filepath)
# recenter positions
poss .-= minimum(poss,dims=1)
poss = poss ./ σ # nondimensionalize to radius

neighbors = readNeighbors_gpt(neighbor_lists)

XYZToSpherical(xyz) = (
    # (x y z) -> (r θ ϕ) == (r polar azimuthal)
    sqrt(xyz⋅xyz) |>
        (r -> [r,
        acos(xyz[3]/r),
        sign(xyz[2])*acos(xyz[1]/sqrt(xyz[[1, 2]]⋅xyz[[1, 2]]))
        ])
)

# Spherical Harmonics

function getBondList(nbr_list,r_list)
    bond_lists = Vector{Array{Float64}}(undef, size(nbr_list,1))
    for i in 1:size(nbr_list,1)
        bond_lists[i] = (
            r_list[nbr_list[i],:] 
                |> (y -> y.-reshape(r_list[i,:],(1,d)))
                |> (y -> 0.5y) # halfway between particle i and particle j
                |> (y -> mapslices(XYZToSpherical, y; dims=2))
        )
    end
    return bond_lists
end

# Calculate bond positions from neighborlist
# bond positions are (1/2) * r_i - r_j for particle r_i
# as per Ernst, Nagel, Grest (1991)
bonds = getBondList(neighbors,poss)

# Sum Ylm
# maybe figure out fastsphericalharmonics?

function calculateQList(bond_list, l_max)
    # r weighting?? Needs to drop off somehow
    #=
    Q_l = \sqrt{
        4 pi / (2l + 1) *
        \Sum_{-l}^l Magnitude(
            (1/ N_bonds) *
            \Sum_{i=1}^N_bonds(
                Y_lm(\theta,\phi)
            )
        )^2
    }
    =#

    N_b = size(bond_list, 1)

    Ylms = [computeYlm(bond[2], bond[3], l_max) for bond in eachrow(bond_list)]
    Q_ls = Array{Float64}(undef, l_max + 1)

    for l in 0:l_max
        sum_m = sum(-l:l) do m
            avg_Ylm = sum(Ylm[(l, m)] for Ylm in Ylms) / N_b
            abs2(avg_Ylm)
        end

        Q_ls[l + 1] = sqrt(4π / (2l + 1) * sum_m)
    end

    return Q_ls
end

Q_ls = [calculateQList(particle,8) for particle in bonds]

Q_averages = mean(stack(Q_ls), dims=2)
histogram([q[7] for q in Q_ls])
xlabel!(L"Q_{l = 6}")
ylabel!("Counts")