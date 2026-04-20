using Plots
using LinearAlgebra
using SpecialFunctions
using LaTeXStrings

directory = """/Users/alec/Desktop/S26/Spaepen+Weitz/\
Project/data/example_data/"""
filepath = directory*"example_slice.xyz"

σ = 1 # particle diameter

σ_c = 5σ # cell size: hard spheres would be >1σ


readSliceXYZ(filepath) = (
    # function to read slice file and return array of particle positions
        read(filepath,String)
        |> (y -> replace(y,"\n" => " ; ")) 
        |> (y -> string("[",y,"]"))
        |> Meta.parse
        |> eval
    )

data = readSliceXYZ(filepath)
# recenter positions
data .-= minimum(data,dims=1)
data = data ./ σ # nondimensionalize to radius

function createNeighborLists(positions, σ_c) 
    #=
    take N x 3 array and convert to N x n array, with double-counting,
    where n is number of neighbors within distance σ_c
    =#

    N = size(data,1)
    neighbors = [Int[] for _ in 1:N]

    for (i,j) in ((i,j) for i in 1:N for j in i+1:N)
        distance = (positions[i,:]) - (positions[j,:])
        if norm(distance) < σ_c
            push!(neighbors[i],j)
            push!(neighbors[j],i)
        end
            # println("Distance: $(string(norm(distance))), σ_c: $(σ_c), $(norm(distance) < σ_c)")
    end
    return neighbors
end

neighborList = createNeighborLists(data, σ_c) # O(N^2) EFFICIENCY (BAD)

neighborString = (
    # convert neighborlist into string for saving
    neighborList
    |> (y -> map(x -> string.(x), y))
    |> (y -> map(x -> x.*" ", y))
    |> (y -> map(x -> push!(vec(x), "\n"), y))
    |> Iterators.flatten
    |> join
)

open(directory*"neighbor_list.txt", "w") do file
    write(file, neighborString)
end


small_length_cutoff = 10σ

data_small = data[map(norm,eachrow(data)) .< small_length_cutoff,:] 

d = 3
N = size(data_small,1)
box_volume = small_length_cutoff^3



# Pair Correlation

differences(positions) = (
    # antisymmetric matrix with all differences
    reshape(positions, (N,1,d)) .- reshape(positions, (1,N,d))
)

distances(positions) = (
    # function to return N x N pair distances
    # currently crashes computer when using full dataset
    differences(positions)
    |> (y -> reshape(y, (N*N,d)))
    |> eachrow
    |> (y -> map(norm,y))
    |> (dist -> reshape(dist, (N,N)))
)

sVolume(r,d) = (π^(d/2)) * (r^d) / gamma((d/2) + 1)

#pair correlation values
discrete_pcf(dists, box_length, buckets, N, d) = (
    (box_length / buckets)
    |> (dr -> (collect(range(0, box_length, step=dr)), dr))
    |> (t -> map(r -> begin
            shell_vol = sVolume(r + t[2], d) - sVolume(r, d)
            count(x -> r < x <= r + t[2], dists) / ((N / box_length^d) * N * shell_vol)
        end, t[1]))
)


plot(discrete_pcf(vec(distances(data_small)), small_length_cutoff, 500, N, d))

xlims!(0,500σ)
xlabel!(L"\textrm{r~}\sigma"* " ($(σ))")
title!(L"$\textrm{Pair~correlation~function~} g(r)\textrm{,~normalized~by~} 4\rho\pi r^2 (\Delta r)$")
plot!(titlefontsize=10,label=[L"$g(r)$"])
plot!(legend=false)
plot!(yformatter=_->"")
plot!(dpi=600)
# histogram(vec(distances(data_small)), bins = 500)
# xlims!(0,2σ_c)



# dataTranspose = copy(data')

# plt = scatter(
#     dataTranspose[1,:],
#     dataTranspose[2,:],
#     dataTranspose[3,:],
#     title="Packing",
#     xlims=(-Inf,5), # focus on small limits
#     ylims=(-Inf,5),
#     zlims=(-Inf,5),
#     )
# display(plt)



# TO DO

# Neighborlists Done!
# PAIR CORRELATION FUNCTION Done!
# 'STEINHARDT' Q6 PARAMETER: SPHERICAL HARMONICS, freud?
# VORONOI VOLUME
# CURVATURE https://graphics.stanford.edu/courses/cs468-10-fall/LectureSlides/05_Diff_Geo.pdf