module TPLib

include("Tropical.jl")
using .Tropical
export Infinite, SemiRing, MaxPlus, MinPlus, isinfinite, isplusinfinite, isminusinfinite, elt, ∞
export compute_ext_rays, compute_ext_rays_polar, compute_halfspaces, compute_tropical_complex, compute_tangent_hypergraph

"""
    number_to_string(x::Number)
Returns a string representing the number in the TPLib format.
"""
function number_to_string(x::Number)
    return typeof(x) <: Rational ? string([x.num, "/", x.den]...) : string(x)
end

"""
    vector_to_string(V::Vector{<:SemiRing})
Returns the vector V as a string in the TPLib format, converting the tropical zero to ±oo.
"""
function vector_to_string(V::Vector{T}) where {T<:SemiRing}
    # set "tropical 0" element depending on the SemiRing supplied
    T <: MaxPlus ? trop_zero = "-oo" : trop_zero = "+oo"
    S = []
    push!(S,"[")
    for (i, v) in enumerate(V)
        i > 1 && push!(S,",")
        isinfinite(v) ? push!(S, trop_zero) : push!(S,number_to_string(elt(v)))
    end
    push!(S,"]")
    return string(S...)
end

"""
    write_matrix(M::Matrix{<:Real}, io::IO)
Print the matrix M to io in the TPLib format, converting the tropical zero to ±oo.
"""
function write_matrix(io::IO, M::Matrix{<:SemiRing})
    for i in 1:size(M)[1]
        println(io,vector_to_string(M[i,:]))
    end
end

"""
    write_halfspaces(io::IO, H::Matrix{<:SemiRing}, A::Vector{Vector{Int64}})
Print the halfspaces given by H with sectors given by A to io in the TPLib format, converting the tropical zero to ±oo.
"""
function write_halfspaces(io::IO, H::Matrix{T}, A::Vector{Vector{Int64}}) where {T<:SemiRing}
    size(H)[1] == length(A) || error("Mismatch in the number of halfspaces (size(H)[1]) and the number of sectors (length(A)).")  
    
    # set "tropical 0" element depending on the SemiRing supplied
    T <: MaxPlus ? trop_zero = "-oo" : trop_zero = "+oo"

    # print each row of the matrix followed by its sectors in TPLib format
    for i in 1:length(A)
        print(io,"([")
        for (j,x) in enumerate(H[i,:])
            j > 1 && print(io,",")
            isinfinite(x) ? print(io, trop_zero) : print(io,elt(x))
        end
        print(io,"],[")
        for (j,elt) in enumerate(A[i])
            j > 1 && print(io,",")
            print(io,elt-1) # indexes in caml start at 0
        end
        print(io,"])\n")
    end
end
   

"""
    string_to_number(s::String, T::Type{<:SemiRing}
Returns a string representing the number given in TPLib format.
"""
function string_to_number(s::AbstractString, T::Type{<:Number})
    return parse(T,s)
end

function string_to_number(s::AbstractString, T::Type{<:SemiRing{S}}) where {S}
    return T(parse(S,s))
end

function string_to_number(s::AbstractString, T::Type{<:SemiRing{S}}) where {S <: Rational{U}} where {U}
    if '/' in s
        num, den = split(s,'/')
        return T(parse(U,num)//parse(U,den))
    end
    return T(parse(U,s))
end

"""
    string_to_vector(s::String, T::Type{<:SemiRing})
Converts a string in the TPLib format of the type "([_,_,_],[_,_,_],...)" to a vector of vectors of type T.
"""
function string_to_vector(s::AbstractString, T::Type{<:Number})
    A = Vector{Vector{T}}()
    s = s[findfirst('[',s):findlast(']',s)]
    s = split(s,r"(?<=]),")
    for vect in s
        V = Vector{T}()
        vect = strip(vect,['[',']'])
        elts = split(vect, ',')
        for elt in elts
            if elt == "-oo" || elt == "+oo"
                push!(V,T(Infinite()))
            else
                push!(V,string_to_number(elt,T))
            end
        end
        push!(A,V)
    end
    return A
end

"""
    string_to_halfspace(s::String, T::Type{<:SemiRing})
Converts a string in the TPLib format describing a halfspace and its sector to a pair of vectors of type T.
"""
function string_to_halfspace(s::AbstractString, T::Type{<:SemiRing{S}}) where {S<:Number}
    s = s[findfirst('[',s):findlast(']',s)]
    s = split(s,r"(?<=]),")
    h = split(strip(s[1],['[',']']), ',')
    H = map(e -> e == "-oo" || e == "+oo" ? T(Infinite()) : string_to_number(e,T), h)
    a = split(strip(s[2],['[',']']), ',')
    A = map(e -> parse(Int64,e)+1, a) # index in caml start at 0
    return H,A
end

"""
    read_matrix(io::IOStream, T::Type{<:SemiRing}, n::Integer)za
Converts an IOStream composed of a list of vectors in the TPLib format of size n to a matrix M with n columns of type T.
"""
function read_matrix(io::IOStream, T::Type{<:SemiRing{S}}, n::Integer) where {S <: Number}
    M = Matrix{T}(undef,0,n)
    for line in eachline(io)
        M = [M; permutedims(string_to_vector(line,T)[1])]
    end
    return M
end

"""
    read_halfspaces(io::IOStream, T::DataType, n::Integer)
Converts an IOStream composed of a list of halfspaces and theirs sectors in TPLib format to a matrix H of type T describing the halfspaces and a vector of vectors A describing the sectors.
"""
function read_halfspaces(io::IOStream, T::Type{<:SemiRing{S}}, n::Integer) where {S<:Number}
    H = Matrix{T}(undef,0,n)
    A = Vector{Vector{Int64}}()
    for line in eachline(io)
            h,a = string_to_halfspace(line,T)
            H = [H; permutedims(h)]
            push!(A, a)
    end
    return (H,A)
end

"""
    read_complex(io::IOStream, T::Type{<:SemiRing}, n::Integer)
Converts an IOStream composed of a list of vertices followed by a list of cells describing a tropical complex in the TPLib format to a matrix V of type T describing the vertices and a vector of vectors A describing the maximal cells.
"""
function read_complex(io::IOStream, T::Type{<:SemiRing{S}}, n::Integer) where {S<:Number}
    V = Matrix{T}(undef,0,n)
    A = Vector{Vector{Int64}}()
    readline(io) # discard the heading
    for line in eachline(io) # first part contains vertices
        length(line) > 0 || break # the two parts are seperated by a blank line
        V = [V; permutedims(string_to_vector(line,T)[1])]
    end
    readline(io) # discard the heading
    for line in eachline(io) # second part contains cells
        push!(A, string_to_vector(line,Int64)[1].+1) # indexes in caml start at 0
     end
     return (V,A)
end

"""
    read_hypergraph_ineq(io::IOStream, T::Type{<:SemiRing}, n::Integer)
Converts an IOStream composed of a list of hyperedges and inequalities describing a tangent hypergraph to a triplet with the number of vertices of the hypergraph, a vector of vectors describing the hyperedges, and a matrix of type T describing the inequalities.
"""
function read_hypergraph_ineq(io::IOStream, T::Type{<:SemiRing{S}}, n::Integer) where {S<:Number}
    V = Any[]
    M = Matrix{T}(undef,0,2*n)
    nbvert = parse(Int64,split(readline(io),' ')[1])
    
    for line in eachline(io)
        length(line) > 0 || break
        i1, i2, i3, i4 = findall(c -> c == '{' || c == '}', line)
        tail = split(line[i1+2:i2-2],' ')
        head = split(line[i3+2:i4-2],' ')
        push!(V,([parse(Int64,e)+1 for e in tail], [parse(Int64,e)+1 for e in head]))

        line = line[findfirst('[',line)+2:findlast(']',line)-2]
        M = [M; permutedims(string_to_vector(line,T)[1])]
    end
    return nbvert, V, M
end

"""
    read_hypergraph_gen(io::IOStream, T::Type{<:SemiRing}, n::Integer)
Converts an IOStream composed of a list of hyperedges and generators describing a tangent hypergraph to a triplet with the number of vertices of the hypergraph, a vector of vectors describing the hyperedges, and a matrix of type T describing the inequalities.
"""
function read_hypergraph_gen(io::IOStream, T::Type{<:SemiRing{S}}, n::Integer) where {S<:Number}
    V = Any[]
    nbvert = parse(Int64,split(readline(io),' ')[1])
    
    for line in eachline(io)
        length(line) > 0 || break
        i1, i2, i3, i4 = findall(c -> c == '{' || c == '}', line)
        tail = split(line[i1+2:i2-2],' ')
        head = split(line[i3+2:i4-2],' ')
        push!(V,([parse(Int64,e)+1 for e in tail], [parse(Int64,e)+1 for e in head]))
    end
    return nbvert, V
end

"""
    read_hypergraph_halfspace(io::IOStream, T::Type{<:SemiRing}, n::Integer)
Converts an IOStream composed of a list of hyperedges and halfspaces describing a tangent hypergraph to a triplet with the number of vertices of the hypergraph, a vector of vectors describing the hyperedges, and a matrix of type T describing the inequalities.
"""
function read_hypergraph_halfspace(io::IOStream, T::Type{<:SemiRing}, n::Integer)
    V = Any[]
    H = Matrix{T}(undef,0,n)
    A = Vector{Vector{Int64}}()
    nbvert = parse(Int64,split(readline(io),' ')[1])
    
    for line in eachline(io)
        length(line) > 0 || break
        i1, i2, i3, i4 = findall(c -> c == '{' || c == '}', line)
        tail = split(line[i1+2:i2-2],' ')
        head = split(line[i3+2:i4-2],' ')
        push!(V,([parse(Int64,e)+1 for e in tail], [parse(Int64,e)+1 for e in head]))

        line = line[findfirst('[',line)+2:findlast(']',line)-2]
        h, a = string_to_halfspace(line,T)
        H = [H; permutedims(h)]
        push!(A,a)
    end
    return nbvert, V, H, A
end

"""
    @tropicalize(func::Symbol)
Defines a new method for the function func(M::Matrix{<:SemiRing}, n::Integer) which can take as argument a matrix M::Matrix{T<:Number} and an extra argument semiring which specifies if the elements of M should be converted to MaxPlus{T} or MinPlus{T}.
"""
macro tropicalize(func::Symbol)
    quote 
        function $(esc(func))(M::Matrix{<:Number}, n::Integer, semiring=:max)
            if semiring === :max
                if !(eltype(M) <: MaxPlus)
                    S = eltype([m for m in M if !isminusinfinite(m)])
                    M = convert(Matrix{MaxPlus{S}},M)
                end
            elseif semiring === :min
                if !(eltype(M) <: MinPlus)
                    S = eltype([m for m in M if !isplusinfinite(m)])
                    M = convert(Matrix{MinPlus{S}},M)
                end
            else
                error("Invalid value for semiring. Must be :max or :min.")
            end
            return $(esc(func))(M,n)
        end
    end
end

"""
    compute_ext_rays(I::Matrix{<:SemiRing}, n::Integer)
    compute_ext_rays(I::Matrix{<:Number}, n::Integer, semiring=:max)
Computes the set of the extreme rays of a tropical cone given by the inequalities `I` in dimension `n`, see [3]. Each row of `I` contains `2n` coefficients. When written ``[a_{i1} \\, \\dots \\, a_{in} \\, b_{i1} \\, \\dots \\, b_{in}]``, they represent the tropical inequality 
```math
\\max(a_{i1} + x_1, \\dots, a_{in} + x_n) \\geqslant \\max(b_{i1} + x_1, \\dots, b_{in} + x_n)
```
If the coefficients of `I` have type `MaxPlus{T}` or `MinPlus{T}`, then the computations are done in the max-plus or min-plus tropical semirings respectively. Otherwise, you can specify wether the coefficients should be converted to the semiring `:max` or `:min`.
Returns a matrix `G::Matrix{<:SemiRing}` in which every line is a tropical point generating the cone.

# References

[1] X. Allamigeon. Static analysis of memory manipulations by abstract interpretation -- Algorithmics of tropical polyhedra, and application to abstract interpretation. PhD thesis. 

[2] X. Allamigeon, S. Gaubert, E. Goubault. Computing the vertices of tropical polyhedra using directed hypergraphs. Discrete & Computational Geometry, 49(2):247–279, 2013. E-print arXiv:0904.3436v4.
"""
function compute_ext_rays(I::Matrix{T}, n::Integer) where {T <: SemiRing{S}} where {S} 
    input, io_input = mktemp()
    output, io_output = mktemp()
    write_matrix(io_input, I)
    seekstart(io_input)
    
    flags = Vector{String}()
    # flag for MaxPlus/MinPlus
    if T <: MinPlus
        push!(flags, "-min-plus")
    end

    # flag for numerical type
    if S<:Float64
        push!(flags, "-numerical-data", "ocaml_float")
    elseif S<:BigInt
        push!(flags, "-numerical-data", "zarith_int")
    elseif S<:Rational
        push!(flags, "-numerical-data", "zarith_rat")
    elseif !(S<:Int64)
        error("Only supports Ints, Floats, BigInts and Rationals.")
    end

    run(pipeline(`compute_ext_rays $flags $n`, stdin=input, stdout=output))

    G = read_matrix(io_output, T, n)
    close(io_input)
    close(io_output)
    return G
end

@tropicalize compute_ext_rays

"""
    compute_ext_rays_polar(M::Matrix{<:SemiRing}, n::Integer)
    compute_ext_rays_polar(M::Matrix{<:Number}, n::Integer, semiring=:max)
Computes the set of the extreme rays of the polar of a tropical cone given by a generating set `M` in dimension `n`, see [3]. Each row of `M` contains `n` coefficients and represents a generator of the tropical cone. If the coefficients of `M` have type `MaxPlus{T}` or `MinPlus{T}`, then the computations are done in the max-plus or min-plus tropical semirings respectively. Otherwise, you can specify wether the coefficients should be converted to the semiring `:max` or `:min`.
Returns a matrix `::Matrix{<:SemiRing}` in which every row is a tropical point generating the polar cone.

# References

[1] X. Allamigeon, S. Gaubert, and R. D. Katz. Tropical polar cones, hypergraph transversals, and mean payoff games. Linear Algebra Appl., 435(7):1549–1574, 2011. E-print arXiv:1004.2778.
"""
function compute_ext_rays_polar(M::Matrix{T}, n::Integer) where {T<:SemiRing{S}} where {S}
    input, io_input = mktemp()
    output, io_output = mktemp()
    write_matrix(io_input, M)
    seekstart(io_input)


    flags = Vector{String}()
    # flag for MaxPlus/MinPlus
    if T <: MinPlus
        push!(flags, "-min-plus")
    end

    # flag for numerical type
    if S<:Float64
        push!(flags, "-numerical-data", "ocaml_float")
    elseif S<:BigInt
        push!(flags, "-numerical-data", "zarith_int")
    elseif S<:Rational
        push!(flags, "-numerical-data", "zarith_rat")
    elseif !(S<:Int64)
        error("Only supports Ints, Floats, BigInts and Rationals.")
    end

    run(pipeline(`compute_ext_rays_polar $flags $n`, stdin=input, stdout=output))

    G = read_matrix(io_output, T, 2*n)
    close(io_input)
    close(io_output)
    return G
end

@tropicalize compute_ext_rays_polar

"""
    compute_halfspaces(M::Matrix{<:SemiRing}, n::Integer)
    compute_halfspaces(M::Matrix{<:Number}, n::Integer, semiring=:max)
Computes a minimal external representation by means of tropical half-spaces of a tropical cone given by a generating set `M` in dimension `n`, see [4]. The set of generators must be nonempty, ie `M` cannot be empty. It currently handles only the case of generating sets in which every vector has finite entries.
Returns a pair `(H::Matrix{<:SemiRing},A::Vector{Vector{Integer}})` where the rows of `H` are the equations for the tropical half-spaces and the vectors of `A` are the sectors of each half-space. 

# References

[1] X. Allamigeon and R.D. Katz. Minimal external representations of tropical polyhedra. Journal of Combinatorial Theory, Series A, 120(4):907–940, 2013.  Eprint arXiv:1205.6314.
"""
function compute_halfspaces(M::Matrix{T}, n::Integer) where {T<:SemiRing{S}} where {S}
    # the function returns an error if M is empty
    isempty(M) && error("The family of generators is empty.")

    # TPLib restricts the computation of halfspaces to points with finite coefficients
    any(m->isinfinite(m), M) && error("The points must have finite coefficients.")

    input, io_input = mktemp()
    output, io_output = mktemp()
    write_matrix(io_input, M)
    seekstart(io_input)

    flags = Vector{String}()
    # flag for MaxPlus/MinPlus
    if T <: MinPlus
        push!(flags, "-min-plus")
    end

    # flag for numerical type
    if S==Float64
        push!(flags, "-numerical-data", "ocaml_float")
    elseif S != Int64
        error("Only supports types Int64 and Float64")
    end

    run(pipeline(`compute_halfspaces $flags $n`, stdin=input, stdout=output))

    I = read_halfspaces(io_output, T, n)
    close(io_input)
    close(io_output)
    return I
end

@tropicalize compute_halfspaces

"""
    compute_tropical_complex(M::Matrix{<:SemiRing}, n::Integer)
    compute_tropical_complex(M::Matrix{<:Number}, n::Integer, semiring=:max)
Computes the tropical complex associated with a tropical cone given by a generating set. Only generating sets which are minimal and in which every vector has finite entries are supported.
Returns a pair `(V::Matrix{<:Number},A::Vector{Vector{Int64}})` where the rows of `V` are the vertices of the complex, and `A` is the collection of maximal cells (defined by adjacency with vertices). Note that each cell is seen as a usual polytope, and not as a tropical one, and that the vertices are not necessarily unique.

# References

[1] M. Develin and B. Sturmfels. Tropical convexity. Doc. Math., 9:1–27 (electronic), 2004. E-print arXiv:math.MG/0308254.
"""
function compute_tropical_complex(M::Matrix{T}, n::Integer) where {T<:SemiRing{S}} where {S}
    # TPLib restricts the computation of halfspaces to points with finite coefficients
    any(m->isinfinite(m), M) && error("The points must have finite coefficients.")

    input, io_input = mktemp()
    output, io_output = mktemp()
    write_matrix(io_input, M)
    seekstart(io_input)

    flags = Vector{String}()
    # flag for MaxPlus/MinPlus
    if T <: MinPlus
        push!(flags, "-min-plus")
    end

    # flag for numerical type
    if S<:Float64
        push!(flags, "-numerical-data", "ocaml_float")
    elseif S<:BigInt
        push!(flags, "-numerical-data", "zarith_int")
    elseif S<:Rational
        push!(flags, "-numerical-data", "zarith_rat")
    elseif S != Int64
        error("Only supports types Int64 and Float64")
    end

    run(pipeline(`compute_tropical_complex $flags $n`, stdin=input, stdout=output))

    I = read_complex(io_output, T, n)
    close(io_input)
    close(io_output)
    return I
end

@tropicalize compute_tropical_complex

"""
    compute_tangent_hypergraph(M::Matrix{<:SemiRing}, P::Vector{<:SemiRing}, n::Integer)
    compute_tangent_hypergraph(H::Matrix{<:SemiRing}, A::Vector{Vector{Int64}},  P::Vector{<:SemiRing}, n::Integer) 
    compute_tangent_hypergraph(M::Matrix{<:Number}, P::Vector{<:Number}, n::Integer, semiring=:max)
    compute_tangent_hypergraph(H::Matrix{<:Number}, A::Vector{Vector{Int64}},  P::Vector{<:Number}, n::Integer, semiring=:max) 
Computes the tangent directed hypergraph of a point `P` in a tropical cone specified by its generators or inequalities `M` (the difference is made by wether `M` is composed of `n` or `2n` columns) or by halfspaces H with their sectors `A`.
Returns the number of vertices of the hypergraph, its hyperedges, and in the case where `M` was specified by inequalities or halfspaces, then the active inequalities or halfspaces associated to the hyperedges. Note that the vertices of the hypergraph are labeled starting at 1, whereas they are labeled starting at 0 in TPLib.

# References
[1] X. Allamigeon. Static analysis of memory manipulations by abstract interpretation -- Algorithmics of tropical polyhedra, and application to abstract interpretation. PhD thesis. 

[2] X. Allamigeon, S. Gaubert, E. Goubault. Computing the vertices of tropical polyhedra using directed hypergraphs. Discrete & Computational Geometry, 49(2):247–279, 2013. E-print arXiv:0904.3436v4.
"""
function compute_tangent_hypergraph(M::Matrix{S}, P::Vector{T}, n::Integer) where {S<:SemiRing{U},T<:SemiRing{V}} where {U,V}
    input, io_input = mktemp()
    output, io_output = mktemp()
    write_matrix(io_input, M)
    seekstart(io_input)
    point = vector_to_string(P)

    flags = Vector{String}()
    # flag for MaxPlus/MinPlus
    if promote_type(S,T) <: MinPlus
        push!(flags, "-min-plus")
    end

    # flag for numerical type
    if promote_type(U,V)<:Float64
        push!(flags, "-numerical-data", "ocaml_float")
    elseif promote_type(U,V)<:BigInt
        push!(flags, "-numerical-data", "zarith_int")
    elseif promote_type(U,V)<:Rational
        push!(flags, "-numerical-data", "zarith_rat")
    elseif promote_type(U,V)!=Int64
        error("Only supports types Int64 and Float64")
    end
    
    if size(M)[2] == n
        run(pipeline(`compute_tangent_hypergraph $flags -generators $n $point`, stdin=input, stdout=output))
        F = read_hypergraph_gen(io_output, promote_type(S,T), n)
    elseif size(M)[2] == 2n
        run(pipeline(`compute_tangent_hypergraph $flags $n $point`, stdin=input, stdout=output))
        F = read_hypergraph_ineq(io_output, promote_type(S,T), n)
    else
        error("M must have n (if the cone is defined by generators) or 2n (if it is defined by inequalities) columns.")
    end

    close(io_input)
    close(io_output)
    return F
end

function compute_tangent_hypergraph(H::Matrix{S}, A::Vector{Vector{Int64}}, P::Vector{T}, n::Integer) where {S<:SemiRing{U},T<:SemiRing{V}} where {U,V}
    input, io_input = mktemp()
    output, io_output = mktemp()
    write_halfspaces(io_input, H, A)
    seekstart(io_input)
    point = vector_to_string(P)

    flags = Vector{String}()
    # flag for MaxPlus/MinPlus
    if promote_type(S,T) <: MinPlus
        push!(flags, "-min-plus")
    end

    # flag for numerical type
    if promote_type(U,V)<:Float64
        push!(flags, "-numerical-data", "ocaml_float")
    elseif promote_type(U,V)<:BigInt
        push!(flags, "-numerical-data", "zarith_int")
    elseif promote_type(U,V)<:Rational
        push!(flags, "-numerical-data", "zarith_rat")
    elseif promote_type(U,V)!=Int64
        error("Only supports types Int64 and Float64")
    end
    
    run(pipeline(`compute_tangent_hypergraph $flags -half-spaces $n $point`, stdin=input, stdout=output))
    F = read_hypergraph_halfspace(io_output, promote_type(S,T), n)

    close(io_input)
    close(io_output)
    return F
end

function compute_tangent_hypergraph(M::Matrix{<:Number}, P::Vector{<:Number}, n::Integer, semiring=:max) 
    if semiring === :max
        if !(eltype(M) <: MaxPlus)
            T = eltype([m for m in M if !isminusinfinite(m)])
            M = convert(Matrix{MaxPlus{T}},M)
        end
        if !(eltype(P) <: MaxPlus)
            S = eltype([p for p in P if !isminusinfinite(p)])
            P = convert(Vector{MaxPlus{S}},P)
        end
    elseif semiring === :min
        if !(eltype(M) <: MinPlus)
            T = eltype([m for m in M if !isplusinfinite(m)])
            M = convert(Matrix{MinPlus{T}},M)
        end
        if !(eltype(P) <: MinPlus)
            S = eltype([p for p in P if !isplusinfinite(p)])
            P = convert(Vector{MinPlus{S}},P)
        end
    else
        error("Invalid value for semiring. Must be :max or :min.")
    end
    return compute_tangent_hypergraph(M,P,n)
end

function compute_tangent_hypergraph(H::Matrix{<:Number}, A::Vector{Vector{Int64}}, P::Vector{<:Number}, n::Integer, semiring=:max) 
    if semiring === :max
        if !(eltype(H) <: MaxPlus)
            T = eltype([h for h in H if !isminusinfinite(h)])
            H = convert(Matrix{MaxPlus{T}},H)
        end
        if !(eltype(P) <: MaxPlus)
            S = eltype([p for p in P if !isminusinfinite(p)])
            P = convert(Vector{MaxPlus{S}},P)
        end
    elseif semiring === :min
        if !(eltype(H) <: MinPlus)
            T = eltype([h for h in H if !isplusinfinite(h)])
            H = convert(Matrix{MinPlus{T}},H)
        end
        if !(eltype(P) <: MinPlus)
            S = eltype([p for p in P if !isplusinfinite(p)])
            P = convert(Vector{MinPlus{S}},P)
        end
    else
        error("Invalid value for semiring. Must be :max or :min.")
    end
    println(typeof(H), "\t", typeof(P))
    return compute_tangent_hypergraph(H,A,P,n)
end

end
