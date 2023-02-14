# ToDo: Deal with NaNs

module Tropical

export Infinite, SemiRing, MaxPlus, MinPlus, isinfinite, isminusinfinite, isplusinfinite, elt, ∞

abstract type SemiRing{T}<:Number end

struct Infinite<:Number end
const ∞ = Infinite()

struct MaxPlus{T}<:SemiRing{T}
    elt::Union{T, Infinite}
    MaxPlus{T}(x::Number) where {T<:Number} = isminusinfinite(x) ? new{T}(Infinite()) : new{T}(T(x))
end

struct MinPlus{T}<:SemiRing{T}
    elt::Union{T, Infinite}
    MinPlus{T}(x::Number) where {T<:Number} = isplusinfinite(x) ? new{T}(Infinite()) : new{T}(T(x))
end

# Pretty printing
Base.show(io::IO, ::Infinite) = print(io,"⋅")
Base.show(io::IO, x::SemiRing) = print(io,elt(x))

# Fundamental functions
"""
    elt(x::SemiRing)
Returns the tropical element.
"""
elt(x::SemiRing) = x.elt

isplusinfinite(x::Number) = x isa Infinite || x == Inf || x == Inf32 || x == Inf16
isminusinfinite(x::Number) = x isa Infinite || x == -Inf || x == -Inf32 || x == -Inf16

"""
    isinfinite(x::SemiRing)
Returns `true` if `x` is the infinite element of the tropical field (ie the tropical 0), `false` otherwise.
"""
isinfinite(x::MaxPlus) = isminusinfinite(elt(x))
isinfinite(x::MinPlus) = isplusinfinite(elt(x))

# Constructors, conversions, promotions
"""
    MaxPlus{T}(x::Number) where {T<:Number}
    MaxPlus(x::Number)
Constructor for the type `MaxPlus`. If a type `T` is given, `x` is converted to that type, otherwise the type is inferred from `x`.
"""
MaxPlus(x::T) where {T<:Number} = T <: MaxPlus ? x : MaxPlus{T}(x)
MaxPlus{T}(x::MaxPlus{S}) where {T<:Number, S<:Number} = MaxPlus{T}(elt(x))

"""
    MinPlus{T}(x::Number) where {T<:Number}
    MinPlus(x::Number)
Constructor for the type `MinPlus`. If a type `T` is given, `x` is converted to that type, otherwise the type is inferred from `x`.
"""
MinPlus(x::T) where {T<:Number} = T <: MinPlus ? x : MinPlus{T}(x)
MinPlus{T}(x::MinPlus{S}) where {T<:Number, S<:Number} = MinPlus{T}(elt(x))

Base.convert(::Type{MaxPlus{T}}, x::MaxPlus) where {T<:Number} = MaxPlus{T}(elt(x))
Base.convert(::Type{MinPlus{T}}, x::MinPlus) where {T<:Number} = MinPlus{T}(elt(x))
Base.convert(::Type{MaxPlus{T}}, x::Number) where {T<:Number} = MaxPlus{T}(x)
Base.convert(::Type{MinPlus{T}}, x::Number) where {T<:Number} = MinPlus{T}(x)
Base.convert(::Type{T}, x::MaxPlus) where {T<:Number} = isinfinite(x) ? T(-Inf) : T(elt(x))
Base.convert(::Type{T}, x::MinPlus) where {T<:Number} = isinfinite(x) ? T(Inf) : T(elt(x))
Base.promote_rule(::Type{MaxPlus{T}}, ::Type{S}) where {T<:Number,S<:Number} = MaxPlus{promote_type(T,S)}
Base.promote_rule(::Type{MaxPlus{T}}, ::Type{MaxPlus{S}}) where {T<:Number,S<:Number} = MaxPlus{promote_type(T,S)}
Base.promote_rule(::Type{MinPlus{T}}, ::Type{S}) where {T<:Number,S<:Number} = MinPlus{promote_type(T,S)}
Base.promote_rule(::Type{MinPlus{T}}, ::Type{MinPlus{S}}) where {T<:Number,S<:Number} = MinPlus{promote_type(T,S)}

# Arithmetic operations
"""
    Base.:(==)(x::S,y::S) where {S<:SemiRing}
Equality of tropical numbers.
"""
Base.:(==)(x::SemiRing, y::SemiRing) = (isinfinite(x) || isinfinite(y)) ? (isinfinite(x) && isinfinite(y)) : elt(x) == elt(y)

"""
    Base.:(*)(x::S,y::S) where {S<:SemiRing}
Multiplication of tropical numbers.

### Example
```jldoctest
julia> using TPLib

julia> MaxPlus(4) * MaxPlus(8)
12

julia> MinPlus(2) * MinPlus(-3.5)
-1.5

julia> MaxPlus(10.0) * -Inf
⋅
```
"""
function Base.:(*)(x::S,y::S) where {S<:SemiRing}
    isinfinite(x) && return x
    isinfinite(y) && return y
    return S(elt(x) + elt(y))
end

"""
    Base.isless(x::S,y::S) where {S<:SemiRing}
Comparison of tropical numbers.
"""
function Base.isless(x::MaxPlus, y::MaxPlus)
    isinfinite(x) && !isinfinite(y) && return true
    !isinfinite(x) && isinfinite(y) && return false
    isinfinite(x) && isinfinite(y) && return false
    !isinfinite(x) && !isinfinite(y) && return elt(x) < elt(y)
end

function Base.isless(x::MinPlus, y::MinPlus)
    isinfinite(x) && !isinfinite(y) && return false
    !isinfinite(x) && isinfinite(y) && return true
    isinfinite(x) && isinfinite(y) && return false
    !isinfinite(x) && !isinfinite(y) && return elt(x) < elt(y)
end

"""
    Base.:(+)(x::S,y::S) where {S<:SemiRing}
Addition of tropical numbers.

### Example
```jldoctest
julia> using TPLib

julia> MaxPlus(3) + MaxPlus(8)
8

julia> MinPlus(4.) + MinPlus{Int64}(∞)
4.0

julia> MaxPlus(6.2) + 4
6.2
```
"""
Base.:(+)(x::MaxPlus,y::MaxPlus) = max(x,y)
Base.:(+)(x::MinPlus,y::MinPlus) = min(x,y)

end
