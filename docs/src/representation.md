# Representation of tropical numbers

```@setup example
using TPLib
```

A tropical number is an element of `MaxPlus{T}` or `MinPlus{T}` where `T` is a numeric type. These types are a subtype of the abstract type `SemiRing{T}`, which is itself a subtype of `Number`. They have a single field defined as the union `Union{T,Infinite}` where `Infinite` is a constant type without field representing the infinite element of the tropical ring,

```julia
 struct Infinite<:Number end
```

The unique element of this type is denoted `∞`, and is displayed as ``\cdot`` in the REPL.

```repl example
∞
typeof(∞)
```

In the case where the numeric type has an infinite element, it will be equal to `∞`.

```repl example
MaxPlus{Float64}(∞) == -Inf
MaxPlus{Float64}(-Inf) == ∞
MinPlus{Float64}(Inf) == ∞
```

In this last example, the comparison of `∞` to an element of type `MaxPlus` forced its conversion to an element of type `MaxPlus`. In general, `∞` can be either of type `MaxPlus` or `MinPlus`, allowing for some flexibility.

```repl example
∞ == -Inf
V = [∞,0,-3,5,∞,-8]
convert(Vector{MaxPlus{Int64}},V)
convert(Vector{MinPlus{Int64}},V)
```

The field of a type `SemiRing{T}` is accessed through the function `elt`.

```repl example
a = MaxPlus(3)
elt(a)
b = MaxPlus{Int64}(∞)
elt(a)
```

To construct an element of `MaxPlus{T}` or `MinPlus{T}`, you can either use the constructor without precising the type, in which case the type will be inferred from the element, or you can precise the type of the semiring, in which case the element will be converted.
