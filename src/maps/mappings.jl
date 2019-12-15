import Base.*
import Base.inv

using OffsetArrays

struct Map{T}
	img :: OffsetVector{T}
end


function Map(x) where T
	v = OffsetVector(x, 0 : length(x) - 1)
	return Map(v)
end

#=
function Map(f :: Family{T}) where T
	n = length(f)
	x = zeros(Int, 2^n)
	for i in 1 : 2^n
		x[i] = f(i - 1)
	end
	return Map(x)
end
=#

function (m :: Map)(x)
	return m.img[x]
end


function len(m :: Map{T}) where T
	return length(m.img)
end


function *(m1 :: Map{T}, m2 :: Map{T}) where T
	n = len(m1)
	x = zeros(T, n)
	for i in 0 : n-1
		# i -> m2(i) -> m1(m2(i))
		x[i+1] = m1.img[m2.img[i]]
	end
	return Map(x)
end


function inv(m :: Map{T}) where T
	n = length(m.img)
	x = zeros(Int, n)
	for i in 0 : n - 1
		x[m.img[i] + 1] = i
	end
	return Map(x)
end


function fixpt(m :: Map{T}) where T
	n = len(m)
	fix = T[]
	for i in 0 : n - 1
		if m(i) == i
			push!(fix, i)
		end
	end
	return fix
end


function comm(m1 :: Map, m2 :: Map)
	return m1 * m2 * inv(m1) * inv(m2)
end
	
