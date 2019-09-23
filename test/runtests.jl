using GibbsSampling
using Test
using CSV


## Tips de RAz para ser bien pro
# Em;iezas un Pkg ocn PkgTemplates.jl
# generas un paquete y le haces `dev MyPkg` para que genere el path correcto
# haces todo tu codigo en runtests.jl e iteras con tests

@testset "GibbsSampling.jl" begin
    # Write your own tests here.
    @test 1 == 1
end

abstract type AbstractVector end
struct PVec{T,S} <: AbstractVector

    x :: T
    y :: S
end

pvecx(v :: PVec{T,S}) where {T,S} = v.x
pvecy(v :: PVec{T,S}) where {T,S} = v.y
# TODO Setters en vez de getters

p1 = PVec(1,1)
p2 = PVec(1,0)
p3 = PVec(true,false)


@test p1.x==1
@test p1.y==1
@test p2.x==1
@test p2.y==0
@test p3.x==true
@test p3.y==false

@test typeof(p3.x) <: Bool && typeof(p3.y) <: Bool


path = "/home/pepefv97/.julia/dev/GibbsSampling/src/"
cd(path)
xa = CSV.read("centroids.csv",header=false)
xa = convert(Matrix, xa)
# Si entradas estan al reves correr metodo distinto para escribir matriz desde numpy
# probalemente cambiar el a.flatten()
mat = [PVec(xa[i],xa[i+1]) for i in 1:2:length(xa)-1]
mat = reshape(mat,192,202)

kronecker(i,j, x) = i == j ? one(x) : zero(x)
Base.one(v :: PVec{T,S}) where {T,S} = one(pvecx(v))
Base.zero(v :: PVec{T,S}) where {T,S} = zero(pvecx(v))

@test kronecker(1,1,1) == 1
@test kronecker(1,1,1.0) == 1.0
@test kronecker(1,2,1) == 0
@test kronecker(1,2,1.0) == 0.0
@test kronecker(1,1,true) == true
@test kronecker(1,2,true) == false

p4 = PVec(4.0,Bool)
@test one(1) == 1.0
@test zero(p4) == 0.0

choice() = rand(1:10)
@test 1 <= choice() <= 10


function energia(i,j,mat, r)
    @assert 2 == length(size(mat))
    res1 :: = interactuar_vecinos(i, j, mat, r)
    res2 = res1 - delta_vecinos(i, j, mat, r)
end

delta_vecinos(i,j,mat, r) = kronecker(pvecx(mat[i,j]), r, mat[i,j])
@test delta_vecinos(6,4,mat, choice()) == 0.0

#@inline function interactuar_vecinos(i,j,x, r) :: typeof(one(x[1,1]))
#    sum(skipmissing((vecino(findvecinoup(i,j,x),r) ,vecino(findvecinodown(i,j,x),r),vecino(findvecinoleft(i,j,x),r),vecino(findvecinoright(i,j,x),r))))
#end
@inline function interactuar_vecinos(i,j,x,r)
    m=vecino(findvecinoup(i,j,x),r)+vecino(findvecinodown(i,j,x),r)+vecino(findvecinoleft(i,j,x),r)+vecino(findvecinoright(i,j,x),r)
end
# Para cambiar a forntera toroidal usar mod1
@inbounds findvecinoup(i,j,x) =  i > 1 ? x[i-1,j] : missing
@inbounds findvecinodown(i,j,x) = i < size(x)[1] ? x[i+1,j] : missing
@inbounds findvecinoleft(i,j,x) = j > 1 ? x[i,j-1] : missing
@inbounds findvecinoright(i,j,x) = j < size(x)[2] ? x[i,j+1] : missing



vecino(v :: PVec{T,S}, r) where {T,S} = (pvecx(v) - r)^2
#vecino(x :: Missing, r) = missing
vecino(x :: Missing, r) = 0

# TODO


@test findvecinoup(1,1,mat) === missing
@test findvecinodown(192,1,mat) === missing
@test findvecinoright(1,202,mat) === missing
@test findvecinoleft(1,1,mat) === missing

@test findvecinoup(1,10,mat) === missing

abstract type AbstractSimulation end
struct CircularSimulation <: AbstractSimulation
    N :: Int
    mat :: Matrix
end
struct MCSimulation <: AbstractSimulation
    N :: Int
    mat :: Matrix
end

# energias(s :: CircularSimulation, mat) = println("tada")
# energias(s :: MCSimulation, mat) = println("asa")

"""
exp(- energia(i,j,mat,r) * log(i + j + k))/ exp( - energia(i,j,mat,pvecx(mat[i,j])) * log(i + j + k))

counter incrementa con cada recorrido entero de la matriz
"""
@inbounds function cociente(i,j,mat,counter,r)
    set_zero_subnormals(true)
    exp(- energia(i,j,mat,r) * log(i + j + counter) + energia(i,j,mat,pvecx(mat[i,j])) * log(i + j + counter))
end

@inbounds function update!(i,j,mat, α, r)
    if α >= 1.0
        mat[i,j] = PVec(convert(typeof(one(mat[i,j])),r),pvecy(mat[i,j])) # TODO Setters
    else
        ρ = rand() # esto vive en [0,1] uniforme
        if ρ <= α
            return
        else
            mat[i,j] = PVec(convert(typeof(one(mat[i,j])),r), pvecy(mat[i,j]))
        end
    end
end

using UnicodePlots
miniplot(mat) = spy(sim.mat .|> pvecx)
@inbounds function simulate(s :: AbstractSimulation)
    # TODO usar Random.seed!(7)
    I,J = [size(s.mat)...]
    N = s.N

    @inbounds for counter in 1:N
        for idxj in 1:J
            for idxi in 1:I

                # Desechar algunos casos
                if pvecy(mat[idxi,idxj]) == one(mat[1,1])
                    continue
                end
                r = choice()
                # En casos de interes, checar propension energetica
                α = cociente(idxi,idxj,mat, counter, r)

                # Reemplazar valores
                update!(idxi,idxj, mat, α, r)
            end
        end
    end

    return mat
    # TODO guardar las ultimas mil
end

### Benchmarking

### Threads
mattest = deepcopy(mat2)
Random.seed!(7)

@btime PVec(1,2)
@btime cociente(1,1,$mattest,1,4)

@test update!(1,1,mattest, )
