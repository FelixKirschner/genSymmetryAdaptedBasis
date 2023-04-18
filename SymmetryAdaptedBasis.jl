include("libs.jl")

function fillmonomials(n, r)
    #this function returns an array of exponent vectors for all monomials in n variables up to degree r
    monlist = []
    for p in Combinatorics.combinations(1:(n+r), r)
        sort!(p)
        c = zeros(Int64, 1, n + 1)
        pos = 1
        lastPos = 0
        for i in p
            pos = pos + (i - lastPos - 1)
            c[pos] += 1
            lastPos = i
        end
        push!(monlist, c[1:n])
    end
    return monlist
end


function flips(gamma)
    n = length(gamma)
    binarr = [[in(a, c) ? 1 : -1 for a = 1:n] for c in combinations(1:n)]
    push!(binarr, [-1 for i = 1:n])
    binarr2 = []
    for el in binarr
        if !in((-1) * el, binarr2)
            push!(binarr2, el)
        end
    end

    list = []
    for el in binarr2
        push!(list, [gamma[i] * el[i] for i = 1:n])
    end
    return unique(list)
end

function symmetrizeBaseProduct(baseElement1::Vector{Any}, baseElement2::Vector{Any})
    #This function returns the symmetrized product of two base elements in the symmetry adapted basis
    n = length(baseElement1[1]) - 1

    perms = collect(permutations([i for i = 1:n]))

    facn = factorial(n)

    pairs = Vector{Int64}[]

    coeffs = []

    for el1 in baseElement1
        for el2 in baseElement2
            tmp = el1[1:n] - el2[1:n]
            tmp2 = el1[n+1] * el2[n+1]

            for p in perms
                tmp3 = [tmp[p[i]] for i = 1:n]
                push!(coeffs, tmp2)
                push!(pairs, tmp3)
            end
        end
    end

    auxlistPairs = unique(pairs)
    auxlistCoeffs = []
    for el in auxlistPairs
        push!(auxlistCoeffs, sum(coeffs[i] for i in findall(x -> x == el, pairs)))
    end
    retListMon = Vector{Int64}[]
    retListCoeff = []
    for i = 1:length(auxlistCoeffs)
        if auxlistCoeffs[i] != 0.0
            push!(retListMon, auxlistPairs[i])
            push!(retListCoeff, auxlistCoeffs[i] / facn)
        end
    end
    return [retListMon, retListCoeff]
end


function checkOrth(base)
    #NOT USED in final code
    #test function, are the elements of the base pairwise orthogonal ? true : false 
    orth = true
    for i = 1:length(base)
        for j = i+1:length(base)
            for l = 1:length(base[i])
                for m = 1:length(base[j])
                    list = symmetrizeBaseProduct(base[i][l], base[j][m])
                    if !isempty(list[1])
                        orth = false
                        display((i, j, k, l))
                    end
                end
            end
        end
    end
    return orth
end

function initializeBasis(n, r)
    #Construct the symmetry adapted bases for n variables up to degree r
    λ = [AbstractAlgebra.Partition(el) for el in collect(partitions(n))]

    βf, μf = init(n, r)

    fbas = [getPartOfBasis(λ[i], μf, βf, r) for i = 1:length(λ)]

    #just to be sure...
    for base in fbas
        unique!(base)
    end
    retbas = []
    for base in fbas
        if length(base) != 0
            push!(retbas, base)
        end
    end
    return retbas
end


function orbitSize(mon)
    #NOT USED in final code
    #returns the orbit size of a monomial under S_n
    perms = collect(permutations([i for i = 1:length(mon)]))
    list = []
    for p in perms
        push!(list, [mon[p[i]] for i = 1:length(mon)])
    end
    return length(unique!(list))
end



function returnProductsDict(n::Integer, r::Integer)

    dictMons = Dict()
    rMons = []
    baseProducts = []
    dictBaseproducts = Dict()

    @info("Generating symmetry adapted basis")
    fbas = initializeBasis(n, r)

    @info("Generating symmetrized dictBaseproducts")

    for i1 = 1:length(fbas)
        for i2 = 1:length(fbas[i1])
            for i3 = 1:length(fbas[i1])
                symmProd = symmetrizeBaseProduct(fbas[i1][i2], fbas[i1][i3])
                push!(baseProducts, symmProd)
                dictBaseproducts[(i1, i2, i3)] = symmProd
                for el in symmProd[1]
                    tmpSum = sum(abs(el[i]) for i = 1:n)
                    if !haskey(dictMons, el)
                        dictMons[el] = [[i1, i2, i3]]
                    else
                        tmp = dictMons[el]
                        dictMons[el] = push!(tmp, [i1, i2, i3])
                    end
                    if tmpSum > r
                        push!(rMons, el)
                    end
                end
            end
        end
    end
    return baseProducts, dictBaseproducts, fbas, unique(rMons), dictMons
end


function init(n, k)
    β = []
    μ = []

    β = fillmonomials(n, k)

    for el in β
        sort!(el)
    end
    unique!(β)

    for b in β
        btmp = unique(b)
        mu = []
        for el in btmp
            cnt = count(x -> x == el, b)
            push!(mu, cnt)
        end
        sort!(mu, rev = true)
        push!(μ, AbstractAlgebra.Partition([mu[i] for i = 1:length(mu)])) #
    end

    return β, μ
end


# Split permutation into transpositions of form (j,j+1)
function permDecomp(p::Perm)
    res = []
    arr = copy(p.d)
    n = length(arr)
    for i = 1:n
        for j = n:-1:(i+1)
            if arr[j-1] > arr[j]
                push!(res, perm(union(1:(j-2), [j, j - 1], (j+1):n)))
                tmp = arr[j-1]
                arr[j-1] = arr[j]
                arr[j] = tmp
            end
        end
    end
    return res
end

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


#set entry p in tableau Y to value c
function setTabEntry(Y, p, c)
    linPos = 0
    for i = 1:(p[1]-1)
        linPos += Y.part[i]
    end
    linPos += p[2]
    Y.fill[linPos] = c
end

#fill tableau with entries p
function permTab(p, Y)
    return YoungTableau(Y.part, [p[i] for i in Y.fill])
end

#returns true if tableau Y is semi standard
function isSemiStandard(Y)
    s = size(Y)
    #check rows
    for row = 1:s[1]
        last = Y[row, 1]
        for col = 2:(Generic.rowlength(Y, row, 1)+1)
            if last > Y[row, col]
                return false
            end
            last = Y[row, col]
        end
    end
    # Check columns
    for col = 1:s[2]
        last = Y[1, col]
        for row = 2:(Generic.collength(Y, 1, col)+1)
            if last >= Y[row, col]
                return false
            end
            last = Y[row, col]
        end
    end
    return true
end

#generates all semi standard tableaux of shape lam and content mu
function generateGenSemStd(lam, mu)
    if lam.n != mu.n
        return nothing
    end

    Y = YoungTableau(lam.part)
    mu = mu.part
    content = [i for i = 1:length(mu) for j = 1:mu[i]]
    posComb = unique(permutations(content) |> collect)
    cnt = 0
    semiStdTab = []

    for el in posComb
        tmpY = permTab(el, Y)
        if isSemiStandard(tmpY)
            cnt += 1
            push!(semiStdTab, tmpY)
        end
    end

    return [cnt, semiStdTab]
end

#get b, important in order to construct the specht polynomials
function getb(beta, mu , r)
    if length(beta) != sum(mu)
        return nothing
    end

    b = [0 for it = 1:length(mu)]
    vec = [0 for it = 1:r+1]
    for j = 1:r+1
        vec[j] = count(x -> x == j - 1, beta)
    end
    for i = 1:length(mu)
        indx = findfirst(x -> x == mu[i], vec)
        b[i] = indx - 1
        vec[indx] = 0
    end
    return b
end

#get the Column Stabilizer Set CStabt
function getCStabt(t)
    s = size(t)
    colStabt = []
    colStabtfin = []
    for i = 1:s[2]
        col = t[:, i]
        z = findall(x -> x == 0, col)
        for k = length(z):-1:1
            deleteat!(col, z[k])
        end
        colPer = unique(collect(permutations(col)))
        if i == 1 && s[2] != 1
            colStabt = colPer
            continue
        elseif i == 1 && s[2] == 1
            colStabt = colPer
            return colStabt
        else
            for el in colStabt
                for el2 in colPer
                    push!(colStabtfin, [el; el2])
                end
            end
        end
        if i != s[2]
            colStabt = colStabtfin
            colStabtfin = []
        elseif i == s[2]
            return colStabtfin
        end
    end
end

#returns the sign of a permutation where the entries are [1:n]
function getsignDan(p)
    res = []
    arr = copy(p)
    n = length(arr)
    for i = 1:n
        for j = n:-1:(i+1)
            if arr[j-1] > arr[j]
                push!(res, perm(union(1:(j-2), [j, j - 1], (j+1):n)))
                tmp = arr[j-1]
                arr[j-1] = arr[j]
                arr[j] = tmp
            end
        end
    end
    return iseven(length(res)) ? 1 : -1
end

#returns the sign of an element of CStabt
function getsign(p, tab)
    s = size(tab)
    cols = [0 for i = 1:s[2]]
    cols2 = []
    for i = 1:s[2]
        colcont = []
        cnt = 1
        for j = 2:s[1]
            if tab[j,i] != 0 && j != s[1]
                cnt += 1
                push!(colcont, tab[j,i])
            elseif tab[j,i] != 0 && j == s[1]
                cols[i] = cnt + 1
                push!(colcont, tab[j,i])
            elseif tab[j,i] == 0
                cols[i] = cnt
                push!(cols2, colcont)
                break
            end
        end
    end
    sign = 1
    ColStab = []
    for i = 1:s[2]
        if i == 1
            ColStab = p[1:cols[1]]
        else
            summe = sum(cols[k] for k=1:i-1)
            ColStab = p[summe+1:summe + cols[i]]
        end
        vect = [0 for k = 1:length(ColStab)]
        for j = 1:length(ColStab)
            vect[j] = findfirst( x -> x == tab[j,i], ColStab)
        end
        sign = sign*getsignDan(vect)
    end
    return sign
end

#returns {T} i.e. the row equivalence class of T
function getEquClaT(Tab)
    equClass = []
    list = unique(getCStabt(conj(Tab)))
    for el in list
        tmpTab = YoungTableau(Tab.part, el)
        push!(equClass, tmpTab)
    end
    return equClass
end

#returns sign(sigma) * sigma (X^(t,T)) i think
function getsubBaseEl(tab, STab, perm, bExp)
    n = length(perm)
    tmp = YoungTableau(conj(tab).part, perm)
    tmp = conj(tmp)
    s = size(tmp)
    x = [0 for i = 1:n+1]
    x[n+1] = getsign(perm, tab)
    for i = 1:s[1]
        for j = 1:s[2]
            if tab[i,j] != 0
                ind = tmp[i, j]
                expo = bExp[STab[i, j]]
                x[ind] = expo
            end
        end
    end
    return x
end

#for given lambda, mu, beta this return the subbasis belonging to those
function getPartOfBasis(lambda, mu, beta,r)
    t = YoungTableau(lambda)
    retBase = []
    for i = 1:length(mu)
        if generateGenSemStd(lambda, mu[i])[1] == 0
            continue
        end
        b = getb(beta[i], mu[i].part,r)
        listssTabs = generateGenSemStd(lambda, mu[i])[2]
        for T in listssTabs
            Slist = getEquClaT(T)
            sigmaCSt = getCStabt(t)
            bas = []
            for sigma in sigmaCSt
                for S in Slist
                    x = getsubBaseEl(t, S, sigma, getb(beta[i], mu[i].part,r))
                    push!(bas, x)
                end
            end
            push!(retBase, bas)
        end
    end
    return retBase
end