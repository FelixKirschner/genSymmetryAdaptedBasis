# Symmetry Adapted Basis

This code produces a symmetry adapted basis of polynomials invariant under the action of the symmetric group S_n. The software is used in this paper https://arxiv.org/abs/2203.05892.

To run the programm set n to be the number of variables, r to the degree. Run initializeBasis(n,r) to obtain a symmetry adapted basis in n variables up to degree r. The program returns a list of subbases. Each subbasis contains lists of arrays. Each of these list correspond to one element of the symmetry adapted basis. The elements of the lists are vectors, whose first n elements indicate the exponent of the monomial they correspond to, while the last index of the array indicates the coefficient. For example, the monomial m(x^a)(y^b)(z^c) appears as [a, b, c, m].
