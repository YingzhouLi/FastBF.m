FastBF
======

Fast Butterfly Factorization for Fourier Integral Operators

This is also known as Interpolative Butterfly Factorization [1].

The current implementation supports any dimension
Fourier integral operators with/without singularity at the origin.

Fourier integral operatiors without singularity could be solved via `fastBF`
in the `src` folder, whereas the Fourier integral operators with singularity at the origin
could be solved via either `fastBF` with polar transform or `fastMBF`.
The former adopts fast butterfly factorization with the idea given in [2]
and the later adopts the idea of multiscale domain decomposition [3] together
with the interpolative butterfly factorization.

1. The example for `fastBF` can be found in `test/test_fastbf_1D` and `test/test_fastbf_2D`,

2. The example for `fastBF` with polar transform can be found in `test/test_fastpbf_2D`,

3. The example for `fastMBF` can be found in `test/test_fastmbf_2D`.

More examples of special function transforms can be found in `test` folder as well.

Reference:

[1] Y. Li and H. Yang. Interpolative Butterfly Factorization. Submitted. [PDF][ibf]

[2] E. Candes, L. Demanet and L. Ying. A fast butterfly algorithm for the computation of Fourier integral operators. SIAM Multiscale Modeling and Simulation 7 (2009). [PDF][pbf]

[3] Y. Li, H. Yang, and L. Ying. Multidimensional butterfly factorization. Submitted. [PDF][mbf]

[ibf]: http://arxiv.org/abs/1605.03616
[pbf]: http://epubs.siam.org/doi/abs/10.1137/080734339
[mbf]: http://arxiv.org/abs/1509.07925
