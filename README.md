# Leopard codec

This is a Rust rewrite of a Go `O(n*log n)` implementation of Reed-Solomon codes,
ported from the https://github.com/klauspost/reedsolomon,
which is a port from the C++ library https://github.com/catid/leopard.


The original implementation is based on the paper:

S.-J. Lin, T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,

"Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"

IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.

## Features support

The leopard algorithm uses either 8-bit or 16-bit Galois fields with Cantor basis.
The 8-bit implementation should be used with up to 256 total shards and 16-bit when more
shards are needed.

- [x] Encoding parity shards using 8-bit leopard algorithm
- [ ] Reconstructing shards using 8-bit leopard algorithm
- [ ] Encoding parity shards using 16-bit leopard algorithm
- [ ] Reconstructing shards using 16-bit leopard algorithm
