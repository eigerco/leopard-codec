# Leopard codec

This is a Rust rewrite of a Go `O(n*log n)` implementation of Reed-Solomon codes,
ported from the [klauspost/reedsolomon](https://github.com/klauspost/reedsolomon),
which is a port of the C++ library [catid/leopard](https://github.com/catid/leopard).


The original implementation is based on the paper:

S.-J. Lin, T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,

"Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"

IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.

## Features support

The leopard algorithm uses either 8-bit or 16-bit Galois fields with Cantor basis.
The 8-bit implementation should be used with up to 256 total shards and 16-bit when more
shards are needed.

- [x] Encoding parity shards using 8-bit leopard algorithm
- [x] Reconstructing shards using 8-bit leopard algorithm
- [ ] Encoding parity shards using 16-bit leopard algorithm
- [ ] Reconstructing shards using 16-bit leopard algorithm

## Contributing

We welcome contributions! Please fork the repository and submit a pull request.

## License

Leopard codec is licensed under the Apache 2.0 License. See the [LICENSE](./LICENSE) file for more details.

## About [Eiger](https://www.eiger.co)

We are engineers. We contribute to various ecosystems by building low level implementations and core components.

Contact us at hello@eiger.co
Follow us on [X/Twitter](https://x.com/eiger_co)
