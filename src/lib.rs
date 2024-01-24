use std::{mem, ops::Range};

use bytes::{Buf, BufMut};
use thiserror::Error;
use utils::alloc_aligned;

// Number of bits per element
pub const BITS: usize = 8;
// Finite field order: Number of elements in the field
pub const ORDER: usize = u8::MAX as usize + 1;
// Modulus for field operations
pub const MODULUS: u8 = u8::MAX;
// LFSR Polynomial that generates the field elements
pub const POLYNOMIAL: usize = 0x11D;
// Basis used for generating logarithm tables
pub const CANTOR_BASIS: [u8; BITS] = [1, 214, 152, 146, 86, 200, 88, 230];

// Size of blocks during encoding
const BLOCK_SIZE: usize = 32 << 10;
// const INVERSION8_BYTES: usize = 256 / 8;

/// lookup tables
pub mod lut;

mod utils;

/// Possible errors that can happen when interacting with Leopard.
#[derive(Debug, Error)]
pub enum LeopardError {
    #[error("Maximum shard number ({}) exceeded: {0}", u16::MAX)]
    MaxShardNumberExceeded(usize),

    #[error("Incorrect amount of shards ({0}), expected ({1})")]
    IncorrectAmountOfShards(usize, usize),

    #[error("Shards contain no data")]
    EmptyShards,

    #[error("Shards of different lengths found")]
    UnequalShardsLengths,

    #[error("Shard size ({0}) should be a multiple of 64")]
    InvalidShardSize(usize),
}

/// A result type with [`LeopardError`].
pub type Result<T, E = LeopardError> = std::result::Result<T, E>;

trait SliceShrink {
    fn shrink(&mut self, range: Range<usize>);
}

impl<T> SliceShrink for &mut [T] {
    fn shrink(&mut self, range: Range<usize>) {
        let (_, mut last) = mem::take(self).split_at_mut(range.start);
        let (mid, _) = mem::take(&mut last).split_at_mut(range.end);
        *self = mid;
    }
}

pub struct LeopardFF8 {
    data_shards: usize,
    parity_shards: usize,
    total_shards: usize,
    // workPool    sync.Pool
    // inversion   map[[inversion8Bytes]byte]leopardGF8cache
    // inversionMu sync.Mutex

    // o options
}

impl LeopardFF8 {
    pub fn new(data_shards: usize, parity_shards: usize) -> Result<Self> {
        let total_shards = data_shards + parity_shards;

        // if opt.inversionCache && (r.totalShards <= 64 || opt.forcedInversionCache) {
        // 	   // Inversion cache is relatively ineffective for big shard counts and takes up potentially lots of memory
        // 	   // r.totalShards is not covering the space, but an estimate.
        // 	   r.inversion = make(map[[inversion8Bytes]byte]leopardGF8cache, r.totalShards)
        // }

        if total_shards > u16::MAX as usize {
            Err(LeopardError::MaxShardNumberExceeded(total_shards))
        } else {
            Ok(Self {
                data_shards,
                parity_shards,
                total_shards,
            })
        }
    }

    pub fn encode(&self, shards: &mut [&mut [u8]]) -> Result<()> {
        if shards.len() != self.total_shards {
            return Err(LeopardError::IncorrectAmountOfShards(
                shards.len(),
                self.total_shards,
            ));
        }

        check_shards(shards, false)?;

        self.encode_inner(shards)
    }

    fn encode_inner(&self, shards: &mut [&mut [u8]]) -> Result<()> {
        let shard_size = shard_size(shards);
        if shard_size % 64 != 0 {
            return Err(LeopardError::InvalidShardSize(shard_size));
        }

        let m = ceil_pow2(self.parity_shards);
        // TODO:
        // var work [][]byte
        // if w, ok := r.workPool.Get().([][]byte); ok {
        // 	work = w
        // }
        // TODO: I don't understand why do we create a work mod after so let's just go like so
        let mut work_mod = alloc_aligned::<BLOCK_SIZE>(m * 2)?;

        // defer r.workPool.Put(work)

        let mtrunc = m.min(self.data_shards);

        // TODO: what? it has only 256 values
        let skew_lut = &lut::FFT_SKEW[m - 1..];

        // TODO: what? we've checked all have the same size
        // Split large shards.
        // More likely on lower shard count.
        let mut off = 0;
        let mut sh = vec![vec![0; shard_size]; shards.len()];

        // let mut work_mod = alloc_aligned::<BLOCK_SIZE>(m * 2)?;
        // work_mod
        //     .iter_mut()
        //     .zip(work.iter())
        //     .for_each(|(wm, w)| wm.copy_from_slice(w));

        let (data_shards, parity_shards) = shards.split_at_mut(self.data_shards);
        let mut data_shards_chunk: Vec<_> = data_shards.iter_mut().collect();
        let mut parity_shards_chunk: Vec<_> = parity_shards.iter_mut().collect();

        while off < shard_size {
            let mut work: Vec<_> = work_mod.iter_mut().map(|shard| &mut shard[..]).collect();
            let mut sh: Vec<_> = sh.iter_mut().map(|shard| &mut shard[..]).collect();
            let mut end = off + BLOCK_SIZE;
            if end > shard_size {
                end = shard_size;
                let sz = shard_size - off;
                work.iter_mut().for_each(|shard| {
                    *shard = &mut mem::take(shard)[..sz];
                })
            }

            sh.iter_mut()
                .zip(data_shards.iter_mut())
                .for_each(|(sh, input)| {
                    *sh = &mut input[off..end];
                });

            work.iter_mut()
                .zip(parity_shards.iter_mut())
                .for_each(|(work, parity)| {
                    *work = &mut parity[off..end];
                });

            ifft_dit_encoder(
                &mut sh[..self.data_shards],
                mtrunc,
                &mut work,
                None, // No xor output
                m,
                skew_lut,
            );

            let last_count = self.data_shards % m;
            let mut skew_lut2 = skew_lut;

            // goto skip_body
            if m < self.data_shards {
                // for sets of m data pieces
                for _ in (m..self.data_shards - m).step_by(m) {
                    // work <- work xor IFFT(data + i, m, m + i)
                    skew_lut2 = &skew_lut2[m..];
                    let (xor_output, temp_work) = work.split_at_mut(m);
                    ifft_dit_encoder(
                        &mut sh[m..],
                        m,
                        temp_work,
                        Some(xor_output), // No xor output
                        m,
                        &skew_lut2[m..],
                    );
                }

                // Handle final partial set of m pieces:
                if last_count != 0 {
                    let (xor_output, temp_work) = work.split_at_mut(m);
                    ifft_dit_encoder(
                        &mut sh[m..],
                        m,
                        temp_work,
                        Some(xor_output), // No xor output
                        m,
                        &skew_lut2[m..],
                    );
                }
            }

            // work <- FFT(work, m, 0)
            fft_dit(work.as_mut_slice(), self.parity_shards, m, &*lut::FFT_SKEW);
            off += BLOCK_SIZE;
        }

        Ok(())
    }
}

fn shard_size(shards: &[impl AsRef<[u8]>]) -> usize {
    shards
        .iter()
        .map(|shard| shard.as_ref().len())
        .find(|&len| len != 0)
        .unwrap_or(0)
}

// checkShards will check if shards are the same size
// or 0, if allowed. An error is returned if this fails.
// An error is also returned if all shards are size 0.
fn check_shards(shards: &[impl AsRef<[u8]>], allow_zero: bool) -> Result<()> {
    let size = shard_size(shards);

    if size == 0 {
        if allow_zero {
            return Ok(());
        } else {
            return Err(LeopardError::EmptyShards);
        }
    }

    // NOTE: for happy case fold would be faster
    let are_all_same_size = shards.iter().all(|shard| shard.as_ref().len() == size);

    if !are_all_same_size {
        return Err(LeopardError::UnequalShardsLengths);
    }

    Ok(())
}

// pub fn ref_mul_add(x: &mut [u8], y: &[u8], log_m: u8) {
//     #[cfg(target_arch = "aarch64")]
//     {
//         let mul8_lut = lut::MUL8[log_m as usize];
//         for (x, y) in x.chunks_exact_mut(64).zip(y.chunks_exact(64)) {
//             for i in 0..64 {
//                 x[i] ^= mul8_lut[y[i] as usize];
//             }
//         }
//         return;
//     }

//     let (prefix, x8, suffix) = unsafe { x.align_to_mut::<u64>() };
//     debug_assert!(prefix.is_empty());
//     debug_assert!(suffix.is_empty());
// }

// z = x + y (mod Modulus)
#[inline]
const fn add_mod(a: u8, b: u8) -> u8 {
    let sum = a as u32 + b as u32;

    // Partial reduction step, allowing for kModulus to be returned
    (sum + (sum >> BITS)) as u8
}

// z = x - y (mod Modulus)
#[inline]
const fn sub_mod(a: u8, b: u8) -> u8 {
    let b = if a < b { b + 1 } else { b };
    // make sure we don't underflow
    let a = a as u32 + ORDER as u32;
    let dif = a - b as u32;

    dif as u8
}

// Note that this operation is not a normal multiplication in a finite
// field because the right operand is already a logarithm.  This is done
// because it moves K table lookups from the Decode() method into the
// initialization step that is less performance critical.  The LogWalsh[]
// table below contains precalculated logarithms so it is easier to do
// all the other multiplies in that form as well.
#[inline]
fn mul_log(a: u8, log_b: u8) -> u8 {
    if a == 0 {
        0
    } else {
        let log_a = lut::log(a);
        lut::exp(add_mod(log_a, log_b))
    }
}

fn mul_add(x: &mut [u8], y: &[u8], log_m: u8) {
    x.iter_mut().zip(y.iter()).for_each(|(x, y)| {
        *x ^= lut::mul(*y, log_m);
    })
}

// Decimation in time (DIT) Fast Walsh-Hadamard Transform
// Unrolls pairs of layers to perform cross-layer operations in registers
// mtrunc: Number of elements that are non-zero at the front of data
fn fwht(data: &mut [u8; ORDER], m: usize, mtrunc: usize) {
    // Decimation in time: Unroll 2 layers at a time
    let mut dist: usize = 1;
    let mut dist4: usize = 4;

    while dist4 <= m {
        for offset in (0..mtrunc).step_by(dist4) {
            let mut offset = offset;

            for _ in 0..dist {
                // TODO: maybe in rust not faster
                // TODO: bound checks
                // fwht4(data[i:], dist) inlined...
                // Reading values appear faster than updating pointers.
                // Casting to uint is not faster.
                let t0 = data[offset];
                let t1 = data[offset + dist];
                let t2 = data[offset + dist * 2];
                let t3 = data[offset + dist * 3];

                let (t0, t1) = fwht2alt(t0, t1);
                let (t2, t3) = fwht2alt(t2, t3);
                let (t0, t2) = fwht2alt(t0, t2);
                let (t1, t3) = fwht2alt(t1, t3);

                data[offset] = t0;
                data[offset + dist] = t1;
                data[offset + dist * 2] = t2;
                data[offset + dist * 3] = t3;

                offset += 1
            }
        }
        dist = dist4;
        dist4 <<= 2;
    }

    // If there is one layer left:
    if dist < m {
        for i in 0..dist {
            let (first, second) = data.split_at_mut(i + 1);
            fwht2(&mut first[i], &mut second[dist]);
        }
    }
}

// {a, b} = {a + b, a - b} (Mod Q)
#[inline]
fn fwht2(a: &mut u8, b: &mut u8) {
    let sum = add_mod(*a, *b);
    let dif = sub_mod(*a, *b);

    *a = sum;
    *b = dif;
}

// fwht2 not in place
#[inline]
fn fwht2alt(a: u8, b: u8) -> (u8, u8) {
    (add_mod(a, b), sub_mod(a, b))
}

#[inline]
const fn ceil_pow2(x: usize) -> usize {
    let bitwidth = usize::BITS;
    1 << (bitwidth - (x - 1).leading_zeros())
}

// Unrolled IFFT for encoder
fn ifft_dit_encoder(
    data: &mut [&mut [u8]],
    mtrunc: usize,
    work: &mut [&mut [u8]],
    xor_output: Option<&mut [&mut [u8]]>,
    m: usize,
    skew_lut: &[u8],
) {
    // NOTE: may not be valid in rust
    // I tried rolling the memcpy/memset into the first layer of the FFT and
    // found that it only yields a 4% performance improvement, which is not
    // worth the extra complexity.
    for i in 0..mtrunc {
        data[i].copy_from_slice(work[i]);
    }
    for row in work[mtrunc..m].iter_mut() {
        row.fill(0);
    }

    // Decimation in time: Unroll 2 layers at a time
    let mut dist = 1;
    let mut dist4 = 4;

    while dist4 <= m {
        for r in (0..mtrunc).step_by(dist4) {
            let iend = r + dist;
            let log_m01 = skew_lut[iend];
            let log_m02 = skew_lut[iend + dist];
            let log_m23 = skew_lut[iend + dist * 2];

            // For each set of dist elements:
            for i in r..iend {
                ifft_dit4(&mut work[i..], dist, log_m01, log_m23, log_m02);
            }
        }

        dist = dist4;
        dist4 <<= 2;
        // orig:
        // I tried alternating sweeps left->right and right->left to reduce cache misses.
        // It provides about 1% performance boost when done for both FFT and IFFT, so it
        // does not seem to be worth the extra complexity.
    }

    // If there is one layer left:
    if dist < m {
        // assume that dist = m / 2
        assert_eq!(dist * 2, m);

        let log_m = lut::fft_skew(dist as u8);

        if log_m == MODULUS {
            let (input, output) = work.split_at_mut(dist);
            slices_xor(&mut output[..dist], input.iter_mut().map(|elem| &**elem));
        } else {
            let (x, y) = work.split_at_mut(dist);
            for i in 0..dist {
                ifft_dit2(x[i], y[i], log_m);
            }
        }
    }

    // orig:
    // I tried unrolling this but it does not provide more than 5% performance
    // improvement for 16-bit finite fields, so it's not worth the complexity.
    if let Some(xor_output) = xor_output {
        slices_xor(
            &mut xor_output[..m],
            work[..m].iter_mut().map(|elem| &**elem),
        );
    }
}

fn ifft_dit4(work: &mut [&mut [u8]], dist: usize, log_m01: u8, log_m23: u8, log_m02: u8) {
    if work[0].is_empty() {
        return;
    }

    // TODO: support AVX2 if enabled
    // No simd version

    let (dist0, dist1) = work.split_at_mut(dist);
    let (dist1, dist2) = dist1.split_at_mut(dist);
    let (dist2, dist3) = dist2.split_at_mut(dist);

    // First layer:
    if log_m01 == MODULUS {
        slice_xor(&*dist0[0], dist1[0]);
    } else {
        ifft_dit2(dist0[0], dist1[0], log_m01);
    }

    if log_m23 == MODULUS {
        slice_xor(&*dist2[0], dist3[0]);
    } else {
        ifft_dit2(dist2[0], dist3[0], log_m23);
    }

    // Second layer:
    if log_m02 == MODULUS {
        slice_xor(&*dist0[0], dist2[0]);
        slice_xor(&*dist1[0], dist2[0]);
    } else {
        ifft_dit2(dist0[0], dist2[0], log_m02);
        ifft_dit2(dist1[0], dist3[0], log_m02);
    }
}

fn ifft_dit2(x: &mut [u8], y: &mut [u8], log_m: u8) {
    slice_xor(&*x, y);
    mul_add(x, y, log_m);
}

// In-place FFT for encoder and decoder
fn fft_dit(work: &mut [&mut [u8]], mtrunc: usize, m: usize, skew_lut: &[u8]) {
    // Decimation in time: Unroll 2 layers at a time
    let mut dist4 = m;
    let mut dist = m >> 2;

    while dist != 0 {
        // For each set of dist*4 elements:
        for r in (0..mtrunc).step_by(dist4) {
            let iend = r + dist;
            let log_m01 = skew_lut[iend - 1];
            let log_m02 = skew_lut[iend + dist - 1];
            let log_m23 = skew_lut[iend + 2 * dist - 1];

            // For each set of dist elements:
            for i in r..iend {
                fft_dit4(&mut work[i..], dist, log_m01, log_m23, log_m02);
            }
        }

        dist4 = dist;
        dist >>= 2;
    }

    // If there is one layer left:
    if dist4 == 2 {
        for r in (0..mtrunc).step_by(2) {
            let log_m = skew_lut[r];
            let (x, y) = work.split_at_mut(r + 1);

            if log_m == MODULUS {
                slice_xor(&*x[r], y[0]);
            } else {
                fft_dit2(x[r], y[0], log_m);
            }
        }
    }
}

// 4-way butterfly
fn fft_dit4(work: &mut [&mut [u8]], dist: usize, log_m01: u8, log_m23: u8, log_m02: u8) {
    if work[0].is_empty() {
        return;
    }

    // TODO: support AVX2 if enabled
    // No simd version

    // First layer:
    let (dist0, dist1) = work.split_at_mut(dist);
    let (dist1, dist2) = dist1.split_at_mut(dist);
    let (dist2, dist3) = dist2.split_at_mut(dist);

    // First layer:
    if log_m02 == MODULUS {
        slice_xor(&*dist0[0], dist2[0]);
        slice_xor(&*dist1[0], dist3[0]);
    } else {
        fft_dit2(dist0[0], dist2[0], log_m02);
        fft_dit2(dist1[0], dist3[0], log_m02);
    }

    // Second layer:
    if log_m01 == MODULUS {
        slice_xor(&*dist0[0], dist1[0]);
    } else {
        fft_dit2(dist0[0], dist1[0], log_m01);
    }

    if log_m23 == MODULUS {
        slice_xor(&*dist2[0], dist3[0]);
    } else {
        ifft_dit2(dist2[0], dist3[0], log_m23);
    }
}

// 2-way butterfly forward
fn fft_dit2(x: &mut [u8], y: &mut [u8], log_m: u8) {
    if x.is_empty() {
        return;
    }

    mul_add(x, y, log_m);
    slice_xor(&*x, y);
}

fn slices_xor(output: &mut [&mut [u8]], input: impl Iterator<Item = impl Buf>) {
    output
        .iter_mut()
        .zip(input)
        .for_each(|(out, inp)| slice_xor(inp, out));
}

fn slice_xor(mut input: impl Buf, mut output: &mut [u8]) {
    // TODO: this unroll is inherited from go code, however it might not be needed
    while output.remaining_mut() >= 32 && input.remaining() >= 32 {
        let mut output_buf = &*output;
        let v0 = output_buf.get_u64_le() ^ input.get_u64_le();
        let v1 = output_buf.get_u64_le() ^ input.get_u64_le();
        let v2 = output_buf.get_u64_le() ^ input.get_u64_le();
        let v3 = output_buf.get_u64_le() ^ input.get_u64_le();

        output.put_u64_le(v0);
        output.put_u64_le(v1);
        output.put_u64_le(v2);
        output.put_u64_le(v3);
    }

    let rest = output.remaining_mut().min(input.remaining());
    for _ in 0..rest {
        let xor = (&*output).get_u8() ^ input.get_u8();
        output.put_u8(xor);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encoding() {
        let leo = LeopardFF8::new(1, 1).unwrap();

        let mut input = include!("../test_data/encode-ff8-1-64-input");
        let expected = include!("../test_data/encode-ff8-1-64-output");

        println!("{input:?}");

        let mut input_ref: Vec<_> = input.iter_mut().map(|shard| shard.as_mut_slice()).collect();

        leo.encode(&mut input_ref).unwrap();

        assert_eq!(input, expected);
    }
}
