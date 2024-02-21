use bytes::{Buf, BufMut};
use thiserror::Error;

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
mod lut;

/// Possible errors that can happen when interacting with Leopard.
#[derive(Debug, Error)]
pub enum LeopardError {
    /// Maximum number of shards exceeded.
    #[error("Maximum shard number ({}) exceeded: {0}", ORDER)]
    MaxShardNumberExceeded(usize),

    /// Some shards contain no data.
    #[error("Shards contain no data")]
    EmptyShards,

    /// Some shards are of different lengths.
    #[error("Shards of different lengths found")]
    UnequalShardsLengths,

    /// Shard size is invalid.
    #[error("Shard size ({0}) should be a multiple of 64")]
    InvalidShardSize(usize),
}

/// A result type with [`LeopardError`].
pub type Result<T, E = LeopardError> = std::result::Result<T, E>;

/// Encode parity data into given shards.
///
/// The first `data_shards` shards will be the treated as data shards
/// and the rest as a parity shards.
///
/// # Errors
///
/// If too many shards provided or shards were of incorrect or different lengths.
pub fn encode(shards: &mut [&mut [u8]], data_shards: usize) -> Result<()> {
    if shards.len() > ORDER {
        return Err(LeopardError::MaxShardNumberExceeded(shards.len()));
    }

    let shard_size = check_shards(shards, false)?;
    if shard_size % 64 != 0 {
        return Err(LeopardError::InvalidShardSize(shard_size));
    }

    encode_inner(shards, data_shards, shard_size);

    Ok(())
}

fn encode_inner(shards: &mut [&mut [u8]], data_shards: usize, shard_size: usize) {
    let parity_shards = shards.len() - data_shards;

    let m = ceil_pow2(parity_shards);
    let mtrunc = m.min(data_shards);

    let mut work_mem = vec![0; 2 * m * shard_size];
    let mut work: Vec<_> = work_mem.chunks_exact_mut(shard_size).collect();

    let mut xor_out_mem = vec![0; 2 * m * shard_size];
    let mut xor_out: Vec<_> = xor_out_mem.chunks_exact_mut(shard_size).collect();

    let skew_lut = &lut::FFT_SKEW[m - 1..];

    // ceil division
    for offset in (0..shard_size).step_by(BLOCK_SIZE) {
        let end = (offset + BLOCK_SIZE).min(shard_size);
        let mut shards_chunk = shards
            .iter_mut()
            .map(|shard| &mut shard[offset..end])
            .collect::<Vec<_>>();
        let mut work_chunk = work
            .iter_mut()
            .map(|shard| &mut shard[offset..end])
            .collect::<Vec<_>>();
        let mut xor_out_chunk = xor_out
            .iter_mut()
            .map(|shard| &mut shard[offset..end])
            .collect::<Vec<_>>();

        // copy the input to the work table
        for (shard, work) in shards_chunk[data_shards..]
            .iter()
            .zip(work_chunk.iter_mut())
        {
            work.copy_from_slice(shard);
        }

        ifft_dit_encoder(
            &shards_chunk[..data_shards],
            mtrunc,
            &mut work_chunk,
            None, // No xor output
            m,
            skew_lut,
        );

        // copy the work back to input
        for (shard, work) in shards_chunk[data_shards..]
            .iter_mut()
            .zip(work_chunk.iter())
        {
            shard.copy_from_slice(work);
        }

        let last_count = data_shards % m;
        let mut skew_lut2 = skew_lut;

        // goto skip_body
        if m < data_shards {
            // for sets of m data pieces
            for _ in (m..data_shards - m).step_by(m) {
                // work <- work xor IFFT(data + i, m, m + i)
                skew_lut2 = &skew_lut2[m..];
                ifft_dit_encoder(
                    &shards_chunk[m..],
                    m,
                    &mut work_chunk[m..],
                    Some(&mut xor_out_chunk),
                    m,
                    &skew_lut2[m..],
                );

                // copy the xor to work
                for (xor, work) in xor_out_chunk[..m].iter().zip(work_chunk.iter_mut()) {
                    work.copy_from_slice(xor);
                }

                // copy the work to input
                for (shard, work) in shards_chunk[m..].iter_mut().zip(work_chunk[m..].iter()) {
                    shard.copy_from_slice(work);
                }
            }

            // Handle final partial set of m pieces:
            if last_count != 0 {
                ifft_dit_encoder(
                    &shards_chunk[m..],
                    m,
                    &mut work_chunk[m..],
                    Some(&mut xor_out_chunk),
                    m,
                    &skew_lut2[m..],
                );

                // copy the xor to work
                for (xor, work) in xor_out_chunk[..m].iter().zip(work_chunk.iter_mut()) {
                    work.copy_from_slice(xor);
                }

                // copy the work to input
                for (shard, work) in shards_chunk[m..].iter_mut().zip(work_chunk[m..].iter()) {
                    shard.copy_from_slice(work);
                }
            }
        }

        // work <- FFT(work, m, 0)
        fft_dit(&mut work_chunk, parity_shards, m, &*lut::FFT_SKEW);

        for (shard, work) in shards_chunk[data_shards..]
            .iter_mut()
            .zip(work_chunk.iter())
        {
            shard.copy_from_slice(work);
        }
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
//
// returns the length of the single shard
fn check_shards(shards: &[impl AsRef<[u8]>], allow_zero: bool) -> Result<usize> {
    let size = shard_size(shards);

    if size == 0 {
        if allow_zero {
            return Ok(0);
        } else {
            return Err(LeopardError::EmptyShards);
        }
    }

    // NOTE: for happy case fold would be faster
    let are_all_same_size = shards.iter().all(|shard| shard.as_ref().len() == size);

    if !are_all_same_size {
        return Err(LeopardError::UnequalShardsLengths);
    }

    Ok(size)
}

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
    data: &[impl AsRef<[u8]>],
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
        work[i].copy_from_slice(data[i].as_ref());
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

        let log_m = skew_lut[dist];

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
        slice_xor(&*dist1[0], dist3[0]);
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
        fft_dit2(dist2[0], dist3[0], log_m23);
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
    use std::fs;

    use super::*;

    #[test]
    fn test_slice_xor() {
        let mut x = [1, 2, 3, 4, 5, 6, 7];
        let y = [7, 6, 5, 4, 3, 2, 1];

        slice_xor(&y[..], &mut x);

        assert_eq!(x, [6, 4, 6, 0, 6, 4, 6]);
    }

    #[test]
    fn test_encoding_small() {
        let mut input = include!("../test_data/encode-ff8-1-64-input");
        let expected = include!("../test_data/encode-ff8-1-64-output");

        let mut input_ref: Vec<_> = input.iter_mut().map(|shard| shard.as_mut_slice()).collect();

        encode(&mut input_ref, 1).unwrap();

        assert_eq!(input, expected);
    }

    #[test]
    fn test_encoding_big() {
        let mut input = include!("../test_data/encode-ff8-100-128-input");
        let expected = include!("../test_data/encode-ff8-100-128-output");

        let mut input_ref: Vec<_> = input.iter_mut().map(|shard| shard.as_mut_slice()).collect();

        encode(&mut input_ref, 100).unwrap();

        assert_eq!(input, expected);
    }

    #[test]
    fn test_ifft_dit_encoder() {
        let data_shards = 100;
        let mtrunc = 100;
        let m = 128;
        let skew_lut = &lut::FFT_SKEW[m - 1..];

        let input = include!("../test_data/ifft-encoder-input");
        let expected = include!("../test_data/ifft-encoder-output");

        let mut work = expected;
        for shard in work.iter_mut() {
            shard.fill(0);
        }
        let mut work_ref = work
            .iter_mut()
            .map(|shard| shard.as_mut_slice())
            .collect::<Vec<_>>();

        ifft_dit_encoder(
            &input[..data_shards],
            mtrunc,
            &mut work_ref,
            None,
            m,
            skew_lut,
        );

        assert_eq!(work, expected);
    }

    #[test]
    fn test_ifft_dit4() {
        let testcases = fs::read_to_string("test_data/ifft-dit-4").unwrap();

        let mut input: Vec<Vec<u8>> = vec![];
        let mut dist = 0usize;
        let mut log_m01 = 0u8;
        let mut log_m23 = 0u8;
        let mut log_m02 = 0u8;

        for (n, line) in testcases.lines().enumerate() {
            match n % 6 {
                0 => {
                    input = serde_json::from_str(line).unwrap();
                }
                1 => {
                    dist = line.parse().unwrap();
                }
                2 => {
                    log_m01 = line.parse().unwrap();
                }
                3 => {
                    log_m23 = line.parse().unwrap();
                }
                4 => {
                    log_m02 = line.parse().unwrap();
                }
                5 => {
                    let expected: Vec<Vec<u8>> = serde_json::from_str(line).unwrap();
                    let mut input_ref = input
                        .iter_mut()
                        .map(|shard| shard.as_mut_slice())
                        .collect::<Vec<_>>();

                    ifft_dit4(&mut input_ref, dist, log_m01, log_m23, log_m02);
                    assert_eq!(input, expected);
                }
                _ => unreachable!(),
            }
        }
    }
}
