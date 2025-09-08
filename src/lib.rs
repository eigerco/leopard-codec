#![doc = include_str!("../README.md")]

use bytes::{Buf, BufMut};
use thiserror::Error;

/// Number of bits per element
pub const BITS: usize = 8;
/// Finite field order: Number of elements in the field
pub const ORDER: usize = u8::MAX as usize + 1;
/// Modulus for field operations
pub const MODULUS: u8 = u8::MAX;
/// LFSR Polynomial that generates the field elements
pub const POLYNOMIAL: usize = 0x11D;
/// Basis used for generating logarithm tables
pub const CANTOR_BASIS: [u8; BITS] = [1, 214, 152, 146, 86, 200, 88, 230];

/// lookup tables
mod lut;

/// Possible errors that can happen when interacting with Leopard.
#[derive(Debug, Error)]
pub enum LeopardError {
    /// Maximum number of shards exceeded.
    #[error("Maximum shard number ({order}) exceeded: {0}", order=ORDER)]
    MaxShardNumberExceeded(usize),

    /// Maximum number of parity shards exceeded.
    #[error("Maximum parity shard number ({0}) exceeded: {1}")]
    MaxParityShardNumberExceeded(usize, usize),

    /// This amount of (data, parity) shards is not supported by Leopard algorithm
    /// and would result in buffer overflow on skew lookup table during encoding.
    /// Please try using different amounts of shards.
    #[error("Unsupported number of data ({0}) and parity ({1}) shards")]
    UnsupportedShardsAmounts(usize, usize),

    /// Some shards contain no data.
    #[error("Shards contain no data")]
    EmptyShards,

    /// Some shards are of different lengths.
    #[error("Shards of different lengths found")]
    UnequalShardsLengths,

    /// Shard size is invalid.
    #[error("Shard size ({0}) should be a multiple of 64")]
    InvalidShardSize(usize),

    /// To few shards to reconstruct data.
    #[error("Too few shards ({0}) to reconstruct data, at least {1} needed")]
    TooFewShards(usize, usize),
}

/// A result type with [`LeopardError`].
pub type Result<T, E = LeopardError> = std::result::Result<T, E>;

/// Encode parity data into given shards.
///
/// The first `data_shards` shards will be the treated as data shards
/// and the rest as parity shards.
///
/// # Errors
///
/// If too many shards provided or shards were of incorrect or different lengths.
pub fn encode(shards: &mut [impl AsMut<[u8]>], data_shards: usize) -> Result<()> {
    if shards.len() > ORDER {
        return Err(LeopardError::MaxShardNumberExceeded(shards.len()));
    }
    let parity_shards = shards.len() - data_shards;
    if parity_shards > data_shards {
        return Err(LeopardError::MaxParityShardNumberExceeded(
            parity_shards,
            data_shards,
        ));
    }
    if is_encode_buf_overflow(data_shards, parity_shards) {
        return Err(LeopardError::UnsupportedShardsAmounts(
            data_shards,
            parity_shards,
        ));
    }

    let mut shards: Vec<&mut [u8]> = shards.iter_mut().map(|shard| shard.as_mut()).collect();
    let shard_size = check_shards(&shards, false)?;

    if shard_size % 64 != 0 {
        return Err(LeopardError::InvalidShardSize(shard_size));
    }

    encode_inner(&mut shards, data_shards, shard_size);

    Ok(())
}

fn encode_inner(shards: &mut [&mut [u8]], data_shards: usize, shard_size: usize) {
    let parity_shards = shards.len() - data_shards;

    let m = ceil_pow2(parity_shards);
    let mtrunc = m.min(data_shards);

    // 'work' is a temporary buffer where the parity shards are computed.
    // the first half of it is where the resulting parity shards will end up.
    let mut work_mem = vec![0; 2 * m * shard_size];
    let mut work: Vec<_> = work_mem.chunks_exact_mut(shard_size).collect();

    let skew_lut = &lut::FFT_SKEW[m - 1..];

    // copy the input to the work table
    for (shard, work) in shards[data_shards..].iter().zip(work.iter_mut()) {
        work.copy_from_slice(shard);
    }

    ifft_dit_encoder(
        &shards[..data_shards],
        mtrunc,
        &mut work,
        None, // No xor output
        m,
        skew_lut,
    );

    let last_count = data_shards % m;

    // goto skip_body
    if m < data_shards {
        let (xor_out, work) = work.split_at_mut(m);
        let mut n = m;

        // for sets of m data pieces
        while n <= data_shards - m {
            // work <- work xor IFFT(data + i, m, m + i)
            ifft_dit_encoder(&shards[n..], m, work, Some(xor_out), m, &skew_lut[n..]);
            n += m;
        }

        // Handle final partial set of m pieces:
        if last_count != 0 {
            ifft_dit_encoder(
                &shards[n..],
                last_count,
                work,
                Some(xor_out),
                m,
                &skew_lut[n..],
            );
        }
    }

    // work <- FFT(work, m, 0)
    fft_dit(&mut work, parity_shards, m, &*lut::FFT_SKEW);

    for (shard, work) in shards[data_shards..].iter_mut().zip(work.iter()) {
        shard.copy_from_slice(work);
    }
}

/// Leopard algorithm is imperfect and can result in a buffer overflow on the 'lut::FFT_SKEW'.
/// Encoding happens in passes where each pass encodes parity for data shards in chunks
/// of the smallest power of 2 bigger (or equal) than parity shards amount.
/// If the last chunk is not a full pass, we can hit the overflow for some pairs
/// of (data_shards, parity_shards).
/// This function detects such inputs.
fn is_encode_buf_overflow(data_shards: usize, parity_shards: usize) -> bool {
    debug_assert!(data_shards >= parity_shards);
    debug_assert!(data_shards + parity_shards <= ORDER);

    let m = ceil_pow2(parity_shards);
    let last_count = data_shards % m;

    // we can finish encoding with only full passes
    if m >= data_shards || last_count == 0 {
        return false;
    }

    let full_passes = data_shards / m;
    // if this is 'true', we would overflow the fft skew table which has the size of `MODULUS`
    (full_passes + 1) * m + 1 > MODULUS as usize
}

/// Reconstructs original shards from the provided slice.
///
/// The shards which are missing should be provided as empty `Vec`s.
///
/// Reconstruction can only happen if the amount of present data and parity shards
/// is equal or greater than the `data_shards`.
///
/// The first `data_shards` shards will be the treated as data shards
/// and the rest as parity shards.
///
/// # Errors
///
/// If too few shards are present to reconstruct original data or shards were of incorrect or different lengths.
pub fn reconstruct(shards: &mut [impl AsMut<Vec<u8>>], data_shards: usize) -> Result<()> {
    if shards.len() > ORDER {
        return Err(LeopardError::MaxShardNumberExceeded(shards.len()));
    }
    let parity_shards = shards.len() - data_shards;
    if parity_shards > data_shards {
        return Err(LeopardError::MaxParityShardNumberExceeded(
            parity_shards,
            data_shards,
        ));
    }

    let mut shards: Vec<_> = shards.iter_mut().map(|shard| shard.as_mut()).collect();
    let shard_size = check_shards(&shards, true)?;

    let present_shards = shards.iter().filter(|shard| !shard.is_empty()).count();
    if present_shards == shards.len() {
        // all shards present, nothing to do
        return Ok(());
    }

    // Check if we have enough to reconstruct.
    if present_shards < data_shards {
        return Err(LeopardError::TooFewShards(present_shards, data_shards));
    }

    if shard_size % 64 != 0 {
        return Err(LeopardError::InvalidShardSize(shard_size));
    }

    reconstruct_inner(&mut shards, data_shards, shard_size);

    Ok(())
}

fn reconstruct_inner(shards: &mut [&mut Vec<u8>], data_shards: usize, shard_size: usize) {
    let parity_shards = shards.len() - data_shards;

    // TODO: errorbitfields for avoiding unnecessary fft steps
    // orig:
    // Use only if we are missing less than 1/4 parity,
    // And we are restoring a significant amount of data.
    // useBits := r.totalShards-numberPresent <= r.parityShards/4 && shardSize*r.totalShards >= 64<<10

    let m = ceil_pow2(parity_shards);
    let n = ceil_pow2(m + data_shards);

    // save the info which shards were empty
    let empty_shards_mask: Vec<_> = shards.iter().map(|shard| shard.is_empty()).collect();
    // and recreate them
    for shard in shards.iter_mut().filter(|shard| shard.is_empty()) {
        shard.resize(shard_size, 0);
    }

    let mut err_locs = [0u8; ORDER];

    for (&is_empty, err_loc) in empty_shards_mask
        .iter()
        .skip(data_shards)
        .zip(err_locs.iter_mut())
    {
        if is_empty {
            *err_loc = 1;
        }
    }

    for err in &mut err_locs[parity_shards..m] {
        *err = 1;
    }

    for (&is_empty, err_loc) in empty_shards_mask
        .iter()
        .take(data_shards)
        .zip(err_locs[m..].iter_mut())
    {
        if is_empty {
            *err_loc = 1;
        }
    }

    // TODO: No inversion...

    // Evaluate error locator polynomial8
    fwht(&mut err_locs, ORDER, m + data_shards);

    for (err, &log_walsh) in err_locs.iter_mut().zip(lut::LOG_WALSH.iter()) {
        let mul = (*err) as usize * log_walsh as usize;
        *err = (mul % MODULUS as usize) as u8;
    }

    fwht(&mut err_locs, ORDER, ORDER);

    let mut work_mem = vec![0u8; shard_size * n];
    let mut work: Vec<_> = work_mem.chunks_exact_mut(shard_size).collect();

    for i in 0..parity_shards {
        if !empty_shards_mask[i + data_shards] {
            mul_gf(work[i], shards[i + data_shards], err_locs[i]);
        } else {
            work[i].fill(0);
        }
    }
    for work in work.iter_mut().take(m).skip(parity_shards) {
        work.fill(0);
    }

    // work <- original data
    for i in 0..data_shards {
        if !empty_shards_mask[i] {
            mul_gf(work[m + i], shards[i], err_locs[m + i])
        } else {
            work[m + i].fill(0);
        }
    }
    for work in work.iter_mut().take(n).skip(m + data_shards) {
        work.fill(0);
    }

    // work <- IFFT(work, n, 0)
    ifft_dit_decoder(m + data_shards, &mut work, n, &lut::FFT_SKEW[..]);

    // work <- FormalDerivative(work, n)
    for i in 1..n {
        let width = ((i ^ (i - 1)) + 1) >> 1;
        let (output, input) = work.split_at_mut(i);
        slices_xor(
            &mut output[i - width..],
            input.iter_mut().map(|elem| &**elem),
        );
    }

    // work <- FFT(work, n, 0) truncated to m + dataShards
    fft_dit(&mut work, m + data_shards, n, &lut::FFT_SKEW[..]);

    // Reveal erasures
    //
    //  Original = -ErrLocator * FFT( Derivative( IFFT( ErrLocator * ReceivedData ) ) )
    //  mul_mem(x, y, log_m, ) equals x[] = y[] * log_m
    //
    // mem layout: [Recovery Data (Power of Two = M)] [Original Data (K)] [Zero Padding out to N]
    for (i, shard) in shards.iter_mut().enumerate() {
        if !empty_shards_mask[i] {
            continue;
        }

        if i >= data_shards {
            // parity shard
            mul_gf(
                shard,
                work[i - data_shards],
                MODULUS - err_locs[i - data_shards],
            );
        } else {
            // data shard
            mul_gf(shard, work[i + m], MODULUS - err_locs[i + m]);
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
    let are_all_same_size = shards.iter().all(|shard| {
        let shard = shard.as_ref();
        if allow_zero && shard.is_empty() {
            true
        } else {
            shard.len() == size
        }
    });

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
    let b = if a < b { b as u32 + 1 } else { b as u32 };
    // make sure we don't underflow
    let a = a as u32 + ORDER as u32;
    let dif = a - b;

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

fn mul_gf(out: &mut [u8], input: &[u8], log_m: u8) {
    let mul_lut = lut::MUL[log_m as usize];
    for (out, &input) in out.iter_mut().zip(input.iter()) {
        *out = mul_lut[input as usize];
    }
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
            // NOTE: this is compatible with klauspost/reedsolomon but diverages from catid/leopard
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
        debug_assert_eq!(dist * 2, m);

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

    // NOTE: this is compatible with klauspost/reedsolomon but diverages from catid/leopard
    if let Some(xor_output) = xor_output {
        slices_xor(
            &mut xor_output[..m],
            work[..m].iter_mut().map(|elem| &**elem),
        );
    }
}

// Basic no-frills version for decoder
fn ifft_dit_decoder(mtrunc: usize, work: &mut [&mut [u8]], m: usize, skew_lut: &[u8]) {
    // Decimation in time: Unroll 2 layers at a time
    let mut dist = 1;
    let mut dist4 = 4;

    while dist4 <= m {
        // For each set of dist*4 elements:
        for r in (0..mtrunc).step_by(dist4) {
            let iend = r + dist;
            let log_m01 = skew_lut[iend - 1];
            let log_m02 = skew_lut[iend + dist - 1];
            let log_m23 = skew_lut[iend + 2 * dist - 1];

            // For each set of dist elements:
            for i in r..iend {
                ifft_dit4(&mut work[i..], dist, log_m01, log_m23, log_m02);
            }
        }

        dist = dist4;
        dist4 <<= 2;
    }

    // If there is one layer left:
    if dist < m {
        // Assuming that dist = m / 2
        debug_assert_eq!(2 * dist, m);

        let log_m = skew_lut[dist - 1];

        if log_m == MODULUS {
            let (input, output) = work.split_at_mut(dist);
            slices_xor(&mut output[..dist], input.iter_mut().map(|elem| &**elem));
        } else {
            let (x, y) = work.split_at_mut(dist);
            for i in 0..dist {
                ifft_dit2(x[i], y[i], log_m)
            }
        }
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
    use std::panic::catch_unwind;

    use rand::{seq::index, Fill, Rng};
    use test_strategy::{proptest, Arbitrary};

    use super::*;

    #[proptest]
    fn go_reedsolomon_encode_compatibility(input: TestCase) {
        let TestCase {
            data_shards,
            parity_shards,
            shard_size,
        } = input;
        let total_shards = data_shards + parity_shards;
        let test_shards = random_shards(total_shards, shard_size);

        let mut shards = test_shards.clone();
        encode(&mut shards, data_shards).unwrap();

        let mut expected = test_shards;
        go_leopard::encode(&mut expected, data_shards, shard_size).unwrap();

        if expected != shards {
            panic!("Go and Rust encoding differ for {input:#?}")
        }
    }

    #[proptest]
    fn encode_reconstruct(input: TestCase) {
        let TestCase {
            data_shards,
            parity_shards,
            shard_size,
        } = input;
        let total_shards = data_shards + parity_shards;
        let mut shards = random_shards(total_shards, shard_size);

        encode(&mut shards, data_shards).unwrap();

        let expected = shards.clone();

        let mut rng = rand::rng();
        let missing_shards = rng.random_range(1..=parity_shards);
        for idx in index::sample(&mut rng, total_shards, missing_shards) {
            shards[idx] = vec![];
        }

        reconstruct(&mut shards, data_shards).unwrap();

        if expected != shards {
            panic!("shares differ after reconstruction");
        }
    }

    #[test]
    fn overflow_detection() {
        for data_shards in 1..MODULUS as usize {
            for parity_shards in 1..data_shards {
                let total_shards = data_shards + parity_shards;

                // too many shards
                if total_shards > ORDER {
                    continue;
                }

                let overflow = is_encode_buf_overflow(data_shards, parity_shards);

                let result = catch_unwind(|| {
                    let mut shards = random_shards(total_shards, 64);
                    let mut shards_ref: Vec<_> = shards
                        .iter_mut()
                        .map(|shard| shard.as_mut_slice())
                        .collect();
                    encode_inner(&mut shards_ref, data_shards, 64);
                });

                assert_eq!(result.is_err(), overflow, "{data_shards} {parity_shards}");
            }
        }
    }

    #[derive(Arbitrary, Debug)]
    #[filter(!is_encode_buf_overflow(#data_shards, #parity_shards))]
    struct TestCase {
        #[strategy(1..ORDER - 1)]
        data_shards: usize,

        #[strategy(1..=(ORDER - #data_shards).min(#data_shards))]
        parity_shards: usize,

        #[strategy(1usize..1024)]
        #[map(|x| x * 64)]
        shard_size: usize,
    }

    fn random_shards(shards: usize, shard_size: usize) -> Vec<Vec<u8>> {
        let mut rng = rand::rng();
        (0..shards)
            .map(|_| {
                let mut shard = vec![0; shard_size];
                Fill::fill(shard.as_mut_slice(), &mut rng);
                shard
            })
            .collect()
    }
}
