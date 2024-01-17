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

const INVERSION8_BYTES: usize = 256 / 8;

/// lookup tables
pub mod lut;

mod utils;

/// Possible errors that can happen when interacting with Leopard.
#[derive(Debug, Error)]
pub enum LeopardError {
    #[error("Maximum shard number ({}) exceeded: {0}", u16::MAX)]
    MaxShardNumberExceeded(usize),
}

/// A result type with [`LeopardError`].
pub type Result<T, E = LeopardError> = std::result::Result<T, E>;

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

// Decimation in time (DIT) Fast Walsh-Hadamard Transform
// Unrolls pairs of layers to perform cross-layer operations in registers
// mtrunc: Number of elements that are non-zero at the front of data
fn fwht(data: &mut [u8; ORDER], m: usize, m_trunc: usize) {
    // Decimation in time: Unroll 2 layers at a time
    let mut dist: usize = 1;
    let mut dist4: usize = 4;

    while dist4 <= m {
        for offset in (0..m_trunc).step_by(dist4) {
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

#[cfg(test)]
mod tests {}
