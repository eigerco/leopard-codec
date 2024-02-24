//! Rust wrapper over Go's reedsolomon module

use std::ffi::c_uchar;

use thiserror::Error;

#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
#[allow(unused)]
mod go {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

/// Possible errors that can happen when interacting with `go-leopard`
#[derive(Debug, Error)]
pub enum Error {
    /// Creating codec failed
    #[error("Creating the codec failed")]
    CreatingCodecFailed,

    /// Encoding failed
    #[error("Encoding failed")]
    EncodingFailed,
}

/// Encode the shards using the Go's reedsolomon module
pub fn encode(
    shards: &mut [&mut [u8]],
    data_shards: usize,
    shard_size: usize,
) -> Result<(), Error> {
    let mut shards_out: Vec<_> = shards
        .iter_mut()
        .map(|shard| shard.as_mut_ptr() as *mut c_uchar)
        .collect();

    let res = unsafe {
        go::encode(
            shards_out.as_mut_ptr(),
            shard_size as u64,
            data_shards as u64,
            shards.len() as u64,
        )
    };

    match res {
        0 => Ok(()),
        -1 => Err(Error::CreatingCodecFailed),
        -2 => Err(Error::EncodingFailed),
        _ => unreachable!(),
    }
}
