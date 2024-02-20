use crate::{LeopardError, Result};

const ALIGN: usize = 64;

#[repr(C, align(64))]
#[derive(Clone)]
struct Aligned([u8; ALIGN]);

impl Aligned {
    fn zeroed() -> Self {
        Self([0; ALIGN])
    }
}

#[allow(unused)]
pub fn is_aligned(shards: &[&[u8]]) -> bool {
    !shards.is_empty()
        && shards
            .iter()
            .all(|shard| shard.as_ptr().align_offset(ALIGN) != 0)
}

pub fn alloc_aligned(shards: usize, shard_size: usize) -> Result<Vec<u8>> {
    if shard_size % ALIGN != 0 {
        return Err(LeopardError::InvalidShardSize(shard_size));
    }

    // create a vec of Aligned to make sure they are aligned properly
    let aligned_mem = vec![Aligned::zeroed(); shards * shard_size / ALIGN];

    // grab the raw parts
    let length = aligned_mem.len() * ALIGN;
    let capacity = aligned_mem.capacity() * ALIGN;
    let mem_ptr = aligned_mem.leak().as_mut_ptr() as *mut u8;

    // SAFETY: this is safe because original vector points to [u8; ALIGN]
    // and SHARD_SIZE is a multiple of align. And both the length and capacity
    // are the same size in bytes.
    // recreate the vec over shards
    Ok(unsafe { Vec::from_raw_parts(mem_ptr, length, capacity) })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn size_of_aligned() {
        assert_eq!(std::mem::size_of::<Aligned>(), ALIGN);
    }

    #[test]
    fn shard_size_must_be_multiple_of_align() {
        alloc_aligned(10, 1).unwrap_err();
        alloc_aligned(10, 10).unwrap_err();
        alloc_aligned(10, 32).unwrap_err();
        alloc_aligned(10, 63).unwrap_err();
        alloc_aligned(10, 127).unwrap_err();
        alloc_aligned(10, 257).unwrap_err();

        alloc_aligned(10, 64).unwrap();
        alloc_aligned(10, 128).unwrap();
        alloc_aligned(10, 192).unwrap();
        alloc_aligned(10, 256).unwrap();
    }

    fn is_properly_aligned(shards: usize, shard_size: usize) {
        let mem = alloc_aligned(shards, shard_size).unwrap();

        for shard in mem.chunks(shard_size) {
            assert_eq!(shard.as_ptr().align_offset(ALIGN), 0);
            assert!(shard.iter().all(|&byte| byte == 0));
        }
    }

    #[test]
    fn check_alignment() {
        for shards in [2, 3, 7, 13, 21, 50, 77, 100] {
            is_properly_aligned(shards, 64);
            is_properly_aligned(shards, 128);
            is_properly_aligned(shards, 192);
            is_properly_aligned(shards, 256);
            is_properly_aligned(shards, 512);
            is_properly_aligned(shards, 1024);
        }
    }
}
