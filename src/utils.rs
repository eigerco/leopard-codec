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

pub fn is_aligned(shards: &[&[u8]]) -> bool {
    !shards.is_empty()
        && shards
            .iter()
            .all(|shard| shard.as_ptr().align_offset(ALIGN) != 0)
}

pub fn alloc_aligned<const SHARD_SIZE: usize>(shards: usize) -> Result<Vec<[u8; SHARD_SIZE]>> {
    if SHARD_SIZE % ALIGN != 0 {
        return Err(LeopardError::InvalidShardSize(SHARD_SIZE));
    }

    let ratio = SHARD_SIZE / ALIGN;

    // create a vec of Aligned to make sure they are aligned properly
    let aligned_mem = vec![Aligned::zeroed(); shards * ratio];

    // grab the raw parts
    let length = aligned_mem.len() / ratio;
    let capacity = aligned_mem.capacity() / ratio;
    let mem_ptr = aligned_mem.leak().as_mut_ptr() as *mut [u8; SHARD_SIZE];

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
        alloc_aligned::<1>(10).unwrap_err();
        alloc_aligned::<10>(10).unwrap_err();
        alloc_aligned::<32>(10).unwrap_err();
        alloc_aligned::<63>(10).unwrap_err();
        alloc_aligned::<127>(10).unwrap_err();
        alloc_aligned::<257>(10).unwrap_err();

        alloc_aligned::<64>(10).unwrap();
        alloc_aligned::<128>(10).unwrap();
        alloc_aligned::<192>(10).unwrap();
        alloc_aligned::<256>(10).unwrap();
    }

    fn is_properly_aligned<const S: usize>(size: usize) {
        let mem = alloc_aligned::<S>(size).unwrap();

        for shard in mem.iter() {
            assert_eq!(shard.as_ptr().align_offset(ALIGN), 0);
            assert!(shard.iter().all(|&byte| byte == 0));
        }
    }

    #[test]
    fn check_alignment() {
        for size in [2, 3, 7, 13, 21, 50, 77, 100] {
            is_properly_aligned::<64>(size);
            is_properly_aligned::<128>(size);
            is_properly_aligned::<192>(size);
            is_properly_aligned::<256>(size);
            is_properly_aligned::<512>(size);
            is_properly_aligned::<1024>(size);
        }
    }
}
