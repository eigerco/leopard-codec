use crate::{LeopardError, Result};

const ALIGN: usize = 64;

pub(crate) fn alloc_zeroed_with_padding(size: usize) -> Vec<u8> {
    // in worst case, we must drop (ALIGN - 1) bytes from both
    // the start and the end of the vec when aligning
    vec![0; size + 2 * ALIGN]
}

fn aligned_mut(bytes: &mut [u8]) -> &mut [u8] {
    let head_off = bytes.as_ptr().align_offset(ALIGN);
    if head_off > bytes.len() {
        return &mut [];
    }
    let (_, aligned_head) = bytes.split_at_mut(head_off);

    let tail_off = (aligned_head.len() / ALIGN) * ALIGN;
    // SAFETY: can't be more than aligned_head.len()
    let (aligned, _) = aligned_head.split_at_mut(tail_off);

    aligned
}

pub(crate) fn shards_aligned_mut(bytes: &mut [u8], shard_size: usize) -> Result<Vec<&mut [u8]>> {
    if shard_size % ALIGN != 0 {
        return Err(LeopardError::InvalidShardSize(shard_size));
    }
    let aligned = aligned_mut(bytes);
    Ok(aligned.chunks_exact_mut(shard_size).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shard_size_must_be_multiple_of_align() {
        let mut mem = vec![0; 1024];
        for shard_size in [1, 10, 32, 63, 127, 257] {
            shards_aligned_mut(&mut mem, shard_size).unwrap_err();
        }

        for shard_size in (1..10).map(|x| x * 64) {
            shards_aligned_mut(&mut mem, shard_size).unwrap();
        }
    }

    #[test]
    fn check_alignment() {
        for shards_len in [2, 3, 7, 13, 21, 50, 77, 100] {
            for shard_size in [64, 128, 512, 1024] {
                let mut mem = alloc_zeroed_with_padding(shards_len * shard_size);
                let shards = shards_aligned_mut(&mut mem, shard_size).unwrap();

                for shard in shards {
                    assert_eq!(shard.as_ptr().align_offset(ALIGN), 0);
                    assert!(shard.iter().all(|&byte| byte == 0));
                }
            }
        }
    }
}
