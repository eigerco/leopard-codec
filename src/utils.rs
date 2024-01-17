use std::ops::Index;

const ALIGN: usize = 64;

#[repr(C, align(64))]
#[derive(Clone)]
struct Aligned([u8; ALIGN]);

impl Aligned {
    fn zeroed() -> Self {
        Self([0; ALIGN])
    }
}

pub struct AlignedShards {
    shards: usize,
    shard_size: usize,
    aligned_shard_size: usize,
    memory: Vec<u8>,
}

impl AlignedShards {
    pub fn new(shards: usize, shard_size: usize) -> Self {
        // round up to the closest multiply of align
        let aligned_shard_size = (shard_size + ALIGN - 1) / ALIGN * ALIGN;
        let total = shards * aligned_shard_size;

        // create a vec of Aligned to make sure they are aligned properly
        let aligned_mem = vec![Aligned::zeroed(); total / ALIGN];

        // grab the raw parts
        let length = aligned_mem.len() * ALIGN;
        let capacity = aligned_mem.capacity() * ALIGN;
        let mem_ptr = aligned_mem.leak().as_mut_ptr() as *mut u8;

        // recreate the vec over bytes
        let memory = unsafe { Vec::from_raw_parts(mem_ptr, length, capacity) };

        Self {
            shards,
            shard_size,
            aligned_shard_size,
            memory,
        }
    }
}

impl Index<usize> for AlignedShards {
    type Output = [u8];

    fn index(&self, index: usize) -> &Self::Output {
        let start = self.aligned_shard_size * index;
        let end = start + self.shard_size;
        &self.memory[start..end]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn size_of_aligned() {
        assert_eq!(std::mem::size_of::<Aligned>(), ALIGN);
    }

    #[test]
    fn is_properly_aligned() {
        let mem = AlignedShards::new(50, 500);

        println!("{}", std::mem::size_of::<Aligned>());
        println!("{}", std::mem::align_of::<Aligned>());

        for n in 0..50 {
            assert_eq!(mem[n].as_ptr().align_offset(ALIGN), 0);
        }
    }
}
