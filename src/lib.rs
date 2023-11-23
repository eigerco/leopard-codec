// Number of bits per element
const K_BITS: usize = 8;
// Finite field order: Number of elements in the field
const K_ORDER: usize = u8::MAX as usize + 1;
// Modulus for field operations
const K_MODULUS: u8 = u8::MAX;
// LFSR Polynomial that generates the field elements
const K_POLYNOMIAL: usize = 0x11D;
// Basis used for generating logarithm tables
const K_CANTOR_BASIS: [u8; K_BITS] = [1, 214, 152, 146, 86, 200, 88, 230];

/// lookup tables
mod lut {
    include!(concat!(env!("OUT_DIR"), "/table.rs"));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_lookup_table() {
        let mut log_lut = lut::LOG.to_vec();

        // this array should have all distinct fields from the group
        assert_eq!(log_lut.len(), K_ORDER);

        log_lut.sort_unstable();
        log_lut.dedup();

        assert_eq!(log_lut.len(), K_ORDER);
    }

    #[test]
    fn exp_lookup_table() {
        let mut exp_lut = lut::EXP.to_vec();

        // this array should have all distinct fields from the group
        assert_eq!(exp_lut.len(), K_ORDER);

        exp_lut.sort_unstable();

        // except for 0, as 1 is an exponent of 0 and 1
        assert_eq!(&exp_lut[..2], &[1, 1]);

        exp_lut.dedup();

        assert_eq!(exp_lut.len(), K_ORDER - 1);
    }
}
