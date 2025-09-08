use std::array;
use std::sync::LazyLock;

use crate::{add_mod, fwht, mul_log, BITS, MODULUS, ORDER};

include!(concat!(env!("OUT_DIR"), "/table.rs"));

// TODO: generate those in build.rs too
pub static FFT_SKEW: LazyLock<[u8; MODULUS as usize]> = LazyLock::new(init_fft_skew_table);
pub static LOG_WALSH: LazyLock<[u8; ORDER]> = LazyLock::new(init_log_walsh_table);

#[inline]
pub fn log(x: u8) -> u8 {
    // SAFETY: it has the size of ORDER, so u8 must be in bounds
    unsafe { *LOG.get_unchecked(x as usize) }
}

#[inline]
pub fn exp(x: u8) -> u8 {
    // SAFETY: it has the size of ORDER, so u8 must be in bounds
    unsafe { *EXP.get_unchecked(x as usize) }
}

#[inline]
pub fn mul(x: u8, y: u8) -> u8 {
    // SAFETY: it has the sizes of ORDER, so u8 must be in bounds
    unsafe { *MUL.get_unchecked(y as usize).get_unchecked(x as usize) }
}

/// # Panics
/// If x equals to [`MODULUS`]
#[allow(unused)]
#[inline]
pub fn fft_skew(x: u8) -> u8 {
    FFT_SKEW[x as usize]
}

#[allow(unused)]
#[inline]
pub fn log_walsh(x: u8) -> u8 {
    // SAFETY: it has the size of ORDER, so u8 must be in bounds
    unsafe { *LOG_WALSH.get_unchecked(x as usize) }
}

fn init_fft_skew_table() -> [u8; MODULUS as usize] {
    let mut temp: [u8; BITS - 1] = array::from_fn(|n| 1 << (n + 1));

    let mut fft_skew = [0u8; MODULUS as usize];

    for n in 0..BITS - 1 {
        let step = 1 << (n + 1);
        fft_skew[(1 << n) - 1] = 0;

        #[allow(clippy::needless_range_loop)]
        for i in n..BITS - 1 {
            let s = 1 << (i + 1);
            let start = (1 << n) - 1;

            for j in (start..s).step_by(step) {
                fft_skew[j + s] = fft_skew[j] ^ temp[i];
            }
        }

        temp[n] = MODULUS - log(mul_log(temp[n], log(temp[n] ^ 1)));

        for i in n + 1..BITS - 1 {
            let sum = add_mod(log(temp[i] ^ 1), temp[n]);
            temp[i] = mul_log(temp[i], sum);
        }
    }

    for skew in fft_skew.iter_mut() {
        *skew = log(*skew);
    }

    fft_skew
}

fn init_log_walsh_table() -> [u8; ORDER] {
    let mut log_walsh = [0u8; ORDER];

    for (log_walsh, &log) in log_walsh.iter_mut().zip(LOG.iter()) {
        *log_walsh = log;
    }

    fwht(&mut log_walsh, ORDER, ORDER);

    log_walsh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_lookup_table() {
        let expected = [
            255, 0, 85, 170, 17, 68, 34, 136, 153, 119, 102, 221, 51, 238, 187, 204, 219, 189, 7,
            112, 111, 246, 28, 193, 183, 123, 224, 14, 131, 56, 222, 237, 196, 23, 76, 113, 54,
            208, 99, 13, 19, 92, 49, 197, 67, 216, 52, 141, 137, 46, 226, 152, 26, 198, 108, 161,
            27, 104, 134, 177, 38, 184, 139, 98, 241, 210, 254, 24, 31, 45, 239, 129, 73, 236, 200,
            8, 148, 206, 140, 128, 199, 75, 251, 96, 180, 124, 6, 191, 35, 32, 37, 179, 2, 50, 59,
            82, 227, 165, 48, 253, 3, 223, 62, 90, 1, 25, 41, 157, 146, 217, 16, 145, 164, 118, 4,
            100, 70, 64, 103, 74, 143, 150, 192, 247, 127, 12, 105, 248, 160, 81, 29, 181, 167, 87,
            144, 135, 10, 21, 209, 91, 122, 117, 9, 120, 242, 86, 126, 110, 53, 147, 88, 175, 101,
            47, 230, 231, 57, 83, 250, 133, 130, 69, 116, 214, 93, 158, 30, 66, 71, 109, 40, 84,
            225, 36, 213, 233, 78, 212, 190, 97, 203, 89, 249, 185, 22, 235, 77, 228, 155, 159,
            149, 188, 65, 162, 107, 58, 15, 33, 79, 174, 240, 18, 244, 234, 20, 42, 182, 163, 11,
            245, 114, 166, 202, 94, 207, 205, 229, 172, 220, 252, 95, 176, 106, 39, 121, 43, 55,
            63, 44, 215, 201, 154, 156, 169, 194, 125, 115, 243, 151, 178, 5, 138, 173, 232, 132,
            60, 186, 61, 211, 171, 195, 72, 142, 218, 168, 80,
        ];
        assert_eq!(LOG.len(), expected.len());

        for (log, expected) in LOG.iter().zip(expected.iter()) {
            assert_eq!(log, expected);
        }
    }

    #[test]
    fn exp_lookup_table() {
        let expected = [
            1, 104, 92, 100, 114, 240, 86, 18, 75, 142, 136, 208, 125, 39, 27, 196, 110, 4, 201,
            40, 204, 137, 184, 33, 67, 105, 52, 56, 22, 130, 166, 68, 89, 197, 6, 88, 173, 90, 60,
            223, 170, 106, 205, 225, 228, 69, 49, 153, 98, 42, 93, 12, 46, 148, 36, 226, 29, 156,
            195, 94, 245, 247, 102, 227, 117, 192, 167, 44, 5, 161, 116, 168, 251, 72, 119, 81, 34,
            186, 176, 198, 255, 129, 95, 157, 171, 2, 145, 133, 150, 181, 103, 139, 41, 164, 213,
            220, 83, 179, 63, 38, 115, 152, 10, 118, 57, 126, 222, 194, 54, 169, 147, 20, 19, 35,
            210, 236, 162, 141, 113, 9, 143, 224, 140, 25, 85, 235, 146, 124, 79, 71, 160, 28, 244,
            159, 58, 135, 7, 48, 241, 62, 78, 47, 252, 120, 134, 111, 108, 149, 76, 190, 121, 238,
            51, 8, 231, 188, 232, 107, 165, 189, 128, 55, 193, 207, 112, 97, 211, 132, 254, 233, 3,
            249, 217, 242, 199, 151, 221, 59, 239, 91, 84, 131, 206, 24, 61, 183, 246, 14, 191, 17,
            178, 87, 122, 23, 234, 250, 32, 43, 53, 80, 74, 230, 212, 180, 15, 215, 77, 214, 37,
            138, 65, 248, 177, 174, 163, 229, 45, 109, 253, 16, 218, 11, 30, 101, 26, 172, 50, 96,
            187, 216, 154, 155, 243, 175, 203, 185, 73, 31, 13, 70, 200, 64, 144, 237, 202, 209,
            21, 123, 127, 182, 158, 82, 219, 99, 66, 1,
        ];
        assert_eq!(EXP.len(), expected.len());

        for (exp, expected) in EXP.iter().zip(expected.iter()) {
            assert_eq!(exp, expected);
        }
    }

    #[test]
    fn fft_skew_lookup_table() {
        let expected = [
            255, 255, 85, 255, 17, 85, 34, 255, 153, 17, 102, 85, 51, 34, 187, 255, 219, 153, 7,
            17, 111, 102, 28, 85, 183, 51, 224, 34, 131, 187, 222, 255, 196, 219, 76, 153, 54, 7,
            99, 17, 19, 111, 49, 102, 67, 28, 52, 85, 137, 183, 226, 51, 26, 224, 108, 34, 27, 131,
            134, 187, 38, 222, 139, 255, 241, 196, 254, 219, 31, 76, 239, 153, 73, 54, 200, 7, 148,
            99, 140, 17, 199, 19, 251, 111, 180, 49, 6, 102, 35, 67, 37, 28, 2, 52, 59, 85, 227,
            137, 48, 183, 3, 226, 62, 51, 1, 26, 41, 224, 146, 108, 16, 34, 164, 27, 4, 131, 70,
            134, 103, 187, 143, 38, 192, 222, 127, 139, 105, 255, 160, 241, 29, 196, 167, 254, 144,
            219, 10, 31, 209, 76, 122, 239, 9, 153, 242, 73, 126, 54, 53, 200, 88, 7, 101, 148,
            230, 99, 57, 140, 250, 17, 130, 199, 116, 19, 93, 251, 30, 111, 71, 180, 40, 49, 225,
            6, 213, 102, 78, 35, 190, 67, 203, 37, 249, 28, 22, 2, 77, 52, 155, 59, 149, 85, 65,
            227, 107, 137, 15, 48, 79, 183, 240, 3, 244, 226, 20, 62, 182, 51, 11, 1, 114, 26, 202,
            41, 207, 224, 229, 146, 220, 108, 95, 16, 106, 34, 121, 164, 55, 27, 44, 4, 201, 131,
            156, 70, 194, 134, 115, 103, 151, 187, 5, 143, 173, 38, 132, 192, 186, 222, 211, 127,
            195, 139, 142, 105, 168,
        ];

        assert_eq!(FFT_SKEW.len(), expected.len());

        for (skew, expected) in FFT_SKEW.iter().zip(expected.iter()) {
            assert_eq!(skew, expected);
        }
    }

    #[test]
    fn walsh_lookup_table() {
        let expected = [
            255, 252, 236, 126, 27, 63, 118, 40, 193, 20, 59, 87, 141, 132, 62, 159, 238, 207, 31,
            133, 198, 151, 127, 66, 224, 162, 96, 10, 76, 171, 157, 78, 17, 39, 206, 88, 38, 154,
            64, 213, 112, 240, 37, 81, 197, 5, 48, 235, 119, 114, 241, 231, 46, 194, 143, 214, 222,
            33, 191, 122, 99, 218, 145, 203, 170, 229, 200, 86, 177, 121, 4, 109, 111, 163, 50,
            144, 67, 61, 223, 149, 187, 94, 129, 57, 108, 243, 248, 93, 123, 107, 199, 91, 23, 195,
            1, 97, 136, 95, 140, 147, 98, 44, 103, 36, 219, 234, 32, 29, 19, 232, 165, 77, 102,
            245, 24, 101, 226, 169, 247, 130, 56, 79, 180, 120, 208, 168, 146, 185, 0, 220, 73,
            110, 104, 55, 164, 84, 28, 42, 82, 167, 52, 60, 90, 155, 51, 205, 45, 250, 26, 178, 12,
            30, 14, 65, 251, 21, 113, 211, 41, 212, 68, 106, 148, 175, 184, 201, 70, 233, 7, 18,
            179, 160, 49, 138, 253, 22, 153, 166, 210, 230, 137, 125, 150, 116, 237, 15, 6, 117,
            13, 142, 16, 89, 85, 172, 8, 242, 134, 43, 100, 71, 246, 182, 2, 135, 216, 186, 3, 188,
            204, 202, 239, 83, 161, 115, 105, 158, 183, 58, 75, 209, 196, 72, 25, 190, 34, 176,
            128, 53, 139, 215, 74, 225, 189, 244, 35, 181, 92, 173, 227, 228, 221, 11, 254, 47,
            152, 156, 192, 69, 131, 174, 124, 9, 54, 80, 217, 249,
        ];
        assert_eq!(LOG_WALSH.len(), expected.len());

        for (walsh, expected) in LOG_WALSH.iter().zip(expected.iter()) {
            assert_eq!(walsh, expected);
        }
    }

    #[test]
    fn mul_lookup_table() {
        let expected: [[u8; ORDER]; ORDER] = include!("../test_data/expected_mul_lut.rs");

        assert_eq!(MUL.len(), expected.len());

        for (mul_row, expected_row) in MUL.iter().zip(expected.iter()) {
            assert_eq!(mul_row.len(), expected.len());
            for (mul, expected) in mul_row.iter().zip(expected_row.iter()) {
                assert_eq!(mul, expected);
            }
        }
    }
}
