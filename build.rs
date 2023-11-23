use std::{env, fs::File, path::Path};

// Number of bits per element
const K_BITS: usize = 8;
// Finite field order: Number of elements in the field
const K_ORDER: usize = 256;
// Modulus for field operations
const K_MODULUS: u8 = 255;
// LFSR Polynomial that generates the field elements
const K_POLYNOMIAL: usize = 0x11D;
// Basis used for generating logarithm tables
const K_CANTOR_BASIS: [u8; K_BITS] = [1, 214, 152, 146, 86, 200, 88, 230];

macro_rules! write_table {
    ($dest:expr, $table:ident) => {{
        let typ = type_name_of_val(&$table);
        let name = stringify!($table).to_ascii_uppercase();
        let table = format!("pub static {name}: {typ} = {:?};\n", $table);
        ::std::io::Write::write_all($dest, table.as_bytes()).expect(&format!(
            "Failed to write lookup table: {}",
            stringify!($table)
        ));
    }};
}

fn main() {
    let (log, exp) = create_lut();

    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("table.rs");
    let mut outfile = File::create(dest_path).expect("Failed to create file for lookup tables");

    write_table!(&mut outfile, log);
    write_table!(&mut outfile, exp);
}

fn create_lut() -> ([u8; K_ORDER], [u8; K_ORDER]) {
    let mut log_table = [0u8; K_ORDER];
    let mut exp_table = [0u8; K_ORDER];
    let mut state: usize = 1;

    // Generate field group with LFSR
    for i in 0..K_MODULUS as usize {
        exp_table[state] = i.try_into().expect("must be valid u8");
        state <<= 1;
        if state >= K_ORDER {
            state ^= K_POLYNOMIAL;
        }
    }
    exp_table[0] = K_MODULUS;

    // Convert to Cantor basis
    log_table[0] = 0;
    for (i, basis) in K_CANTOR_BASIS.iter().enumerate() {
        let width = 1 << i;
        for j in 0..width {
            log_table[j + width] = log_table[j] ^ basis;
        }
    }
    log_table = std::array::from_fn(|i| exp_table[log_table[i] as usize]);

    // generate exponent table from logarithm table
    for i in 0..K_ORDER {
        exp_table[log_table[i] as usize] = i.try_into().expect("must be valid u8");
    }

    // Note: Handles modulus wrap around with LUT
    exp_table[K_MODULUS as usize] = exp_table[0];

    (log_table, exp_table)
}

fn type_name_of_val<T: ?Sized>(_val: &T) -> &'static str {
    std::any::type_name::<T>()
}
