//! build go leopard library with bindings

use std::{env, path::PathBuf, process::Command};

use bindgen::Builder;

fn main() {
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    // build go as a static lib and a header with cgo
    Command::new("go")
        .arg("build")
        .arg("-buildmode=c-archive")
        .arg("-o")
        .arg(out_path.join("libgoleopard.a"))
        .arg("leopard.go")
        .current_dir("go")
        .status()
        .expect("Building libgoleopard.a failed");

    // generate bindings
    Builder::default()
        .header(out_path.join("libgoleopard.h").to_str().unwrap())
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");

    println!("cargo:rerun-if-changed=go");
    println!(
        "cargo:rustc-link-search=native={}",
        out_path.to_str().unwrap()
    );
    println!("cargo:rustc-link-lib=static=goleopard");
}
