#!/bin/sh
git submodule init
git submodule update --recursive
srcdir=$(pwd)
mkdir build
cd build
cmake ../osqp -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCTRLC=OFF -DDFLOAT=OFF -DDLONG=ON
c2rust transpile --overwrite-existing --fail-on-error --translate-const-macros -b osqp_demo -b example -o osqp  compile_commands.json -- -D_POSIX_C_SOURCE
pwd
rsync -Pai osqp/. ../.
cd "$srcdir"
sed -i -e '1,3d' Cargo.toml
sed -i -e 's/"osqp"/"osqp-rust-sys"/g' Cargo.toml
sed -i -e '3,$s/"osqp-rust-sys"/"osqp_rust_sys"/g' Cargo.toml
sed -i -e 's/0.0.0/0.6.2/g' Cargo.toml
sed -i -e 's/C2Rust/Michael Stone <michael.r.stone@gmail.com>/g' Cargo.toml
sed -i -e 's/use ::osqp::\*;/use ::osqp_rust_sys::*;/' src/examples/osqp_demo.rs
sed -i -e 's/use ::osqp::\*;/use ::osqp_rust_sys::*;/' src/lin_sys/direct/qdldl/qdldl_sources/examples/example.rs
(cd build/osqp; for f in $(find . -type f -print); do echo git add "$srcdir/$f"; git add "$srcdir/$f"; done)
rm -rf build
