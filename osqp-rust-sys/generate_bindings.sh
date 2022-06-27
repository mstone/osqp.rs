#!/bin/sh
git submodule init
git submodule update --recursive
srcdir=$(pwd)

mkdir build32
cd build32
cmake ../osqp -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCTRLC=OFF -DDFLOAT=OFF -DDLONG=OFF
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
find src -name '*.rs' -execdir sed -i -e 's/libc::/::std::os::raw::/g' -e '/::libc/d' -e 's/::std::os::raw::intptr_t/isize/g' {} +
sed -i -e '/libc/d' Cargo.toml
sed -i -e '/libc/d' lib.rs
git rm -rf src32
mv src src32

mkdir build64
cd build64
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
find src -name '*.rs' -execdir sed -i -e 's/libc::/::std::os::raw::/g' -e '/::libc/d' -e 's/::std::os::raw::intptr_t/isize/g' {} +
sed -i -e '/libc/d' Cargo.toml
sed -i -e '/libc/d' lib.rs
git rm -rf src64
mv src src64

awk '/pub mod src/ && !x {print "#[cfg_attr(target_pointer_width=\"32\", path=\"src32\")]\n#[cfg_attr(target_pointer_width=\"64\", path=\"src64\")]"; x=1} 1' lib.rs > lib.rs.tmp
mv lib.rs.tmp lib.rs

sed -i -e 's/src/src64/' Cargo.toml

rm -rf build32 build64

find . -name '*.rs' -execdir git add {} +
git add Cargo.toml
