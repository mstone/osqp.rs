#![allow(dead_code)]
#![allow(mutable_transmutes)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_assignments)]
#![allow(unused_mut)]
#![feature(extern_types)]
#![feature(linkage)]
#![feature(register_tool)]
#![register_tool(c2rust)]


extern crate libc;
pub mod src {
pub mod lin_sys {
pub mod direct {
pub mod pardiso {
pub mod pardiso_interface;
pub mod pardiso_loader;
} // mod pardiso
pub mod qdldl {
pub mod amd {
pub mod src {
pub mod SuiteSparse_config;
pub mod amd_1;
pub mod amd_2;
pub mod amd_aat;
pub mod amd_control;
pub mod amd_defaults;
pub mod amd_info;
pub mod amd_order;
pub mod amd_post_tree;
pub mod amd_postorder;
pub mod amd_preprocess;
pub mod amd_valid;
} // mod src
} // mod amd
pub mod qdldl_interface;
pub mod qdldl_sources {
pub mod src {
pub mod qdldl;
} // mod src
} // mod qdldl_sources
} // mod qdldl
} // mod direct
pub mod lib_handler;
} // mod lin_sys
pub mod src {
pub mod auxil;
pub mod cs;
pub mod error;
pub mod kkt;
pub mod lin_alg;
pub mod lin_sys;
pub mod osqp;
pub mod polish;
pub mod proj;
pub mod scaling;
pub mod util;
} // mod src
} // mod src
