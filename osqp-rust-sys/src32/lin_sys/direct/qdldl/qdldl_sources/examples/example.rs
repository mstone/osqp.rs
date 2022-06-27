#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![register_tool(c2rust)]
#![feature(register_tool)]
use ::osqp_rust_sys::*;
extern "C" {
    fn QDLDL_etree(
        n: QDLDL_int,
        Ap_0: *const QDLDL_int,
        Ai_0: *const QDLDL_int,
        work: *mut QDLDL_int,
        Lnz: *mut QDLDL_int,
        etree: *mut QDLDL_int,
    ) -> QDLDL_int;
    fn QDLDL_factor(
        n: QDLDL_int,
        Ap_0: *const QDLDL_int,
        Ai_0: *const QDLDL_int,
        Ax_0: *const QDLDL_float,
        Lp: *mut QDLDL_int,
        Li: *mut QDLDL_int,
        Lx: *mut QDLDL_float,
        D: *mut QDLDL_float,
        Dinv: *mut QDLDL_float,
        Lnz: *const QDLDL_int,
        etree: *const QDLDL_int,
        bwork: *mut QDLDL_bool,
        iwork: *mut QDLDL_int,
        fwork: *mut QDLDL_float,
    ) -> QDLDL_int;
    fn QDLDL_solve(
        n: QDLDL_int,
        Lp: *const QDLDL_int,
        Li: *const QDLDL_int,
        Lx: *const QDLDL_float,
        Dinv: *const QDLDL_float,
        x: *mut QDLDL_float,
    );
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn free(_: *mut ::std::os::raw::c_void);
    fn printf(_: *const ::std::os::raw::c_char, _: ...) -> ::std::os::raw::c_int;
}
pub type QDLDL_int = ::std::os::raw::c_int;
pub type QDLDL_float = ::std::os::raw::c_double;
pub type QDLDL_bool = ::std::os::raw::c_uchar;
#[no_mangle]
pub static mut An: QDLDL_int = 10 as ::std::os::raw::c_int;
#[no_mangle]
pub static mut Ap: [QDLDL_int; 11] = [
    0 as ::std::os::raw::c_int,
    1 as ::std::os::raw::c_int,
    2 as ::std::os::raw::c_int,
    4 as ::std::os::raw::c_int,
    5 as ::std::os::raw::c_int,
    6 as ::std::os::raw::c_int,
    8 as ::std::os::raw::c_int,
    10 as ::std::os::raw::c_int,
    12 as ::std::os::raw::c_int,
    14 as ::std::os::raw::c_int,
    17 as ::std::os::raw::c_int,
];
#[no_mangle]
pub static mut Ai: [QDLDL_int; 17] = [
    0 as ::std::os::raw::c_int,
    1 as ::std::os::raw::c_int,
    1 as ::std::os::raw::c_int,
    2 as ::std::os::raw::c_int,
    3 as ::std::os::raw::c_int,
    4 as ::std::os::raw::c_int,
    1 as ::std::os::raw::c_int,
    5 as ::std::os::raw::c_int,
    0 as ::std::os::raw::c_int,
    6 as ::std::os::raw::c_int,
    3 as ::std::os::raw::c_int,
    7 as ::std::os::raw::c_int,
    6 as ::std::os::raw::c_int,
    8 as ::std::os::raw::c_int,
    1 as ::std::os::raw::c_int,
    2 as ::std::os::raw::c_int,
    9 as ::std::os::raw::c_int,
];
#[no_mangle]
pub static mut Ax: [QDLDL_float; 17] = [
    1.0f64,
    0.460641f64,
    -0.121189f64,
    0.417928f64,
    0.177828f64,
    0.1f64,
    -0.0290058f64,
    -1.0f64,
    0.350321f64,
    -0.441092f64,
    -0.0845395f64,
    -0.316228f64,
    0.178663f64,
    -0.299077f64,
    0.182452f64,
    -1.56506f64,
    -0.1f64,
];
#[no_mangle]
pub static mut b: [QDLDL_float; 10] = [
    1 as ::std::os::raw::c_int as QDLDL_float,
    2 as ::std::os::raw::c_int as QDLDL_float,
    3 as ::std::os::raw::c_int as QDLDL_float,
    4 as ::std::os::raw::c_int as QDLDL_float,
    5 as ::std::os::raw::c_int as QDLDL_float,
    6 as ::std::os::raw::c_int as QDLDL_float,
    7 as ::std::os::raw::c_int as QDLDL_float,
    8 as ::std::os::raw::c_int as QDLDL_float,
    9 as ::std::os::raw::c_int as QDLDL_float,
    10 as ::std::os::raw::c_int as QDLDL_float,
];
unsafe fn main_0() -> ::std::os::raw::c_int {
    let mut i: QDLDL_int = 0;
    let mut Ln: QDLDL_int = An;
    let mut Lp: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut Li: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut Lx: *mut QDLDL_float = 0 as *mut QDLDL_float;
    let mut D: *mut QDLDL_float = 0 as *mut QDLDL_float;
    let mut Dinv: *mut QDLDL_float = 0 as *mut QDLDL_float;
    let mut etree: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut Lnz: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut sumLnz: QDLDL_int = 0;
    let mut iwork: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut bwork: *mut QDLDL_bool = 0 as *mut QDLDL_bool;
    let mut fwork: *mut QDLDL_float = 0 as *mut QDLDL_float;
    let mut x: *mut QDLDL_float = 0 as *mut QDLDL_float;
    etree = malloc(
        (::std::mem::size_of::<QDLDL_int>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_int;
    Lnz = malloc(
        (::std::mem::size_of::<QDLDL_int>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_int;
    Lp = malloc(
        (::std::mem::size_of::<QDLDL_int>() as ::std::os::raw::c_ulong)
            .wrapping_mul((An + 1 as ::std::os::raw::c_int) as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_int;
    D = malloc(
        (::std::mem::size_of::<QDLDL_float>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_float;
    Dinv = malloc(
        (::std::mem::size_of::<QDLDL_float>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_float;
    iwork = malloc(
        (::std::mem::size_of::<QDLDL_int>() as ::std::os::raw::c_ulong)
            .wrapping_mul((3 as ::std::os::raw::c_int * An) as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_int;
    bwork = malloc(
        (::std::mem::size_of::<QDLDL_bool>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_bool;
    fwork = malloc(
        (::std::mem::size_of::<QDLDL_float>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_float;
    sumLnz = QDLDL_etree(An, Ap.as_ptr(), Ai.as_ptr(), iwork, Lnz, etree);
    Li = malloc(
        (::std::mem::size_of::<QDLDL_int>() as ::std::os::raw::c_ulong)
            .wrapping_mul(sumLnz as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_int;
    Lx = malloc(
        (::std::mem::size_of::<QDLDL_float>() as ::std::os::raw::c_ulong)
            .wrapping_mul(sumLnz as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_float;
    QDLDL_factor(
        An,
        Ap.as_ptr(),
        Ai.as_ptr(),
        Ax.as_ptr(),
        Lp,
        Li,
        Lx,
        D,
        Dinv,
        Lnz,
        etree,
        bwork,
        iwork,
        fwork,
    );
    x = malloc(
        (::std::mem::size_of::<QDLDL_float>() as ::std::os::raw::c_ulong)
            .wrapping_mul(An as ::std::os::raw::c_ulong),
    ) as *mut QDLDL_float;
    i = 0 as ::std::os::raw::c_int;
    while i < Ln {
        *x.offset(i as isize) = b[i as usize];
        i += 1;
    }
    QDLDL_solve(Ln, Lp, Li, Lx, Dinv, x);
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(b"A (CSC format):\n\0" as *const u8 as *const ::std::os::raw::c_char);
    print_line();
    print_arrayi(
        Ap.as_ptr(),
        An + 1 as ::std::os::raw::c_int,
        b"A.p\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayi(
        Ai.as_ptr(),
        Ap[An as usize],
        b"A.i\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayf(
        Ax.as_ptr(),
        Ap[An as usize],
        b"A.x\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    printf(b"\n\n\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(b"elimination tree:\n\0" as *const u8 as *const ::std::os::raw::c_char);
    print_line();
    print_arrayi(
        etree,
        Ln,
        b"etree\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayi(
        Lnz,
        Ln,
        b"Lnz\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    printf(b"\n\n\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(b"L (CSC format):\n\0" as *const u8 as *const ::std::os::raw::c_char);
    print_line();
    print_arrayi(
        Lp,
        Ln + 1 as ::std::os::raw::c_int,
        b"L.p\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayi(
        Li,
        *Lp.offset(Ln as isize),
        b"L.i\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayf(
        Lx,
        *Lp.offset(Ln as isize),
        b"L.x\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    printf(b"\n\n\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(b"D:\n\0" as *const u8 as *const ::std::os::raw::c_char);
    print_line();
    print_arrayf(
        D,
        An,
        b"diag(D)     \0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayf(
        Dinv,
        An,
        b"diag(D^{-1})\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    printf(b"\n\n\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(b"solve results:\n\0" as *const u8 as *const ::std::os::raw::c_char);
    print_line();
    print_arrayf(
        b.as_ptr(),
        An,
        b"b\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    print_arrayf(
        x,
        An,
        b"A\\b\0" as *const u8 as *const ::std::os::raw::c_char as *mut ::std::os::raw::c_char,
    );
    printf(b"\n\n\0" as *const u8 as *const ::std::os::raw::c_char);
    free(Lp as *mut ::std::os::raw::c_void);
    free(Li as *mut ::std::os::raw::c_void);
    free(Lx as *mut ::std::os::raw::c_void);
    free(D as *mut ::std::os::raw::c_void);
    free(Dinv as *mut ::std::os::raw::c_void);
    free(etree as *mut ::std::os::raw::c_void);
    free(Lnz as *mut ::std::os::raw::c_void);
    free(iwork as *mut ::std::os::raw::c_void);
    free(bwork as *mut ::std::os::raw::c_void);
    free(fwork as *mut ::std::os::raw::c_void);
    free(x as *mut ::std::os::raw::c_void);
    return 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn print_line() {
    printf(b"--------------------------\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn print_arrayi(
    mut data: *const QDLDL_int,
    mut n: QDLDL_int,
    mut varName: *mut ::std::os::raw::c_char,
) {
    let mut i: QDLDL_int = 0;
    printf(b"%s = [\0" as *const u8 as *const ::std::os::raw::c_char, varName);
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        printf(b"%i,\0" as *const u8 as *const ::std::os::raw::c_char, *data.offset(i as isize));
        i += 1;
    }
    printf(b"]\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn print_arrayf(
    mut data: *const QDLDL_float,
    mut n: QDLDL_int,
    mut varName: *mut ::std::os::raw::c_char,
) {
    let mut i: QDLDL_int = 0;
    printf(b"%s = [\0" as *const u8 as *const ::std::os::raw::c_char, varName);
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        printf(b"%.3g,\0" as *const u8 as *const ::std::os::raw::c_char, *data.offset(i as isize));
        i += 1;
    }
    printf(b"]\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
pub fn main() {
    unsafe { ::std::process::exit(main_0() as i32) }
}
