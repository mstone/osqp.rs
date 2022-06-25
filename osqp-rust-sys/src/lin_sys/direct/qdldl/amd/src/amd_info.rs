use ::libc;
extern "C" {
    static mut SuiteSparse_config: SuiteSparse_config_struct;
}
pub type __darwin_size_t = libc::c_ulong;
pub type size_t = __darwin_size_t;
pub type c_float = libc::c_double;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct SuiteSparse_config_struct {
    pub malloc_func: Option::<unsafe extern "C" fn(size_t) -> *mut libc::c_void>,
    pub realloc_func: Option::<
        unsafe extern "C" fn(*mut libc::c_void, size_t) -> *mut libc::c_void,
    >,
    pub free_func: Option::<unsafe extern "C" fn(*mut libc::c_void) -> ()>,
    pub printf_func: Option::<
        unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int,
    >,
    pub hypot_func: Option::<unsafe extern "C" fn(c_float, c_float) -> c_float>,
    pub divcomplex_func: Option::<
        unsafe extern "C" fn(
            c_float,
            c_float,
            c_float,
            c_float,
            *mut c_float,
            *mut c_float,
        ) -> libc::c_int,
    >,
}
pub const NULL: libc::c_int = 0 as libc::c_int;
pub const AMD_OK_BUT_JUMBLED: libc::c_int = 1 as libc::c_int;
pub const AMD_STATUS: libc::c_int = 0 as libc::c_int;
pub const AMD_INVALID: libc::c_int = -(2 as libc::c_int);
pub const AMD_OUT_OF_MEMORY: libc::c_int = -(1 as libc::c_int);
pub const AMD_OK: libc::c_int = 0 as libc::c_int;
pub const AMD_LNZ: libc::c_int = 9 as libc::c_int;
pub const AMD_NMULTSUBS_LU: libc::c_int = 12 as libc::c_int;
pub const AMD_NMULTSUBS_LDL: libc::c_int = 11 as libc::c_int;
pub const AMD_NDIV: libc::c_int = 10 as libc::c_int;
pub const AMD_N: libc::c_int = 1 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_info(mut Info: *mut c_float) {
    let mut n: c_float = 0.;
    let mut ndiv: c_float = 0.;
    let mut nmultsubs_ldl: c_float = 0.;
    let mut nmultsubs_lu: c_float = 0.;
    let mut lnz: c_float = 0.;
    let mut lnzd: c_float = 0.;
    if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(
            b"\nAMD version %d.%d.%d, %s, results:\n\0" as *const u8
                as *const libc::c_char,
            2 as libc::c_int,
            4 as libc::c_int,
            6 as libc::c_int,
            b"May 4, 2016\0" as *const u8 as *const libc::c_char,
        );
    }
    if Info.is_null() {
        return;
    }
    n = *Info.offset(AMD_N as isize);
    ndiv = *Info.offset(AMD_NDIV as isize);
    nmultsubs_ldl = *Info.offset(AMD_NMULTSUBS_LDL as isize);
    nmultsubs_lu = *Info.offset(AMD_NMULTSUBS_LU as isize);
    lnz = *Info.offset(AMD_LNZ as isize);
    lnzd = if n >= 0 as libc::c_int as libc::c_double
        && lnz >= 0 as libc::c_int as libc::c_double
    {
        n + lnz
    } else {
        -(1 as libc::c_int) as libc::c_double
    };
    if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(b"    status: \0" as *const u8 as *const libc::c_char);
    }
    if *Info.offset(AMD_STATUS as isize) == AMD_OK as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(b"OK\n\0" as *const u8 as *const libc::c_char);
        }
    } else if *Info.offset(AMD_STATUS as isize) == AMD_OUT_OF_MEMORY as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(b"out of memory\n\0" as *const u8 as *const libc::c_char);
        }
    } else if *Info.offset(AMD_STATUS as isize) == AMD_INVALID as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(b"invalid matrix\n\0" as *const u8 as *const libc::c_char);
        }
    } else if *Info.offset(AMD_STATUS as isize) == AMD_OK_BUT_JUMBLED as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(b"OK, but jumbled\n\0" as *const u8 as *const libc::c_char);
        }
    } else if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(b"unknown\n\0" as *const u8 as *const libc::c_char);
    }
    if n >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    n, dimension of A:                                  %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                n,
            );
        }
    }
    if *Info.offset(2 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    nz, number of nonzeros in A:                        %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(2 as libc::c_int as isize),
            );
        }
    }
    if *Info.offset(3 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    symmetry of A:                                      %.4f\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(3 as libc::c_int as isize),
            );
        }
    }
    if *Info.offset(4 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    number of nonzeros on diagonal:                     %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(4 as libc::c_int as isize),
            );
        }
    }
    if *Info.offset(5 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    nonzeros in pattern of A+A' (excl. diagonal):       %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(5 as libc::c_int as isize),
            );
        }
    }
    if *Info.offset(6 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    # dense rows/columns of A+A':                       %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(6 as libc::c_int as isize),
            );
        }
    }
    if *Info.offset(7 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    memory used, in bytes:                              %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(7 as libc::c_int as isize),
            );
        }
    }
    if *Info.offset(8 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    # of memory compactions:                            %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(8 as libc::c_int as isize),
            );
        }
    }
    if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(
            b"\n    The following approximate statistics are for a subsequent\n    factorization of A(P,P) + A(P,P)'.  They are slight upper\n    bounds if there are no dense rows/columns in A+A', and become\n    looser if dense rows/columns exist.\n\n\0"
                as *const u8 as *const libc::c_char,
        );
    }
    if lnz >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    nonzeros in L (excluding diagonal):                 %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                lnz,
            );
        }
    }
    if lnzd >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    nonzeros in L (including diagonal):                 %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                lnzd,
            );
        }
    }
    if ndiv >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    # divide operations for LDL' or LU:                 %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                ndiv,
            );
        }
    }
    if nmultsubs_ldl >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    # multiply-subtract operations for LDL':            %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                nmultsubs_ldl,
            );
        }
    }
    if nmultsubs_lu >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    # multiply-subtract operations for LU:              %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                nmultsubs_lu,
            );
        }
    }
    if *Info.offset(13 as libc::c_int as isize) >= 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    max nz. in any column of L (incl. diagonal):        %.20g\n\0"
                    as *const u8 as *const libc::c_char,
                *Info.offset(13 as libc::c_int as isize),
            );
        }
    }
    if n >= 0 as libc::c_int as libc::c_double
        && ndiv >= 0 as libc::c_int as libc::c_double
        && nmultsubs_ldl >= 0 as libc::c_int as libc::c_double
        && nmultsubs_lu >= 0 as libc::c_int as libc::c_double
    {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"\n    chol flop count for real A, sqrt counted as 1 flop: %.20g\n    LDL' flop count for real A:                         %.20g\n    LDL' flop count for complex A:                      %.20g\n    LU flop count for real A (with no pivoting):        %.20g\n    LU flop count for complex A (with no pivoting):     %.20g\n\n\0"
                    as *const u8 as *const libc::c_char,
                n + ndiv + 2 as libc::c_int as libc::c_double * nmultsubs_ldl,
                ndiv + 2 as libc::c_int as libc::c_double * nmultsubs_ldl,
                9 as libc::c_int as libc::c_double * ndiv
                    + 8 as libc::c_int as libc::c_double * nmultsubs_ldl,
                ndiv + 2 as libc::c_int as libc::c_double * nmultsubs_lu,
                9 as libc::c_int as libc::c_double * ndiv
                    + 8 as libc::c_int as libc::c_double * nmultsubs_lu,
            );
        }
    }
}
