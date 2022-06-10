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
pub const AMD_DEFAULT_AGGRESSIVE: libc::c_int = 1 as libc::c_int;
pub const AMD_DEFAULT_DENSE: libc::c_double = 10.0f64;
pub const AMD_AGGRESSIVE: libc::c_int = 1 as libc::c_int;
pub const AMD_DENSE: libc::c_int = 0 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_control(mut Control: *mut c_float) {
    let mut alpha: c_float = 0.;
    let mut aggressive: libc::c_longlong = 0;
    if !Control.is_null() {
        alpha = *Control.offset(AMD_DENSE as isize);
        aggressive = (*Control.offset(AMD_AGGRESSIVE as isize)
            != 0 as libc::c_int as libc::c_double) as libc::c_int as libc::c_longlong;
    } else {
        alpha = AMD_DEFAULT_DENSE;
        aggressive = AMD_DEFAULT_AGGRESSIVE as libc::c_longlong;
    }
    if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(
            b"\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n    dense row parameter: %g\n\0"
                as *const u8 as *const libc::c_char,
            2 as libc::c_int,
            4 as libc::c_int,
            6 as libc::c_int,
            b"May 4, 2016\0" as *const u8 as *const libc::c_char,
            alpha,
        );
    }
    if alpha < 0 as libc::c_int as libc::c_double {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    no rows treated as dense\n\0" as *const u8 as *const libc::c_char,
            );
        }
    } else if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(
            b"    (rows with more than max (%g * sqrt (n), 16) entries are\n    considered \"dense\", and placed last in output permutation)\n\0"
                as *const u8 as *const libc::c_char,
            alpha,
        );
    }
    if aggressive != 0 {
        if (SuiteSparse_config.printf_func).is_some() {
            (SuiteSparse_config.printf_func)
                .expect(
                    "non-null function pointer",
                )(
                b"    aggressive absorption:  yes\n\0" as *const u8
                    as *const libc::c_char,
            );
        }
    } else if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(b"    aggressive absorption:  no\n\0" as *const u8 as *const libc::c_char);
    }
    if (SuiteSparse_config.printf_func).is_some() {
        (SuiteSparse_config.printf_func)
            .expect(
                "non-null function pointer",
            )(
            b"    size of AMD integer: %d\n\n\0" as *const u8 as *const libc::c_char,
            ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
        );
    }
}
