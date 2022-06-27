extern "C" {
    fn fabs(_: ::std::os::raw::c_double) -> ::std::os::raw::c_double;
    fn sqrt(_: ::std::os::raw::c_double) -> ::std::os::raw::c_double;
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn realloc(_: *mut ::std::os::raw::c_void, _: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn free(_: *mut ::std::os::raw::c_void);
}
pub type __darwin_size_t = ::std::os::raw::c_ulong;
pub type size_t = __darwin_size_t;
pub type c_float = ::std::os::raw::c_double;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct SuiteSparse_config_struct {
    pub malloc_func: Option::<unsafe extern "C" fn(size_t) -> *mut ::std::os::raw::c_void>,
    pub realloc_func: Option::<
        unsafe extern "C" fn(*mut ::std::os::raw::c_void, size_t) -> *mut ::std::os::raw::c_void,
    >,
    pub free_func: Option::<unsafe extern "C" fn(*mut ::std::os::raw::c_void) -> ()>,
    pub printf_func: Option::<
        unsafe extern "C" fn(*const ::std::os::raw::c_char, ...) -> ::std::os::raw::c_int,
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
        ) -> ::std::os::raw::c_int,
    >,
}
pub const c_malloc: unsafe extern "C" fn(::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void = malloc;
pub const c_realloc: unsafe extern "C" fn(
    *mut ::std::os::raw::c_void,
    ::std::os::raw::c_ulong,
) -> *mut ::std::os::raw::c_void = realloc;
pub const c_free: unsafe extern "C" fn(*mut ::std::os::raw::c_void) -> () = free;
pub const NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const __DARWIN_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const SUITESPARSE_MAIN_VERSION: ::std::os::raw::c_int = 4 as ::std::os::raw::c_int;
pub const SUITESPARSE_SUB_VERSION: ::std::os::raw::c_int = 5 as ::std::os::raw::c_int;
pub const SUITESPARSE_SUBSUB_VERSION: ::std::os::raw::c_int = 3 as ::std::os::raw::c_int;
#[no_mangle]
pub static mut SuiteSparse_config: SuiteSparse_config_struct = unsafe {
    {
        let mut init = SuiteSparse_config_struct {
            malloc_func: Some(c_malloc),
            realloc_func: Some(c_realloc),
            free_func: Some(c_free),
            printf_func: ::std::mem::transmute::<
                isize,
                Option::<unsafe extern "C" fn(*const ::std::os::raw::c_char, ...) -> ::std::os::raw::c_int>,
            >(NULL as isize),
            hypot_func: Some(
                SuiteSparse_hypot as unsafe extern "C" fn(c_float, c_float) -> c_float,
            ),
            divcomplex_func: Some(
                SuiteSparse_divcomplex
                    as unsafe extern "C" fn(
                        c_float,
                        c_float,
                        c_float,
                        c_float,
                        *mut c_float,
                        *mut c_float,
                    ) -> ::std::os::raw::c_int,
            ),
        };
        init
    }
};
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_malloc(
    mut nitems: size_t,
    mut size_of_item: size_t,
) -> *mut ::std::os::raw::c_void {
    let mut p: *mut ::std::os::raw::c_void = 0 as *mut ::std::os::raw::c_void;
    let mut size: size_t = 0;
    if nitems < 1 as ::std::os::raw::c_int as ::std::os::raw::c_ulong {
        nitems = 1 as ::std::os::raw::c_int as size_t;
    }
    if size_of_item < 1 as ::std::os::raw::c_int as ::std::os::raw::c_ulong {
        size_of_item = 1 as ::std::os::raw::c_int as size_t;
    }
    size = nitems.wrapping_mul(size_of_item);
    if size as ::std::os::raw::c_double != nitems as c_float * size_of_item as ::std::os::raw::c_double {
        p = NULL as *mut ::std::os::raw::c_void;
    } else {
        p = (SuiteSparse_config.malloc_func).expect("non-null function pointer")(size);
    }
    return p;
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_realloc(
    mut nitems_new: size_t,
    mut nitems_old: size_t,
    mut size_of_item: size_t,
    mut p: *mut ::std::os::raw::c_void,
    mut ok: *mut ::std::os::raw::c_int,
) -> *mut ::std::os::raw::c_void {
    let mut size: size_t = 0;
    if nitems_old < 1 as ::std::os::raw::c_int as ::std::os::raw::c_ulong {
        nitems_old = 1 as ::std::os::raw::c_int as size_t;
    }
    if nitems_new < 1 as ::std::os::raw::c_int as ::std::os::raw::c_ulong {
        nitems_new = 1 as ::std::os::raw::c_int as size_t;
    }
    if size_of_item < 1 as ::std::os::raw::c_int as ::std::os::raw::c_ulong {
        size_of_item = 1 as ::std::os::raw::c_int as size_t;
    }
    size = nitems_new.wrapping_mul(size_of_item);
    if size as ::std::os::raw::c_double != nitems_new as c_float * size_of_item as ::std::os::raw::c_double {
        *ok = 0 as ::std::os::raw::c_int;
    } else if p.is_null() {
        p = SuiteSparse_malloc(nitems_new, size_of_item);
        *ok = (p != NULL as *mut ::std::os::raw::c_void) as ::std::os::raw::c_int;
    } else if nitems_old == nitems_new {
        *ok = 1 as ::std::os::raw::c_int;
    } else {
        let mut pnew: *mut ::std::os::raw::c_void = 0 as *mut ::std::os::raw::c_void;
        pnew = realloc(p, size);
        if pnew.is_null() {
            if nitems_new < nitems_old {
                *ok = 1 as ::std::os::raw::c_int;
            } else {
                *ok = 0 as ::std::os::raw::c_int;
            }
        } else {
            p = pnew;
            *ok = 1 as ::std::os::raw::c_int;
        }
    }
    return p;
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_free(
    mut p: *mut ::std::os::raw::c_void,
) -> *mut ::std::os::raw::c_void {
    if !p.is_null() {
        (SuiteSparse_config.free_func).expect("non-null function pointer")(p);
    }
    return 0 as *mut ::std::os::raw::c_void;
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_tic(mut tic: *mut c_float) {
    *tic.offset(0 as ::std::os::raw::c_int as isize) = 0 as ::std::os::raw::c_int as c_float;
    *tic.offset(1 as ::std::os::raw::c_int as isize) = 0 as ::std::os::raw::c_int as c_float;
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_toc(mut tic: *mut c_float) -> c_float {
    let mut toc: [c_float; 2] = [0.; 2];
    SuiteSparse_tic(toc.as_mut_ptr());
    return toc[0 as ::std::os::raw::c_int as usize] - *tic.offset(0 as ::std::os::raw::c_int as isize)
        + 1e-9f64
            * (toc[1 as ::std::os::raw::c_int as usize] - *tic.offset(1 as ::std::os::raw::c_int as isize));
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_time() -> c_float {
    let mut toc: [c_float; 2] = [0.; 2];
    SuiteSparse_tic(toc.as_mut_ptr());
    return toc[0 as ::std::os::raw::c_int as usize] + 1e-9f64 * toc[1 as ::std::os::raw::c_int as usize];
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_version(
    mut version: *mut ::std::os::raw::c_int,
) -> ::std::os::raw::c_int {
    if !version.is_null() {
        *version.offset(0 as ::std::os::raw::c_int as isize) = SUITESPARSE_MAIN_VERSION;
        *version.offset(1 as ::std::os::raw::c_int as isize) = SUITESPARSE_SUB_VERSION;
        *version.offset(2 as ::std::os::raw::c_int as isize) = SUITESPARSE_SUBSUB_VERSION;
    }
    return 4 as ::std::os::raw::c_int * 1000 as ::std::os::raw::c_int + 5 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_hypot(mut x: c_float, mut y: c_float) -> c_float {
    let mut s: c_float = 0.;
    let mut r: c_float = 0.;
    x = fabs(x);
    y = fabs(y);
    if x >= y {
        if x + y == x {
            s = x;
        } else {
            r = y / x;
            s = x * sqrt(1.0f64 + r * r);
        }
    } else if y + x == y {
        s = y;
    } else {
        r = x / y;
        s = y * sqrt(1.0f64 + r * r);
    }
    return s;
}
#[no_mangle]
pub unsafe extern "C" fn SuiteSparse_divcomplex(
    mut ar: c_float,
    mut ai: c_float,
    mut br: c_float,
    mut bi: c_float,
    mut cr: *mut c_float,
    mut ci: *mut c_float,
) -> ::std::os::raw::c_int {
    let mut tr: c_float = 0.;
    let mut ti: c_float = 0.;
    let mut r: c_float = 0.;
    let mut den: c_float = 0.;
    if fabs(br) >= fabs(bi) {
        r = bi / br;
        den = br + r * bi;
        tr = (ar + ai * r) / den;
        ti = (ai - ar * r) / den;
    } else {
        r = br / bi;
        den = r * br + bi;
        tr = (ar * r + ai) / den;
        ti = (ai * r - ar) / den;
    }
    *cr = tr;
    *ci = ti;
    return (den == 0.0f64) as ::std::os::raw::c_int;
}
