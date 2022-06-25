use ::libc;
extern "C" {
    fn dlsym(
        __handle: *mut libc::c_void,
        __symbol: *const libc::c_char,
    ) -> *mut libc::c_void;
    fn dlopen(__path: *const libc::c_char, __mode: libc::c_int) -> *mut libc::c_void;
    fn dlerror() -> *mut libc::c_char;
    fn dlclose(__handle: *mut libc::c_void) -> libc::c_int;
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
    fn __toupper(_: __darwin_ct_rune_t) -> __darwin_ct_rune_t;
    fn __tolower(_: __darwin_ct_rune_t) -> __darwin_ct_rune_t;
    fn c_strcpy(dest: *mut libc::c_char, source: *const libc::c_char);
}
pub type __darwin_ct_rune_t = libc::c_int;
pub type __darwin_size_t = libc::c_ulong;
pub type size_t = __darwin_size_t;
pub type c_int = libc::c_longlong;
pub type soHandle_t = *mut libc::c_void;
pub type symtype = *mut libc::c_void;
pub const RTLD_LAZY: libc::c_int = 0x1 as libc::c_int;
pub const c_print: unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int = printf;
#[no_mangle]
#[inline]
#[linkage = "external"]
pub unsafe extern "C" fn tolower(mut _c: libc::c_int) -> libc::c_int {
    return __tolower(_c);
}
#[no_mangle]
#[inline]
#[linkage = "external"]
pub unsafe extern "C" fn toupper(mut _c: libc::c_int) -> libc::c_int {
    return __toupper(_c);
}
pub const OSQP_NULL: libc::c_int = 0 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn lh_load_lib(mut libName: *const libc::c_char) -> soHandle_t {
    let mut h: soHandle_t = OSQP_NULL as soHandle_t;
    if libName.is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<&[u8; 12], &[libc::c_char; 12]>(b"lh_load_lib\0"))
                .as_ptr(),
        );
        printf(b"no library name given\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return OSQP_NULL as soHandle_t;
    }
    h = dlopen(libName, RTLD_LAZY);
    if h.is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<&[u8; 12], &[libc::c_char; 12]>(b"lh_load_lib\0"))
                .as_ptr(),
        );
        printf(
            b"Error while loading dynamic library %s: %s\0" as *const u8
                as *const libc::c_char,
            libName,
            dlerror(),
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
    }
    return h;
}
#[no_mangle]
pub unsafe extern "C" fn lh_unload_lib(mut h: soHandle_t) -> c_int {
    let mut rc: c_int = 1 as libc::c_int as c_int;
    rc = dlclose(h) as c_int;
    return rc;
}
#[no_mangle]
pub unsafe extern "C" fn lh_load_sym(
    mut h: soHandle_t,
    mut symName: *const libc::c_char,
) -> symtype {
    let mut s: symtype = 0 as *mut libc::c_void;
    let mut from: *const libc::c_char = 0 as *const libc::c_char;
    let mut to: *mut libc::c_char = 0 as *mut libc::c_char;
    let mut tripSym: *const libc::c_char = 0 as *const libc::c_char;
    let mut err: *mut libc::c_char = 0 as *mut libc::c_char;
    let mut lcbuf: [libc::c_char; 257] = [0; 257];
    let mut ucbuf: [libc::c_char; 257] = [0; 257];
    let mut ocbuf: [libc::c_char; 257] = [0; 257];
    let mut symLen: size_t = 0;
    let mut trip: libc::c_int = 0;
    s = OSQP_NULL as symtype;
    err = OSQP_NULL as *mut libc::c_char;
    symLen = 0 as libc::c_int as size_t;
    trip = 1 as libc::c_int;
    while trip <= 6 as libc::c_int {
        match trip {
            1 => {
                tripSym = symName;
            }
            2 => {
                from = symName;
                to = lcbuf.as_mut_ptr();
                while *from != 0 {
                    *to = tolower(*from as libc::c_int) as libc::c_char;
                    from = from.offset(1);
                    to = to.offset(1);
                }
                symLen = from.offset_from(symName) as libc::c_long as size_t;
                let fresh0 = to;
                to = to.offset(1);
                *fresh0 = '_' as i32 as libc::c_char;
                *to = '\u{0}' as i32 as libc::c_char;
                tripSym = lcbuf.as_mut_ptr();
            }
            3 => {
                from = symName;
                to = ucbuf.as_mut_ptr();
                while *from != 0 {
                    *to = toupper(*from as libc::c_int) as libc::c_char;
                    from = from.offset(1);
                    to = to.offset(1);
                }
                let fresh1 = to;
                to = to.offset(1);
                *fresh1 = '_' as i32 as libc::c_char;
                *to = '\u{0}' as i32 as libc::c_char;
                tripSym = ucbuf.as_mut_ptr();
            }
            4 => {
                c_strcpy(ocbuf.as_mut_ptr(), symName);
                ocbuf[symLen as usize] = '_' as i32 as libc::c_char;
                ocbuf[symLen.wrapping_add(1 as libc::c_int as libc::c_ulong)
                    as usize] = '\u{0}' as i32 as libc::c_char;
                tripSym = ocbuf.as_mut_ptr();
            }
            5 => {
                lcbuf[symLen as usize] = '\u{0}' as i32 as libc::c_char;
                tripSym = lcbuf.as_mut_ptr();
            }
            6 => {
                ucbuf[symLen as usize] = '\u{0}' as i32 as libc::c_char;
                tripSym = ucbuf.as_mut_ptr();
            }
            _ => {
                tripSym = symName;
            }
        }
        s = dlsym(h, tripSym);
        err = dlerror();
        if !err.is_null() {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 12],
                    &[libc::c_char; 12],
                >(b"lh_load_sym\0"))
                    .as_ptr(),
            );
            printf(
                b"Cannot find symbol %s in dynamic library, error = %s\0" as *const u8
                    as *const libc::c_char,
                symName,
                err,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
        } else {
            return s
        }
        trip += 1;
    }
    return OSQP_NULL as symtype;
}
