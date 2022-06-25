use ::libc;
extern "C" {
    fn amd_l_valid(
        n_row: libc::c_longlong,
        n_col: libc::c_longlong,
        Ap: *const libc::c_longlong,
        Ai: *const libc::c_longlong,
    ) -> libc::c_longlong;
    fn amd_l_preprocess(
        n: libc::c_longlong,
        Ap: *const libc::c_longlong,
        Ai: *const libc::c_longlong,
        Rp: *mut libc::c_longlong,
        Ri: *mut libc::c_longlong,
        W: *mut libc::c_longlong,
        Flag: *mut libc::c_longlong,
    );
    fn amd_l_aat(
        n: libc::c_longlong,
        Ap: *const libc::c_longlong,
        Ai: *const libc::c_longlong,
        Len: *mut libc::c_longlong,
        Tp: *mut libc::c_longlong,
        Info: *mut c_float,
    ) -> size_t;
    fn SuiteSparse_malloc(nitems: size_t, size_of_item: size_t) -> *mut libc::c_void;
    fn SuiteSparse_free(p: *mut libc::c_void) -> *mut libc::c_void;
    fn amd_l1(
        n: libc::c_longlong,
        Ap: *const libc::c_longlong,
        Ai: *const libc::c_longlong,
        P: *mut libc::c_longlong,
        Pinv: *mut libc::c_longlong,
        Len: *mut libc::c_longlong,
        slen: libc::c_longlong,
        S: *mut libc::c_longlong,
        Control: *mut c_float,
        Info: *mut c_float,
    );
}
pub type __darwin_size_t = libc::c_ulong;
pub type size_t = __darwin_size_t;
pub type c_float = libc::c_double;
pub const AMD_INFO: libc::c_int = 20 as libc::c_int;
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
pub const AMD_N: libc::c_int = 1 as libc::c_int;
pub const AMD_OK: libc::c_int = 0 as libc::c_int;
pub const AMD_NZ: libc::c_int = 2 as libc::c_int;
pub const AMD_valid: unsafe extern "C" fn(
    libc::c_longlong,
    libc::c_longlong,
    *const libc::c_longlong,
    *const libc::c_longlong,
) -> libc::c_longlong = amd_l_valid;
pub const AMD_INVALID: libc::c_int = -(2 as libc::c_int);
pub const AMD_OK_BUT_JUMBLED: libc::c_int = 1 as libc::c_int;
pub const AMD_preprocess: unsafe extern "C" fn(
    libc::c_longlong,
    *const libc::c_longlong,
    *const libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
) -> () = amd_l_preprocess;
pub const AMD_aat: unsafe extern "C" fn(
    libc::c_longlong,
    *const libc::c_longlong,
    *const libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut c_float,
) -> size_t = amd_l_aat;
pub const NULL: libc::c_int = 0 as libc::c_int;
pub const SIZE_T_MAX: libc::c_int = -(1 as libc::c_int);
pub const Int_MAX: libc::c_long = 0x7fffffffffffffff as libc::c_long;
pub const SuiteSparse_long_max: libc::c_long = 0x7fffffffffffffff as libc::c_long;
pub const LONG_MAX: libc::c_long = 0x7fffffffffffffff as libc::c_long;
pub const AMD_OUT_OF_MEMORY: libc::c_int = -(1 as libc::c_int);
pub const AMD_MEMORY: libc::c_int = 7 as libc::c_int;
pub const AMD_1: unsafe extern "C" fn(
    libc::c_longlong,
    *const libc::c_longlong,
    *const libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    libc::c_longlong,
    *mut libc::c_longlong,
    *mut c_float,
    *mut c_float,
) -> () = amd_l1;
pub const AMD_STATUS: libc::c_int = 0 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_order(
    mut n: libc::c_longlong,
    mut Ap: *const libc::c_longlong,
    mut Ai: *const libc::c_longlong,
    mut P: *mut libc::c_longlong,
    mut Control: *mut c_float,
    mut Info: *mut c_float,
) -> libc::c_longlong {
    let mut Len: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut S: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut nz: libc::c_longlong = 0;
    let mut i: libc::c_longlong = 0;
    let mut Pinv: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut info: libc::c_longlong = 0;
    let mut status: libc::c_longlong = 0;
    let mut Rp: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Ri: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Cp: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Ci: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut ok: libc::c_longlong = 0;
    let mut nzaat: size_t = 0;
    let mut slen: size_t = 0;
    let mut mem: c_float = 0 as libc::c_int as c_float;
    info = (Info != 0 as *mut c_float) as libc::c_int as libc::c_longlong;
    if info != 0 {
        i = 0 as libc::c_int as libc::c_longlong;
        while i < AMD_INFO as libc::c_longlong {
            *Info.offset(i as isize) = EMPTY as c_float;
            i += 1;
        }
        *Info.offset(AMD_N as isize) = n as c_float;
        *Info.offset(AMD_STATUS as isize) = AMD_OK as c_float;
    }
    if Ai.is_null() || Ap.is_null() || P.is_null()
        || n < 0 as libc::c_int as libc::c_longlong
    {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_INVALID as c_float;
        }
        return -(2 as libc::c_int) as libc::c_longlong;
    }
    if n == 0 as libc::c_int as libc::c_longlong {
        return 0 as libc::c_int as libc::c_longlong;
    }
    nz = *Ap.offset(n as isize);
    if info != 0 {
        *Info.offset(AMD_NZ as isize) = nz as c_float;
    }
    if nz < 0 as libc::c_int as libc::c_longlong {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_INVALID as c_float;
        }
        return -(2 as libc::c_int) as libc::c_longlong;
    }
    if n as size_t
        >= (SIZE_T_MAX as size_t)
            .wrapping_div(::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong)
        || nz as size_t
            >= (SIZE_T_MAX as size_t)
                .wrapping_div(::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong)
    {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
        }
        return -(1 as libc::c_int) as libc::c_longlong;
    }
    status = amd_l_valid(n, n, Ap, Ai);
    if status == AMD_INVALID as libc::c_longlong {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_INVALID as c_float;
        }
        return -(2 as libc::c_int) as libc::c_longlong;
    }
    Len = SuiteSparse_malloc(
        n as size_t,
        ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
    ) as *mut libc::c_longlong;
    Pinv = SuiteSparse_malloc(
        n as size_t,
        ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
    ) as *mut libc::c_longlong;
    mem += n as libc::c_double;
    mem += n as libc::c_double;
    if Len.is_null() || Pinv.is_null() {
        SuiteSparse_free(Len as *mut libc::c_void);
        SuiteSparse_free(Pinv as *mut libc::c_void);
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
        }
        return -(1 as libc::c_int) as libc::c_longlong;
    }
    if status == AMD_OK_BUT_JUMBLED as libc::c_longlong {
        Rp = SuiteSparse_malloc(
            (n + 1 as libc::c_int as libc::c_longlong) as size_t,
            ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
        ) as *mut libc::c_longlong;
        Ri = SuiteSparse_malloc(
            nz as size_t,
            ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
        ) as *mut libc::c_longlong;
        mem += (n + 1 as libc::c_int as libc::c_longlong) as libc::c_double;
        mem
            += (if nz > 1 as libc::c_int as libc::c_longlong {
                nz
            } else {
                1 as libc::c_int as libc::c_longlong
            }) as libc::c_double;
        if Rp.is_null() || Ri.is_null() {
            SuiteSparse_free(Rp as *mut libc::c_void);
            SuiteSparse_free(Ri as *mut libc::c_void);
            SuiteSparse_free(Len as *mut libc::c_void);
            SuiteSparse_free(Pinv as *mut libc::c_void);
            if info != 0 {
                *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
            }
            return -(1 as libc::c_int) as libc::c_longlong;
        }
        amd_l_preprocess(n, Ap, Ai, Rp, Ri, Len, Pinv);
        Cp = Rp;
        Ci = Ri;
    } else {
        Rp = NULL as *mut libc::c_longlong;
        Ri = NULL as *mut libc::c_longlong;
        Cp = Ap as *mut libc::c_longlong;
        Ci = Ai as *mut libc::c_longlong;
    }
    nzaat = amd_l_aat(
        n,
        Cp as *const libc::c_longlong,
        Ci as *const libc::c_longlong,
        Len,
        P,
        Info,
    );
    S = NULL as *mut libc::c_longlong;
    slen = nzaat;
    ok = (slen.wrapping_add(nzaat.wrapping_div(5 as libc::c_int as libc::c_ulong))
        >= slen) as libc::c_int as libc::c_longlong;
    slen = (slen as libc::c_ulong)
        .wrapping_add(nzaat.wrapping_div(5 as libc::c_int as libc::c_ulong)) as size_t
        as size_t;
    i = 0 as libc::c_int as libc::c_longlong;
    while ok != 0 && i < 7 as libc::c_int as libc::c_longlong {
        ok = ((slen as libc::c_ulonglong).wrapping_add(n as libc::c_ulonglong)
            > slen as libc::c_ulonglong) as libc::c_int as libc::c_longlong;
        slen = (slen as libc::c_ulonglong).wrapping_add(n as libc::c_ulonglong) as size_t
            as size_t;
        i += 1;
    }
    mem += slen as libc::c_double;
    ok = (ok != 0
        && slen
            < (SIZE_T_MAX as size_t)
                .wrapping_div(
                    ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
                )) as libc::c_int as libc::c_longlong;
    ok = (ok != 0 && slen < Int_MAX as libc::c_ulong) as libc::c_int as libc::c_longlong;
    if ok != 0 {
        S = SuiteSparse_malloc(
            slen,
            ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong,
        ) as *mut libc::c_longlong;
    }
    if S.is_null() {
        SuiteSparse_free(Rp as *mut libc::c_void);
        SuiteSparse_free(Ri as *mut libc::c_void);
        SuiteSparse_free(Len as *mut libc::c_void);
        SuiteSparse_free(Pinv as *mut libc::c_void);
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
        }
        return -(1 as libc::c_int) as libc::c_longlong;
    }
    if info != 0 {
        *Info
            .offset(
                AMD_MEMORY as isize,
            ) = mem
            * ::std::mem::size_of::<libc::c_longlong>() as libc::c_ulong
                as libc::c_double;
    }
    amd_l1(
        n,
        Cp as *const libc::c_longlong,
        Ci as *const libc::c_longlong,
        P,
        Pinv,
        Len,
        slen as libc::c_longlong,
        S,
        Control,
        Info,
    );
    SuiteSparse_free(Rp as *mut libc::c_void);
    SuiteSparse_free(Ri as *mut libc::c_void);
    SuiteSparse_free(Len as *mut libc::c_void);
    SuiteSparse_free(Pinv as *mut libc::c_void);
    SuiteSparse_free(S as *mut libc::c_void);
    if info != 0 {
        *Info.offset(AMD_STATUS as isize) = status as c_float;
    }
    return status;
}
