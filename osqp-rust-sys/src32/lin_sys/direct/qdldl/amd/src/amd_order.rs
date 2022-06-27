extern "C" {
    fn SuiteSparse_malloc(nitems: size_t, size_of_item: size_t) -> *mut ::std::os::raw::c_void;
    fn SuiteSparse_free(p: *mut ::std::os::raw::c_void) -> *mut ::std::os::raw::c_void;
    fn amd_1(
        n: ::std::os::raw::c_int,
        Ap: *const ::std::os::raw::c_int,
        Ai: *const ::std::os::raw::c_int,
        P: *mut ::std::os::raw::c_int,
        Pinv: *mut ::std::os::raw::c_int,
        Len: *mut ::std::os::raw::c_int,
        slen: ::std::os::raw::c_int,
        S: *mut ::std::os::raw::c_int,
        Control: *mut c_float,
        Info: *mut c_float,
    );
    fn amd_aat(
        n: ::std::os::raw::c_int,
        Ap: *const ::std::os::raw::c_int,
        Ai: *const ::std::os::raw::c_int,
        Len: *mut ::std::os::raw::c_int,
        Tp: *mut ::std::os::raw::c_int,
        Info: *mut c_float,
    ) -> size_t;
    fn amd_preprocess(
        n: ::std::os::raw::c_int,
        Ap: *const ::std::os::raw::c_int,
        Ai: *const ::std::os::raw::c_int,
        Rp: *mut ::std::os::raw::c_int,
        Ri: *mut ::std::os::raw::c_int,
        W: *mut ::std::os::raw::c_int,
        Flag: *mut ::std::os::raw::c_int,
    );
    fn amd_valid(
        n_row: ::std::os::raw::c_int,
        n_col: ::std::os::raw::c_int,
        Ap: *const ::std::os::raw::c_int,
        Ai: *const ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}
pub type __darwin_size_t = ::std::os::raw::c_ulong;
pub type size_t = __darwin_size_t;
pub type c_float = ::std::os::raw::c_double;
pub const AMD_STATUS: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_1: unsafe extern "C" fn(
    ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut c_float,
    *mut c_float,
) -> () = amd_1;
pub const AMD_MEMORY: ::std::os::raw::c_int = 7 as ::std::os::raw::c_int;
pub const AMD_OUT_OF_MEMORY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const INT_MAX: ::std::os::raw::c_int = 2147483647 as ::std::os::raw::c_int;
pub const Int_MAX: ::std::os::raw::c_int = 2147483647 as ::std::os::raw::c_int;
pub const SIZE_T_MAX: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_aat: unsafe extern "C" fn(
    ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut c_float,
) -> size_t = amd_aat;
pub const AMD_preprocess: unsafe extern "C" fn(
    ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
) -> () = amd_preprocess;
pub const AMD_OK_BUT_JUMBLED: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const AMD_INVALID: ::std::os::raw::c_int = -(2 as ::std::os::raw::c_int);
pub const AMD_valid: unsafe extern "C" fn(
    ::std::os::raw::c_int,
    ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
) -> ::std::os::raw::c_int = amd_valid;
pub const AMD_NZ: ::std::os::raw::c_int = 2 as ::std::os::raw::c_int;
pub const AMD_OK: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_N: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const AMD_INFO: ::std::os::raw::c_int = 20 as ::std::os::raw::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_order(
    mut n: ::std::os::raw::c_int,
    mut Ap: *const ::std::os::raw::c_int,
    mut Ai: *const ::std::os::raw::c_int,
    mut P: *mut ::std::os::raw::c_int,
    mut Control: *mut c_float,
    mut Info: *mut c_float,
) -> ::std::os::raw::c_int {
    let mut Len: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut S: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut nz: ::std::os::raw::c_int = 0;
    let mut i: ::std::os::raw::c_int = 0;
    let mut Pinv: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut info: ::std::os::raw::c_int = 0;
    let mut status: ::std::os::raw::c_int = 0;
    let mut Rp: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut Ri: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut Cp: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut Ci: *mut ::std::os::raw::c_int = 0 as *mut ::std::os::raw::c_int;
    let mut ok: ::std::os::raw::c_int = 0;
    let mut nzaat: size_t = 0;
    let mut slen: size_t = 0;
    let mut mem: c_float = 0 as ::std::os::raw::c_int as c_float;
    info = (Info != 0 as *mut c_float) as ::std::os::raw::c_int;
    if info != 0 {
        i = 0 as ::std::os::raw::c_int;
        while i < AMD_INFO {
            *Info.offset(i as isize) = EMPTY as c_float;
            i += 1;
        }
        *Info.offset(AMD_N as isize) = n as c_float;
        *Info.offset(AMD_STATUS as isize) = AMD_OK as c_float;
    }
    if Ai.is_null() || Ap.is_null() || P.is_null() || n < 0 as ::std::os::raw::c_int {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_INVALID as c_float;
        }
        return -(2 as ::std::os::raw::c_int);
    }
    if n == 0 as ::std::os::raw::c_int {
        return 0 as ::std::os::raw::c_int;
    }
    nz = *Ap.offset(n as isize);
    if info != 0 {
        *Info.offset(AMD_NZ as isize) = nz as c_float;
    }
    if nz < 0 as ::std::os::raw::c_int {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_INVALID as c_float;
        }
        return -(2 as ::std::os::raw::c_int);
    }
    if n as size_t
        >= (SIZE_T_MAX as size_t)
            .wrapping_div(::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong)
        || nz as size_t
            >= (SIZE_T_MAX as size_t)
                .wrapping_div(::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong)
    {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
        }
        return -(1 as ::std::os::raw::c_int);
    }
    status = amd_valid(n, n, Ap, Ai);
    if status == AMD_INVALID {
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_INVALID as c_float;
        }
        return -(2 as ::std::os::raw::c_int);
    }
    Len = SuiteSparse_malloc(
        n as size_t,
        ::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong,
    ) as *mut ::std::os::raw::c_int;
    Pinv = SuiteSparse_malloc(
        n as size_t,
        ::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong,
    ) as *mut ::std::os::raw::c_int;
    mem += n as ::std::os::raw::c_double;
    mem += n as ::std::os::raw::c_double;
    if Len.is_null() || Pinv.is_null() {
        SuiteSparse_free(Len as *mut ::std::os::raw::c_void);
        SuiteSparse_free(Pinv as *mut ::std::os::raw::c_void);
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
        }
        return -(1 as ::std::os::raw::c_int);
    }
    if status == AMD_OK_BUT_JUMBLED {
        Rp = SuiteSparse_malloc(
            (n + 1 as ::std::os::raw::c_int) as size_t,
            ::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong,
        ) as *mut ::std::os::raw::c_int;
        Ri = SuiteSparse_malloc(
            nz as size_t,
            ::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong,
        ) as *mut ::std::os::raw::c_int;
        mem += (n + 1 as ::std::os::raw::c_int) as ::std::os::raw::c_double;
        mem
            += (if nz > 1 as ::std::os::raw::c_int { nz } else { 1 as ::std::os::raw::c_int })
                as ::std::os::raw::c_double;
        if Rp.is_null() || Ri.is_null() {
            SuiteSparse_free(Rp as *mut ::std::os::raw::c_void);
            SuiteSparse_free(Ri as *mut ::std::os::raw::c_void);
            SuiteSparse_free(Len as *mut ::std::os::raw::c_void);
            SuiteSparse_free(Pinv as *mut ::std::os::raw::c_void);
            if info != 0 {
                *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
            }
            return -(1 as ::std::os::raw::c_int);
        }
        amd_preprocess(n, Ap, Ai, Rp, Ri, Len, Pinv);
        Cp = Rp;
        Ci = Ri;
    } else {
        Rp = NULL as *mut ::std::os::raw::c_int;
        Ri = NULL as *mut ::std::os::raw::c_int;
        Cp = Ap as *mut ::std::os::raw::c_int;
        Ci = Ai as *mut ::std::os::raw::c_int;
    }
    nzaat = amd_aat(n, Cp as *const ::std::os::raw::c_int, Ci as *const ::std::os::raw::c_int, Len, P, Info);
    S = NULL as *mut ::std::os::raw::c_int;
    slen = nzaat;
    ok = (slen.wrapping_add(nzaat.wrapping_div(5 as ::std::os::raw::c_int as ::std::os::raw::c_ulong))
        >= slen) as ::std::os::raw::c_int;
    slen = (slen as ::std::os::raw::c_ulong)
        .wrapping_add(nzaat.wrapping_div(5 as ::std::os::raw::c_int as ::std::os::raw::c_ulong)) as size_t
        as size_t;
    i = 0 as ::std::os::raw::c_int;
    while ok != 0 && i < 7 as ::std::os::raw::c_int {
        ok = (slen.wrapping_add(n as ::std::os::raw::c_ulong) > slen) as ::std::os::raw::c_int;
        slen = (slen as ::std::os::raw::c_ulong).wrapping_add(n as ::std::os::raw::c_ulong) as size_t
            as size_t;
        i += 1;
    }
    mem += slen as ::std::os::raw::c_double;
    ok = (ok != 0
        && slen
            < (SIZE_T_MAX as size_t)
                .wrapping_div(::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong))
        as ::std::os::raw::c_int;
    ok = (ok != 0 && slen < Int_MAX as ::std::os::raw::c_ulong) as ::std::os::raw::c_int;
    if ok != 0 {
        S = SuiteSparse_malloc(
            slen,
            ::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong,
        ) as *mut ::std::os::raw::c_int;
    }
    if S.is_null() {
        SuiteSparse_free(Rp as *mut ::std::os::raw::c_void);
        SuiteSparse_free(Ri as *mut ::std::os::raw::c_void);
        SuiteSparse_free(Len as *mut ::std::os::raw::c_void);
        SuiteSparse_free(Pinv as *mut ::std::os::raw::c_void);
        if info != 0 {
            *Info.offset(AMD_STATUS as isize) = AMD_OUT_OF_MEMORY as c_float;
        }
        return -(1 as ::std::os::raw::c_int);
    }
    if info != 0 {
        *Info
            .offset(
                AMD_MEMORY as isize,
            ) = mem
            * ::std::mem::size_of::<::std::os::raw::c_int>() as ::std::os::raw::c_ulong as ::std::os::raw::c_double;
    }
    amd_1(
        n,
        Cp as *const ::std::os::raw::c_int,
        Ci as *const ::std::os::raw::c_int,
        P,
        Pinv,
        Len,
        slen as ::std::os::raw::c_int,
        S,
        Control,
        Info,
    );
    SuiteSparse_free(Rp as *mut ::std::os::raw::c_void);
    SuiteSparse_free(Ri as *mut ::std::os::raw::c_void);
    SuiteSparse_free(Len as *mut ::std::os::raw::c_void);
    SuiteSparse_free(Pinv as *mut ::std::os::raw::c_void);
    SuiteSparse_free(S as *mut ::std::os::raw::c_void);
    if info != 0 {
        *Info.offset(AMD_STATUS as isize) = status as c_float;
    }
    return status;
}
