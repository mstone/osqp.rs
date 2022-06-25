use ::libc;
extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
    fn prea_vec_copy(a: *const c_float, b: *mut c_float, n: c_int);
    fn prea_int_vec_copy(a: *const c_int, b: *mut c_int, n: c_int);
}
pub type c_int = libc::c_longlong;
pub type c_float = libc::c_double;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct csc {
    pub nzmax: c_int,
    pub m: c_int,
    pub n: c_int,
    pub p: *mut c_int,
    pub i: *mut c_int,
    pub x: *mut c_float,
    pub nz: c_int,
}
pub const c_print: unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int = printf;
pub const c_free: unsafe extern "C" fn(*mut libc::c_void) -> () = free;
pub const c_calloc: unsafe extern "C" fn(
    libc::c_ulong,
    libc::c_ulong,
) -> *mut libc::c_void = calloc;
pub const OSQP_NULL: libc::c_int = 0 as libc::c_int;
pub const c_malloc: unsafe extern "C" fn(libc::c_ulong) -> *mut libc::c_void = malloc;
unsafe extern "C" fn csc_malloc(mut n: c_int, mut size: c_int) -> *mut libc::c_void {
    return malloc((n * size) as libc::c_ulong);
}
unsafe extern "C" fn csc_calloc(mut n: c_int, mut size: c_int) -> *mut libc::c_void {
    return calloc(n as libc::c_ulong, size as libc::c_ulong);
}
#[no_mangle]
pub unsafe extern "C" fn csc_matrix(
    mut m: c_int,
    mut n: c_int,
    mut nzmax: c_int,
    mut x: *mut c_float,
    mut i: *mut c_int,
    mut p: *mut c_int,
) -> *mut csc {
    let mut M: *mut csc = malloc(::std::mem::size_of::<csc>() as libc::c_ulong)
        as *mut csc;
    if M.is_null() {
        return OSQP_NULL as *mut csc;
    }
    (*M).m = m;
    (*M).n = n;
    (*M).nz = -(1 as libc::c_int) as c_int;
    (*M).nzmax = nzmax;
    let ref mut fresh0 = (*M).x;
    *fresh0 = x;
    let ref mut fresh1 = (*M).i;
    *fresh1 = i;
    let ref mut fresh2 = (*M).p;
    *fresh2 = p;
    return M;
}
#[no_mangle]
pub unsafe extern "C" fn csc_spalloc(
    mut m: c_int,
    mut n: c_int,
    mut nzmax: c_int,
    mut values: c_int,
    mut triplet: c_int,
) -> *mut csc {
    let mut A: *mut csc = csc_calloc(
        1 as libc::c_int as c_int,
        ::std::mem::size_of::<csc>() as libc::c_ulong as c_int,
    ) as *mut csc;
    if A.is_null() {
        return OSQP_NULL as *mut csc;
    }
    (*A).m = m;
    (*A).n = n;
    nzmax = if nzmax > 1 as libc::c_int as libc::c_longlong {
        nzmax
    } else {
        1 as libc::c_int as libc::c_longlong
    };
    (*A).nzmax = nzmax;
    (*A)
        .nz = (if triplet != 0 { 0 as libc::c_int } else { -(1 as libc::c_int) })
        as c_int;
    let ref mut fresh3 = (*A).p;
    *fresh3 = csc_malloc(
        if triplet != 0 { nzmax } else { n + 1 as libc::c_int as libc::c_longlong },
        ::std::mem::size_of::<c_int>() as libc::c_ulong as c_int,
    ) as *mut c_int;
    let ref mut fresh4 = (*A).i;
    *fresh4 = csc_malloc(nzmax, ::std::mem::size_of::<c_int>() as libc::c_ulong as c_int)
        as *mut c_int;
    let ref mut fresh5 = (*A).x;
    *fresh5 = (if values != 0 {
        csc_malloc(nzmax, ::std::mem::size_of::<c_float>() as libc::c_ulong as c_int)
    } else {
        OSQP_NULL as *mut libc::c_void
    }) as *mut c_float;
    if ((*A).p).is_null() || ((*A).i).is_null() || values != 0 && ((*A).x).is_null() {
        csc_spfree(A);
        return OSQP_NULL as *mut csc;
    } else {
        return A
    };
}
#[no_mangle]
pub unsafe extern "C" fn csc_spfree(mut A: *mut csc) {
    if !A.is_null() {
        if !((*A).p).is_null() {
            free((*A).p as *mut libc::c_void);
        }
        if !((*A).i).is_null() {
            free((*A).i as *mut libc::c_void);
        }
        if !((*A).x).is_null() {
            free((*A).x as *mut libc::c_void);
        }
        free(A as *mut libc::c_void);
    }
}
#[no_mangle]
pub unsafe extern "C" fn triplet_to_csc(
    mut T: *const csc,
    mut TtoC: *mut c_int,
) -> *mut csc {
    let mut m: c_int = 0;
    let mut n: c_int = 0;
    let mut nz: c_int = 0;
    let mut p: c_int = 0;
    let mut k: c_int = 0;
    let mut Cp: *mut c_int = 0 as *mut c_int;
    let mut Ci: *mut c_int = 0 as *mut c_int;
    let mut w: *mut c_int = 0 as *mut c_int;
    let mut Ti: *mut c_int = 0 as *mut c_int;
    let mut Tj: *mut c_int = 0 as *mut c_int;
    let mut Cx: *mut c_float = 0 as *mut c_float;
    let mut Tx: *mut c_float = 0 as *mut c_float;
    let mut C: *mut csc = 0 as *mut csc;
    m = (*T).m;
    n = (*T).n;
    Ti = (*T).i;
    Tj = (*T).p;
    Tx = (*T).x;
    nz = (*T).nz;
    C = csc_spalloc(
        m,
        n,
        nz,
        (Tx != OSQP_NULL as *mut c_float) as libc::c_int as c_int,
        0 as libc::c_int as c_int,
    );
    w = csc_calloc(n, ::std::mem::size_of::<c_int>() as libc::c_ulong as c_int)
        as *mut c_int;
    if C.is_null() || w.is_null() {
        return csc_done(
            C,
            w as *mut libc::c_void,
            OSQP_NULL as *mut libc::c_void,
            0 as libc::c_int as c_int,
        );
    }
    Cp = (*C).p;
    Ci = (*C).i;
    Cx = (*C).x;
    k = 0 as libc::c_int as c_int;
    while k < nz {
        let ref mut fresh6 = *w.offset(*Tj.offset(k as isize) as isize);
        *fresh6 += 1;
        k += 1;
    }
    csc_cumsum(Cp, w, n);
    k = 0 as libc::c_int as c_int;
    while k < nz {
        let ref mut fresh7 = *w.offset(*Tj.offset(k as isize) as isize);
        let fresh8 = *fresh7;
        *fresh7 = *fresh7 + 1;
        p = fresh8;
        *Ci.offset(p as isize) = *Ti.offset(k as isize);
        if !Cx.is_null() {
            *Cx.offset(p as isize) = *Tx.offset(k as isize);
            if !TtoC.is_null() {
                *TtoC.offset(k as isize) = p;
            }
        }
        k += 1;
    }
    return csc_done(
        C,
        w as *mut libc::c_void,
        OSQP_NULL as *mut libc::c_void,
        1 as libc::c_int as c_int,
    );
}
#[no_mangle]
pub unsafe extern "C" fn triplet_to_csr(
    mut T: *const csc,
    mut TtoC: *mut c_int,
) -> *mut csc {
    let mut m: c_int = 0;
    let mut n: c_int = 0;
    let mut nz: c_int = 0;
    let mut p: c_int = 0;
    let mut k: c_int = 0;
    let mut Cp: *mut c_int = 0 as *mut c_int;
    let mut Cj: *mut c_int = 0 as *mut c_int;
    let mut w: *mut c_int = 0 as *mut c_int;
    let mut Ti: *mut c_int = 0 as *mut c_int;
    let mut Tj: *mut c_int = 0 as *mut c_int;
    let mut Cx: *mut c_float = 0 as *mut c_float;
    let mut Tx: *mut c_float = 0 as *mut c_float;
    let mut C: *mut csc = 0 as *mut csc;
    m = (*T).m;
    n = (*T).n;
    Ti = (*T).i;
    Tj = (*T).p;
    Tx = (*T).x;
    nz = (*T).nz;
    C = csc_spalloc(
        m,
        n,
        nz,
        (Tx != OSQP_NULL as *mut c_float) as libc::c_int as c_int,
        0 as libc::c_int as c_int,
    );
    w = csc_calloc(m, ::std::mem::size_of::<c_int>() as libc::c_ulong as c_int)
        as *mut c_int;
    if C.is_null() || w.is_null() {
        return csc_done(
            C,
            w as *mut libc::c_void,
            OSQP_NULL as *mut libc::c_void,
            0 as libc::c_int as c_int,
        );
    }
    Cp = (*C).p;
    Cj = (*C).i;
    Cx = (*C).x;
    k = 0 as libc::c_int as c_int;
    while k < nz {
        let ref mut fresh9 = *w.offset(*Ti.offset(k as isize) as isize);
        *fresh9 += 1;
        k += 1;
    }
    csc_cumsum(Cp, w, m);
    k = 0 as libc::c_int as c_int;
    while k < nz {
        let ref mut fresh10 = *w.offset(*Ti.offset(k as isize) as isize);
        let fresh11 = *fresh10;
        *fresh10 = *fresh10 + 1;
        p = fresh11;
        *Cj.offset(p as isize) = *Tj.offset(k as isize);
        if !Cx.is_null() {
            *Cx.offset(p as isize) = *Tx.offset(k as isize);
            if !TtoC.is_null() {
                *TtoC.offset(k as isize) = p;
            }
        }
        k += 1;
    }
    return csc_done(
        C,
        w as *mut libc::c_void,
        OSQP_NULL as *mut libc::c_void,
        1 as libc::c_int as c_int,
    );
}
#[no_mangle]
pub unsafe extern "C" fn csc_cumsum(
    mut p: *mut c_int,
    mut c: *mut c_int,
    mut n: c_int,
) -> c_int {
    let mut i: c_int = 0;
    let mut nz: c_int = 0 as libc::c_int as c_int;
    if p.is_null() || c.is_null() {
        return -(1 as libc::c_int) as c_int;
    }
    i = 0 as libc::c_int as c_int;
    while i < n {
        *p.offset(i as isize) = nz;
        nz += *c.offset(i as isize);
        *c.offset(i as isize) = *p.offset(i as isize);
        i += 1;
    }
    *p.offset(n as isize) = nz;
    return nz;
}
#[no_mangle]
pub unsafe extern "C" fn csc_pinv(mut p: *const c_int, mut n: c_int) -> *mut c_int {
    let mut k: c_int = 0;
    let mut pinv: *mut c_int = 0 as *mut c_int;
    if p.is_null() {
        return OSQP_NULL as *mut c_int;
    }
    pinv = csc_malloc(n, ::std::mem::size_of::<c_int>() as libc::c_ulong as c_int)
        as *mut c_int;
    if pinv.is_null() {
        return OSQP_NULL as *mut c_int;
    }
    k = 0 as libc::c_int as c_int;
    while k < n {
        *pinv.offset(*p.offset(k as isize) as isize) = k;
        k += 1;
    }
    return pinv;
}
#[no_mangle]
pub unsafe extern "C" fn csc_symperm(
    mut A: *const csc,
    mut pinv: *const c_int,
    mut AtoC: *mut c_int,
    mut values: c_int,
) -> *mut csc {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut p: c_int = 0;
    let mut q: c_int = 0;
    let mut i2: c_int = 0;
    let mut j2: c_int = 0;
    let mut n: c_int = 0;
    let mut Ap: *mut c_int = 0 as *mut c_int;
    let mut Ai: *mut c_int = 0 as *mut c_int;
    let mut Cp: *mut c_int = 0 as *mut c_int;
    let mut Ci: *mut c_int = 0 as *mut c_int;
    let mut w: *mut c_int = 0 as *mut c_int;
    let mut Cx: *mut c_float = 0 as *mut c_float;
    let mut Ax: *mut c_float = 0 as *mut c_float;
    let mut C: *mut csc = 0 as *mut csc;
    n = (*A).n;
    Ap = (*A).p;
    Ai = (*A).i;
    Ax = (*A).x;
    C = csc_spalloc(
        n,
        n,
        *Ap.offset(n as isize),
        (values != 0 && !Ax.is_null()) as libc::c_int as c_int,
        0 as libc::c_int as c_int,
    );
    w = csc_calloc(n, ::std::mem::size_of::<c_int>() as libc::c_ulong as c_int)
        as *mut c_int;
    if C.is_null() || w.is_null() {
        return csc_done(
            C,
            w as *mut libc::c_void,
            OSQP_NULL as *mut libc::c_void,
            0 as libc::c_int as c_int,
        );
    }
    Cp = (*C).p;
    Ci = (*C).i;
    Cx = (*C).x;
    j = 0 as libc::c_int as c_int;
    while j < n {
        j2 = if !pinv.is_null() { *pinv.offset(j as isize) } else { j };
        p = *Ap.offset(j as isize);
        while p < *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize) {
            i = *Ai.offset(p as isize);
            if !(i > j) {
                i2 = if !pinv.is_null() { *pinv.offset(i as isize) } else { i };
                let ref mut fresh12 = *w
                    .offset((if i2 > j2 { i2 } else { j2 }) as isize);
                *fresh12 += 1;
            }
            p += 1;
        }
        j += 1;
    }
    csc_cumsum(Cp, w, n);
    j = 0 as libc::c_int as c_int;
    while j < n {
        j2 = if !pinv.is_null() { *pinv.offset(j as isize) } else { j };
        p = *Ap.offset(j as isize);
        while p < *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize) {
            i = *Ai.offset(p as isize);
            if !(i > j) {
                i2 = if !pinv.is_null() { *pinv.offset(i as isize) } else { i };
                let ref mut fresh13 = *w
                    .offset((if i2 > j2 { i2 } else { j2 }) as isize);
                let fresh14 = *fresh13;
                *fresh13 = *fresh13 + 1;
                q = fresh14;
                *Ci.offset(q as isize) = if i2 < j2 { i2 } else { j2 };
                if !Cx.is_null() {
                    *Cx.offset(q as isize) = *Ax.offset(p as isize);
                }
                if !AtoC.is_null() {
                    *AtoC.offset(p as isize) = q;
                }
            }
            p += 1;
        }
        j += 1;
    }
    return csc_done(
        C,
        w as *mut libc::c_void,
        OSQP_NULL as *mut libc::c_void,
        1 as libc::c_int as c_int,
    );
}
#[no_mangle]
pub unsafe extern "C" fn copy_csc_mat(mut A: *const csc) -> *mut csc {
    let mut B: *mut csc = csc_spalloc(
        (*A).m,
        (*A).n,
        *((*A).p).offset((*A).n as isize),
        1 as libc::c_int as c_int,
        0 as libc::c_int as c_int,
    );
    if B.is_null() {
        return OSQP_NULL as *mut csc;
    }
    prea_int_vec_copy((*A).p, (*B).p, (*A).n + 1 as libc::c_int as libc::c_longlong);
    prea_int_vec_copy((*A).i, (*B).i, *((*A).p).offset((*A).n as isize));
    prea_vec_copy((*A).x, (*B).x, *((*A).p).offset((*A).n as isize));
    return B;
}
#[no_mangle]
pub unsafe extern "C" fn prea_copy_csc_mat(mut A: *const csc, mut B: *mut csc) {
    prea_int_vec_copy((*A).p, (*B).p, (*A).n + 1 as libc::c_int as libc::c_longlong);
    prea_int_vec_copy((*A).i, (*B).i, *((*A).p).offset((*A).n as isize));
    prea_vec_copy((*A).x, (*B).x, *((*A).p).offset((*A).n as isize));
    (*B).nzmax = (*A).nzmax;
}
#[no_mangle]
pub unsafe extern "C" fn csc_done(
    mut C: *mut csc,
    mut w: *mut libc::c_void,
    mut x: *mut libc::c_void,
    mut ok: c_int,
) -> *mut csc {
    free(w);
    free(x);
    if ok != 0 {
        return C
    } else {
        csc_spfree(C);
        return OSQP_NULL as *mut csc;
    };
}
#[no_mangle]
pub unsafe extern "C" fn csc_to_triu(mut M: *mut csc) -> *mut csc {
    let mut M_trip: *mut csc = 0 as *mut csc;
    let mut M_triu: *mut csc = 0 as *mut csc;
    let mut nnzorigM: c_int = 0;
    let mut nnzmaxM: c_int = 0;
    let mut n: c_int = 0;
    let mut ptr: c_int = 0;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut z_M: c_int = 0 as libc::c_int as c_int;
    if (*M).m != (*M).n {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<&[u8; 12], &[libc::c_char; 12]>(b"csc_to_triu\0"))
                .as_ptr(),
        );
        printf(b"Matrix M not square\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return OSQP_NULL as *mut csc;
    }
    n = (*M).n;
    nnzorigM = *((*M).p).offset(n as isize);
    nnzmaxM = nnzorigM + n;
    M_trip = csc_spalloc(
        n,
        n,
        nnzmaxM,
        1 as libc::c_int as c_int,
        1 as libc::c_int as c_int,
    );
    if M_trip.is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<&[u8; 12], &[libc::c_char; 12]>(b"csc_to_triu\0"))
                .as_ptr(),
        );
        printf(
            b"Upper triangular matrix extraction failed (out of memory)\0" as *const u8
                as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return OSQP_NULL as *mut csc;
    }
    j = 0 as libc::c_int as c_int;
    while j < n {
        ptr = *((*M).p).offset(j as isize);
        while ptr < *((*M).p).offset((j + 1 as libc::c_int as libc::c_longlong) as isize)
        {
            i = *((*M).i).offset(ptr as isize);
            if i <= j {
                *((*M_trip).i).offset(z_M as isize) = i;
                *((*M_trip).p).offset(z_M as isize) = j;
                *((*M_trip).x).offset(z_M as isize) = *((*M).x).offset(ptr as isize);
                z_M += 1;
            }
            ptr += 1;
        }
        j += 1;
    }
    (*M_trip).nz = z_M;
    M_triu = triplet_to_csc(M_trip, OSQP_NULL as *mut c_int);
    (*M_triu).nzmax = nnzmaxM;
    csc_spfree(M_trip);
    return M_triu;
}
