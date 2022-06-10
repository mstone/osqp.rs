use ::libc;
pub type __darwin_size_t = libc::c_ulong;
pub type size_t = __darwin_size_t;
pub type c_float = libc::c_double;
pub const AMD_INFO: libc::c_int = 20 as libc::c_int;
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
pub const NULL: libc::c_int = 0 as libc::c_int;
pub const AMD_STATUS: libc::c_int = 0 as libc::c_int;
pub const AMD_OK: libc::c_int = 0 as libc::c_int;
pub const AMD_N: libc::c_int = 1 as libc::c_int;
pub const AMD_NZ: libc::c_int = 2 as libc::c_int;
pub const AMD_SYMMETRY: libc::c_int = 3 as libc::c_int;
pub const AMD_NZDIAG: libc::c_int = 4 as libc::c_int;
pub const AMD_NZ_A_PLUS_AT: libc::c_int = 5 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_aat(
    mut n: libc::c_longlong,
    mut Ap: *const libc::c_longlong,
    mut Ai: *const libc::c_longlong,
    mut Len: *mut libc::c_longlong,
    mut Tp: *mut libc::c_longlong,
    mut Info: *mut c_float,
) -> size_t {
    let mut p1: libc::c_longlong = 0;
    let mut p2: libc::c_longlong = 0;
    let mut p: libc::c_longlong = 0;
    let mut i: libc::c_longlong = 0;
    let mut j: libc::c_longlong = 0;
    let mut pj: libc::c_longlong = 0;
    let mut pj2: libc::c_longlong = 0;
    let mut k: libc::c_longlong = 0;
    let mut nzdiag: libc::c_longlong = 0;
    let mut nzboth: libc::c_longlong = 0;
    let mut nz: libc::c_longlong = 0;
    let mut sym: c_float = 0.;
    let mut nzaat: size_t = 0;
    if !Info.is_null() {
        i = 0 as libc::c_int as libc::c_longlong;
        while i < AMD_INFO as libc::c_longlong {
            *Info.offset(i as isize) = EMPTY as c_float;
            i += 1;
        }
        *Info.offset(AMD_STATUS as isize) = AMD_OK as c_float;
    }
    k = 0 as libc::c_int as libc::c_longlong;
    while k < n {
        *Len.offset(k as isize) = 0 as libc::c_int as libc::c_longlong;
        k += 1;
    }
    nzdiag = 0 as libc::c_int as libc::c_longlong;
    nzboth = 0 as libc::c_int as libc::c_longlong;
    nz = *Ap.offset(n as isize);
    k = 0 as libc::c_int as libc::c_longlong;
    while k < n {
        p1 = *Ap.offset(k as isize);
        p2 = *Ap.offset((k + 1 as libc::c_int as libc::c_longlong) as isize);
        p = p1;
        while p < p2 {
            j = *Ai.offset(p as isize);
            if j < k {
                let ref mut fresh0 = *Len.offset(j as isize);
                *fresh0 += 1;
                let ref mut fresh1 = *Len.offset(k as isize);
                *fresh1 += 1;
                p += 1;
                pj2 = *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize);
                pj = *Tp.offset(j as isize);
                while pj < pj2 {
                    i = *Ai.offset(pj as isize);
                    if i < k {
                        let ref mut fresh2 = *Len.offset(i as isize);
                        *fresh2 += 1;
                        let ref mut fresh3 = *Len.offset(j as isize);
                        *fresh3 += 1;
                        pj += 1;
                    } else if i == k {
                        pj += 1;
                        nzboth += 1;
                        break;
                    } else {
                        break;
                    }
                }
                *Tp.offset(j as isize) = pj;
            } else if j == k {
                p += 1;
                nzdiag += 1;
                break;
            } else {
                break;
            }
        }
        *Tp.offset(k as isize) = p;
        k += 1;
    }
    j = 0 as libc::c_int as libc::c_longlong;
    while j < n {
        pj = *Tp.offset(j as isize);
        while pj < *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize) {
            i = *Ai.offset(pj as isize);
            let ref mut fresh4 = *Len.offset(i as isize);
            *fresh4 += 1;
            let ref mut fresh5 = *Len.offset(j as isize);
            *fresh5 += 1;
            pj += 1;
        }
        j += 1;
    }
    if nz == nzdiag {
        sym = 1 as libc::c_int as c_float;
    } else {
        sym = 2 as libc::c_int as libc::c_double * nzboth as c_float
            / (nz - nzdiag) as c_float;
    }
    nzaat = 0 as libc::c_int as size_t;
    k = 0 as libc::c_int as libc::c_longlong;
    while k < n {
        nzaat = (nzaat as libc::c_ulonglong)
            .wrapping_add(*Len.offset(k as isize) as libc::c_ulonglong) as size_t
            as size_t;
        k += 1;
    }
    if !Info.is_null() {
        *Info.offset(AMD_STATUS as isize) = AMD_OK as c_float;
        *Info.offset(AMD_N as isize) = n as c_float;
        *Info.offset(AMD_NZ as isize) = nz as c_float;
        *Info.offset(AMD_SYMMETRY as isize) = sym;
        *Info.offset(AMD_NZDIAG as isize) = nzdiag as c_float;
        *Info.offset(AMD_NZ_A_PLUS_AT as isize) = nzaat as c_float;
    }
    return nzaat;
}
