pub type QDLDL_int = ::std::os::raw::c_int;
pub type QDLDL_float = ::std::os::raw::c_double;
pub type QDLDL_bool = ::std::os::raw::c_uchar;
pub const INT_MAX: ::std::os::raw::c_int = 2147483647 as ::std::os::raw::c_int;
pub const QDLDL_INT_MAX: ::std::os::raw::c_int = INT_MAX;
pub const QDLDL_UNKNOWN: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const QDLDL_USED: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const QDLDL_UNUSED: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
#[no_mangle]
pub unsafe extern "C" fn QDLDL_etree(
    n: QDLDL_int,
    mut Ap: *const QDLDL_int,
    mut Ai: *const QDLDL_int,
    mut work: *mut QDLDL_int,
    mut Lnz: *mut QDLDL_int,
    mut etree: *mut QDLDL_int,
) -> QDLDL_int {
    let mut sumLnz: QDLDL_int = 0;
    let mut i: QDLDL_int = 0;
    let mut j: QDLDL_int = 0;
    let mut p: QDLDL_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *work.offset(i as isize) = 0 as ::std::os::raw::c_int;
        *Lnz.offset(i as isize) = 0 as ::std::os::raw::c_int;
        *etree.offset(i as isize) = QDLDL_UNKNOWN;
        if *Ap.offset(i as isize) == *Ap.offset((i + 1 as ::std::os::raw::c_int) as isize) {
            return -(1 as ::std::os::raw::c_int);
        }
        i += 1;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < n {
        *work.offset(j as isize) = j;
        p = *Ap.offset(j as isize);
        while p < *Ap.offset((j + 1 as ::std::os::raw::c_int) as isize) {
            i = *Ai.offset(p as isize);
            if i > j {
                return -(1 as ::std::os::raw::c_int);
            }
            while *work.offset(i as isize) != j {
                if *etree.offset(i as isize) == QDLDL_UNKNOWN {
                    *etree.offset(i as isize) = j;
                }
                let ref mut fresh0 = *Lnz.offset(i as isize);
                *fresh0 += 1;
                *work.offset(i as isize) = j;
                i = *etree.offset(i as isize);
            }
            p += 1;
        }
        j += 1;
    }
    sumLnz = 0 as ::std::os::raw::c_int;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        if sumLnz > QDLDL_INT_MAX - *Lnz.offset(i as isize) {
            sumLnz = -(2 as ::std::os::raw::c_int);
            break;
        } else {
            sumLnz += *Lnz.offset(i as isize);
            i += 1;
        }
    }
    return sumLnz;
}
#[no_mangle]
pub unsafe extern "C" fn QDLDL_factor(
    n: QDLDL_int,
    mut Ap: *const QDLDL_int,
    mut Ai: *const QDLDL_int,
    mut Ax: *const QDLDL_float,
    mut Lp: *mut QDLDL_int,
    mut Li: *mut QDLDL_int,
    mut Lx: *mut QDLDL_float,
    mut D: *mut QDLDL_float,
    mut Dinv: *mut QDLDL_float,
    mut Lnz: *const QDLDL_int,
    mut etree: *const QDLDL_int,
    mut bwork: *mut QDLDL_bool,
    mut iwork: *mut QDLDL_int,
    mut fwork: *mut QDLDL_float,
) -> QDLDL_int {
    let mut i: QDLDL_int = 0;
    let mut j: QDLDL_int = 0;
    let mut k: QDLDL_int = 0;
    let mut nnzY: QDLDL_int = 0;
    let mut bidx: QDLDL_int = 0;
    let mut cidx: QDLDL_int = 0;
    let mut nextIdx: QDLDL_int = 0;
    let mut nnzE: QDLDL_int = 0;
    let mut tmpIdx: QDLDL_int = 0;
    let mut yIdx: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut elimBuffer: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut LNextSpaceInCol: *mut QDLDL_int = 0 as *mut QDLDL_int;
    let mut yVals: *mut QDLDL_float = 0 as *mut QDLDL_float;
    let mut yVals_cidx: QDLDL_float = 0.;
    let mut yMarkers: *mut QDLDL_bool = 0 as *mut QDLDL_bool;
    let mut positiveValuesInD: QDLDL_int = 0 as ::std::os::raw::c_int;
    yMarkers = bwork;
    yIdx = iwork;
    elimBuffer = iwork.offset(n as isize);
    LNextSpaceInCol = iwork.offset((n * 2 as ::std::os::raw::c_int) as isize);
    yVals = fwork;
    *Lp.offset(0 as ::std::os::raw::c_int as isize) = 0 as ::std::os::raw::c_int;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *Lp
            .offset(
                (i + 1 as ::std::os::raw::c_int) as isize,
            ) = *Lp.offset(i as isize) + *Lnz.offset(i as isize);
        *yMarkers.offset(i as isize) = QDLDL_UNUSED as QDLDL_bool;
        *yVals.offset(i as isize) = 0.0f64;
        *D.offset(i as isize) = 0.0f64;
        *LNextSpaceInCol.offset(i as isize) = *Lp.offset(i as isize);
        i += 1;
    }
    *D.offset(0 as ::std::os::raw::c_int as isize) = *Ax.offset(0 as ::std::os::raw::c_int as isize);
    if *D.offset(0 as ::std::os::raw::c_int as isize) == 0.0f64 {
        return -(1 as ::std::os::raw::c_int);
    }
    if *D.offset(0 as ::std::os::raw::c_int as isize) > 0.0f64 {
        positiveValuesInD += 1;
    }
    *Dinv
        .offset(
            0 as ::std::os::raw::c_int as isize,
        ) = 1 as ::std::os::raw::c_int as ::std::os::raw::c_double / *D.offset(0 as ::std::os::raw::c_int as isize);
    k = 1 as ::std::os::raw::c_int;
    while k < n {
        nnzY = 0 as ::std::os::raw::c_int;
        tmpIdx = *Ap.offset((k + 1 as ::std::os::raw::c_int) as isize);
        i = *Ap.offset(k as isize);
        while i < tmpIdx {
            bidx = *Ai.offset(i as isize);
            if bidx == k {
                *D.offset(k as isize) = *Ax.offset(i as isize);
            } else {
                *yVals.offset(bidx as isize) = *Ax.offset(i as isize);
                nextIdx = bidx;
                if *yMarkers.offset(nextIdx as isize) as ::std::os::raw::c_int == QDLDL_UNUSED {
                    *yMarkers.offset(nextIdx as isize) = QDLDL_USED as QDLDL_bool;
                    *elimBuffer.offset(0 as ::std::os::raw::c_int as isize) = nextIdx;
                    nnzE = 1 as ::std::os::raw::c_int;
                    nextIdx = *etree.offset(bidx as isize);
                    while nextIdx != QDLDL_UNKNOWN && nextIdx < k {
                        if *yMarkers.offset(nextIdx as isize) as ::std::os::raw::c_int
                            == QDLDL_USED
                        {
                            break;
                        }
                        *yMarkers.offset(nextIdx as isize) = QDLDL_USED as QDLDL_bool;
                        *elimBuffer.offset(nnzE as isize) = nextIdx;
                        nnzE += 1;
                        nextIdx = *etree.offset(nextIdx as isize);
                    }
                    while nnzE != 0 {
                        nnzE -= 1;
                        let fresh1 = nnzY;
                        nnzY = nnzY + 1;
                        *yIdx
                            .offset(fresh1 as isize) = *elimBuffer.offset(nnzE as isize);
                    }
                }
            }
            i += 1;
        }
        i = nnzY - 1 as ::std::os::raw::c_int;
        while i >= 0 as ::std::os::raw::c_int {
            cidx = *yIdx.offset(i as isize);
            tmpIdx = *LNextSpaceInCol.offset(cidx as isize);
            yVals_cidx = *yVals.offset(cidx as isize);
            j = *Lp.offset(cidx as isize);
            while j < tmpIdx {
                let ref mut fresh2 = *yVals.offset(*Li.offset(j as isize) as isize);
                *fresh2 -= *Lx.offset(j as isize) * yVals_cidx;
                j += 1;
            }
            *Li.offset(tmpIdx as isize) = k;
            *Lx.offset(tmpIdx as isize) = yVals_cidx * *Dinv.offset(cidx as isize);
            let ref mut fresh3 = *D.offset(k as isize);
            *fresh3 -= yVals_cidx * *Lx.offset(tmpIdx as isize);
            let ref mut fresh4 = *LNextSpaceInCol.offset(cidx as isize);
            *fresh4 += 1;
            *yVals.offset(cidx as isize) = 0.0f64;
            *yMarkers.offset(cidx as isize) = QDLDL_UNUSED as QDLDL_bool;
            i -= 1;
        }
        if *D.offset(k as isize) == 0.0f64 {
            return -(1 as ::std::os::raw::c_int);
        }
        if *D.offset(k as isize) > 0.0f64 {
            positiveValuesInD += 1;
        }
        *Dinv
            .offset(
                k as isize,
            ) = 1 as ::std::os::raw::c_int as ::std::os::raw::c_double / *D.offset(k as isize);
        k += 1;
    }
    return positiveValuesInD;
}
#[no_mangle]
pub unsafe extern "C" fn QDLDL_Lsolve(
    n: QDLDL_int,
    mut Lp: *const QDLDL_int,
    mut Li: *const QDLDL_int,
    mut Lx: *const QDLDL_float,
    mut x: *mut QDLDL_float,
) {
    let mut i: QDLDL_int = 0;
    let mut j: QDLDL_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        let mut val: QDLDL_float = *x.offset(i as isize);
        j = *Lp.offset(i as isize);
        while j < *Lp.offset((i + 1 as ::std::os::raw::c_int) as isize) {
            let ref mut fresh5 = *x.offset(*Li.offset(j as isize) as isize);
            *fresh5 -= *Lx.offset(j as isize) * val;
            j += 1;
        }
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn QDLDL_Ltsolve(
    n: QDLDL_int,
    mut Lp: *const QDLDL_int,
    mut Li: *const QDLDL_int,
    mut Lx: *const QDLDL_float,
    mut x: *mut QDLDL_float,
) {
    let mut i: QDLDL_int = 0;
    let mut j: QDLDL_int = 0;
    i = n - 1 as ::std::os::raw::c_int;
    while i >= 0 as ::std::os::raw::c_int {
        let mut val: QDLDL_float = *x.offset(i as isize);
        j = *Lp.offset(i as isize);
        while j < *Lp.offset((i + 1 as ::std::os::raw::c_int) as isize) {
            val -= *Lx.offset(j as isize) * *x.offset(*Li.offset(j as isize) as isize);
            j += 1;
        }
        *x.offset(i as isize) = val;
        i -= 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn QDLDL_solve(
    n: QDLDL_int,
    mut Lp: *const QDLDL_int,
    mut Li: *const QDLDL_int,
    mut Lx: *const QDLDL_float,
    mut Dinv: *const QDLDL_float,
    mut x: *mut QDLDL_float,
) {
    let mut i: QDLDL_int = 0;
    QDLDL_Lsolve(n, Lp, Li, Lx, x);
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        let ref mut fresh6 = *x.offset(i as isize);
        *fresh6 *= *Dinv.offset(i as isize);
        i += 1;
    }
    QDLDL_Ltsolve(n, Lp, Li, Lx, x);
}
