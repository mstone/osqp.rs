use ::libc;
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
#[no_mangle]
pub unsafe extern "C" fn amd_l_preprocess(
    mut n: libc::c_longlong,
    mut Ap: *const libc::c_longlong,
    mut Ai: *const libc::c_longlong,
    mut Rp: *mut libc::c_longlong,
    mut Ri: *mut libc::c_longlong,
    mut W: *mut libc::c_longlong,
    mut Flag: *mut libc::c_longlong,
) {
    let mut i: libc::c_longlong = 0;
    let mut j: libc::c_longlong = 0;
    let mut p: libc::c_longlong = 0;
    let mut p2: libc::c_longlong = 0;
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        *W.offset(i as isize) = 0 as libc::c_int as libc::c_longlong;
        *Flag.offset(i as isize) = EMPTY as libc::c_longlong;
        i += 1;
    }
    j = 0 as libc::c_int as libc::c_longlong;
    while j < n {
        p2 = *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize);
        p = *Ap.offset(j as isize);
        while p < p2 {
            i = *Ai.offset(p as isize);
            if *Flag.offset(i as isize) != j {
                let ref mut fresh0 = *W.offset(i as isize);
                *fresh0 += 1;
                *Flag.offset(i as isize) = j;
            }
            p += 1;
        }
        j += 1;
    }
    *Rp.offset(0 as libc::c_int as isize) = 0 as libc::c_int as libc::c_longlong;
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        *Rp
            .offset(
                (i + 1 as libc::c_int as libc::c_longlong) as isize,
            ) = *Rp.offset(i as isize) + *W.offset(i as isize);
        i += 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        *W.offset(i as isize) = *Rp.offset(i as isize);
        *Flag.offset(i as isize) = EMPTY as libc::c_longlong;
        i += 1;
    }
    j = 0 as libc::c_int as libc::c_longlong;
    while j < n {
        p2 = *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize);
        p = *Ap.offset(j as isize);
        while p < p2 {
            i = *Ai.offset(p as isize);
            if *Flag.offset(i as isize) != j {
                let ref mut fresh1 = *W.offset(i as isize);
                let fresh2 = *fresh1;
                *fresh1 = *fresh1 + 1;
                *Ri.offset(fresh2 as isize) = j;
                *Flag.offset(i as isize) = j;
            }
            p += 1;
        }
        j += 1;
    }
}
