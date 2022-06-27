pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
#[no_mangle]
pub unsafe extern "C" fn amd_preprocess(
    mut n: ::std::os::raw::c_int,
    mut Ap: *const ::std::os::raw::c_int,
    mut Ai: *const ::std::os::raw::c_int,
    mut Rp: *mut ::std::os::raw::c_int,
    mut Ri: *mut ::std::os::raw::c_int,
    mut W: *mut ::std::os::raw::c_int,
    mut Flag: *mut ::std::os::raw::c_int,
) {
    let mut i: ::std::os::raw::c_int = 0;
    let mut j: ::std::os::raw::c_int = 0;
    let mut p: ::std::os::raw::c_int = 0;
    let mut p2: ::std::os::raw::c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *W.offset(i as isize) = 0 as ::std::os::raw::c_int;
        *Flag.offset(i as isize) = EMPTY;
        i += 1;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < n {
        p2 = *Ap.offset((j + 1 as ::std::os::raw::c_int) as isize);
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
    *Rp.offset(0 as ::std::os::raw::c_int as isize) = 0 as ::std::os::raw::c_int;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *Rp
            .offset(
                (i + 1 as ::std::os::raw::c_int) as isize,
            ) = *Rp.offset(i as isize) + *W.offset(i as isize);
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *W.offset(i as isize) = *Rp.offset(i as isize);
        *Flag.offset(i as isize) = EMPTY;
        i += 1;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < n {
        p2 = *Ap.offset((j + 1 as ::std::os::raw::c_int) as isize);
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
