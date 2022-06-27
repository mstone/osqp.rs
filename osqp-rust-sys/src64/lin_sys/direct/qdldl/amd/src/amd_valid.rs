pub const AMD_OK: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_OK_BUT_JUMBLED: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_valid(
    mut n_row: ::std::os::raw::c_longlong,
    mut n_col: ::std::os::raw::c_longlong,
    mut Ap: *const ::std::os::raw::c_longlong,
    mut Ai: *const ::std::os::raw::c_longlong,
) -> ::std::os::raw::c_longlong {
    let mut nz: ::std::os::raw::c_longlong = 0;
    let mut j: ::std::os::raw::c_longlong = 0;
    let mut p1: ::std::os::raw::c_longlong = 0;
    let mut p2: ::std::os::raw::c_longlong = 0;
    let mut ilast: ::std::os::raw::c_longlong = 0;
    let mut i: ::std::os::raw::c_longlong = 0;
    let mut p: ::std::os::raw::c_longlong = 0;
    let mut result: ::std::os::raw::c_longlong = AMD_OK as ::std::os::raw::c_longlong;
    if n_row < 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
        || n_col < 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong || Ap.is_null() || Ai.is_null()
    {
        return -(2 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong;
    }
    nz = *Ap.offset(n_col as isize);
    if *Ap.offset(0 as ::std::os::raw::c_int as isize) != 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
        || nz < 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
    {
        return -(2 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong;
    }
    j = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    while j < n_col {
        p1 = *Ap.offset(j as isize);
        p2 = *Ap.offset((j + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong) as isize);
        if p1 > p2 {
            return -(2 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong;
        }
        ilast = EMPTY as ::std::os::raw::c_longlong;
        p = p1;
        while p < p2 {
            i = *Ai.offset(p as isize);
            if i < 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong || i >= n_row {
                return -(2 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong;
            }
            if i <= ilast {
                result = AMD_OK_BUT_JUMBLED as ::std::os::raw::c_longlong;
            }
            ilast = i;
            p += 1;
        }
        j += 1;
    }
    return result;
}
