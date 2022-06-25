use ::libc;
pub const AMD_OK: libc::c_int = 0 as libc::c_int;
pub const AMD_OK_BUT_JUMBLED: libc::c_int = 1 as libc::c_int;
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
pub const NULL: libc::c_int = 0 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_valid(
    mut n_row: libc::c_longlong,
    mut n_col: libc::c_longlong,
    mut Ap: *const libc::c_longlong,
    mut Ai: *const libc::c_longlong,
) -> libc::c_longlong {
    let mut nz: libc::c_longlong = 0;
    let mut j: libc::c_longlong = 0;
    let mut p1: libc::c_longlong = 0;
    let mut p2: libc::c_longlong = 0;
    let mut ilast: libc::c_longlong = 0;
    let mut i: libc::c_longlong = 0;
    let mut p: libc::c_longlong = 0;
    let mut result: libc::c_longlong = AMD_OK as libc::c_longlong;
    if n_row < 0 as libc::c_int as libc::c_longlong
        || n_col < 0 as libc::c_int as libc::c_longlong || Ap.is_null() || Ai.is_null()
    {
        return -(2 as libc::c_int) as libc::c_longlong;
    }
    nz = *Ap.offset(n_col as isize);
    if *Ap.offset(0 as libc::c_int as isize) != 0 as libc::c_int as libc::c_longlong
        || nz < 0 as libc::c_int as libc::c_longlong
    {
        return -(2 as libc::c_int) as libc::c_longlong;
    }
    j = 0 as libc::c_int as libc::c_longlong;
    while j < n_col {
        p1 = *Ap.offset(j as isize);
        p2 = *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize);
        if p1 > p2 {
            return -(2 as libc::c_int) as libc::c_longlong;
        }
        ilast = EMPTY as libc::c_longlong;
        p = p1;
        while p < p2 {
            i = *Ai.offset(p as isize);
            if i < 0 as libc::c_int as libc::c_longlong || i >= n_row {
                return -(2 as libc::c_int) as libc::c_longlong;
            }
            if i <= ilast {
                result = AMD_OK_BUT_JUMBLED as libc::c_longlong;
            }
            ilast = i;
            p += 1;
        }
        j += 1;
    }
    return result;
}
