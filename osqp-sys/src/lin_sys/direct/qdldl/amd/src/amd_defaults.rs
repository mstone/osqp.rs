use ::libc;
pub type c_float = libc::c_double;
pub const NULL: libc::c_int = 0 as libc::c_int;
pub const AMD_CONTROL: libc::c_int = 5 as libc::c_int;
pub const AMD_DENSE: libc::c_int = 0 as libc::c_int;
pub const AMD_DEFAULT_DENSE: libc::c_double = 10.0f64;
pub const AMD_AGGRESSIVE: libc::c_int = 1 as libc::c_int;
pub const AMD_DEFAULT_AGGRESSIVE: libc::c_int = 1 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_l_defaults(mut Control: *mut c_float) {
    let mut i: libc::c_longlong = 0;
    if !Control.is_null() {
        i = 0 as libc::c_int as libc::c_longlong;
        while i < AMD_CONTROL as libc::c_longlong {
            *Control.offset(i as isize) = 0 as libc::c_int as c_float;
            i += 1;
        }
        *Control.offset(AMD_DENSE as isize) = AMD_DEFAULT_DENSE;
        *Control.offset(AMD_AGGRESSIVE as isize) = AMD_DEFAULT_AGGRESSIVE as c_float;
    }
}
