pub type c_float = ::std::os::raw::c_double;
pub const NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_CONTROL: ::std::os::raw::c_int = 5 as ::std::os::raw::c_int;
pub const AMD_DENSE: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_DEFAULT_DENSE: ::std::os::raw::c_double = 10.0f64;
pub const AMD_AGGRESSIVE: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const AMD_DEFAULT_AGGRESSIVE: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
#[no_mangle]
pub unsafe extern "C" fn amd_defaults(mut Control: *mut c_float) {
    let mut i: ::std::os::raw::c_int = 0;
    if !Control.is_null() {
        i = 0 as ::std::os::raw::c_int;
        while i < AMD_CONTROL {
            *Control.offset(i as isize) = 0 as ::std::os::raw::c_int as c_float;
            i += 1;
        }
        *Control.offset(AMD_DENSE as isize) = AMD_DEFAULT_DENSE;
        *Control.offset(AMD_AGGRESSIVE as isize) = AMD_DEFAULT_AGGRESSIVE as c_float;
    }
}
