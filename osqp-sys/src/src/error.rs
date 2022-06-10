use ::libc;
extern "C" {
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
}
pub type c_int = libc::c_longlong;
pub type osqp_error_type = libc::c_uint;
pub const OSQP_WORKSPACE_NOT_INIT_ERROR: osqp_error_type = 7;
pub const OSQP_MEM_ALLOC_ERROR: osqp_error_type = 6;
pub const OSQP_NONCVX_ERROR: osqp_error_type = 5;
pub const OSQP_LINSYS_SOLVER_INIT_ERROR: osqp_error_type = 4;
pub const OSQP_LINSYS_SOLVER_LOAD_ERROR: osqp_error_type = 3;
pub const OSQP_SETTINGS_VALIDATION_ERROR: osqp_error_type = 2;
pub const OSQP_DATA_VALIDATION_ERROR: osqp_error_type = 1;
pub const c_print: unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int = printf;
#[no_mangle]
pub static mut OSQP_ERROR_MESSAGE: [*const libc::c_char; 7] = [
    b"Problem data validation.\0" as *const u8 as *const libc::c_char,
    b"Solver settings validation.\0" as *const u8 as *const libc::c_char,
    b"Linear system solver not available.\nTried to obtain it from shared library.\0"
        as *const u8 as *const libc::c_char,
    b"Linear system solver initialization.\0" as *const u8 as *const libc::c_char,
    b"KKT matrix factorization.\nThe problem seems to be non-convex.\0" as *const u8
        as *const libc::c_char,
    b"Memory allocation.\0" as *const u8 as *const libc::c_char,
    b"Solver workspace not initialized.\0" as *const u8 as *const libc::c_char,
];
#[no_mangle]
pub unsafe extern "C" fn _osqp_error(
    mut error_code: osqp_error_type,
    mut function_name: *const libc::c_char,
) -> c_int {
    printf(
        b"ERROR in %s: %s\n\0" as *const u8 as *const libc::c_char,
        function_name,
        OSQP_ERROR_MESSAGE[(error_code as libc::c_uint)
            .wrapping_sub(1 as libc::c_int as libc::c_uint) as usize],
    );
    return error_code as c_int;
}
