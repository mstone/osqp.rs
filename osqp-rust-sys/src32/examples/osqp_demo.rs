#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![register_tool(c2rust)]
#![feature(register_tool)]
use ::osqp_rust_sys::*;
extern "C" {
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn free(_: *mut ::std::os::raw::c_void);
    fn csc_matrix(
        m: c_int,
        n: c_int,
        nzmax: c_int,
        x: *mut c_float,
        i: *mut c_int,
        p: *mut c_int,
    ) -> *mut csc;
    fn osqp_set_default_settings(settings: *mut OSQPSettings);
    fn osqp_setup(
        workp: *mut *mut OSQPWorkspace,
        data: *const OSQPData,
        settings: *const OSQPSettings,
    ) -> c_int;
    fn osqp_solve(work: *mut OSQPWorkspace) -> c_int;
    fn osqp_cleanup(work: *mut OSQPWorkspace) -> c_int;
}
pub type c_int = ::std::os::raw::c_int;
pub type c_float = ::std::os::raw::c_double;
pub type linsys_solver_type = ::std::os::raw::c_uint;
pub const MKL_PARDISO_SOLVER: linsys_solver_type = 1;
pub const QDLDL_SOLVER: linsys_solver_type = 0;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct linsys_solver {
    pub type_0: linsys_solver_type,
    pub solve: Option::<unsafe extern "C" fn(*mut LinSysSolver, *mut c_float) -> c_int>,
    pub free: Option::<unsafe extern "C" fn(*mut LinSysSolver) -> ()>,
    pub update_matrices: Option::<
        unsafe extern "C" fn(*mut LinSysSolver, *const csc, *const csc) -> c_int,
    >,
    pub update_rho_vec: Option::<
        unsafe extern "C" fn(*mut LinSysSolver, *const c_float) -> c_int,
    >,
    pub nthreads: c_int,
}
pub type LinSysSolver = linsys_solver;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQP_TIMER {
    pub tic: uint64_t,
    pub toc: uint64_t,
    pub tinfo: mach_timebase_info_data_t,
}
pub type mach_timebase_info_data_t = mach_timebase_info;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct mach_timebase_info {
    pub numer: uint32_t,
    pub denom: uint32_t,
}
pub type uint32_t = ::std::os::raw::c_uint;
pub type uint64_t = ::std::os::raw::c_ulonglong;
pub type OSQPTimer = OSQP_TIMER;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPScaling {
    pub c: c_float,
    pub D: *mut c_float,
    pub E: *mut c_float,
    pub cinv: c_float,
    pub Dinv: *mut c_float,
    pub Einv: *mut c_float,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPSolution {
    pub x: *mut c_float,
    pub y: *mut c_float,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPInfo {
    pub iter: c_int,
    pub status: [::std::os::raw::c_char; 32],
    pub status_val: c_int,
    pub status_polish: c_int,
    pub obj_val: c_float,
    pub pri_res: c_float,
    pub dua_res: c_float,
    pub setup_time: c_float,
    pub solve_time: c_float,
    pub update_time: c_float,
    pub polish_time: c_float,
    pub run_time: c_float,
    pub rho_updates: c_int,
    pub rho_estimate: c_float,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPPolish {
    pub Ared: *mut csc,
    pub n_low: c_int,
    pub n_upp: c_int,
    pub A_to_Alow: *mut c_int,
    pub A_to_Aupp: *mut c_int,
    pub Alow_to_A: *mut c_int,
    pub Aupp_to_A: *mut c_int,
    pub x: *mut c_float,
    pub z: *mut c_float,
    pub y: *mut c_float,
    pub obj_val: c_float,
    pub pri_res: c_float,
    pub dua_res: c_float,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPData {
    pub n: c_int,
    pub m: c_int,
    pub P: *mut csc,
    pub A: *mut csc,
    pub q: *mut c_float,
    pub l: *mut c_float,
    pub u: *mut c_float,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPSettings {
    pub rho: c_float,
    pub sigma: c_float,
    pub scaling: c_int,
    pub adaptive_rho: c_int,
    pub adaptive_rho_interval: c_int,
    pub adaptive_rho_tolerance: c_float,
    pub adaptive_rho_fraction: c_float,
    pub max_iter: c_int,
    pub eps_abs: c_float,
    pub eps_rel: c_float,
    pub eps_prim_inf: c_float,
    pub eps_dual_inf: c_float,
    pub alpha: c_float,
    pub linsys_solver: linsys_solver_type,
    pub delta: c_float,
    pub polish: c_int,
    pub polish_refine_iter: c_int,
    pub verbose: c_int,
    pub scaled_termination: c_int,
    pub check_termination: c_int,
    pub warm_start: c_int,
    pub time_limit: c_float,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct OSQPWorkspace {
    pub data: *mut OSQPData,
    pub linsys_solver: *mut LinSysSolver,
    pub pol: *mut OSQPPolish,
    pub rho_vec: *mut c_float,
    pub rho_inv_vec: *mut c_float,
    pub constr_type: *mut c_int,
    pub x: *mut c_float,
    pub y: *mut c_float,
    pub z: *mut c_float,
    pub xz_tilde: *mut c_float,
    pub x_prev: *mut c_float,
    pub z_prev: *mut c_float,
    pub Ax: *mut c_float,
    pub Px: *mut c_float,
    pub Aty: *mut c_float,
    pub delta_y: *mut c_float,
    pub Atdelta_y: *mut c_float,
    pub delta_x: *mut c_float,
    pub Pdelta_x: *mut c_float,
    pub Adelta_x: *mut c_float,
    pub D_temp: *mut c_float,
    pub D_temp_A: *mut c_float,
    pub E_temp: *mut c_float,
    pub settings: *mut OSQPSettings,
    pub scaling: *mut OSQPScaling,
    pub solution: *mut OSQPSolution,
    pub info: *mut OSQPInfo,
    pub timer: *mut OSQPTimer,
    pub first_run: c_int,
    pub clear_update_time: c_int,
    pub rho_update_from_solve: c_int,
    pub summary_printed: c_int,
}
pub const c_free: unsafe extern "C" fn(*mut ::std::os::raw::c_void) -> () = free;
pub const c_malloc: unsafe extern "C" fn(::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void = malloc;
unsafe fn main_0(
    mut argc: ::std::os::raw::c_int,
    mut argv: *mut *mut ::std::os::raw::c_char,
) -> ::std::os::raw::c_int {
    let mut P_x: [c_float; 3] = [4.0f64, 1.0f64, 2.0f64];
    let mut P_nnz: c_int = 3 as ::std::os::raw::c_int;
    let mut P_i: [c_int; 3] = [0 as ::std::os::raw::c_int, 0 as ::std::os::raw::c_int, 1 as ::std::os::raw::c_int];
    let mut P_p: [c_int; 3] = [0 as ::std::os::raw::c_int, 1 as ::std::os::raw::c_int, 3 as ::std::os::raw::c_int];
    let mut q: [c_float; 2] = [1.0f64, 1.0f64];
    let mut A_x: [c_float; 4] = [1.0f64, 1.0f64, 1.0f64, 1.0f64];
    let mut A_nnz: c_int = 4 as ::std::os::raw::c_int;
    let mut A_i: [c_int; 4] = [
        0 as ::std::os::raw::c_int,
        1 as ::std::os::raw::c_int,
        0 as ::std::os::raw::c_int,
        2 as ::std::os::raw::c_int,
    ];
    let mut A_p: [c_int; 3] = [0 as ::std::os::raw::c_int, 2 as ::std::os::raw::c_int, 4 as ::std::os::raw::c_int];
    let mut l: [c_float; 3] = [1.0f64, 0.0f64, 0.0f64];
    let mut u: [c_float; 3] = [1.0f64, 0.7f64, 0.7f64];
    let mut n: c_int = 2 as ::std::os::raw::c_int;
    let mut m: c_int = 3 as ::std::os::raw::c_int;
    let mut exitflag: c_int = 0 as ::std::os::raw::c_int;
    let mut work: *mut OSQPWorkspace = 0 as *mut OSQPWorkspace;
    let mut settings: *mut OSQPSettings = malloc(
        ::std::mem::size_of::<OSQPSettings>() as ::std::os::raw::c_ulong,
    ) as *mut OSQPSettings;
    let mut data: *mut OSQPData = malloc(
        ::std::mem::size_of::<OSQPData>() as ::std::os::raw::c_ulong,
    ) as *mut OSQPData;
    if !data.is_null() {
        (*data).n = n;
        (*data).m = m;
        let ref mut fresh0 = (*data).P;
        *fresh0 = csc_matrix(
            (*data).n,
            (*data).n,
            P_nnz,
            P_x.as_mut_ptr(),
            P_i.as_mut_ptr(),
            P_p.as_mut_ptr(),
        );
        let ref mut fresh1 = (*data).q;
        *fresh1 = q.as_mut_ptr();
        let ref mut fresh2 = (*data).A;
        *fresh2 = csc_matrix(
            (*data).m,
            (*data).n,
            A_nnz,
            A_x.as_mut_ptr(),
            A_i.as_mut_ptr(),
            A_p.as_mut_ptr(),
        );
        let ref mut fresh3 = (*data).l;
        *fresh3 = l.as_mut_ptr();
        let ref mut fresh4 = (*data).u;
        *fresh4 = u.as_mut_ptr();
    }
    if !settings.is_null() {
        osqp_set_default_settings(settings);
    }
    exitflag = osqp_setup(&mut work, data, settings);
    osqp_solve(work);
    osqp_cleanup(work);
    if !data.is_null() {
        if !((*data).A).is_null() {
            free((*data).A as *mut ::std::os::raw::c_void);
        }
        if !((*data).P).is_null() {
            free((*data).P as *mut ::std::os::raw::c_void);
        }
        free(data as *mut ::std::os::raw::c_void);
    }
    if !settings.is_null() {
        free(settings as *mut ::std::os::raw::c_void);
    }
    return exitflag;
}
pub fn main() {
    let mut args: Vec::<*mut ::std::os::raw::c_char> = Vec::new();
    for arg in ::std::env::args() {
        args.push(
            (::std::ffi::CString::new(arg))
                .expect("Failed to convert argument into CString.")
                .into_raw(),
        );
    }
    args.push(::std::ptr::null_mut());
    unsafe {
        ::std::process::exit(
            main_0(
                (args.len() - 1) as ::std::os::raw::c_int,
                args.as_mut_ptr() as *mut *mut ::std::os::raw::c_char,
            ) as i32,
        )
    }
}
