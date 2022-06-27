extern "C" {
    pub type OSQP_TIMER;
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
#[no_mangle]
pub unsafe extern "C" fn project(mut work: *mut OSQPWorkspace, mut z: *mut c_float) {
    let mut i: c_int = 0;
    let mut m: c_int = 0;
    m = (*(*work).data).m;
    i = 0 as ::std::os::raw::c_int;
    while i < m {
        *z
            .offset(
                i as isize,
            ) = if (if *z.offset(i as isize) > *((*(*work).data).l).offset(i as isize) {
            *z.offset(i as isize)
        } else {
            *((*(*work).data).l).offset(i as isize)
        }) < *((*(*work).data).u).offset(i as isize)
        {
            if *z.offset(i as isize) > *((*(*work).data).l).offset(i as isize) {
                *z.offset(i as isize)
            } else {
                *((*(*work).data).l).offset(i as isize)
            }
        } else {
            *((*(*work).data).u).offset(i as isize)
        };
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn project_normalcone(
    mut work: *mut OSQPWorkspace,
    mut z: *mut c_float,
    mut y: *mut c_float,
) {
    let mut i: c_int = 0;
    let mut m: c_int = 0;
    m = (*(*work).data).m;
    i = 0 as ::std::os::raw::c_int;
    while i < m {
        *((*work).z_prev)
            .offset(i as isize) = *z.offset(i as isize) + *y.offset(i as isize);
        *z
            .offset(
                i as isize,
            ) = if (if *((*work).z_prev).offset(i as isize)
            > *((*(*work).data).l).offset(i as isize)
        {
            *((*work).z_prev).offset(i as isize)
        } else {
            *((*(*work).data).l).offset(i as isize)
        }) < *((*(*work).data).u).offset(i as isize)
        {
            if *((*work).z_prev).offset(i as isize)
                > *((*(*work).data).l).offset(i as isize)
            {
                *((*work).z_prev).offset(i as isize)
            } else {
                *((*(*work).data).l).offset(i as isize)
            }
        } else {
            *((*(*work).data).u).offset(i as isize)
        };
        *y
            .offset(
                i as isize,
            ) = *((*work).z_prev).offset(i as isize) - *z.offset(i as isize);
        i += 1;
    }
}
