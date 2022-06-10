use ::libc;
extern "C" {
    pub type OSQP_TIMER;
    fn vec_set_scalar(a: *mut c_float, sc: c_float, n: c_int);
    fn vec_mult_scalar(a: *mut c_float, sc: c_float, n: c_int);
    fn vec_norm_inf(v: *const c_float, l: c_int) -> c_float;
    fn vec_mean(a: *const c_float, n: c_int) -> c_float;
    fn vec_ew_recipr(a: *const c_float, b: *mut c_float, n: c_int);
    fn vec_ew_prod(a: *const c_float, b: *const c_float, c: *mut c_float, n: c_int);
    fn vec_ew_sqrt(a: *mut c_float, n: c_int);
    fn vec_ew_max_vec(a: *const c_float, b: *const c_float, c: *mut c_float, n: c_int);
    fn mat_mult_scalar(A: *mut csc, sc: c_float);
    fn mat_premult_diag(A: *mut csc, d: *const c_float);
    fn mat_postmult_diag(A: *mut csc, d: *const c_float);
    fn mat_inf_norm_cols(M: *const csc, E: *mut c_float);
    fn mat_inf_norm_rows(M: *const csc, E: *mut c_float);
    fn mat_inf_norm_cols_sym_triu(M: *const csc, E: *mut c_float);
}
pub type c_int = libc::c_longlong;
pub type c_float = libc::c_double;
pub type linsys_solver_type = libc::c_uint;
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
    pub status: [libc::c_char; 32],
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
pub const MIN_SCALING: libc::c_double = 1e-04f64;
pub const MAX_SCALING: libc::c_double = 1e+04f64;
#[no_mangle]
pub unsafe extern "C" fn limit_scaling(mut D: *mut c_float, mut n: c_int) {
    let mut i: c_int = 0;
    i = 0 as libc::c_int as c_int;
    while i < n {
        *D
            .offset(
                i as isize,
            ) = if *D.offset(i as isize) < MIN_SCALING {
            1.0f64
        } else {
            *D.offset(i as isize)
        };
        *D
            .offset(
                i as isize,
            ) = if *D.offset(i as isize) > MAX_SCALING {
            MAX_SCALING
        } else {
            *D.offset(i as isize)
        };
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn compute_inf_norm_cols_KKT(
    mut P: *const csc,
    mut A: *const csc,
    mut D: *mut c_float,
    mut D_temp_A: *mut c_float,
    mut E: *mut c_float,
    mut n: c_int,
) {
    mat_inf_norm_cols_sym_triu(P, D);
    mat_inf_norm_cols(A, D_temp_A);
    vec_ew_max_vec(D, D_temp_A, D, n);
    mat_inf_norm_rows(A, E);
}
#[no_mangle]
pub unsafe extern "C" fn scale_data(mut work: *mut OSQPWorkspace) -> c_int {
    let mut i: c_int = 0;
    let mut n: c_int = 0;
    let mut m: c_int = 0;
    let mut c_temp: c_float = 0.;
    let mut inf_norm_q: c_float = 0.;
    n = (*(*work).data).n;
    m = (*(*work).data).m;
    (*(*work).scaling).c = 1.0f64;
    vec_set_scalar((*(*work).scaling).D, 1.0f64, (*(*work).data).n);
    vec_set_scalar((*(*work).scaling).Dinv, 1.0f64, (*(*work).data).n);
    vec_set_scalar((*(*work).scaling).E, 1.0f64, (*(*work).data).m);
    vec_set_scalar((*(*work).scaling).Einv, 1.0f64, (*(*work).data).m);
    i = 0 as libc::c_int as c_int;
    while i < (*(*work).settings).scaling {
        compute_inf_norm_cols_KKT(
            (*(*work).data).P,
            (*(*work).data).A,
            (*work).D_temp,
            (*work).D_temp_A,
            (*work).E_temp,
            n,
        );
        limit_scaling((*work).D_temp, n);
        limit_scaling((*work).E_temp, m);
        vec_ew_sqrt((*work).D_temp, n);
        vec_ew_sqrt((*work).E_temp, m);
        vec_ew_recipr((*work).D_temp, (*work).D_temp, n);
        vec_ew_recipr((*work).E_temp, (*work).E_temp, m);
        mat_premult_diag((*(*work).data).P, (*work).D_temp);
        mat_postmult_diag((*(*work).data).P, (*work).D_temp);
        mat_premult_diag((*(*work).data).A, (*work).E_temp);
        mat_postmult_diag((*(*work).data).A, (*work).D_temp);
        vec_ew_prod((*work).D_temp, (*(*work).data).q, (*(*work).data).q, n);
        vec_ew_prod((*(*work).scaling).D, (*work).D_temp, (*(*work).scaling).D, n);
        vec_ew_prod((*(*work).scaling).E, (*work).E_temp, (*(*work).scaling).E, m);
        mat_inf_norm_cols_sym_triu((*(*work).data).P, (*work).D_temp);
        c_temp = vec_mean((*work).D_temp, n);
        inf_norm_q = vec_norm_inf((*(*work).data).q, n);
        limit_scaling(&mut inf_norm_q, 1 as libc::c_int as c_int);
        c_temp = if c_temp > inf_norm_q { c_temp } else { inf_norm_q };
        limit_scaling(&mut c_temp, 1 as libc::c_int as c_int);
        c_temp = 1.0f64 / c_temp;
        mat_mult_scalar((*(*work).data).P, c_temp);
        vec_mult_scalar((*(*work).data).q, c_temp, n);
        let ref mut fresh0 = (*(*work).scaling).c;
        *fresh0 *= c_temp;
        i += 1;
    }
    (*(*work).scaling).cinv = 1.0f64 / (*(*work).scaling).c;
    vec_ew_recipr((*(*work).scaling).D, (*(*work).scaling).Dinv, (*(*work).data).n);
    vec_ew_recipr((*(*work).scaling).E, (*(*work).scaling).Einv, (*(*work).data).m);
    vec_ew_prod(
        (*(*work).scaling).E,
        (*(*work).data).l,
        (*(*work).data).l,
        (*(*work).data).m,
    );
    vec_ew_prod(
        (*(*work).scaling).E,
        (*(*work).data).u,
        (*(*work).data).u,
        (*(*work).data).m,
    );
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn unscale_data(mut work: *mut OSQPWorkspace) -> c_int {
    mat_mult_scalar((*(*work).data).P, (*(*work).scaling).cinv);
    mat_premult_diag((*(*work).data).P, (*(*work).scaling).Dinv);
    mat_postmult_diag((*(*work).data).P, (*(*work).scaling).Dinv);
    vec_mult_scalar((*(*work).data).q, (*(*work).scaling).cinv, (*(*work).data).n);
    vec_ew_prod(
        (*(*work).scaling).Dinv,
        (*(*work).data).q,
        (*(*work).data).q,
        (*(*work).data).n,
    );
    mat_premult_diag((*(*work).data).A, (*(*work).scaling).Einv);
    mat_postmult_diag((*(*work).data).A, (*(*work).scaling).Dinv);
    vec_ew_prod(
        (*(*work).scaling).Einv,
        (*(*work).data).l,
        (*(*work).data).l,
        (*(*work).data).m,
    );
    vec_ew_prod(
        (*(*work).scaling).Einv,
        (*(*work).data).u,
        (*(*work).data).u,
        (*(*work).data).m,
    );
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn unscale_solution(mut work: *mut OSQPWorkspace) -> c_int {
    vec_ew_prod(
        (*(*work).scaling).D,
        (*(*work).solution).x,
        (*(*work).solution).x,
        (*(*work).data).n,
    );
    vec_ew_prod(
        (*(*work).scaling).E,
        (*(*work).solution).y,
        (*(*work).solution).y,
        (*(*work).data).m,
    );
    vec_mult_scalar((*(*work).solution).y, (*(*work).scaling).cinv, (*(*work).data).m);
    return 0 as libc::c_int as c_int;
}
