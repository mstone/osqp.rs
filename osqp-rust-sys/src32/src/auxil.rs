extern "C" {
    fn c_strcpy(dest: *mut ::std::os::raw::c_char, source: *const ::std::os::raw::c_char);
    fn sqrt(_: ::std::os::raw::c_double) -> ::std::os::raw::c_double;
    fn printf(_: *const ::std::os::raw::c_char, _: ...) -> ::std::os::raw::c_int;
    fn osqp_toc(t: *mut OSQPTimer) -> c_float;
    fn prea_vec_copy(a: *const c_float, b: *mut c_float, n: c_int);
    fn vec_set_scalar(a: *mut c_float, sc: c_float, n: c_int);
    fn vec_mult_scalar(a: *mut c_float, sc: c_float, n: c_int);
    fn vec_add_scaled(
        c: *mut c_float,
        a: *const c_float,
        b: *const c_float,
        n: c_int,
        sc: c_float,
    );
    fn vec_norm_inf(v: *const c_float, l: c_int) -> c_float;
    fn vec_scaled_norm_inf(S: *const c_float, v: *const c_float, l: c_int) -> c_float;
    fn vec_prod(a: *const c_float, b: *const c_float, n: c_int) -> c_float;
    fn vec_ew_prod(a: *const c_float, b: *const c_float, c: *mut c_float, n: c_int);
    fn mat_vec(A: *const csc, x: *const c_float, y: *mut c_float, plus_eq: c_int);
    fn mat_tpose_vec(
        A: *const csc,
        x: *const c_float,
        y: *mut c_float,
        plus_eq: c_int,
        skip_diag: c_int,
    );
    fn quad_form(P: *const csc, x: *const c_float) -> c_float;
    fn osqp_update_rho(work: *mut OSQPWorkspace, rho_new: c_float) -> c_int;
    fn project(work: *mut OSQPWorkspace, z: *mut c_float);
    fn unscale_solution(work: *mut OSQPWorkspace) -> c_int;
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
pub const c_print: unsafe extern "C" fn(*const ::std::os::raw::c_char, ...) -> ::std::os::raw::c_int = printf;
pub const OSQP_SOLVED: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const OSQP_SOLVED_INACCURATE: ::std::os::raw::c_int = 2 as ::std::os::raw::c_int;
pub const OSQP_MAX_ITER_REACHED: ::std::os::raw::c_int = -(2 as ::std::os::raw::c_int);
pub const OSQP_TIME_LIMIT_REACHED: ::std::os::raw::c_int = -(6 as ::std::os::raw::c_int);
pub const OSQP_SIGINT: ::std::os::raw::c_int = -(5 as ::std::os::raw::c_int);
pub const OSQP_UNSOLVED: ::std::os::raw::c_int = -(10 as ::std::os::raw::c_int);
pub const OSQP_NAN: ::std::os::raw::c_ulong = 0x7fc00000 as ::std::os::raw::c_ulong;
pub const OSQP_PRIMAL_INFEASIBLE: ::std::os::raw::c_int = -(3 as ::std::os::raw::c_int);
pub const OSQP_PRIMAL_INFEASIBLE_INACCURATE: ::std::os::raw::c_int = 3 as ::std::os::raw::c_int;
pub const OSQP_DUAL_INFEASIBLE: ::std::os::raw::c_int = -(4 as ::std::os::raw::c_int);
pub const OSQP_DUAL_INFEASIBLE_INACCURATE: ::std::os::raw::c_int = 4 as ::std::os::raw::c_int;
pub const OSQP_NON_CVX: ::std::os::raw::c_int = -(7 as ::std::os::raw::c_int);
pub const c_sqrt: unsafe extern "C" fn(::std::os::raw::c_double) -> ::std::os::raw::c_double = sqrt;
pub const OSQP_INFTY: ::std::os::raw::c_double = 1e30f64;
pub const OSQP_DIVISION_TOL: ::std::os::raw::c_double = 1.0f64 / OSQP_INFTY;
pub const MIN_SCALING: ::std::os::raw::c_double = 1e-04f64;
pub const RHO_EQ_OVER_RHO_INEQ: ::std::os::raw::c_double = 1e03f64;
pub const RHO_TOL: ::std::os::raw::c_double = 1e-04f64;
pub const RHO_MIN: ::std::os::raw::c_double = 1e-06f64;
#[no_mangle]
pub unsafe extern "C" fn compute_rho_estimate(mut work: *mut OSQPWorkspace) -> c_float {
    let mut n: c_int = 0;
    let mut m: c_int = 0;
    let mut pri_res: c_float = 0.;
    let mut dua_res: c_float = 0.;
    let mut pri_res_norm: c_float = 0.;
    let mut dua_res_norm: c_float = 0.;
    let mut temp_res_norm: c_float = 0.;
    let mut rho_estimate: c_float = 0.;
    n = (*(*work).data).n;
    m = (*(*work).data).m;
    pri_res = vec_norm_inf((*work).z_prev, m);
    dua_res = vec_norm_inf((*work).x_prev, n);
    pri_res_norm = vec_norm_inf((*work).z, m);
    temp_res_norm = vec_norm_inf((*work).Ax, m);
    pri_res_norm = if pri_res_norm > temp_res_norm {
        pri_res_norm
    } else {
        temp_res_norm
    };
    pri_res /= pri_res_norm + OSQP_DIVISION_TOL;
    dua_res_norm = vec_norm_inf((*(*work).data).q, n);
    temp_res_norm = vec_norm_inf((*work).Aty, n);
    dua_res_norm = if dua_res_norm > temp_res_norm {
        dua_res_norm
    } else {
        temp_res_norm
    };
    temp_res_norm = vec_norm_inf((*work).Px, n);
    dua_res_norm = if dua_res_norm > temp_res_norm {
        dua_res_norm
    } else {
        temp_res_norm
    };
    dua_res /= dua_res_norm + OSQP_DIVISION_TOL;
    rho_estimate = (*(*work).settings).rho * sqrt(pri_res / dua_res);
    rho_estimate = if (if rho_estimate > 1e-06f64 { rho_estimate } else { 1e-06f64 })
        < 1e06f64
    {
        if rho_estimate > 1e-06f64 { rho_estimate } else { 1e-06f64 }
    } else {
        1e06f64
    };
    return rho_estimate;
}
#[no_mangle]
pub unsafe extern "C" fn adapt_rho(mut work: *mut OSQPWorkspace) -> c_int {
    let mut exitflag: c_int = 0;
    let mut rho_new: c_float = 0.;
    exitflag = 0 as ::std::os::raw::c_int;
    rho_new = compute_rho_estimate(work);
    (*(*work).info).rho_estimate = rho_new;
    if rho_new > (*(*work).settings).rho * (*(*work).settings).adaptive_rho_tolerance
        || rho_new < (*(*work).settings).rho / (*(*work).settings).adaptive_rho_tolerance
    {
        exitflag = osqp_update_rho(work, rho_new);
        let ref mut fresh0 = (*(*work).info).rho_updates;
        *fresh0 += 1 as ::std::os::raw::c_int;
    }
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn set_rho_vec(mut work: *mut OSQPWorkspace) {
    let mut i: c_int = 0;
    (*(*work).settings)
        .rho = if (if (*(*work).settings).rho > 1e-06f64 {
        (*(*work).settings).rho
    } else {
        1e-06f64
    }) < 1e06f64
    {
        if (*(*work).settings).rho > 1e-06f64 {
            (*(*work).settings).rho
        } else {
            1e-06f64
        }
    } else {
        1e06f64
    };
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).m {
        if *((*(*work).data).l).offset(i as isize) < -OSQP_INFTY * MIN_SCALING
            && *((*(*work).data).u).offset(i as isize) > OSQP_INFTY * MIN_SCALING
        {
            *((*work).constr_type).offset(i as isize) = -(1 as ::std::os::raw::c_int);
            *((*work).rho_vec).offset(i as isize) = RHO_MIN;
        } else if *((*(*work).data).u).offset(i as isize)
                - *((*(*work).data).l).offset(i as isize) < RHO_TOL
            {
            *((*work).constr_type).offset(i as isize) = 1 as ::std::os::raw::c_int;
            *((*work).rho_vec)
                .offset(i as isize) = RHO_EQ_OVER_RHO_INEQ * (*(*work).settings).rho;
        } else {
            *((*work).constr_type).offset(i as isize) = 0 as ::std::os::raw::c_int;
            *((*work).rho_vec).offset(i as isize) = (*(*work).settings).rho;
        }
        *((*work).rho_inv_vec)
            .offset(i as isize) = 1.0f64 / *((*work).rho_vec).offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn update_rho_vec(mut work: *mut OSQPWorkspace) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0;
    let mut constr_type_changed: c_int = 0;
    exitflag = 0 as ::std::os::raw::c_int;
    constr_type_changed = 0 as ::std::os::raw::c_int;
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).m {
        if *((*(*work).data).l).offset(i as isize) < -OSQP_INFTY * MIN_SCALING
            && *((*(*work).data).u).offset(i as isize) > OSQP_INFTY * MIN_SCALING
        {
            if *((*work).constr_type).offset(i as isize) != -(1 as ::std::os::raw::c_int) {
                *((*work).constr_type).offset(i as isize) = -(1 as ::std::os::raw::c_int);
                *((*work).rho_vec).offset(i as isize) = RHO_MIN;
                *((*work).rho_inv_vec).offset(i as isize) = 1.0f64 / RHO_MIN;
                constr_type_changed = 1 as ::std::os::raw::c_int;
            }
        } else if *((*(*work).data).u).offset(i as isize)
                - *((*(*work).data).l).offset(i as isize) < RHO_TOL
            {
            if *((*work).constr_type).offset(i as isize) != 1 as ::std::os::raw::c_int {
                *((*work).constr_type).offset(i as isize) = 1 as ::std::os::raw::c_int;
                *((*work).rho_vec)
                    .offset(i as isize) = RHO_EQ_OVER_RHO_INEQ * (*(*work).settings).rho;
                *((*work).rho_inv_vec)
                    .offset(i as isize) = 1.0f64 / *((*work).rho_vec).offset(i as isize);
                constr_type_changed = 1 as ::std::os::raw::c_int;
            }
        } else if *((*work).constr_type).offset(i as isize) != 0 as ::std::os::raw::c_int {
            *((*work).constr_type).offset(i as isize) = 0 as ::std::os::raw::c_int;
            *((*work).rho_vec).offset(i as isize) = (*(*work).settings).rho;
            *((*work).rho_inv_vec).offset(i as isize) = 1.0f64 / (*(*work).settings).rho;
            constr_type_changed = 1 as ::std::os::raw::c_int;
        }
        i += 1;
    }
    if constr_type_changed == 1 as ::std::os::raw::c_int {
        exitflag = ((*(*work).linsys_solver).update_rho_vec)
            .expect("non-null function pointer")((*work).linsys_solver, (*work).rho_vec);
    }
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn swap_vectors(
    mut a: *mut *mut c_float,
    mut b: *mut *mut c_float,
) {
    let mut temp: *mut c_float = 0 as *mut c_float;
    temp = *b;
    *b = *a;
    *a = temp;
}
#[no_mangle]
pub unsafe extern "C" fn cold_start(mut work: *mut OSQPWorkspace) {
    vec_set_scalar((*work).x, 0.0f64, (*(*work).data).n);
    vec_set_scalar((*work).z, 0.0f64, (*(*work).data).m);
    vec_set_scalar((*work).y, 0.0f64, (*(*work).data).m);
}
unsafe extern "C" fn compute_rhs(mut work: *mut OSQPWorkspace) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).n {
        *((*work).xz_tilde)
            .offset(
                i as isize,
            ) = (*(*work).settings).sigma * *((*work).x_prev).offset(i as isize)
            - *((*(*work).data).q).offset(i as isize);
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).m {
        *((*work).xz_tilde)
            .offset(
                (i + (*(*work).data).n) as isize,
            ) = *((*work).z_prev).offset(i as isize)
            - *((*work).rho_inv_vec).offset(i as isize)
                * *((*work).y).offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn update_xz_tilde(mut work: *mut OSQPWorkspace) {
    compute_rhs(work);
    ((*(*work).linsys_solver).solve)
        .expect("non-null function pointer")((*work).linsys_solver, (*work).xz_tilde);
}
#[no_mangle]
pub unsafe extern "C" fn update_x(mut work: *mut OSQPWorkspace) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).n {
        *((*work).x)
            .offset(
                i as isize,
            ) = (*(*work).settings).alpha * *((*work).xz_tilde).offset(i as isize)
            + (1.0f64 - (*(*work).settings).alpha)
                * *((*work).x_prev).offset(i as isize);
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).n {
        *((*work).delta_x)
            .offset(
                i as isize,
            ) = *((*work).x).offset(i as isize) - *((*work).x_prev).offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn update_z(mut work: *mut OSQPWorkspace) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).m {
        *((*work).z)
            .offset(
                i as isize,
            ) = (*(*work).settings).alpha
            * *((*work).xz_tilde).offset((i + (*(*work).data).n) as isize)
            + (1.0f64 - (*(*work).settings).alpha) * *((*work).z_prev).offset(i as isize)
            + *((*work).rho_inv_vec).offset(i as isize)
                * *((*work).y).offset(i as isize);
        i += 1;
    }
    project(work, (*work).z);
}
#[no_mangle]
pub unsafe extern "C" fn update_y(mut work: *mut OSQPWorkspace) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).m {
        *((*work).delta_y)
            .offset(
                i as isize,
            ) = *((*work).rho_vec).offset(i as isize)
            * ((*(*work).settings).alpha
                * *((*work).xz_tilde).offset((i + (*(*work).data).n) as isize)
                + (1.0f64 - (*(*work).settings).alpha)
                    * *((*work).z_prev).offset(i as isize)
                - *((*work).z).offset(i as isize));
        let ref mut fresh1 = *((*work).y).offset(i as isize);
        *fresh1 += *((*work).delta_y).offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn compute_obj_val(
    mut work: *mut OSQPWorkspace,
    mut x: *mut c_float,
) -> c_float {
    let mut obj_val: c_float = 0.;
    obj_val = quad_form((*(*work).data).P, x)
        + vec_prod((*(*work).data).q, x, (*(*work).data).n);
    if (*(*work).settings).scaling != 0 {
        obj_val *= (*(*work).scaling).cinv;
    }
    return obj_val;
}
#[no_mangle]
pub unsafe extern "C" fn compute_pri_res(
    mut work: *mut OSQPWorkspace,
    mut x: *mut c_float,
    mut z: *mut c_float,
) -> c_float {
    mat_vec((*(*work).data).A, x, (*work).Ax, 0 as ::std::os::raw::c_int);
    vec_add_scaled(
        (*work).z_prev,
        (*work).Ax,
        z,
        (*(*work).data).m,
        -(1 as ::std::os::raw::c_int) as c_float,
    );
    if (*(*work).settings).scaling != 0 && (*(*work).settings).scaled_termination == 0 {
        return vec_scaled_norm_inf(
            (*(*work).scaling).Einv,
            (*work).z_prev,
            (*(*work).data).m,
        );
    }
    return vec_norm_inf((*work).z_prev, (*(*work).data).m);
}
#[no_mangle]
pub unsafe extern "C" fn compute_pri_tol(
    mut work: *mut OSQPWorkspace,
    mut eps_abs: c_float,
    mut eps_rel: c_float,
) -> c_float {
    let mut max_rel_eps: c_float = 0.;
    let mut temp_rel_eps: c_float = 0.;
    if (*(*work).settings).scaling != 0 && (*(*work).settings).scaled_termination == 0 {
        max_rel_eps = vec_scaled_norm_inf(
            (*(*work).scaling).Einv,
            (*work).z,
            (*(*work).data).m,
        );
        temp_rel_eps = vec_scaled_norm_inf(
            (*(*work).scaling).Einv,
            (*work).Ax,
            (*(*work).data).m,
        );
        max_rel_eps = if max_rel_eps > temp_rel_eps {
            max_rel_eps
        } else {
            temp_rel_eps
        };
    } else {
        max_rel_eps = vec_norm_inf((*work).z, (*(*work).data).m);
        temp_rel_eps = vec_norm_inf((*work).Ax, (*(*work).data).m);
        max_rel_eps = if max_rel_eps > temp_rel_eps {
            max_rel_eps
        } else {
            temp_rel_eps
        };
    }
    return eps_abs + eps_rel * max_rel_eps;
}
#[no_mangle]
pub unsafe extern "C" fn compute_dua_res(
    mut work: *mut OSQPWorkspace,
    mut x: *mut c_float,
    mut y: *mut c_float,
) -> c_float {
    prea_vec_copy((*(*work).data).q, (*work).x_prev, (*(*work).data).n);
    mat_vec((*(*work).data).P, x, (*work).Px, 0 as ::std::os::raw::c_int);
    mat_tpose_vec((*(*work).data).P, x, (*work).Px, 1 as ::std::os::raw::c_int, 1 as ::std::os::raw::c_int);
    vec_add_scaled(
        (*work).x_prev,
        (*work).x_prev,
        (*work).Px,
        (*(*work).data).n,
        1 as ::std::os::raw::c_int as c_float,
    );
    if (*(*work).data).m > 0 as ::std::os::raw::c_int {
        mat_tpose_vec(
            (*(*work).data).A,
            y,
            (*work).Aty,
            0 as ::std::os::raw::c_int,
            0 as ::std::os::raw::c_int,
        );
        vec_add_scaled(
            (*work).x_prev,
            (*work).x_prev,
            (*work).Aty,
            (*(*work).data).n,
            1 as ::std::os::raw::c_int as c_float,
        );
    }
    if (*(*work).settings).scaling != 0 && (*(*work).settings).scaled_termination == 0 {
        return (*(*work).scaling).cinv
            * vec_scaled_norm_inf(
                (*(*work).scaling).Dinv,
                (*work).x_prev,
                (*(*work).data).n,
            );
    }
    return vec_norm_inf((*work).x_prev, (*(*work).data).n);
}
#[no_mangle]
pub unsafe extern "C" fn compute_dua_tol(
    mut work: *mut OSQPWorkspace,
    mut eps_abs: c_float,
    mut eps_rel: c_float,
) -> c_float {
    let mut max_rel_eps: c_float = 0.;
    let mut temp_rel_eps: c_float = 0.;
    if (*(*work).settings).scaling != 0 && (*(*work).settings).scaled_termination == 0 {
        max_rel_eps = vec_scaled_norm_inf(
            (*(*work).scaling).Dinv,
            (*(*work).data).q,
            (*(*work).data).n,
        );
        temp_rel_eps = vec_scaled_norm_inf(
            (*(*work).scaling).Dinv,
            (*work).Aty,
            (*(*work).data).n,
        );
        max_rel_eps = if max_rel_eps > temp_rel_eps {
            max_rel_eps
        } else {
            temp_rel_eps
        };
        temp_rel_eps = vec_scaled_norm_inf(
            (*(*work).scaling).Dinv,
            (*work).Px,
            (*(*work).data).n,
        );
        max_rel_eps = if max_rel_eps > temp_rel_eps {
            max_rel_eps
        } else {
            temp_rel_eps
        };
        max_rel_eps *= (*(*work).scaling).cinv;
    } else {
        max_rel_eps = vec_norm_inf((*(*work).data).q, (*(*work).data).n);
        temp_rel_eps = vec_norm_inf((*work).Aty, (*(*work).data).n);
        max_rel_eps = if max_rel_eps > temp_rel_eps {
            max_rel_eps
        } else {
            temp_rel_eps
        };
        temp_rel_eps = vec_norm_inf((*work).Px, (*(*work).data).n);
        max_rel_eps = if max_rel_eps > temp_rel_eps {
            max_rel_eps
        } else {
            temp_rel_eps
        };
    }
    return eps_abs + eps_rel * max_rel_eps;
}
#[no_mangle]
pub unsafe extern "C" fn is_primal_infeasible(
    mut work: *mut OSQPWorkspace,
    mut eps_prim_inf: c_float,
) -> c_int {
    let mut i: c_int = 0;
    let mut norm_delta_y: c_float = 0.;
    let mut ineq_lhs: c_float = 0.0f64;
    i = 0 as ::std::os::raw::c_int;
    while i < (*(*work).data).m {
        if *((*(*work).data).u).offset(i as isize) > OSQP_INFTY * MIN_SCALING {
            if *((*(*work).data).l).offset(i as isize) < -OSQP_INFTY * MIN_SCALING {
                *((*work).delta_y).offset(i as isize) = 0.0f64;
            } else {
                *((*work).delta_y)
                    .offset(
                        i as isize,
                    ) = if *((*work).delta_y).offset(i as isize) < 0.0f64 {
                    *((*work).delta_y).offset(i as isize)
                } else {
                    0.0f64
                };
            }
        } else if *((*(*work).data).l).offset(i as isize) < -OSQP_INFTY * MIN_SCALING {
            *((*work).delta_y)
                .offset(
                    i as isize,
                ) = if *((*work).delta_y).offset(i as isize) > 0.0f64 {
                *((*work).delta_y).offset(i as isize)
            } else {
                0.0f64
            };
        }
        i += 1;
    }
    if (*(*work).settings).scaling != 0 && (*(*work).settings).scaled_termination == 0 {
        vec_ew_prod(
            (*(*work).scaling).E,
            (*work).delta_y,
            (*work).Adelta_x,
            (*(*work).data).m,
        );
        norm_delta_y = vec_norm_inf((*work).Adelta_x, (*(*work).data).m);
    } else {
        norm_delta_y = vec_norm_inf((*work).delta_y, (*(*work).data).m);
    }
    if norm_delta_y > OSQP_DIVISION_TOL {
        i = 0 as ::std::os::raw::c_int;
        while i < (*(*work).data).m {
            ineq_lhs
                += *((*(*work).data).u).offset(i as isize)
                    * (if *((*work).delta_y).offset(i as isize)
                        > 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
                    {
                        *((*work).delta_y).offset(i as isize)
                    } else {
                        0 as ::std::os::raw::c_int as ::std::os::raw::c_double
                    })
                    + *((*(*work).data).l).offset(i as isize)
                        * (if *((*work).delta_y).offset(i as isize)
                            < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
                        {
                            *((*work).delta_y).offset(i as isize)
                        } else {
                            0 as ::std::os::raw::c_int as ::std::os::raw::c_double
                        });
            i += 1;
        }
        if ineq_lhs < eps_prim_inf * norm_delta_y {
            mat_tpose_vec(
                (*(*work).data).A,
                (*work).delta_y,
                (*work).Atdelta_y,
                0 as ::std::os::raw::c_int,
                0 as ::std::os::raw::c_int,
            );
            if (*(*work).settings).scaling != 0
                && (*(*work).settings).scaled_termination == 0
            {
                vec_ew_prod(
                    (*(*work).scaling).Dinv,
                    (*work).Atdelta_y,
                    (*work).Atdelta_y,
                    (*(*work).data).n,
                );
            }
            return (vec_norm_inf((*work).Atdelta_y, (*(*work).data).n)
                < eps_prim_inf * norm_delta_y) as ::std::os::raw::c_int;
        }
    }
    return 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn is_dual_infeasible(
    mut work: *mut OSQPWorkspace,
    mut eps_dual_inf: c_float,
) -> c_int {
    let mut i: c_int = 0;
    let mut norm_delta_x: c_float = 0.;
    let mut cost_scaling: c_float = 0.;
    if (*(*work).settings).scaling != 0 && (*(*work).settings).scaled_termination == 0 {
        norm_delta_x = vec_scaled_norm_inf(
            (*(*work).scaling).D,
            (*work).delta_x,
            (*(*work).data).n,
        );
        cost_scaling = (*(*work).scaling).c;
    } else {
        norm_delta_x = vec_norm_inf((*work).delta_x, (*(*work).data).n);
        cost_scaling = 1.0f64;
    }
    if norm_delta_x > OSQP_DIVISION_TOL {
        if vec_prod((*(*work).data).q, (*work).delta_x, (*(*work).data).n)
            < cost_scaling * eps_dual_inf * norm_delta_x
        {
            mat_vec(
                (*(*work).data).P,
                (*work).delta_x,
                (*work).Pdelta_x,
                0 as ::std::os::raw::c_int,
            );
            mat_tpose_vec(
                (*(*work).data).P,
                (*work).delta_x,
                (*work).Pdelta_x,
                1 as ::std::os::raw::c_int,
                1 as ::std::os::raw::c_int,
            );
            if (*(*work).settings).scaling != 0
                && (*(*work).settings).scaled_termination == 0
            {
                vec_ew_prod(
                    (*(*work).scaling).Dinv,
                    (*work).Pdelta_x,
                    (*work).Pdelta_x,
                    (*(*work).data).n,
                );
            }
            if vec_norm_inf((*work).Pdelta_x, (*(*work).data).n)
                < cost_scaling * eps_dual_inf * norm_delta_x
            {
                mat_vec(
                    (*(*work).data).A,
                    (*work).delta_x,
                    (*work).Adelta_x,
                    0 as ::std::os::raw::c_int,
                );
                if (*(*work).settings).scaling != 0
                    && (*(*work).settings).scaled_termination == 0
                {
                    vec_ew_prod(
                        (*(*work).scaling).Einv,
                        (*work).Adelta_x,
                        (*work).Adelta_x,
                        (*(*work).data).m,
                    );
                }
                i = 0 as ::std::os::raw::c_int;
                while i < (*(*work).data).m {
                    if *((*(*work).data).u).offset(i as isize) < OSQP_INFTY * MIN_SCALING
                        && *((*work).Adelta_x).offset(i as isize)
                            > eps_dual_inf * norm_delta_x
                        || *((*(*work).data).l).offset(i as isize)
                            > -OSQP_INFTY * MIN_SCALING
                            && *((*work).Adelta_x).offset(i as isize)
                                < -eps_dual_inf * norm_delta_x
                    {
                        return 0 as ::std::os::raw::c_int;
                    }
                    i += 1;
                }
                return 1 as ::std::os::raw::c_int;
            }
        }
    }
    return 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn has_solution(mut info: *mut OSQPInfo) -> c_int {
    return ((*info).status_val != OSQP_PRIMAL_INFEASIBLE
        && (*info).status_val != OSQP_PRIMAL_INFEASIBLE_INACCURATE
        && (*info).status_val != OSQP_DUAL_INFEASIBLE
        && (*info).status_val != OSQP_DUAL_INFEASIBLE_INACCURATE
        && (*info).status_val != OSQP_NON_CVX) as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn store_solution(mut work: *mut OSQPWorkspace) {
    let mut norm_vec: c_float = 0.;
    if has_solution((*work).info) != 0 {
        prea_vec_copy((*work).x, (*(*work).solution).x, (*(*work).data).n);
        prea_vec_copy((*work).y, (*(*work).solution).y, (*(*work).data).m);
        if (*(*work).settings).scaling != 0 {
            unscale_solution(work);
        }
    } else {
        vec_set_scalar((*(*work).solution).x, OSQP_NAN as c_float, (*(*work).data).n);
        vec_set_scalar((*(*work).solution).y, OSQP_NAN as c_float, (*(*work).data).m);
        if (*(*work).info).status_val == OSQP_PRIMAL_INFEASIBLE
            || (*(*work).info).status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE
        {
            norm_vec = vec_norm_inf((*work).delta_y, (*(*work).data).m);
            vec_mult_scalar((*work).delta_y, 1.0f64 / norm_vec, (*(*work).data).m);
        }
        if (*(*work).info).status_val == OSQP_DUAL_INFEASIBLE
            || (*(*work).info).status_val == OSQP_DUAL_INFEASIBLE_INACCURATE
        {
            norm_vec = vec_norm_inf((*work).delta_x, (*(*work).data).n);
            vec_mult_scalar((*work).delta_x, 1.0f64 / norm_vec, (*(*work).data).n);
        }
        cold_start(work);
    };
}
#[no_mangle]
pub unsafe extern "C" fn update_info(
    mut work: *mut OSQPWorkspace,
    mut iter: c_int,
    mut compute_objective: c_int,
    mut polish: c_int,
) {
    let mut x: *mut c_float = 0 as *mut c_float;
    let mut z: *mut c_float = 0 as *mut c_float;
    let mut y: *mut c_float = 0 as *mut c_float;
    let mut obj_val: *mut c_float = 0 as *mut c_float;
    let mut pri_res: *mut c_float = 0 as *mut c_float;
    let mut dua_res: *mut c_float = 0 as *mut c_float;
    let mut run_time: *mut c_float = 0 as *mut c_float;
    if polish != 0 {
        x = (*(*work).pol).x;
        y = (*(*work).pol).y;
        z = (*(*work).pol).z;
        obj_val = &mut (*(*work).pol).obj_val;
        pri_res = &mut (*(*work).pol).pri_res;
        dua_res = &mut (*(*work).pol).dua_res;
        run_time = &mut (*(*work).info).polish_time;
    } else {
        x = (*work).x;
        y = (*work).y;
        z = (*work).z;
        obj_val = &mut (*(*work).info).obj_val;
        pri_res = &mut (*(*work).info).pri_res;
        dua_res = &mut (*(*work).info).dua_res;
        (*(*work).info).iter = iter;
        run_time = &mut (*(*work).info).solve_time;
    }
    if compute_objective != 0 {
        *obj_val = compute_obj_val(work, x);
    }
    if (*(*work).data).m == 0 as ::std::os::raw::c_int {
        *pri_res = 0.0f64;
    } else {
        *pri_res = compute_pri_res(work, x, z);
    }
    *dua_res = compute_dua_res(work, x, y);
    *run_time = osqp_toc((*work).timer);
    (*work).summary_printed = 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn reset_info(mut info: *mut OSQPInfo) {
    (*info).solve_time = 0.0f64;
    (*info).polish_time = 0.0f64;
    update_status(info, OSQP_UNSOLVED);
    (*info).rho_updates = 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn update_status(mut info: *mut OSQPInfo, mut status_val: c_int) {
    (*info).status_val = status_val;
    if status_val == OSQP_SOLVED {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"solved\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    }
    if status_val == OSQP_SOLVED_INACCURATE {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"solved inaccurate\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_PRIMAL_INFEASIBLE {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"primal infeasible\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"primal infeasible inaccurate\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_UNSOLVED {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"unsolved\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_DUAL_INFEASIBLE {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"dual infeasible\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_DUAL_INFEASIBLE_INACCURATE {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"dual infeasible inaccurate\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_MAX_ITER_REACHED {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"maximum iterations reached\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_TIME_LIMIT_REACHED {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"run time limit reached\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_SIGINT {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"interrupted\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    } else if status_val == OSQP_NON_CVX {
        c_strcpy(
            ((*info).status).as_mut_ptr(),
            b"problem non convex\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    }
}
#[no_mangle]
pub unsafe extern "C" fn check_termination(
    mut work: *mut OSQPWorkspace,
    mut approximate: c_int,
) -> c_int {
    let mut eps_prim: c_float = 0.;
    let mut eps_dual: c_float = 0.;
    let mut eps_prim_inf: c_float = 0.;
    let mut eps_dual_inf: c_float = 0.;
    let mut exitflag: c_int = 0;
    let mut prim_res_check: c_int = 0;
    let mut dual_res_check: c_int = 0;
    let mut prim_inf_check: c_int = 0;
    let mut dual_inf_check: c_int = 0;
    let mut eps_abs: c_float = 0.;
    let mut eps_rel: c_float = 0.;
    exitflag = 0 as ::std::os::raw::c_int;
    prim_res_check = 0 as ::std::os::raw::c_int;
    dual_res_check = 0 as ::std::os::raw::c_int;
    prim_inf_check = 0 as ::std::os::raw::c_int;
    dual_inf_check = 0 as ::std::os::raw::c_int;
    eps_abs = (*(*work).settings).eps_abs;
    eps_rel = (*(*work).settings).eps_rel;
    eps_prim_inf = (*(*work).settings).eps_prim_inf;
    eps_dual_inf = (*(*work).settings).eps_dual_inf;
    if (*(*work).info).pri_res > OSQP_INFTY || (*(*work).info).dua_res > OSQP_INFTY {
        update_status((*work).info, OSQP_NON_CVX);
        (*(*work).info).obj_val = OSQP_NAN as c_float;
        return 1 as ::std::os::raw::c_int;
    }
    if approximate != 0 {
        eps_abs *= 10 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        eps_rel *= 10 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        eps_prim_inf *= 10 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        eps_dual_inf *= 10 as ::std::os::raw::c_int as ::std::os::raw::c_double;
    }
    if (*(*work).data).m == 0 as ::std::os::raw::c_int {
        prim_res_check = 1 as ::std::os::raw::c_int;
    } else {
        eps_prim = compute_pri_tol(work, eps_abs, eps_rel);
        if (*(*work).info).pri_res < eps_prim {
            prim_res_check = 1 as ::std::os::raw::c_int;
        } else {
            prim_inf_check = is_primal_infeasible(work, eps_prim_inf);
        }
    }
    eps_dual = compute_dua_tol(work, eps_abs, eps_rel);
    if (*(*work).info).dua_res < eps_dual {
        dual_res_check = 1 as ::std::os::raw::c_int;
    } else {
        dual_inf_check = is_dual_infeasible(work, eps_dual_inf);
    }
    if prim_res_check != 0 && dual_res_check != 0 {
        if approximate != 0 {
            update_status((*work).info, OSQP_SOLVED_INACCURATE);
        } else {
            update_status((*work).info, OSQP_SOLVED);
        }
        exitflag = 1 as ::std::os::raw::c_int;
    } else if prim_inf_check != 0 {
        if approximate != 0 {
            update_status((*work).info, OSQP_PRIMAL_INFEASIBLE_INACCURATE);
        } else {
            update_status((*work).info, OSQP_PRIMAL_INFEASIBLE);
        }
        if (*(*work).settings).scaling != 0
            && (*(*work).settings).scaled_termination == 0
        {
            vec_ew_prod(
                (*(*work).scaling).E,
                (*work).delta_y,
                (*work).delta_y,
                (*(*work).data).m,
            );
        }
        (*(*work).info).obj_val = OSQP_INFTY;
        exitflag = 1 as ::std::os::raw::c_int;
    } else if dual_inf_check != 0 {
        if approximate != 0 {
            update_status((*work).info, OSQP_DUAL_INFEASIBLE_INACCURATE);
        } else {
            update_status((*work).info, OSQP_DUAL_INFEASIBLE);
        }
        if (*(*work).settings).scaling != 0
            && (*(*work).settings).scaled_termination == 0
        {
            vec_ew_prod(
                (*(*work).scaling).D,
                (*work).delta_x,
                (*work).delta_x,
                (*(*work).data).n,
            );
        }
        (*(*work).info).obj_val = -OSQP_INFTY;
        exitflag = 1 as ::std::os::raw::c_int;
    }
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn validate_data(mut data: *const OSQPData) -> c_int {
    let mut j: c_int = 0;
    let mut ptr: c_int = 0;
    if data.is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(b"Missing data\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if ((*data).P).is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(b"Missing matrix P\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if ((*data).A).is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(b"Missing matrix A\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if ((*data).q).is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(b"Missing vector q\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*data).n <= 0 as ::std::os::raw::c_int || (*data).m < 0 as ::std::os::raw::c_int {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(
            b"n must be positive and m nonnegative; n = %i, m = %i\0" as *const u8
                as *const ::std::os::raw::c_char,
            (*data).n,
            (*data).m,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*(*data).P).m != (*data).n {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(
            b"P does not have dimension n x n with n = %i\0" as *const u8
                as *const ::std::os::raw::c_char,
            (*data).n,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*(*data).P).m != (*(*data).P).n {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(b"P is not square\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < (*data).n {
        ptr = *((*(*data).P).p).offset(j as isize);
        while ptr < *((*(*data).P).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            if *((*(*data).P).i).offset(ptr as isize) > j {
                printf(
                    b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
                    (*::std::mem::transmute::<
                        &[u8; 14],
                        &[::std::os::raw::c_char; 14],
                    >(b"validate_data\0"))
                        .as_ptr(),
                );
                printf(
                    b"P is not upper triangular\0" as *const u8 as *const ::std::os::raw::c_char,
                );
                printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
                return 1 as ::std::os::raw::c_int;
            }
            ptr += 1;
        }
        j += 1;
    }
    if (*(*data).A).m != (*data).m || (*(*data).A).n != (*data).n {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[::std::os::raw::c_char; 14],
            >(b"validate_data\0"))
                .as_ptr(),
        );
        printf(
            b"A does not have dimension %i x %i\0" as *const u8 as *const ::std::os::raw::c_char,
            (*data).m,
            (*data).n,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < (*data).m {
        if *((*data).l).offset(j as isize) > *((*data).u).offset(j as isize) {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
                (*::std::mem::transmute::<
                    &[u8; 14],
                    &[::std::os::raw::c_char; 14],
                >(b"validate_data\0"))
                    .as_ptr(),
            );
            printf(
                b"Lower bound at index %d is greater than upper bound: %.4e > %.4e\0"
                    as *const u8 as *const ::std::os::raw::c_char,
                j,
                *((*data).l).offset(j as isize),
                *((*data).u).offset(j as isize),
            );
            printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
            return 1 as ::std::os::raw::c_int;
        }
        j += 1;
    }
    return 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn validate_linsys_solver(mut linsys_solver: c_int) -> c_int {
    if linsys_solver != QDLDL_SOLVER as ::std::os::raw::c_int
        && linsys_solver != MKL_PARDISO_SOLVER as ::std::os::raw::c_int
    {
        return 1 as ::std::os::raw::c_int;
    }
    return 0 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn validate_settings(mut settings: *const OSQPSettings) -> c_int {
    if settings.is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"Missing settings!\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).scaling < 0 as ::std::os::raw::c_int {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"scaling must be nonnegative\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).adaptive_rho != 0 as ::std::os::raw::c_int
        && (*settings).adaptive_rho != 1 as ::std::os::raw::c_int
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"adaptive_rho must be either 0 or 1\0" as *const u8 as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).adaptive_rho_interval < 0 as ::std::os::raw::c_int {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"adaptive_rho_interval must be nonnegative\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).adaptive_rho_fraction <= 0 as ::std::os::raw::c_int as ::std::os::raw::c_double {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"adaptive_rho_fraction must be positive\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).adaptive_rho_tolerance < 1.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"adaptive_rho_tolerance must be >= 1\0" as *const u8 as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).polish_refine_iter < 0 as ::std::os::raw::c_int {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"polish_refine_iter must be nonnegative\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).rho <= 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"rho must be positive\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).sigma <= 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"sigma must be positive\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).delta <= 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"delta must be positive\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).max_iter <= 0 as ::std::os::raw::c_int {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"max_iter must be positive\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).eps_abs < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"eps_abs must be nonnegative\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).eps_rel < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"eps_rel must be nonnegative\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).eps_rel == 0.0f64 && (*settings).eps_abs == 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"at least one of eps_abs and eps_rel must be positive\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).eps_prim_inf <= 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"eps_prim_inf must be positive\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).eps_dual_inf <= 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"eps_dual_inf must be positive\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).alpha <= 0.0f64 || (*settings).alpha >= 2.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"alpha must be strictly between 0 and 2\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if validate_linsys_solver((*settings).linsys_solver as c_int) != 0 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"linsys_solver not recognized\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).verbose != 0 as ::std::os::raw::c_int && (*settings).verbose != 1 as ::std::os::raw::c_int
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(b"verbose must be either 0 or 1\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).scaled_termination != 0 as ::std::os::raw::c_int
        && (*settings).scaled_termination != 1 as ::std::os::raw::c_int
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"scaled_termination must be either 0 or 1\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).check_termination < 0 as ::std::os::raw::c_int {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"check_termination must be nonnegative\0" as *const u8
                as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).warm_start != 0 as ::std::os::raw::c_int
        && (*settings).warm_start != 1 as ::std::os::raw::c_int
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"warm_start must be either 0 or 1\0" as *const u8 as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    if (*settings).time_limit < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[::std::os::raw::c_char; 18],
            >(b"validate_settings\0"))
                .as_ptr(),
        );
        printf(
            b"time_limit must be nonnegative\n\0" as *const u8 as *const ::std::os::raw::c_char,
        );
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
        return 1 as ::std::os::raw::c_int;
    }
    return 0 as ::std::os::raw::c_int;
}
