use ::libc;
extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
    fn copy_csc_mat(A: *const csc) -> *mut csc;
    fn csc_spfree(A: *mut csc);
    fn mat_vec(A: *const csc, x: *const c_float, y: *mut c_float, plus_eq: c_int);
    fn vec_ew_prod(a: *const c_float, b: *const c_float, c: *mut c_float, n: c_int);
    fn vec_mult_scalar(a: *mut c_float, sc: c_float, n: c_int);
    fn prea_vec_copy(a: *const c_float, b: *mut c_float, n: c_int);
    fn vec_copy(a: *mut c_float, n: c_int) -> *mut c_float;
    fn osqp_toc(t: *mut OSQPTimer) -> c_float;
    fn osqp_tic(t: *mut OSQPTimer);
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
    fn copy_settings(settings: *const OSQPSettings) -> *mut OSQPSettings;
    fn print_setup_header(work: *const OSQPWorkspace);
    fn print_header();
    fn print_summary(work: *mut OSQPWorkspace);
    fn print_footer(info: *mut OSQPInfo, polish_0: c_int);
    fn fmod(_: libc::c_double, _: libc::c_double) -> libc::c_double;
    fn compute_rho_estimate(work: *mut OSQPWorkspace) -> c_float;
    fn adapt_rho(work: *mut OSQPWorkspace) -> c_int;
    fn set_rho_vec(work: *mut OSQPWorkspace);
    fn update_rho_vec(work: *mut OSQPWorkspace) -> c_int;
    fn swap_vectors(a: *mut *mut c_float, b: *mut *mut c_float);
    fn cold_start(work: *mut OSQPWorkspace);
    fn update_xz_tilde(work: *mut OSQPWorkspace);
    fn update_x(work: *mut OSQPWorkspace);
    fn update_z(work: *mut OSQPWorkspace);
    fn update_y(work: *mut OSQPWorkspace);
    fn compute_obj_val(work: *mut OSQPWorkspace, x: *mut c_float) -> c_float;
    fn has_solution(info: *mut OSQPInfo) -> c_int;
    fn store_solution(work: *mut OSQPWorkspace);
    fn update_info(
        work: *mut OSQPWorkspace,
        iter: c_int,
        compute_objective: c_int,
        polish_0: c_int,
    );
    fn reset_info(info: *mut OSQPInfo);
    fn update_status(info: *mut OSQPInfo, status_val: c_int);
    fn check_termination(work: *mut OSQPWorkspace, approximate: c_int) -> c_int;
    fn validate_data(data: *const OSQPData) -> c_int;
    fn validate_settings(settings: *const OSQPSettings) -> c_int;
    fn scale_data(work: *mut OSQPWorkspace) -> c_int;
    fn unscale_data(work: *mut OSQPWorkspace) -> c_int;
    fn _osqp_error(
        error_code: osqp_error_type,
        function_name: *const libc::c_char,
    ) -> c_int;
    fn polish(work: *mut OSQPWorkspace) -> c_int;
    fn load_linsys_solver(linsys_solver: linsys_solver_type) -> c_int;
    fn unload_linsys_solver(linsys_solver: linsys_solver_type) -> c_int;
    fn init_linsys_solver(
        s: *mut *mut LinSysSolver,
        P: *const csc,
        A: *const csc,
        sigma: c_float,
        rho_vec: *const c_float,
        linsys_solver: linsys_solver_type,
        polish_0: c_int,
    ) -> c_int;
}
pub type c_int = libc::c_longlong;
pub type c_float = libc::c_double;
pub type linsys_solver_type = libc::c_uint;
pub const MKL_PARDISO_SOLVER: linsys_solver_type = 1;
pub const QDLDL_SOLVER: linsys_solver_type = 0;
pub type osqp_error_type = libc::c_uint;
pub const OSQP_WORKSPACE_NOT_INIT_ERROR: osqp_error_type = 7;
pub const OSQP_MEM_ALLOC_ERROR: osqp_error_type = 6;
pub const OSQP_NONCVX_ERROR: osqp_error_type = 5;
pub const OSQP_LINSYS_SOLVER_INIT_ERROR: osqp_error_type = 4;
pub const OSQP_LINSYS_SOLVER_LOAD_ERROR: osqp_error_type = 3;
pub const OSQP_SETTINGS_VALIDATION_ERROR: osqp_error_type = 2;
pub const OSQP_DATA_VALIDATION_ERROR: osqp_error_type = 1;
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
pub type uint32_t = libc::c_uint;
pub type uint64_t = libc::c_ulonglong;
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
pub const RHO_EQ_OVER_RHO_INEQ: libc::c_double = 1e03f64;
pub const c_free: unsafe extern "C" fn(*mut libc::c_void) -> () = free;
pub const PRINT_INTERVAL: libc::c_int = 200 as libc::c_int;
pub const c_fmod: unsafe extern "C" fn(
    libc::c_double,
    libc::c_double,
) -> libc::c_double = fmod;
pub const c_print: unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int = printf;
pub const OSQP_MAX_ITER_REACHED: libc::c_int = -(2 as libc::c_int);
pub const OSQP_TIME_LIMIT_REACHED: libc::c_int = -(6 as libc::c_int);
pub const OSQP_SOLVED: libc::c_int = 1 as libc::c_int;
pub const OSQP_NULL: libc::c_int = 0 as libc::c_int;
pub const c_malloc: unsafe extern "C" fn(libc::c_ulong) -> *mut libc::c_void = malloc;
pub const c_calloc: unsafe extern "C" fn(
    libc::c_ulong,
    libc::c_ulong,
) -> *mut libc::c_void = calloc;
pub const OSQP_UNSOLVED: libc::c_int = -(10 as libc::c_int);
pub const RHO: libc::c_double = 0.1f64;
pub const SIGMA: libc::c_double = 1E-06f64;
pub const SCALING: libc::c_int = 10 as libc::c_int;
pub const ADAPTIVE_RHO: libc::c_int = 1 as libc::c_int;
pub const ADAPTIVE_RHO_INTERVAL: libc::c_int = 0 as libc::c_int;
pub const ADAPTIVE_RHO_TOLERANCE: libc::c_int = 5 as libc::c_int;
pub const ADAPTIVE_RHO_FRACTION: libc::c_double = 0.4f64;
pub const MAX_ITER: libc::c_int = 4000 as libc::c_int;
pub const EPS_ABS: libc::c_double = 1E-3f64;
pub const EPS_REL: libc::c_double = 1E-3f64;
pub const EPS_PRIM_INF: libc::c_double = 1E-4f64;
pub const EPS_DUAL_INF: libc::c_double = 1E-4f64;
pub const ALPHA: libc::c_double = 1.6f64;
pub const LINSYS_SOLVER: libc::c_int = QDLDL_SOLVER as libc::c_int;
pub const DELTA: libc::c_double = 1E-6f64;
pub const POLISH: libc::c_int = 0 as libc::c_int;
pub const POLISH_REFINE_ITER: libc::c_int = 3 as libc::c_int;
pub const VERBOSE: libc::c_int = 1 as libc::c_int;
pub const SCALED_TERMINATION: libc::c_int = 0 as libc::c_int;
pub const CHECK_TERMINATION: libc::c_int = 25 as libc::c_int;
pub const WARM_START: libc::c_int = 1 as libc::c_int;
pub const TIME_LIMIT: libc::c_int = 0 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn osqp_set_default_settings(mut settings: *mut OSQPSettings) {
    (*settings).rho = RHO;
    (*settings).sigma = SIGMA;
    (*settings).scaling = SCALING as c_int;
    (*settings).adaptive_rho = ADAPTIVE_RHO as c_int;
    (*settings).adaptive_rho_interval = ADAPTIVE_RHO_INTERVAL as c_int;
    (*settings).adaptive_rho_tolerance = ADAPTIVE_RHO_TOLERANCE as c_float;
    (*settings).adaptive_rho_fraction = ADAPTIVE_RHO_FRACTION;
    (*settings).max_iter = MAX_ITER as c_int;
    (*settings).eps_abs = EPS_ABS;
    (*settings).eps_rel = EPS_REL;
    (*settings).eps_prim_inf = EPS_PRIM_INF;
    (*settings).eps_dual_inf = EPS_DUAL_INF;
    (*settings).alpha = ALPHA;
    (*settings).linsys_solver = LINSYS_SOLVER as linsys_solver_type;
    (*settings).delta = DELTA;
    (*settings).polish = POLISH as c_int;
    (*settings).polish_refine_iter = POLISH_REFINE_ITER as c_int;
    (*settings).verbose = VERBOSE as c_int;
    (*settings).scaled_termination = SCALED_TERMINATION as c_int;
    (*settings).check_termination = CHECK_TERMINATION as c_int;
    (*settings).warm_start = WARM_START as c_int;
    (*settings).time_limit = TIME_LIMIT as c_float;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_setup(
    mut workp: *mut *mut OSQPWorkspace,
    mut data: *const OSQPData,
    mut settings: *const OSQPSettings,
) -> c_int {
    let mut exitflag: c_int = 0;
    let mut work: *mut OSQPWorkspace = 0 as *mut OSQPWorkspace;
    if validate_data(data) != 0 {
        return _osqp_error(
            OSQP_DATA_VALIDATION_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    if validate_settings(settings) != 0 {
        return _osqp_error(
            OSQP_SETTINGS_VALIDATION_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    work = calloc(
        1 as libc::c_int as libc::c_ulong,
        ::std::mem::size_of::<OSQPWorkspace>() as libc::c_ulong,
    ) as *mut OSQPWorkspace;
    if work.is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    *workp = work;
    let ref mut fresh0 = (*work).timer;
    *fresh0 = malloc(::std::mem::size_of::<OSQPTimer>() as libc::c_ulong)
        as *mut OSQPTimer;
    if ((*work).timer).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    osqp_tic((*work).timer);
    let ref mut fresh1 = (*work).data;
    *fresh1 = malloc(::std::mem::size_of::<OSQPData>() as libc::c_ulong)
        as *mut OSQPData;
    if ((*work).data).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    (*(*work).data).n = (*data).n;
    (*(*work).data).m = (*data).m;
    let ref mut fresh2 = (*(*work).data).P;
    *fresh2 = copy_csc_mat((*data).P);
    let ref mut fresh3 = (*(*work).data).q;
    *fresh3 = vec_copy((*data).q, (*data).n);
    if ((*(*work).data).P).is_null() || ((*(*work).data).q).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh4 = (*(*work).data).A;
    *fresh4 = copy_csc_mat((*data).A);
    if ((*(*work).data).A).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh5 = (*(*work).data).l;
    *fresh5 = vec_copy((*data).l, (*data).m);
    let ref mut fresh6 = (*(*work).data).u;
    *fresh6 = vec_copy((*data).u, (*data).m);
    if (*data).m != 0 && (((*(*work).data).l).is_null() || ((*(*work).data).u).is_null())
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh7 = (*work).rho_vec;
    *fresh7 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh8 = (*work).rho_inv_vec;
    *fresh8 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    if (*data).m != 0 && (((*work).rho_vec).is_null() || ((*work).rho_inv_vec).is_null())
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh9 = (*work).constr_type;
    *fresh9 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_int>() as libc::c_ulong,
    ) as *mut c_int;
    if (*data).m != 0 && ((*work).constr_type).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh10 = (*work).x;
    *fresh10 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh11 = (*work).z;
    *fresh11 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh12 = (*work).xz_tilde;
    *fresh12 = calloc(
        ((*data).n + (*data).m) as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh13 = (*work).x_prev;
    *fresh13 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh14 = (*work).z_prev;
    *fresh14 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh15 = (*work).y;
    *fresh15 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    if ((*work).x).is_null() || ((*work).xz_tilde).is_null()
        || ((*work).x_prev).is_null()
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    if (*data).m != 0
        && (((*work).z).is_null() || ((*work).z_prev).is_null() || ((*work).y).is_null())
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    cold_start(work);
    let ref mut fresh16 = (*work).Ax;
    *fresh16 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh17 = (*work).Px;
    *fresh17 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh18 = (*work).Aty;
    *fresh18 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh19 = (*work).delta_y;
    *fresh19 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh20 = (*work).Atdelta_y;
    *fresh20 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh21 = (*work).delta_x;
    *fresh21 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh22 = (*work).Pdelta_x;
    *fresh22 = calloc(
        (*data).n as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh23 = (*work).Adelta_x;
    *fresh23 = calloc(
        (*data).m as libc::c_ulong,
        ::std::mem::size_of::<c_float>() as libc::c_ulong,
    ) as *mut c_float;
    if ((*work).Px).is_null() || ((*work).Aty).is_null() || ((*work).Atdelta_y).is_null()
        || ((*work).delta_x).is_null() || ((*work).Pdelta_x).is_null()
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    if (*data).m != 0
        && (((*work).Ax).is_null() || ((*work).delta_y).is_null()
            || ((*work).Adelta_x).is_null())
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh24 = (*work).settings;
    *fresh24 = copy_settings(settings);
    if ((*work).settings).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    if (*settings).scaling != 0 {
        let ref mut fresh25 = (*work).scaling;
        *fresh25 = malloc(::std::mem::size_of::<OSQPScaling>() as libc::c_ulong)
            as *mut OSQPScaling;
        if ((*work).scaling).is_null() {
            return _osqp_error(
                OSQP_MEM_ALLOC_ERROR,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"osqp_setup\0"))
                    .as_ptr(),
            );
        }
        let ref mut fresh26 = (*(*work).scaling).D;
        *fresh26 = malloc(
            ((*data).n as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        let ref mut fresh27 = (*(*work).scaling).Dinv;
        *fresh27 = malloc(
            ((*data).n as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        let ref mut fresh28 = (*(*work).scaling).E;
        *fresh28 = malloc(
            ((*data).m as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        let ref mut fresh29 = (*(*work).scaling).Einv;
        *fresh29 = malloc(
            ((*data).m as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        if ((*(*work).scaling).D).is_null() || ((*(*work).scaling).Dinv).is_null() {
            return _osqp_error(
                OSQP_MEM_ALLOC_ERROR,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"osqp_setup\0"))
                    .as_ptr(),
            );
        }
        if (*data).m != 0
            && (((*(*work).scaling).E).is_null() || ((*(*work).scaling).Einv).is_null())
        {
            return _osqp_error(
                OSQP_MEM_ALLOC_ERROR,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"osqp_setup\0"))
                    .as_ptr(),
            );
        }
        let ref mut fresh30 = (*work).D_temp;
        *fresh30 = malloc(
            ((*data).n as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        let ref mut fresh31 = (*work).D_temp_A;
        *fresh31 = malloc(
            ((*data).n as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        let ref mut fresh32 = (*work).E_temp;
        *fresh32 = malloc(
            ((*data).m as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_float>() as libc::c_ulong
                        as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_float;
        if ((*work).D_temp).is_null() || ((*work).D_temp_A).is_null() {
            return _osqp_error(
                OSQP_MEM_ALLOC_ERROR,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"osqp_setup\0"))
                    .as_ptr(),
            );
        }
        if (*data).m != 0 && ((*work).E_temp).is_null() {
            return _osqp_error(
                OSQP_MEM_ALLOC_ERROR,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"osqp_setup\0"))
                    .as_ptr(),
            );
        }
        scale_data(work);
    } else {
        let ref mut fresh33 = (*work).scaling;
        *fresh33 = OSQP_NULL as *mut OSQPScaling;
        let ref mut fresh34 = (*work).D_temp;
        *fresh34 = OSQP_NULL as *mut c_float;
        let ref mut fresh35 = (*work).D_temp_A;
        *fresh35 = OSQP_NULL as *mut c_float;
        let ref mut fresh36 = (*work).E_temp;
        *fresh36 = OSQP_NULL as *mut c_float;
    }
    set_rho_vec(work);
    if load_linsys_solver((*(*work).settings).linsys_solver) != 0 {
        return _osqp_error(
            OSQP_LINSYS_SOLVER_LOAD_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    exitflag = init_linsys_solver(
        &mut (*work).linsys_solver,
        (*(*work).data).P,
        (*(*work).data).A,
        (*(*work).settings).sigma,
        (*work).rho_vec,
        (*(*work).settings).linsys_solver,
        0 as libc::c_int as c_int,
    );
    if exitflag != 0 {
        return _osqp_error(
            exitflag as osqp_error_type,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh37 = (*work).pol;
    *fresh37 = malloc(::std::mem::size_of::<OSQPPolish>() as libc::c_ulong)
        as *mut OSQPPolish;
    if ((*work).pol).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh38 = (*(*work).pol).Alow_to_A;
    *fresh38 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_int;
    let ref mut fresh39 = (*(*work).pol).Aupp_to_A;
    *fresh39 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_int;
    let ref mut fresh40 = (*(*work).pol).A_to_Alow;
    *fresh40 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_int;
    let ref mut fresh41 = (*(*work).pol).A_to_Aupp;
    *fresh41 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_int;
    let ref mut fresh42 = (*(*work).pol).x;
    *fresh42 = malloc(
        ((*data).n as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh43 = (*(*work).pol).z;
    *fresh43 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh44 = (*(*work).pol).y;
    *fresh44 = malloc(
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    if ((*(*work).pol).x).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    if (*data).m != 0
        && (((*(*work).pol).Alow_to_A).is_null() || ((*(*work).pol).Aupp_to_A).is_null()
            || ((*(*work).pol).A_to_Alow).is_null()
            || ((*(*work).pol).A_to_Aupp).is_null() || ((*(*work).pol).z).is_null()
            || ((*(*work).pol).y).is_null())
    {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh45 = (*work).solution;
    *fresh45 = calloc(
        1 as libc::c_int as libc::c_ulong,
        ::std::mem::size_of::<OSQPSolution>() as libc::c_ulong,
    ) as *mut OSQPSolution;
    if ((*work).solution).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh46 = (*(*work).solution).x;
    *fresh46 = calloc(
        1 as libc::c_int as libc::c_ulong,
        ((*data).n as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh47 = (*(*work).solution).y;
    *fresh47 = calloc(
        1 as libc::c_int as libc::c_ulong,
        ((*data).m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_float;
    if ((*(*work).solution).x).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    if (*data).m != 0 && ((*(*work).solution).y).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    let ref mut fresh48 = (*work).info;
    *fresh48 = calloc(
        1 as libc::c_int as libc::c_ulong,
        ::std::mem::size_of::<OSQPInfo>() as libc::c_ulong,
    ) as *mut OSQPInfo;
    if ((*work).info).is_null() {
        return _osqp_error(
            OSQP_MEM_ALLOC_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_setup\0"))
                .as_ptr(),
        );
    }
    (*(*work).info).status_polish = 0 as libc::c_int as c_int;
    update_status((*work).info, OSQP_UNSOLVED as c_int);
    (*(*work).info).solve_time = 0.0f64;
    (*(*work).info).update_time = 0.0f64;
    (*(*work).info).polish_time = 0.0f64;
    (*(*work).info).run_time = 0.0f64;
    (*(*work).info).setup_time = osqp_toc((*work).timer);
    (*work).first_run = 1 as libc::c_int as c_int;
    (*work).clear_update_time = 0 as libc::c_int as c_int;
    (*work).rho_update_from_solve = 0 as libc::c_int as c_int;
    (*(*work).info).rho_updates = 0 as libc::c_int as c_int;
    (*(*work).info).rho_estimate = (*(*work).settings).rho;
    if (*(*work).settings).verbose != 0 {
        print_setup_header(work);
    }
    (*work).summary_printed = 0 as libc::c_int as c_int;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_solve(mut work: *mut OSQPWorkspace) -> c_int {
    let mut current_block: u64;
    let mut exitflag: c_int = 0;
    let mut iter: c_int = 0;
    let mut compute_cost_function: c_int = 0;
    let mut can_check_termination: c_int = 0;
    let mut temp_run_time: c_float = 0.;
    let mut can_print: c_int = 0;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"osqp_solve\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*(*work).info).update_time = 0.0f64;
    }
    (*work).rho_update_from_solve = 1 as libc::c_int as c_int;
    exitflag = 0 as libc::c_int as c_int;
    can_check_termination = 0 as libc::c_int as c_int;
    can_print = (*(*work).settings).verbose;
    compute_cost_function = (*(*work).settings).verbose;
    osqp_tic((*work).timer);
    if (*(*work).settings).verbose != 0 {
        print_header();
    }
    if (*(*work).settings).warm_start == 0 {
        cold_start(work);
    }
    iter = 1 as libc::c_int as c_int;
    loop {
        if !(iter <= (*(*work).settings).max_iter) {
            current_block = 13125627826496529465;
            break;
        }
        swap_vectors(&mut (*work).x, &mut (*work).x_prev);
        swap_vectors(&mut (*work).z, &mut (*work).z_prev);
        update_xz_tilde(work);
        update_x(work);
        update_z(work);
        update_y(work);
        if (*work).first_run != 0 {
            temp_run_time = (*(*work).info).setup_time + osqp_toc((*work).timer);
        } else {
            temp_run_time = (*(*work).info).update_time + osqp_toc((*work).timer);
        }
        if (*(*work).settings).time_limit != 0.
            && temp_run_time >= (*(*work).settings).time_limit
        {
            update_status((*work).info, OSQP_TIME_LIMIT_REACHED as c_int);
            if (*(*work).settings).verbose != 0 {
                printf(
                    b"run time limit reached\n\0" as *const u8 as *const libc::c_char,
                );
            }
            can_print = 0 as libc::c_int as c_int;
            current_block = 13125627826496529465;
            break;
        } else {
            can_check_termination = ((*(*work).settings).check_termination != 0
                && iter % (*(*work).settings).check_termination
                    == 0 as libc::c_int as libc::c_longlong) as libc::c_int as c_int;
            can_print = ((*(*work).settings).verbose != 0
                && (iter % PRINT_INTERVAL as libc::c_longlong
                    == 0 as libc::c_int as libc::c_longlong
                    || iter == 1 as libc::c_int as libc::c_longlong)) as libc::c_int
                as c_int;
            if can_check_termination != 0 || can_print != 0 {
                update_info(
                    work,
                    iter,
                    compute_cost_function,
                    0 as libc::c_int as c_int,
                );
                if can_print != 0 {
                    print_summary(work);
                }
                if can_check_termination != 0 {
                    if check_termination(work, 0 as libc::c_int as c_int) != 0 {
                        current_block = 13125627826496529465;
                        break;
                    }
                }
            }
            if (*(*work).settings).adaptive_rho != 0
                && (*(*work).settings).adaptive_rho_interval == 0
            {
                if osqp_toc((*work).timer)
                    > (*(*work).settings).adaptive_rho_fraction
                        * (*(*work).info).setup_time
                {
                    if (*(*work).settings).check_termination != 0 {
                        (*(*work).settings)
                            .adaptive_rho_interval = (iter as libc::c_double
                            + 0.5f64
                                * (*(*work).settings).check_termination as libc::c_double
                            - fmod(
                                iter as libc::c_double
                                    + 0.5f64
                                        * (*(*work).settings).check_termination as libc::c_double,
                                (*(*work).settings).check_termination as libc::c_double,
                            )) as c_int;
                    } else {
                        (*(*work).settings)
                            .adaptive_rho_interval = (iter as libc::c_double
                            + 0.5f64 * 25 as libc::c_int as libc::c_double
                            - fmod(
                                iter as libc::c_double
                                    + 0.5f64 * 25 as libc::c_int as libc::c_double,
                                25 as libc::c_int as libc::c_double,
                            )) as c_int;
                    }
                    (*(*work).settings)
                        .adaptive_rho_interval = if (*(*work).settings)
                        .adaptive_rho_interval > (*(*work).settings).check_termination
                    {
                        (*(*work).settings).adaptive_rho_interval
                    } else {
                        (*(*work).settings).check_termination
                    };
                }
            }
            if (*(*work).settings).adaptive_rho != 0
                && (*(*work).settings).adaptive_rho_interval != 0
                && iter % (*(*work).settings).adaptive_rho_interval
                    == 0 as libc::c_int as libc::c_longlong
            {
                if can_check_termination == 0 && can_print == 0 {
                    update_info(
                        work,
                        iter,
                        compute_cost_function,
                        0 as libc::c_int as c_int,
                    );
                }
                if adapt_rho(work) != 0 {
                    printf(
                        b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                        (*::std::mem::transmute::<
                            &[u8; 11],
                            &[libc::c_char; 11],
                        >(b"osqp_solve\0"))
                            .as_ptr(),
                    );
                    printf(b"Failed rho update\0" as *const u8 as *const libc::c_char);
                    printf(b"\n\0" as *const u8 as *const libc::c_char);
                    exitflag = 1 as libc::c_int as c_int;
                    current_block = 4936170274040873418;
                    break;
                }
            }
            iter += 1;
        }
    }
    match current_block {
        13125627826496529465 => {
            if can_check_termination == 0 {
                if can_print == 0 {
                    update_info(
                        work,
                        iter - 1 as libc::c_int as libc::c_longlong,
                        compute_cost_function,
                        0 as libc::c_int as c_int,
                    );
                }
                if (*(*work).settings).verbose != 0 && (*work).summary_printed == 0 {
                    print_summary(work);
                }
                check_termination(work, 0 as libc::c_int as c_int);
            }
            if compute_cost_function == 0 && has_solution((*work).info) != 0 {
                (*(*work).info).obj_val = compute_obj_val(work, (*work).x);
            }
            if (*(*work).settings).verbose != 0 && (*work).summary_printed == 0 {
                print_summary(work);
            }
            if (*(*work).info).status_val == OSQP_UNSOLVED as libc::c_longlong {
                if check_termination(work, 1 as libc::c_int as c_int) == 0 {
                    update_status((*work).info, OSQP_MAX_ITER_REACHED as c_int);
                }
            }
            if (*(*work).info).status_val == OSQP_TIME_LIMIT_REACHED as libc::c_longlong
            {
                if check_termination(work, 1 as libc::c_int as c_int) == 0 {
                    update_status((*work).info, OSQP_TIME_LIMIT_REACHED as c_int);
                }
            }
            (*(*work).info).rho_estimate = compute_rho_estimate(work);
            (*(*work).info).solve_time = osqp_toc((*work).timer);
            if (*(*work).settings).polish != 0
                && (*(*work).info).status_val == OSQP_SOLVED as libc::c_longlong
            {
                polish(work);
            }
            if (*work).first_run != 0 {
                (*(*work).info)
                    .run_time = (*(*work).info).setup_time + (*(*work).info).solve_time
                    + (*(*work).info).polish_time;
            } else {
                (*(*work).info)
                    .run_time = (*(*work).info).update_time + (*(*work).info).solve_time
                    + (*(*work).info).polish_time;
            }
            if (*work).first_run != 0 {
                (*work).first_run = 0 as libc::c_int as c_int;
            }
            (*work).clear_update_time = 1 as libc::c_int as c_int;
            (*work).rho_update_from_solve = 0 as libc::c_int as c_int;
            if (*(*work).settings).verbose != 0 {
                print_footer((*work).info, (*(*work).settings).polish);
            }
            store_solution(work);
        }
        _ => {}
    }
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_cleanup(mut work: *mut OSQPWorkspace) -> c_int {
    let mut exitflag: c_int = 0 as libc::c_int as c_int;
    if !work.is_null() {
        if !((*work).data).is_null() {
            if !((*(*work).data).P).is_null() {
                csc_spfree((*(*work).data).P);
            }
            if !((*(*work).data).A).is_null() {
                csc_spfree((*(*work).data).A);
            }
            if !((*(*work).data).q).is_null() {
                free((*(*work).data).q as *mut libc::c_void);
            }
            if !((*(*work).data).l).is_null() {
                free((*(*work).data).l as *mut libc::c_void);
            }
            if !((*(*work).data).u).is_null() {
                free((*(*work).data).u as *mut libc::c_void);
            }
            free((*work).data as *mut libc::c_void);
        }
        if !((*work).scaling).is_null() {
            if !((*(*work).scaling).D).is_null() {
                free((*(*work).scaling).D as *mut libc::c_void);
            }
            if !((*(*work).scaling).Dinv).is_null() {
                free((*(*work).scaling).Dinv as *mut libc::c_void);
            }
            if !((*(*work).scaling).E).is_null() {
                free((*(*work).scaling).E as *mut libc::c_void);
            }
            if !((*(*work).scaling).Einv).is_null() {
                free((*(*work).scaling).Einv as *mut libc::c_void);
            }
            free((*work).scaling as *mut libc::c_void);
        }
        if !((*work).D_temp).is_null() {
            free((*work).D_temp as *mut libc::c_void);
        }
        if !((*work).D_temp_A).is_null() {
            free((*work).D_temp_A as *mut libc::c_void);
        }
        if !((*work).E_temp).is_null() {
            free((*work).E_temp as *mut libc::c_void);
        }
        if !((*work).linsys_solver).is_null() {
            if ((*(*work).linsys_solver).free).is_some() {
                ((*(*work).linsys_solver).free)
                    .expect("non-null function pointer")((*work).linsys_solver);
            }
        }
        if !((*work).settings).is_null() {
            exitflag = unload_linsys_solver((*(*work).settings).linsys_solver);
        }
        if !((*work).pol).is_null() {
            if !((*(*work).pol).Alow_to_A).is_null() {
                free((*(*work).pol).Alow_to_A as *mut libc::c_void);
            }
            if !((*(*work).pol).Aupp_to_A).is_null() {
                free((*(*work).pol).Aupp_to_A as *mut libc::c_void);
            }
            if !((*(*work).pol).A_to_Alow).is_null() {
                free((*(*work).pol).A_to_Alow as *mut libc::c_void);
            }
            if !((*(*work).pol).A_to_Aupp).is_null() {
                free((*(*work).pol).A_to_Aupp as *mut libc::c_void);
            }
            if !((*(*work).pol).x).is_null() {
                free((*(*work).pol).x as *mut libc::c_void);
            }
            if !((*(*work).pol).z).is_null() {
                free((*(*work).pol).z as *mut libc::c_void);
            }
            if !((*(*work).pol).y).is_null() {
                free((*(*work).pol).y as *mut libc::c_void);
            }
            free((*work).pol as *mut libc::c_void);
        }
        if !((*work).rho_vec).is_null() {
            free((*work).rho_vec as *mut libc::c_void);
        }
        if !((*work).rho_inv_vec).is_null() {
            free((*work).rho_inv_vec as *mut libc::c_void);
        }
        if !((*work).constr_type).is_null() {
            free((*work).constr_type as *mut libc::c_void);
        }
        if !((*work).x).is_null() {
            free((*work).x as *mut libc::c_void);
        }
        if !((*work).z).is_null() {
            free((*work).z as *mut libc::c_void);
        }
        if !((*work).xz_tilde).is_null() {
            free((*work).xz_tilde as *mut libc::c_void);
        }
        if !((*work).x_prev).is_null() {
            free((*work).x_prev as *mut libc::c_void);
        }
        if !((*work).z_prev).is_null() {
            free((*work).z_prev as *mut libc::c_void);
        }
        if !((*work).y).is_null() {
            free((*work).y as *mut libc::c_void);
        }
        if !((*work).Ax).is_null() {
            free((*work).Ax as *mut libc::c_void);
        }
        if !((*work).Px).is_null() {
            free((*work).Px as *mut libc::c_void);
        }
        if !((*work).Aty).is_null() {
            free((*work).Aty as *mut libc::c_void);
        }
        if !((*work).delta_y).is_null() {
            free((*work).delta_y as *mut libc::c_void);
        }
        if !((*work).Atdelta_y).is_null() {
            free((*work).Atdelta_y as *mut libc::c_void);
        }
        if !((*work).delta_x).is_null() {
            free((*work).delta_x as *mut libc::c_void);
        }
        if !((*work).Pdelta_x).is_null() {
            free((*work).Pdelta_x as *mut libc::c_void);
        }
        if !((*work).Adelta_x).is_null() {
            free((*work).Adelta_x as *mut libc::c_void);
        }
        if !((*work).settings).is_null() {
            free((*work).settings as *mut libc::c_void);
        }
        if !((*work).solution).is_null() {
            if !((*(*work).solution).x).is_null() {
                free((*(*work).solution).x as *mut libc::c_void);
            }
            if !((*(*work).solution).y).is_null() {
                free((*(*work).solution).y as *mut libc::c_void);
            }
            free((*work).solution as *mut libc::c_void);
        }
        if !((*work).info).is_null() {
            free((*work).info as *mut libc::c_void);
        }
        if !((*work).timer).is_null() {
            free((*work).timer as *mut libc::c_void);
        }
        free(work as *mut libc::c_void);
    }
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_lin_cost(
    mut work: *mut OSQPWorkspace,
    mut q_new: *const c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 21],
                &[libc::c_char; 21],
            >(b"osqp_update_lin_cost\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    prea_vec_copy(q_new, (*(*work).data).q, (*(*work).data).n);
    if (*(*work).settings).scaling != 0 {
        vec_ew_prod(
            (*(*work).scaling).D,
            (*(*work).data).q,
            (*(*work).data).q,
            (*(*work).data).n,
        );
        vec_mult_scalar((*(*work).data).q, (*(*work).scaling).c, (*(*work).data).n);
    }
    reset_info((*work).info);
    let ref mut fresh49 = (*(*work).info).update_time;
    *fresh49 += osqp_toc((*work).timer);
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_bounds(
    mut work: *mut OSQPWorkspace,
    mut l_new: *const c_float,
    mut u_new: *const c_float,
) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0 as libc::c_int as c_int;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 19],
                &[libc::c_char; 19],
            >(b"osqp_update_bounds\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    i = 0 as libc::c_int as c_int;
    while i < (*(*work).data).m {
        if *l_new.offset(i as isize) > *u_new.offset(i as isize) {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 19],
                    &[libc::c_char; 19],
                >(b"osqp_update_bounds\0"))
                    .as_ptr(),
            );
            printf(
                b"lower bound must be lower than or equal to upper bound\0" as *const u8
                    as *const libc::c_char,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 1 as libc::c_int as c_int;
        }
        i += 1;
    }
    prea_vec_copy(l_new, (*(*work).data).l, (*(*work).data).m);
    prea_vec_copy(u_new, (*(*work).data).u, (*(*work).data).m);
    if (*(*work).settings).scaling != 0 {
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
    }
    reset_info((*work).info);
    exitflag = update_rho_vec(work);
    let ref mut fresh50 = (*(*work).info).update_time;
    *fresh50 += osqp_toc((*work).timer);
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_lower_bound(
    mut work: *mut OSQPWorkspace,
    mut l_new: *const c_float,
) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0 as libc::c_int as c_int;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 24],
                &[libc::c_char; 24],
            >(b"osqp_update_lower_bound\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    prea_vec_copy(l_new, (*(*work).data).l, (*(*work).data).m);
    if (*(*work).settings).scaling != 0 {
        vec_ew_prod(
            (*(*work).scaling).E,
            (*(*work).data).l,
            (*(*work).data).l,
            (*(*work).data).m,
        );
    }
    i = 0 as libc::c_int as c_int;
    while i < (*(*work).data).m {
        if *((*(*work).data).l).offset(i as isize)
            > *((*(*work).data).u).offset(i as isize)
        {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 24],
                    &[libc::c_char; 24],
                >(b"osqp_update_lower_bound\0"))
                    .as_ptr(),
            );
            printf(
                b"upper bound must be greater than or equal to lower bound\0"
                    as *const u8 as *const libc::c_char,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 1 as libc::c_int as c_int;
        }
        i += 1;
    }
    reset_info((*work).info);
    exitflag = update_rho_vec(work);
    let ref mut fresh51 = (*(*work).info).update_time;
    *fresh51 += osqp_toc((*work).timer);
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_upper_bound(
    mut work: *mut OSQPWorkspace,
    mut u_new: *const c_float,
) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0 as libc::c_int as c_int;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 24],
                &[libc::c_char; 24],
            >(b"osqp_update_upper_bound\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    prea_vec_copy(u_new, (*(*work).data).u, (*(*work).data).m);
    if (*(*work).settings).scaling != 0 {
        vec_ew_prod(
            (*(*work).scaling).E,
            (*(*work).data).u,
            (*(*work).data).u,
            (*(*work).data).m,
        );
    }
    i = 0 as libc::c_int as c_int;
    while i < (*(*work).data).m {
        if *((*(*work).data).u).offset(i as isize)
            < *((*(*work).data).l).offset(i as isize)
        {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 24],
                    &[libc::c_char; 24],
                >(b"osqp_update_upper_bound\0"))
                    .as_ptr(),
            );
            printf(
                b"lower bound must be lower than or equal to upper bound\0" as *const u8
                    as *const libc::c_char,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 1 as libc::c_int as c_int;
        }
        i += 1;
    }
    reset_info((*work).info);
    exitflag = update_rho_vec(work);
    let ref mut fresh52 = (*(*work).info).update_time;
    *fresh52 += osqp_toc((*work).timer);
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_warm_start(
    mut work: *mut OSQPWorkspace,
    mut x: *const c_float,
    mut y: *const c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 16],
                &[libc::c_char; 16],
            >(b"osqp_warm_start\0"))
                .as_ptr(),
        );
    }
    if (*(*work).settings).warm_start == 0 {
        (*(*work).settings).warm_start = 1 as libc::c_int as c_int;
    }
    prea_vec_copy(x, (*work).x, (*(*work).data).n);
    prea_vec_copy(y, (*work).y, (*(*work).data).m);
    if (*(*work).settings).scaling != 0 {
        vec_ew_prod((*(*work).scaling).Dinv, (*work).x, (*work).x, (*(*work).data).n);
        vec_ew_prod((*(*work).scaling).Einv, (*work).y, (*work).y, (*(*work).data).m);
        vec_mult_scalar((*work).y, (*(*work).scaling).c, (*(*work).data).m);
    }
    mat_vec((*(*work).data).A, (*work).x, (*work).z, 0 as libc::c_int as c_int);
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_warm_start_x(
    mut work: *mut OSQPWorkspace,
    mut x: *const c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[libc::c_char; 18],
            >(b"osqp_warm_start_x\0"))
                .as_ptr(),
        );
    }
    if (*(*work).settings).warm_start == 0 {
        (*(*work).settings).warm_start = 1 as libc::c_int as c_int;
    }
    prea_vec_copy(x, (*work).x, (*(*work).data).n);
    if (*(*work).settings).scaling != 0 {
        vec_ew_prod((*(*work).scaling).Dinv, (*work).x, (*work).x, (*(*work).data).n);
    }
    mat_vec((*(*work).data).A, (*work).x, (*work).z, 0 as libc::c_int as c_int);
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_warm_start_y(
    mut work: *mut OSQPWorkspace,
    mut y: *const c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[libc::c_char; 18],
            >(b"osqp_warm_start_y\0"))
                .as_ptr(),
        );
    }
    if (*(*work).settings).warm_start == 0 {
        (*(*work).settings).warm_start = 1 as libc::c_int as c_int;
    }
    prea_vec_copy(y, (*work).y, (*(*work).data).m);
    if (*(*work).settings).scaling != 0 {
        vec_ew_prod((*(*work).scaling).Einv, (*work).y, (*work).y, (*(*work).data).m);
        vec_mult_scalar((*work).y, (*(*work).scaling).c, (*(*work).data).m);
    }
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_P(
    mut work: *mut OSQPWorkspace,
    mut Px_new: *const c_float,
    mut Px_new_idx: *const c_int,
    mut P_new_n: c_int,
) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0;
    let mut nnzP: c_int = 0;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[libc::c_char; 14],
            >(b"osqp_update_P\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    nnzP = *((*(*(*work).data).P).p).offset((*(*(*work).data).P).n as isize);
    if !Px_new_idx.is_null() {
        if P_new_n > nnzP {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 14],
                    &[libc::c_char; 14],
                >(b"osqp_update_P\0"))
                    .as_ptr(),
            );
            printf(
                b"new number of elements (%i) greater than elements in P (%i)\0"
                    as *const u8 as *const libc::c_char,
                P_new_n as libc::c_int,
                nnzP as libc::c_int,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 1 as libc::c_int as c_int;
        }
    }
    if (*(*work).settings).scaling != 0 {
        unscale_data(work);
    }
    if !Px_new_idx.is_null() {
        i = 0 as libc::c_int as c_int;
        while i < P_new_n {
            *((*(*(*work).data).P).x)
                .offset(
                    *Px_new_idx.offset(i as isize) as isize,
                ) = *Px_new.offset(i as isize);
            i += 1;
        }
    } else {
        i = 0 as libc::c_int as c_int;
        while i < nnzP {
            *((*(*(*work).data).P).x).offset(i as isize) = *Px_new.offset(i as isize);
            i += 1;
        }
    }
    if (*(*work).settings).scaling != 0 {
        scale_data(work);
    }
    exitflag = ((*(*work).linsys_solver).update_matrices)
        .expect(
            "non-null function pointer",
        )((*work).linsys_solver, (*(*work).data).P, (*(*work).data).A);
    reset_info((*work).info);
    if exitflag < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[libc::c_char; 14],
            >(b"osqp_update_P\0"))
                .as_ptr(),
        );
        printf(
            b"new KKT matrix is not quasidefinite\0" as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
    }
    let ref mut fresh53 = (*(*work).info).update_time;
    *fresh53 += osqp_toc((*work).timer);
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_A(
    mut work: *mut OSQPWorkspace,
    mut Ax_new: *const c_float,
    mut Ax_new_idx: *const c_int,
    mut A_new_n: c_int,
) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0;
    let mut nnzA: c_int = 0;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[libc::c_char; 14],
            >(b"osqp_update_A\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    nnzA = *((*(*(*work).data).A).p).offset((*(*(*work).data).A).n as isize);
    if !Ax_new_idx.is_null() {
        if A_new_n > nnzA {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 14],
                    &[libc::c_char; 14],
                >(b"osqp_update_A\0"))
                    .as_ptr(),
            );
            printf(
                b"new number of elements (%i) greater than elements in A (%i)\0"
                    as *const u8 as *const libc::c_char,
                A_new_n as libc::c_int,
                nnzA as libc::c_int,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 1 as libc::c_int as c_int;
        }
    }
    if (*(*work).settings).scaling != 0 {
        unscale_data(work);
    }
    if !Ax_new_idx.is_null() {
        i = 0 as libc::c_int as c_int;
        while i < A_new_n {
            *((*(*(*work).data).A).x)
                .offset(
                    *Ax_new_idx.offset(i as isize) as isize,
                ) = *Ax_new.offset(i as isize);
            i += 1;
        }
    } else {
        i = 0 as libc::c_int as c_int;
        while i < nnzA {
            *((*(*(*work).data).A).x).offset(i as isize) = *Ax_new.offset(i as isize);
            i += 1;
        }
    }
    if (*(*work).settings).scaling != 0 {
        scale_data(work);
    }
    exitflag = ((*(*work).linsys_solver).update_matrices)
        .expect(
            "non-null function pointer",
        )((*work).linsys_solver, (*(*work).data).P, (*(*work).data).A);
    reset_info((*work).info);
    if exitflag < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 14],
                &[libc::c_char; 14],
            >(b"osqp_update_A\0"))
                .as_ptr(),
        );
        printf(
            b"new KKT matrix is not quasidefinite\0" as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
    }
    let ref mut fresh54 = (*(*work).info).update_time;
    *fresh54 += osqp_toc((*work).timer);
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_P_A(
    mut work: *mut OSQPWorkspace,
    mut Px_new: *const c_float,
    mut Px_new_idx: *const c_int,
    mut P_new_n: c_int,
    mut Ax_new: *const c_float,
    mut Ax_new_idx: *const c_int,
    mut A_new_n: c_int,
) -> c_int {
    let mut i: c_int = 0;
    let mut exitflag: c_int = 0;
    let mut nnzP: c_int = 0;
    let mut nnzA: c_int = 0;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 16],
                &[libc::c_char; 16],
            >(b"osqp_update_P_A\0"))
                .as_ptr(),
        );
    }
    if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
        (*work).clear_update_time = 0 as libc::c_int as c_int;
        (*(*work).info).update_time = 0.0f64;
    }
    osqp_tic((*work).timer);
    nnzP = *((*(*(*work).data).P).p).offset((*(*(*work).data).P).n as isize);
    nnzA = *((*(*(*work).data).A).p).offset((*(*(*work).data).A).n as isize);
    if !Px_new_idx.is_null() {
        if P_new_n > nnzP {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 16],
                    &[libc::c_char; 16],
                >(b"osqp_update_P_A\0"))
                    .as_ptr(),
            );
            printf(
                b"new number of elements (%i) greater than elements in P (%i)\0"
                    as *const u8 as *const libc::c_char,
                P_new_n as libc::c_int,
                nnzP as libc::c_int,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 1 as libc::c_int as c_int;
        }
    }
    if !Ax_new_idx.is_null() {
        if A_new_n > nnzA {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 16],
                    &[libc::c_char; 16],
                >(b"osqp_update_P_A\0"))
                    .as_ptr(),
            );
            printf(
                b"new number of elements (%i) greater than elements in A (%i)\0"
                    as *const u8 as *const libc::c_char,
                A_new_n as libc::c_int,
                nnzA as libc::c_int,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return 2 as libc::c_int as c_int;
        }
    }
    if (*(*work).settings).scaling != 0 {
        unscale_data(work);
    }
    if !Px_new_idx.is_null() {
        i = 0 as libc::c_int as c_int;
        while i < P_new_n {
            *((*(*(*work).data).P).x)
                .offset(
                    *Px_new_idx.offset(i as isize) as isize,
                ) = *Px_new.offset(i as isize);
            i += 1;
        }
    } else {
        i = 0 as libc::c_int as c_int;
        while i < nnzP {
            *((*(*(*work).data).P).x).offset(i as isize) = *Px_new.offset(i as isize);
            i += 1;
        }
    }
    if !Ax_new_idx.is_null() {
        i = 0 as libc::c_int as c_int;
        while i < A_new_n {
            *((*(*(*work).data).A).x)
                .offset(
                    *Ax_new_idx.offset(i as isize) as isize,
                ) = *Ax_new.offset(i as isize);
            i += 1;
        }
    } else {
        i = 0 as libc::c_int as c_int;
        while i < nnzA {
            *((*(*(*work).data).A).x).offset(i as isize) = *Ax_new.offset(i as isize);
            i += 1;
        }
    }
    if (*(*work).settings).scaling != 0 {
        scale_data(work);
    }
    exitflag = ((*(*work).linsys_solver).update_matrices)
        .expect(
            "non-null function pointer",
        )((*work).linsys_solver, (*(*work).data).P, (*(*work).data).A);
    reset_info((*work).info);
    if exitflag < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 16],
                &[libc::c_char; 16],
            >(b"osqp_update_P_A\0"))
                .as_ptr(),
        );
        printf(
            b"new KKT matrix is not quasidefinite\0" as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
    }
    let ref mut fresh55 = (*(*work).info).update_time;
    *fresh55 += osqp_toc((*work).timer);
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_rho(
    mut work: *mut OSQPWorkspace,
    mut rho_new: c_float,
) -> c_int {
    let mut exitflag: c_int = 0;
    let mut i: c_int = 0;
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 16],
                &[libc::c_char; 16],
            >(b"osqp_update_rho\0"))
                .as_ptr(),
        );
    }
    if rho_new <= 0 as libc::c_int as libc::c_double {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 16],
                &[libc::c_char; 16],
            >(b"osqp_update_rho\0"))
                .as_ptr(),
        );
        printf(b"rho must be positive\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    if (*work).rho_update_from_solve == 0 as libc::c_int as libc::c_longlong {
        if (*work).clear_update_time == 1 as libc::c_int as libc::c_longlong {
            (*work).clear_update_time = 0 as libc::c_int as c_int;
            (*(*work).info).update_time = 0.0f64;
        }
        osqp_tic((*work).timer);
    }
    (*(*work).settings)
        .rho = if (if rho_new > 1e-06f64 { rho_new } else { 1e-06f64 }) < 1e06f64 {
        if rho_new > 1e-06f64 { rho_new } else { 1e-06f64 }
    } else {
        1e06f64
    };
    i = 0 as libc::c_int as c_int;
    while i < (*(*work).data).m {
        if *((*work).constr_type).offset(i as isize)
            == 0 as libc::c_int as libc::c_longlong
        {
            *((*work).rho_vec).offset(i as isize) = (*(*work).settings).rho;
            *((*work).rho_inv_vec).offset(i as isize) = 1.0f64 / (*(*work).settings).rho;
        } else if *((*work).constr_type).offset(i as isize)
                == 1 as libc::c_int as libc::c_longlong
            {
            *((*work).rho_vec)
                .offset(i as isize) = RHO_EQ_OVER_RHO_INEQ * (*(*work).settings).rho;
            *((*work).rho_inv_vec)
                .offset(i as isize) = 1.0f64 / *((*work).rho_vec).offset(i as isize);
        }
        i += 1;
    }
    exitflag = ((*(*work).linsys_solver).update_rho_vec)
        .expect("non-null function pointer")((*work).linsys_solver, (*work).rho_vec);
    if (*work).rho_update_from_solve == 0 as libc::c_int as libc::c_longlong {
        let ref mut fresh56 = (*(*work).info).update_time;
        *fresh56 += osqp_toc((*work).timer);
    }
    return exitflag;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_max_iter(
    mut work: *mut OSQPWorkspace,
    mut max_iter_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 21],
                &[libc::c_char; 21],
            >(b"osqp_update_max_iter\0"))
                .as_ptr(),
        );
    }
    if max_iter_new <= 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 21],
                &[libc::c_char; 21],
            >(b"osqp_update_max_iter\0"))
                .as_ptr(),
        );
        printf(b"max_iter must be positive\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).max_iter = max_iter_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_eps_abs(
    mut work: *mut OSQPWorkspace,
    mut eps_abs_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 20],
                &[libc::c_char; 20],
            >(b"osqp_update_eps_abs\0"))
                .as_ptr(),
        );
    }
    if eps_abs_new < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 20],
                &[libc::c_char; 20],
            >(b"osqp_update_eps_abs\0"))
                .as_ptr(),
        );
        printf(b"eps_abs must be nonnegative\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).eps_abs = eps_abs_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_eps_rel(
    mut work: *mut OSQPWorkspace,
    mut eps_rel_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 20],
                &[libc::c_char; 20],
            >(b"osqp_update_eps_rel\0"))
                .as_ptr(),
        );
    }
    if eps_rel_new < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 20],
                &[libc::c_char; 20],
            >(b"osqp_update_eps_rel\0"))
                .as_ptr(),
        );
        printf(b"eps_rel must be nonnegative\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).eps_rel = eps_rel_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_eps_prim_inf(
    mut work: *mut OSQPWorkspace,
    mut eps_prim_inf_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 25],
                &[libc::c_char; 25],
            >(b"osqp_update_eps_prim_inf\0"))
                .as_ptr(),
        );
    }
    if eps_prim_inf_new < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 25],
                &[libc::c_char; 25],
            >(b"osqp_update_eps_prim_inf\0"))
                .as_ptr(),
        );
        printf(
            b"eps_prim_inf must be nonnegative\0" as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).eps_prim_inf = eps_prim_inf_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_eps_dual_inf(
    mut work: *mut OSQPWorkspace,
    mut eps_dual_inf_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 25],
                &[libc::c_char; 25],
            >(b"osqp_update_eps_dual_inf\0"))
                .as_ptr(),
        );
    }
    if eps_dual_inf_new < 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 25],
                &[libc::c_char; 25],
            >(b"osqp_update_eps_dual_inf\0"))
                .as_ptr(),
        );
        printf(
            b"eps_dual_inf must be nonnegative\0" as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).eps_dual_inf = eps_dual_inf_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_alpha(
    mut work: *mut OSQPWorkspace,
    mut alpha_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[libc::c_char; 18],
            >(b"osqp_update_alpha\0"))
                .as_ptr(),
        );
    }
    if alpha_new <= 0.0f64 || alpha_new >= 2.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[libc::c_char; 18],
            >(b"osqp_update_alpha\0"))
                .as_ptr(),
        );
        printf(b"alpha must be between 0 and 2\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).alpha = alpha_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_warm_start(
    mut work: *mut OSQPWorkspace,
    mut warm_start_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 23],
                &[libc::c_char; 23],
            >(b"osqp_update_warm_start\0"))
                .as_ptr(),
        );
    }
    if warm_start_new != 0 as libc::c_int as libc::c_longlong
        && warm_start_new != 1 as libc::c_int as libc::c_longlong
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 23],
                &[libc::c_char; 23],
            >(b"osqp_update_warm_start\0"))
                .as_ptr(),
        );
        printf(
            b"warm_start should be either 0 or 1\0" as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).warm_start = warm_start_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_scaled_termination(
    mut work: *mut OSQPWorkspace,
    mut scaled_termination_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 31],
                &[libc::c_char; 31],
            >(b"osqp_update_scaled_termination\0"))
                .as_ptr(),
        );
    }
    if scaled_termination_new != 0 as libc::c_int as libc::c_longlong
        && scaled_termination_new != 1 as libc::c_int as libc::c_longlong
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 31],
                &[libc::c_char; 31],
            >(b"osqp_update_scaled_termination\0"))
                .as_ptr(),
        );
        printf(
            b"scaled_termination should be either 0 or 1\0" as *const u8
                as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).scaled_termination = scaled_termination_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_check_termination(
    mut work: *mut OSQPWorkspace,
    mut check_termination_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 30],
                &[libc::c_char; 30],
            >(b"osqp_update_check_termination\0"))
                .as_ptr(),
        );
    }
    if check_termination_new < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 30],
                &[libc::c_char; 30],
            >(b"osqp_update_check_termination\0"))
                .as_ptr(),
        );
        printf(
            b"check_termination should be nonnegative\0" as *const u8
                as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).check_termination = check_termination_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_delta(
    mut work: *mut OSQPWorkspace,
    mut delta_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[libc::c_char; 18],
            >(b"osqp_update_delta\0"))
                .as_ptr(),
        );
    }
    if delta_new <= 0.0f64 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 18],
                &[libc::c_char; 18],
            >(b"osqp_update_delta\0"))
                .as_ptr(),
        );
        printf(b"delta must be positive\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).delta = delta_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_polish(
    mut work: *mut OSQPWorkspace,
    mut polish_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 19],
                &[libc::c_char; 19],
            >(b"osqp_update_polish\0"))
                .as_ptr(),
        );
    }
    if polish_new != 0 as libc::c_int as libc::c_longlong
        && polish_new != 1 as libc::c_int as libc::c_longlong
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 19],
                &[libc::c_char; 19],
            >(b"osqp_update_polish\0"))
                .as_ptr(),
        );
        printf(b"polish should be either 0 or 1\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).polish = polish_new;
    (*(*work).info).polish_time = 0.0f64;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_polish_refine_iter(
    mut work: *mut OSQPWorkspace,
    mut polish_refine_iter_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 31],
                &[libc::c_char; 31],
            >(b"osqp_update_polish_refine_iter\0"))
                .as_ptr(),
        );
    }
    if polish_refine_iter_new < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 31],
                &[libc::c_char; 31],
            >(b"osqp_update_polish_refine_iter\0"))
                .as_ptr(),
        );
        printf(
            b"polish_refine_iter must be nonnegative\0" as *const u8
                as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).polish_refine_iter = polish_refine_iter_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_verbose(
    mut work: *mut OSQPWorkspace,
    mut verbose_new: c_int,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 20],
                &[libc::c_char; 20],
            >(b"osqp_update_verbose\0"))
                .as_ptr(),
        );
    }
    if verbose_new != 0 as libc::c_int as libc::c_longlong
        && verbose_new != 1 as libc::c_int as libc::c_longlong
    {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 20],
                &[libc::c_char; 20],
            >(b"osqp_update_verbose\0"))
                .as_ptr(),
        );
        printf(b"verbose should be either 0 or 1\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).verbose = verbose_new;
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_update_time_limit(
    mut work: *mut OSQPWorkspace,
    mut time_limit_new: c_float,
) -> c_int {
    if work.is_null() {
        return _osqp_error(
            OSQP_WORKSPACE_NOT_INIT_ERROR,
            (*::std::mem::transmute::<
                &[u8; 23],
                &[libc::c_char; 23],
            >(b"osqp_update_time_limit\0"))
                .as_ptr(),
        );
    }
    if time_limit_new < 0.0f64 {
        printf(
            b"time_limit must be nonnegative\n\0" as *const u8 as *const libc::c_char,
        );
        return 1 as libc::c_int as c_int;
    }
    (*(*work).settings).time_limit = time_limit_new;
    return 0 as libc::c_int as c_int;
}
