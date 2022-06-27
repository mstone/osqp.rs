extern "C" {
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn printf(_: *const ::std::os::raw::c_char, _: ...) -> ::std::os::raw::c_int;
    static mut LINSYS_SOLVER_NAME: [*const ::std::os::raw::c_char; 0];
    fn mach_absolute_time() -> uint64_t;
    fn mach_timebase_info(info: mach_timebase_info_t) -> kern_return_t;
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
pub type kern_return_t = ::std::os::raw::c_int;
pub type mach_timebase_info_t = *mut mach_timebase_info;
pub const OSQP_VERSION: [::std::os::raw::c_char; 6] = unsafe {
    *::std::mem::transmute::<&[u8; 6], &[::std::os::raw::c_char; 6]>(b"0.6.2\0")
};
pub const c_malloc: unsafe extern "C" fn(::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void = malloc;
pub const OSQP_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const c_print: unsafe extern "C" fn(*const ::std::os::raw::c_char, ...) -> ::std::os::raw::c_int = printf;
pub const OSQP_SOLVED_INACCURATE: ::std::os::raw::c_int = 2 as ::std::os::raw::c_int;
pub const OSQP_SOLVED: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
#[no_mangle]
pub unsafe extern "C" fn osqp_version() -> *const ::std::os::raw::c_char {
    return OSQP_VERSION.as_ptr();
}
pub const HEADER_LINE_LEN: ::std::os::raw::c_int = 65 as ::std::os::raw::c_int;
#[no_mangle]
pub unsafe extern "C" fn c_strcpy(
    mut dest: *mut ::std::os::raw::c_char,
    mut source: *const ::std::os::raw::c_char,
) {
    let mut i: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
    loop {
        *dest.offset(i as isize) = *source.offset(i as isize);
        if *dest.offset(i as isize) as ::std::os::raw::c_int == '\u{0}' as i32 {
            break;
        }
        i += 1;
    };
}
unsafe extern "C" fn print_line() {
    let mut the_line: [::std::os::raw::c_char; 66] = [0; 66];
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < HEADER_LINE_LEN {
        the_line[i as usize] = '-' as i32 as ::std::os::raw::c_char;
        i += 1;
    }
    the_line[HEADER_LINE_LEN as usize] = '\u{0}' as i32 as ::std::os::raw::c_char;
    printf(b"%s\n\0" as *const u8 as *const ::std::os::raw::c_char, the_line.as_mut_ptr());
}
#[no_mangle]
pub unsafe extern "C" fn print_header() {
    printf(b"iter   \0" as *const u8 as *const ::std::os::raw::c_char);
    printf(
        b"objective    pri res    dua res    rho\0" as *const u8 as *const ::std::os::raw::c_char,
    );
    printf(b"        time\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn print_setup_header(mut work: *const OSQPWorkspace) {
    let mut data: *mut OSQPData = 0 as *mut OSQPData;
    let mut settings: *mut OSQPSettings = 0 as *mut OSQPSettings;
    let mut nnz: c_int = 0;
    data = (*work).data;
    settings = (*work).settings;
    nnz = *((*(*data).P).p).offset((*(*data).P).n as isize)
        + *((*(*data).A).p).offset((*(*data).A).n as isize);
    print_line();
    printf(
        b"           OSQP v%s  -  Operator Splitting QP Solver\n              (c) Bartolomeo Stellato,  Goran Banjac\n        University of Oxford  -  Stanford University 2021\n\0"
            as *const u8 as *const ::std::os::raw::c_char,
        OSQP_VERSION.as_ptr(),
    );
    print_line();
    printf(b"problem:  \0" as *const u8 as *const ::std::os::raw::c_char);
    printf(
        b"variables n = %i, constraints m = %i\n          \0" as *const u8
            as *const ::std::os::raw::c_char,
        (*data).n,
        (*data).m,
    );
    printf(b"nnz(P) + nnz(A) = %i\n\0" as *const u8 as *const ::std::os::raw::c_char, nnz);
    printf(b"settings: \0" as *const u8 as *const ::std::os::raw::c_char);
    printf(
        b"linear system solver = %s\0" as *const u8 as *const ::std::os::raw::c_char,
        *LINSYS_SOLVER_NAME.as_mut_ptr().offset((*settings).linsys_solver as isize),
    );
    if (*(*work).linsys_solver).nthreads != 1 as ::std::os::raw::c_int {
        printf(
            b" (%d threads)\0" as *const u8 as *const ::std::os::raw::c_char,
            (*(*work).linsys_solver).nthreads,
        );
    }
    printf(b",\n          \0" as *const u8 as *const ::std::os::raw::c_char);
    printf(
        b"eps_abs = %.1e, eps_rel = %.1e,\n          \0" as *const u8
            as *const ::std::os::raw::c_char,
        (*settings).eps_abs,
        (*settings).eps_rel,
    );
    printf(
        b"eps_prim_inf = %.1e, eps_dual_inf = %.1e,\n          \0" as *const u8
            as *const ::std::os::raw::c_char,
        (*settings).eps_prim_inf,
        (*settings).eps_dual_inf,
    );
    printf(b"rho = %.2e \0" as *const u8 as *const ::std::os::raw::c_char, (*settings).rho);
    if (*settings).adaptive_rho != 0 {
        printf(b"(adaptive)\0" as *const u8 as *const ::std::os::raw::c_char);
    }
    printf(b",\n          \0" as *const u8 as *const ::std::os::raw::c_char);
    printf(
        b"sigma = %.2e, alpha = %.2f, \0" as *const u8 as *const ::std::os::raw::c_char,
        (*settings).sigma,
        (*settings).alpha,
    );
    printf(
        b"max_iter = %i\n\0" as *const u8 as *const ::std::os::raw::c_char,
        (*settings).max_iter,
    );
    if (*settings).check_termination != 0 {
        printf(
            b"          check_termination: on (interval %i),\n\0" as *const u8
                as *const ::std::os::raw::c_char,
            (*settings).check_termination,
        );
    } else {
        printf(
            b"          check_termination: off,\n\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    }
    if (*settings).time_limit != 0. {
        printf(
            b"          time_limit: %.2e sec,\n\0" as *const u8 as *const ::std::os::raw::c_char,
            (*settings).time_limit,
        );
    }
    if (*settings).scaling != 0 {
        printf(b"          scaling: on, \0" as *const u8 as *const ::std::os::raw::c_char);
    } else {
        printf(b"          scaling: off, \0" as *const u8 as *const ::std::os::raw::c_char);
    }
    if (*settings).scaled_termination != 0 {
        printf(b"scaled_termination: on\n\0" as *const u8 as *const ::std::os::raw::c_char);
    } else {
        printf(b"scaled_termination: off\n\0" as *const u8 as *const ::std::os::raw::c_char);
    }
    if (*settings).warm_start != 0 {
        printf(b"          warm start: on, \0" as *const u8 as *const ::std::os::raw::c_char);
    } else {
        printf(b"          warm start: off, \0" as *const u8 as *const ::std::os::raw::c_char);
    }
    if (*settings).polish != 0 {
        printf(b"polish: on, \0" as *const u8 as *const ::std::os::raw::c_char);
    } else {
        printf(b"polish: off, \0" as *const u8 as *const ::std::os::raw::c_char);
    }
    if (*settings).time_limit != 0. {
        printf(
            b"time_limit: %.2e sec\n\0" as *const u8 as *const ::std::os::raw::c_char,
            (*settings).time_limit,
        );
    } else {
        printf(b"time_limit: off\n\0" as *const u8 as *const ::std::os::raw::c_char);
    }
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn print_summary(mut work: *mut OSQPWorkspace) {
    let mut info: *mut OSQPInfo = 0 as *mut OSQPInfo;
    info = (*work).info;
    printf(b"%4i\0" as *const u8 as *const ::std::os::raw::c_char, (*info).iter);
    printf(b" %12.4e\0" as *const u8 as *const ::std::os::raw::c_char, (*info).obj_val);
    printf(b"  %9.2e\0" as *const u8 as *const ::std::os::raw::c_char, (*info).pri_res);
    printf(b"  %9.2e\0" as *const u8 as *const ::std::os::raw::c_char, (*info).dua_res);
    printf(b"  %9.2e\0" as *const u8 as *const ::std::os::raw::c_char, (*(*work).settings).rho);
    if (*work).first_run != 0 {
        printf(
            b"  %9.2es\0" as *const u8 as *const ::std::os::raw::c_char,
            (*info).setup_time + (*info).solve_time,
        );
    } else {
        printf(
            b"  %9.2es\0" as *const u8 as *const ::std::os::raw::c_char,
            (*info).update_time + (*info).solve_time,
        );
    }
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
    (*work).summary_printed = 1 as ::std::os::raw::c_int;
}
#[no_mangle]
pub unsafe extern "C" fn print_polish(mut work: *mut OSQPWorkspace) {
    let mut info: *mut OSQPInfo = 0 as *mut OSQPInfo;
    info = (*work).info;
    printf(
        b"%4s\0" as *const u8 as *const ::std::os::raw::c_char,
        b"plsh\0" as *const u8 as *const ::std::os::raw::c_char,
    );
    printf(b" %12.4e\0" as *const u8 as *const ::std::os::raw::c_char, (*info).obj_val);
    printf(b"  %9.2e\0" as *const u8 as *const ::std::os::raw::c_char, (*info).pri_res);
    printf(b"  %9.2e\0" as *const u8 as *const ::std::os::raw::c_char, (*info).dua_res);
    printf(b"   --------\0" as *const u8 as *const ::std::os::raw::c_char);
    if (*work).first_run != 0 {
        printf(
            b"  %9.2es\0" as *const u8 as *const ::std::os::raw::c_char,
            (*info).setup_time + (*info).solve_time + (*info).polish_time,
        );
    } else {
        printf(
            b"  %9.2es\0" as *const u8 as *const ::std::os::raw::c_char,
            (*info).update_time + (*info).solve_time + (*info).polish_time,
        );
    }
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn print_footer(mut info: *mut OSQPInfo, mut polish: c_int) {
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
    printf(
        b"status:               %s\n\0" as *const u8 as *const ::std::os::raw::c_char,
        ((*info).status).as_mut_ptr(),
    );
    if polish != 0 && (*info).status_val == OSQP_SOLVED {
        if (*info).status_polish == 1 as ::std::os::raw::c_int {
            printf(
                b"solution polish:      successful\n\0" as *const u8
                    as *const ::std::os::raw::c_char,
            );
        } else if (*info).status_polish < 0 as ::std::os::raw::c_int {
            printf(
                b"solution polish:      unsuccessful\n\0" as *const u8
                    as *const ::std::os::raw::c_char,
            );
        }
    }
    printf(
        b"number of iterations: %i\n\0" as *const u8 as *const ::std::os::raw::c_char,
        (*info).iter,
    );
    if (*info).status_val == OSQP_SOLVED || (*info).status_val == OSQP_SOLVED_INACCURATE
    {
        printf(
            b"optimal objective:    %.4f\n\0" as *const u8 as *const ::std::os::raw::c_char,
            (*info).obj_val,
        );
    }
    printf(
        b"run time:             %.2es\n\0" as *const u8 as *const ::std::os::raw::c_char,
        (*info).run_time,
    );
    printf(
        b"optimal rho estimate: %.2e\n\0" as *const u8 as *const ::std::os::raw::c_char,
        (*info).rho_estimate,
    );
    printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
}
#[no_mangle]
pub unsafe extern "C" fn copy_settings(
    mut settings: *const OSQPSettings,
) -> *mut OSQPSettings {
    let mut new: *mut OSQPSettings = malloc(
        ::std::mem::size_of::<OSQPSettings>() as ::std::os::raw::c_ulong,
    ) as *mut OSQPSettings;
    if new.is_null() {
        return OSQP_NULL as *mut OSQPSettings;
    }
    (*new).rho = (*settings).rho;
    (*new).sigma = (*settings).sigma;
    (*new).scaling = (*settings).scaling;
    (*new).adaptive_rho = (*settings).adaptive_rho;
    (*new).adaptive_rho_interval = (*settings).adaptive_rho_interval;
    (*new).adaptive_rho_tolerance = (*settings).adaptive_rho_tolerance;
    (*new).adaptive_rho_fraction = (*settings).adaptive_rho_fraction;
    (*new).max_iter = (*settings).max_iter;
    (*new).eps_abs = (*settings).eps_abs;
    (*new).eps_rel = (*settings).eps_rel;
    (*new).eps_prim_inf = (*settings).eps_prim_inf;
    (*new).eps_dual_inf = (*settings).eps_dual_inf;
    (*new).alpha = (*settings).alpha;
    (*new).linsys_solver = (*settings).linsys_solver;
    (*new).delta = (*settings).delta;
    (*new).polish = (*settings).polish;
    (*new).polish_refine_iter = (*settings).polish_refine_iter;
    (*new).verbose = (*settings).verbose;
    (*new).scaled_termination = (*settings).scaled_termination;
    (*new).check_termination = (*settings).check_termination;
    (*new).warm_start = (*settings).warm_start;
    (*new).time_limit = (*settings).time_limit;
    return new;
}
#[no_mangle]
pub unsafe extern "C" fn osqp_tic(mut t: *mut OSQPTimer) {
    (*t).tic = mach_absolute_time();
}
#[no_mangle]
pub unsafe extern "C" fn osqp_toc(mut t: *mut OSQPTimer) -> c_float {
    let mut duration: uint64_t = 0;
    (*t).toc = mach_absolute_time();
    duration = ((*t).toc).wrapping_sub((*t).tic);
    mach_timebase_info(&mut (*t).tinfo);
    duration = (duration as ::std::os::raw::c_ulonglong)
        .wrapping_mul((*t).tinfo.numer as ::std::os::raw::c_ulonglong) as uint64_t as uint64_t;
    duration = (duration as ::std::os::raw::c_ulonglong)
        .wrapping_div((*t).tinfo.denom as ::std::os::raw::c_ulonglong) as uint64_t as uint64_t;
    return duration as c_float / 1e9f64;
}
