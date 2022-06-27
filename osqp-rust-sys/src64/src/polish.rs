extern "C" {
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn free(_: *mut ::std::os::raw::c_void);
    fn vec_copy(a: *mut c_float, n: c_int) -> *mut c_float;
    fn prea_vec_copy(a: *const c_float, b: *mut c_float, n: c_int);
    fn vec_set_scalar(a: *mut c_float, sc: c_float, n: c_int);
    fn int_vec_set_scalar(a: *mut c_int, sc: c_int, n: c_int);
    fn mat_vec(A: *const csc, x: *const c_float, y: *mut c_float, plus_eq: c_int);
    fn mat_tpose_vec(
        A: *const csc,
        x: *const c_float,
        y: *mut c_float,
        plus_eq: c_int,
        skip_diag: c_int,
    );
    fn print_polish(work: *mut OSQPWorkspace);
    fn osqp_tic(t: *mut OSQPTimer);
    fn update_info(
        work: *mut OSQPWorkspace,
        iter: c_int,
        compute_objective: c_int,
        polish_0: c_int,
    );
    fn init_linsys_solver(
        s: *mut *mut LinSysSolver,
        P: *const csc,
        A: *const csc,
        sigma: c_float,
        rho_vec: *const c_float,
        linsys_solver: linsys_solver_type,
        polish_0: c_int,
    ) -> c_int;
    fn csc_spalloc(
        m: c_int,
        n: c_int,
        nzmax: c_int,
        values: c_int,
        triplet: c_int,
    ) -> *mut csc;
    fn csc_spfree(A: *mut csc);
    fn project_normalcone(work: *mut OSQPWorkspace, z: *mut c_float, y: *mut c_float);
    fn _osqp_error(
        error_code: osqp_error_type,
        function_name: *const ::std::os::raw::c_char,
    ) -> c_int;
}
pub type c_int = ::std::os::raw::c_longlong;
pub type c_float = ::std::os::raw::c_double;
pub type linsys_solver_type = ::std::os::raw::c_uint;
pub const MKL_PARDISO_SOLVER: linsys_solver_type = 1;
pub const QDLDL_SOLVER: linsys_solver_type = 0;
pub type osqp_error_type = ::std::os::raw::c_uint;
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
pub const OSQP_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const c_malloc: unsafe extern "C" fn(::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void = malloc;
pub const c_free: unsafe extern "C" fn(*mut ::std::os::raw::c_void) -> () = free;
unsafe extern "C" fn form_Ared(mut work: *mut OSQPWorkspace) -> c_int {
    let mut j: c_int = 0;
    let mut ptr: c_int = 0;
    let mut Ared_nnz: c_int = 0 as ::std::os::raw::c_int as c_int;
    (*(*work).pol).n_low = 0 as ::std::os::raw::c_int as c_int;
    (*(*work).pol).n_upp = 0 as ::std::os::raw::c_int as c_int;
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).data).m {
        if *((*work).z).offset(j as isize) - *((*(*work).data).l).offset(j as isize)
            < -*((*work).y).offset(j as isize)
        {
            *((*(*work).pol).Alow_to_A).offset((*(*work).pol).n_low as isize) = j;
            let ref mut fresh0 = (*(*work).pol).n_low;
            let fresh1 = *fresh0;
            *fresh0 = *fresh0 + 1;
            *((*(*work).pol).A_to_Alow).offset(j as isize) = fresh1;
        } else {
            *((*(*work).pol).A_to_Alow)
                .offset(j as isize) = -(1 as ::std::os::raw::c_int) as c_int;
        }
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).data).m {
        if *((*(*work).data).u).offset(j as isize) - *((*work).z).offset(j as isize)
            < *((*work).y).offset(j as isize)
        {
            *((*(*work).pol).Aupp_to_A).offset((*(*work).pol).n_upp as isize) = j;
            let ref mut fresh2 = (*(*work).pol).n_upp;
            let fresh3 = *fresh2;
            *fresh2 = *fresh2 + 1;
            *((*(*work).pol).A_to_Aupp).offset(j as isize) = fresh3;
        } else {
            *((*(*work).pol).A_to_Aupp)
                .offset(j as isize) = -(1 as ::std::os::raw::c_int) as c_int;
        }
        j += 1;
    }
    if (*(*work).pol).n_low + (*(*work).pol).n_upp
        == 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
    {
        let ref mut fresh4 = (*(*work).pol).Ared;
        *fresh4 = csc_spalloc(
            0 as ::std::os::raw::c_int as c_int,
            (*(*work).data).n,
            0 as ::std::os::raw::c_int as c_int,
            1 as ::std::os::raw::c_int as c_int,
            0 as ::std::os::raw::c_int as c_int,
        );
        if ((*(*work).pol).Ared).is_null() {
            return -(1 as ::std::os::raw::c_int) as c_int;
        }
        int_vec_set_scalar(
            (*(*(*work).pol).Ared).p,
            0 as ::std::os::raw::c_int as c_int,
            (*(*work).data).n + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong,
        );
        return 0 as ::std::os::raw::c_int as c_int;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < *((*(*(*work).data).A).p).offset((*(*(*work).data).A).n as isize) {
        if *((*(*work).pol).A_to_Alow)
            .offset(*((*(*(*work).data).A).i).offset(j as isize) as isize)
            != -(1 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong
            || *((*(*work).pol).A_to_Aupp)
                .offset(*((*(*(*work).data).A).i).offset(j as isize) as isize)
                != -(1 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong
        {
            Ared_nnz += 1;
        }
        j += 1;
    }
    let ref mut fresh5 = (*(*work).pol).Ared;
    *fresh5 = csc_spalloc(
        (*(*work).pol).n_low + (*(*work).pol).n_upp,
        (*(*work).data).n,
        Ared_nnz,
        1 as ::std::os::raw::c_int as c_int,
        0 as ::std::os::raw::c_int as c_int,
    );
    if ((*(*work).pol).Ared).is_null() {
        return -(1 as ::std::os::raw::c_int) as c_int;
    }
    Ared_nnz = 0 as ::std::os::raw::c_int as c_int;
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).data).n {
        *((*(*(*work).pol).Ared).p).offset(j as isize) = Ared_nnz;
        ptr = *((*(*(*work).data).A).p).offset(j as isize);
        while ptr
            < *((*(*(*work).data).A).p)
                .offset((j + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong) as isize)
        {
            if *((*(*work).pol).A_to_Alow)
                .offset(*((*(*(*work).data).A).i).offset(ptr as isize) as isize)
                != -(1 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong
            {
                *((*(*(*work).pol).Ared).i)
                    .offset(
                        Ared_nnz as isize,
                    ) = *((*(*work).pol).A_to_Alow)
                    .offset(*((*(*(*work).data).A).i).offset(ptr as isize) as isize);
                let fresh6 = Ared_nnz;
                Ared_nnz = Ared_nnz + 1;
                *((*(*(*work).pol).Ared).x)
                    .offset(
                        fresh6 as isize,
                    ) = *((*(*(*work).data).A).x).offset(ptr as isize);
            } else if *((*(*work).pol).A_to_Aupp)
                    .offset(*((*(*(*work).data).A).i).offset(ptr as isize) as isize)
                    != -(1 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong
                {
                *((*(*(*work).pol).Ared).i)
                    .offset(
                        Ared_nnz as isize,
                    ) = *((*(*work).pol).A_to_Aupp)
                    .offset(*((*(*(*work).data).A).i).offset(ptr as isize) as isize)
                    + (*(*work).pol).n_low;
                let fresh7 = Ared_nnz;
                Ared_nnz = Ared_nnz + 1;
                *((*(*(*work).pol).Ared).x)
                    .offset(
                        fresh7 as isize,
                    ) = *((*(*(*work).data).A).x).offset(ptr as isize);
            }
            ptr += 1;
        }
        j += 1;
    }
    *((*(*(*work).pol).Ared).p).offset((*(*work).data).n as isize) = Ared_nnz;
    return (*(*work).pol).n_low + (*(*work).pol).n_upp;
}
unsafe extern "C" fn form_rhs_red(mut work: *mut OSQPWorkspace, mut rhs: *mut c_float) {
    let mut j: c_int = 0;
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).data).n {
        *rhs.offset(j as isize) = -*((*(*work).data).q).offset(j as isize);
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).pol).n_low {
        *rhs
            .offset(
                ((*(*work).data).n + j) as isize,
            ) = *((*(*work).data).l)
            .offset(*((*(*work).pol).Alow_to_A).offset(j as isize) as isize);
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).pol).n_upp {
        *rhs
            .offset(
                ((*(*work).data).n + (*(*work).pol).n_low + j) as isize,
            ) = *((*(*work).data).u)
            .offset(*((*(*work).pol).Aupp_to_A).offset(j as isize) as isize);
        j += 1;
    }
}
unsafe extern "C" fn iterative_refinement(
    mut work: *mut OSQPWorkspace,
    mut p: *mut LinSysSolver,
    mut z: *mut c_float,
    mut b: *mut c_float,
) -> c_int {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut n: c_int = 0;
    let mut rhs: *mut c_float = 0 as *mut c_float;
    if (*(*work).settings).polish_refine_iter > 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
        n = (*(*work).data).n + (*(*(*work).pol).Ared).m;
        rhs = malloc(
            (::std::mem::size_of::<c_float>() as ::std::os::raw::c_ulong as ::std::os::raw::c_ulonglong)
                .wrapping_mul(n as ::std::os::raw::c_ulonglong) as ::std::os::raw::c_ulong,
        ) as *mut c_float;
        if rhs.is_null() {
            return _osqp_error(
                OSQP_MEM_ALLOC_ERROR,
                (*::std::mem::transmute::<
                    &[u8; 21],
                    &[::std::os::raw::c_char; 21],
                >(b"iterative_refinement\0"))
                    .as_ptr(),
            )
        } else {
            i = 0 as ::std::os::raw::c_int as c_int;
            while i < (*(*work).settings).polish_refine_iter {
                prea_vec_copy(b, rhs, n);
                mat_vec((*(*work).data).P, z, rhs, -(1 as ::std::os::raw::c_int) as c_int);
                mat_tpose_vec(
                    (*(*work).data).P,
                    z,
                    rhs,
                    -(1 as ::std::os::raw::c_int) as c_int,
                    1 as ::std::os::raw::c_int as c_int,
                );
                mat_tpose_vec(
                    (*(*work).pol).Ared,
                    z.offset((*(*work).data).n as isize),
                    rhs,
                    -(1 as ::std::os::raw::c_int) as c_int,
                    0 as ::std::os::raw::c_int as c_int,
                );
                mat_vec(
                    (*(*work).pol).Ared,
                    z,
                    rhs.offset((*(*work).data).n as isize),
                    -(1 as ::std::os::raw::c_int) as c_int,
                );
                ((*p).solve).expect("non-null function pointer")(p, rhs);
                j = 0 as ::std::os::raw::c_int as c_int;
                while j < n {
                    let ref mut fresh8 = *z.offset(j as isize);
                    *fresh8 += *rhs.offset(j as isize);
                    j += 1;
                }
                i += 1;
            }
        }
        if !rhs.is_null() {
            free(rhs as *mut ::std::os::raw::c_void);
        }
    }
    return 0 as ::std::os::raw::c_int as c_int;
}
unsafe extern "C" fn get_ypol_from_yred(
    mut work: *mut OSQPWorkspace,
    mut yred: *mut c_float,
) {
    let mut j: c_int = 0;
    if (*(*work).pol).n_low + (*(*work).pol).n_upp
        == 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
    {
        vec_set_scalar((*(*work).pol).y, 0.0f64, (*(*work).data).m);
        return;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*(*work).data).m {
        if *((*(*work).pol).A_to_Alow).offset(j as isize)
            != -(1 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong
        {
            *((*(*work).pol).y)
                .offset(
                    j as isize,
                ) = *yred
                .offset(*((*(*work).pol).A_to_Alow).offset(j as isize) as isize);
        } else if *((*(*work).pol).A_to_Aupp).offset(j as isize)
                != -(1 as ::std::os::raw::c_int) as ::std::os::raw::c_longlong
            {
            *((*(*work).pol).y)
                .offset(
                    j as isize,
                ) = *yred
                .offset(
                    (*((*(*work).pol).A_to_Aupp).offset(j as isize)
                        + (*(*work).pol).n_low) as isize,
                );
        } else {
            *((*(*work).pol).y).offset(j as isize) = 0.0f64;
        }
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn polish(mut work: *mut OSQPWorkspace) -> c_int {
    let mut mred: c_int = 0;
    let mut polish_successful: c_int = 0;
    let mut exitflag: c_int = 0;
    let mut rhs_red: *mut c_float = 0 as *mut c_float;
    let mut plsh: *mut LinSysSolver = 0 as *mut LinSysSolver;
    let mut pol_sol: *mut c_float = 0 as *mut c_float;
    osqp_tic((*work).timer);
    mred = form_Ared(work);
    if mred < 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
        (*(*work).info).status_polish = -(1 as ::std::os::raw::c_int) as c_int;
        return -(1 as ::std::os::raw::c_int) as c_int;
    }
    exitflag = init_linsys_solver(
        &mut plsh,
        (*(*work).data).P,
        (*(*work).pol).Ared,
        (*(*work).settings).delta,
        OSQP_NULL as *const c_float,
        (*(*work).settings).linsys_solver,
        1 as ::std::os::raw::c_int as c_int,
    );
    if exitflag != 0 {
        (*(*work).info).status_polish = -(1 as ::std::os::raw::c_int) as c_int;
        if !((*(*work).pol).Ared).is_null() {
            csc_spfree((*(*work).pol).Ared);
        }
        return 1 as ::std::os::raw::c_int as c_int;
    }
    rhs_red = malloc(
        (::std::mem::size_of::<c_float>() as ::std::os::raw::c_ulong as ::std::os::raw::c_ulonglong)
            .wrapping_mul(((*(*work).data).n + mred) as ::std::os::raw::c_ulonglong)
            as ::std::os::raw::c_ulong,
    ) as *mut c_float;
    if rhs_red.is_null() {
        (*(*work).info).status_polish = -(1 as ::std::os::raw::c_int) as c_int;
        csc_spfree((*(*work).pol).Ared);
        return -(1 as ::std::os::raw::c_int) as c_int;
    }
    form_rhs_red(work, rhs_red);
    pol_sol = vec_copy(rhs_red, (*(*work).data).n + mred);
    if pol_sol.is_null() {
        (*(*work).info).status_polish = -(1 as ::std::os::raw::c_int) as c_int;
        csc_spfree((*(*work).pol).Ared);
        free(rhs_red as *mut ::std::os::raw::c_void);
        return -(1 as ::std::os::raw::c_int) as c_int;
    }
    ((*plsh).solve).expect("non-null function pointer")(plsh, pol_sol);
    exitflag = iterative_refinement(work, plsh, pol_sol, rhs_red);
    if exitflag != 0 {
        (*(*work).info).status_polish = -(1 as ::std::os::raw::c_int) as c_int;
        csc_spfree((*(*work).pol).Ared);
        free(rhs_red as *mut ::std::os::raw::c_void);
        free(pol_sol as *mut ::std::os::raw::c_void);
        return -(1 as ::std::os::raw::c_int) as c_int;
    }
    prea_vec_copy(pol_sol, (*(*work).pol).x, (*(*work).data).n);
    mat_vec(
        (*(*work).data).A,
        (*(*work).pol).x,
        (*(*work).pol).z,
        0 as ::std::os::raw::c_int as c_int,
    );
    get_ypol_from_yred(work, pol_sol.offset((*(*work).data).n as isize));
    project_normalcone(work, (*(*work).pol).z, (*(*work).pol).y);
    update_info(
        work,
        0 as ::std::os::raw::c_int as c_int,
        1 as ::std::os::raw::c_int as c_int,
        1 as ::std::os::raw::c_int as c_int,
    );
    polish_successful = ((*(*work).pol).pri_res < (*(*work).info).pri_res
        && (*(*work).pol).dua_res < (*(*work).info).dua_res
        || (*(*work).pol).pri_res < (*(*work).info).pri_res
            && (*(*work).info).dua_res < 1e-10f64
        || (*(*work).pol).dua_res < (*(*work).info).dua_res
            && (*(*work).info).pri_res < 1e-10f64) as ::std::os::raw::c_int as c_int;
    if polish_successful != 0 {
        (*(*work).info).obj_val = (*(*work).pol).obj_val;
        (*(*work).info).pri_res = (*(*work).pol).pri_res;
        (*(*work).info).dua_res = (*(*work).pol).dua_res;
        (*(*work).info).status_polish = 1 as ::std::os::raw::c_int as c_int;
        prea_vec_copy((*(*work).pol).x, (*work).x, (*(*work).data).n);
        prea_vec_copy((*(*work).pol).z, (*work).z, (*(*work).data).m);
        prea_vec_copy((*(*work).pol).y, (*work).y, (*(*work).data).m);
        if (*(*work).settings).verbose != 0 {
            print_polish(work);
        }
    } else {
        (*(*work).info).status_polish = -(1 as ::std::os::raw::c_int) as c_int;
    }
    ((*plsh).free).expect("non-null function pointer")(plsh);
    csc_spfree((*(*work).pol).Ared);
    free(rhs_red as *mut ::std::os::raw::c_void);
    free(pol_sol as *mut ::std::os::raw::c_void);
    return 0 as ::std::os::raw::c_int as c_int;
}
