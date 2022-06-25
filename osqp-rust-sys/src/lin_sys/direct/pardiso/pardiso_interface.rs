use ::libc;
extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
    fn update_KKT_param2(
        KKT: *mut csc,
        param2: *const c_float,
        param2toKKT: *const c_int,
        m: c_int,
    );
    fn update_KKT_A(KKT: *mut csc, A: *const csc, AtoKKT: *const c_int);
    fn update_KKT_P(
        KKT: *mut csc,
        P: *const csc,
        PtoKKT: *const c_int,
        param1: c_float,
        Pdiag_idx: *const c_int,
        Pdiag_n: c_int,
    );
    fn form_KKT(
        P: *const csc,
        A: *const csc,
        format: c_int,
        param1: c_float,
        param2: *mut c_float,
        PtoKKT: *mut c_int,
        AtoKKT: *mut c_int,
        Pdiag_idx: *mut *mut c_int,
        Pdiag_n: *mut c_int,
        param2toKKT: *mut c_int,
    ) -> *mut csc;
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
    fn csc_spfree(A: *mut csc);
    fn pardiso(
        _: *mut *mut libc::c_void,
        _: *const c_int,
        _: *const c_int,
        _: *const c_int,
        _: *const c_int,
        _: *const c_int,
        _: *const c_float,
        _: *const c_int,
        _: *const c_int,
        _: *mut c_int,
        _: *const c_int,
        _: *mut c_int,
        _: *const c_int,
        _: *mut c_float,
        _: *mut c_float,
        _: *mut c_int,
    );
    fn mkl_set_interface_layer(_: c_int) -> c_int;
    fn mkl_get_max_threads() -> c_int;
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
pub struct pardiso {
    pub type_0: linsys_solver_type,
    pub solve: Option::<unsafe extern "C" fn(*mut pardiso, *mut c_float) -> c_int>,
    pub free: Option::<unsafe extern "C" fn(*mut pardiso) -> ()>,
    pub update_matrices: Option::<
        unsafe extern "C" fn(*mut pardiso, *const csc, *const csc) -> c_int,
    >,
    pub update_rho_vec: Option::<
        unsafe extern "C" fn(*mut pardiso, *const c_float) -> c_int,
    >,
    pub nthreads: c_int,
    pub KKT: *mut csc,
    pub KKT_i: *mut c_int,
    pub KKT_p: *mut c_int,
    pub bp: *mut c_float,
    pub sol: *mut c_float,
    pub rho_inv_vec: *mut c_float,
    pub sigma: c_float,
    pub polish: c_int,
    pub n: c_int,
    pub m: c_int,
    pub pt: [*mut libc::c_void; 64],
    pub iparm: [c_int; 64],
    pub nKKT: c_int,
    pub mtype: c_int,
    pub nrhs: c_int,
    pub maxfct: c_int,
    pub mnum: c_int,
    pub phase: c_int,
    pub error: c_int,
    pub msglvl: c_int,
    pub idum: c_int,
    pub fdum: c_float,
    pub Pdiag_idx: *mut c_int,
    pub Pdiag_n: c_int,
    pub PtoKKT: *mut c_int,
    pub AtoKKT: *mut c_int,
    pub rhotoKKT: *mut c_int,
}
pub type pardiso_solver = pardiso;
pub const c_calloc: unsafe extern "C" fn(
    libc::c_ulong,
    libc::c_ulong,
) -> *mut libc::c_void = calloc;
pub const c_malloc: unsafe extern "C" fn(libc::c_ulong) -> *mut libc::c_void = malloc;
pub const c_print: unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int = printf;
pub const c_free: unsafe extern "C" fn(*mut libc::c_void) -> () = free;
pub const OSQP_NULL: libc::c_int = 0 as libc::c_int;
pub const MKL_INTERFACE_ILP64: libc::c_int = 0x1 as libc::c_int;
pub const PARDISO_SYMBOLIC: libc::c_int = 11 as libc::c_int;
pub const PARDISO_NUMERIC: libc::c_int = 22 as libc::c_int;
pub const PARDISO_SOLVE: libc::c_int = 33 as libc::c_int;
pub const PARDISO_CLEANUP: libc::c_int = -(1 as libc::c_int);
#[no_mangle]
pub unsafe extern "C" fn free_linsys_solver_pardiso(mut s: *mut pardiso_solver) {
    if !s.is_null() {
        (*s).phase = PARDISO_CLEANUP as c_int;
        pardiso(
            ((*s).pt).as_mut_ptr(),
            &mut (*s).maxfct,
            &mut (*s).mnum,
            &mut (*s).mtype,
            &mut (*s).phase,
            &mut (*s).nKKT,
            &mut (*s).fdum,
            (*s).KKT_p,
            (*s).KKT_i,
            &mut (*s).idum,
            &mut (*s).nrhs,
            ((*s).iparm).as_mut_ptr(),
            &mut (*s).msglvl,
            &mut (*s).fdum,
            &mut (*s).fdum,
            &mut (*s).error,
        );
        if (*s).error != 0 as libc::c_int as libc::c_longlong {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 27],
                    &[libc::c_char; 27],
                >(b"free_linsys_solver_pardiso\0"))
                    .as_ptr(),
            );
            printf(
                b"Error during MKL Pardiso cleanup: %d\0" as *const u8
                    as *const libc::c_char,
                (*s).error as libc::c_int,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
        }
        if !((*s).KKT).is_null() {
            csc_spfree((*s).KKT);
        }
        if !((*s).KKT_i).is_null() {
            free((*s).KKT_i as *mut libc::c_void);
        }
        if !((*s).KKT_p).is_null() {
            free((*s).KKT_p as *mut libc::c_void);
        }
        if !((*s).bp).is_null() {
            free((*s).bp as *mut libc::c_void);
        }
        if !((*s).sol).is_null() {
            free((*s).sol as *mut libc::c_void);
        }
        if !((*s).rho_inv_vec).is_null() {
            free((*s).rho_inv_vec as *mut libc::c_void);
        }
        if !((*s).Pdiag_idx).is_null() {
            free((*s).Pdiag_idx as *mut libc::c_void);
        }
        if !((*s).PtoKKT).is_null() {
            free((*s).PtoKKT as *mut libc::c_void);
        }
        if !((*s).AtoKKT).is_null() {
            free((*s).AtoKKT as *mut libc::c_void);
        }
        if !((*s).rhotoKKT).is_null() {
            free((*s).rhotoKKT as *mut libc::c_void);
        }
        free(s as *mut libc::c_void);
    }
}
#[no_mangle]
pub unsafe extern "C" fn init_linsys_solver_pardiso(
    mut sp: *mut *mut pardiso_solver,
    mut P: *const csc,
    mut A: *const csc,
    mut sigma: c_float,
    mut rho_vec: *const c_float,
    mut polish: c_int,
) -> c_int {
    let mut i: c_int = 0;
    let mut nnzKKT: c_int = 0;
    let mut n_plus_m: c_int = 0;
    let mut s: *mut pardiso_solver = 0 as *mut pardiso_solver;
    s = calloc(
        1 as libc::c_int as libc::c_ulong,
        ::std::mem::size_of::<pardiso_solver>() as libc::c_ulong,
    ) as *mut pardiso_solver;
    *sp = s;
    (*s).n = (*P).n;
    (*s).m = (*A).m;
    n_plus_m = (*s).n + (*s).m;
    (*s).nKKT = n_plus_m;
    (*s).sigma = sigma;
    (*s).polish = polish;
    let ref mut fresh0 = (*s).solve;
    *fresh0 = Some(
        solve_linsys_pardiso
            as unsafe extern "C" fn(*mut pardiso_solver, *mut c_float) -> c_int,
    );
    let ref mut fresh1 = (*s).free;
    *fresh1 = Some(
        free_linsys_solver_pardiso as unsafe extern "C" fn(*mut pardiso_solver) -> (),
    );
    let ref mut fresh2 = (*s).update_matrices;
    *fresh2 = Some(
        update_linsys_solver_matrices_pardiso
            as unsafe extern "C" fn(*mut pardiso_solver, *const csc, *const csc) -> c_int,
    );
    let ref mut fresh3 = (*s).update_rho_vec;
    *fresh3 = Some(
        update_linsys_solver_rho_vec_pardiso
            as unsafe extern "C" fn(*mut pardiso_solver, *const c_float) -> c_int,
    );
    (*s).type_0 = MKL_PARDISO_SOLVER;
    let ref mut fresh4 = (*s).bp;
    *fresh4 = malloc(
        (::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh5 = (*s).sol;
    *fresh5 = malloc(
        (::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh6 = (*s).rho_inv_vec;
    *fresh6 = malloc(
        (::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut c_float;
    if polish != 0 {
        i = 0 as libc::c_int as c_int;
        while i < (*A).m {
            *((*s).rho_inv_vec).offset(i as isize) = sigma;
            i += 1;
        }
        let ref mut fresh7 = (*s).KKT;
        *fresh7 = form_KKT(
            P,
            A,
            1 as libc::c_int as c_int,
            sigma,
            (*s).rho_inv_vec,
            OSQP_NULL as *mut c_int,
            OSQP_NULL as *mut c_int,
            OSQP_NULL as *mut *mut c_int,
            OSQP_NULL as *mut c_int,
            OSQP_NULL as *mut c_int,
        );
    } else {
        let ref mut fresh8 = (*s).PtoKKT;
        *fresh8 = malloc(
            (*((*P).p).offset((*P).n as isize) as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        let ref mut fresh9 = (*s).AtoKKT;
        *fresh9 = malloc(
            (*((*A).p).offset((*A).n as isize) as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        let ref mut fresh10 = (*s).rhotoKKT;
        *fresh10 = malloc(
            ((*A).m as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        i = 0 as libc::c_int as c_int;
        while i < (*A).m {
            *((*s).rho_inv_vec)
                .offset(i as isize) = 1.0f64 / *rho_vec.offset(i as isize);
            i += 1;
        }
        let ref mut fresh11 = (*s).KKT;
        *fresh11 = form_KKT(
            P,
            A,
            1 as libc::c_int as c_int,
            sigma,
            (*s).rho_inv_vec,
            (*s).PtoKKT,
            (*s).AtoKKT,
            &mut (*s).Pdiag_idx,
            &mut (*s).Pdiag_n,
            (*s).rhotoKKT,
        );
    }
    if ((*s).KKT).is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 27],
                &[libc::c_char; 27],
            >(b"init_linsys_solver_pardiso\0"))
                .as_ptr(),
        );
        printf(b"Error in forming KKT matrix\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        free_linsys_solver_pardiso(s);
        return OSQP_LINSYS_SOLVER_INIT_ERROR as libc::c_int as c_int;
    } else {
        nnzKKT = *((*(*s).KKT).p).offset((*(*s).KKT).m as isize);
        let ref mut fresh12 = (*s).KKT_i;
        *fresh12 = malloc(
            (nnzKKT as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        let ref mut fresh13 = (*s).KKT_p;
        *fresh13 = malloc(
            (((*(*s).KKT).m + 1 as libc::c_int as libc::c_longlong) as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        i = 0 as libc::c_int as c_int;
        while i < nnzKKT {
            *((*s).KKT_i)
                .offset(
                    i as isize,
                ) = *((*(*s).KKT).i).offset(i as isize)
                + 1 as libc::c_int as libc::c_longlong;
            i += 1;
        }
        i = 0 as libc::c_int as c_int;
        while i < n_plus_m + 1 as libc::c_int as libc::c_longlong {
            *((*s).KKT_p)
                .offset(
                    i as isize,
                ) = *((*(*s).KKT).p).offset(i as isize)
                + 1 as libc::c_int as libc::c_longlong;
            i += 1;
        }
    }
    mkl_set_interface_layer(MKL_INTERFACE_ILP64 as c_int);
    (*s).mtype = -(2 as libc::c_int) as c_int;
    (*s).nrhs = 1 as libc::c_int as c_int;
    (*s).maxfct = 1 as libc::c_int as c_int;
    (*s).mnum = 1 as libc::c_int as c_int;
    (*s).msglvl = 0 as libc::c_int as c_int;
    (*s).error = 0 as libc::c_int as c_int;
    i = 0 as libc::c_int as c_int;
    while i < 64 as libc::c_int as libc::c_longlong {
        (*s).iparm[i as usize] = 0 as libc::c_int as c_int;
        let ref mut fresh14 = (*s).pt[i as usize];
        *fresh14 = 0 as *mut libc::c_void;
        i += 1;
    }
    (*s).iparm[0 as libc::c_int as usize] = 1 as libc::c_int as c_int;
    (*s).iparm[1 as libc::c_int as usize] = 3 as libc::c_int as c_int;
    if polish != 0 {
        (*s).iparm[5 as libc::c_int as usize] = 1 as libc::c_int as c_int;
    } else {
        (*s).iparm[5 as libc::c_int as usize] = 0 as libc::c_int as c_int;
    }
    (*s).iparm[7 as libc::c_int as usize] = 0 as libc::c_int as c_int;
    (*s).iparm[9 as libc::c_int as usize] = 13 as libc::c_int as c_int;
    (*s).iparm[34 as libc::c_int as usize] = 0 as libc::c_int as c_int;
    (*s).nthreads = mkl_get_max_threads();
    (*s).phase = PARDISO_SYMBOLIC as c_int;
    pardiso(
        ((*s).pt).as_mut_ptr(),
        &mut (*s).maxfct,
        &mut (*s).mnum,
        &mut (*s).mtype,
        &mut (*s).phase,
        &mut (*s).nKKT,
        (*(*s).KKT).x,
        (*s).KKT_p,
        (*s).KKT_i,
        &mut (*s).idum,
        &mut (*s).nrhs,
        ((*s).iparm).as_mut_ptr(),
        &mut (*s).msglvl,
        &mut (*s).fdum,
        &mut (*s).fdum,
        &mut (*s).error,
    );
    if (*s).error != 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 27],
                &[libc::c_char; 27],
            >(b"init_linsys_solver_pardiso\0"))
                .as_ptr(),
        );
        printf(
            b"Error during symbolic factorization: %d\0" as *const u8
                as *const libc::c_char,
            (*s).error as libc::c_int,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        free_linsys_solver_pardiso(s);
        *sp = OSQP_NULL as *mut pardiso_solver;
        return OSQP_LINSYS_SOLVER_INIT_ERROR as libc::c_int as c_int;
    }
    (*s).phase = PARDISO_NUMERIC as c_int;
    pardiso(
        ((*s).pt).as_mut_ptr(),
        &mut (*s).maxfct,
        &mut (*s).mnum,
        &mut (*s).mtype,
        &mut (*s).phase,
        &mut (*s).nKKT,
        (*(*s).KKT).x,
        (*s).KKT_p,
        (*s).KKT_i,
        &mut (*s).idum,
        &mut (*s).nrhs,
        ((*s).iparm).as_mut_ptr(),
        &mut (*s).msglvl,
        &mut (*s).fdum,
        &mut (*s).fdum,
        &mut (*s).error,
    );
    if (*s).error != 0 {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 27],
                &[libc::c_char; 27],
            >(b"init_linsys_solver_pardiso\0"))
                .as_ptr(),
        );
        printf(
            b"Error during numerical factorization: %d\0" as *const u8
                as *const libc::c_char,
            (*s).error as libc::c_int,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        free_linsys_solver_pardiso(s);
        *sp = OSQP_NULL as *mut pardiso_solver;
        return OSQP_LINSYS_SOLVER_INIT_ERROR as libc::c_int as c_int;
    }
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn solve_linsys_pardiso(
    mut s: *mut pardiso_solver,
    mut b: *mut c_float,
) -> c_int {
    let mut j: c_int = 0;
    (*s).phase = PARDISO_SOLVE as c_int;
    pardiso(
        ((*s).pt).as_mut_ptr(),
        &mut (*s).maxfct,
        &mut (*s).mnum,
        &mut (*s).mtype,
        &mut (*s).phase,
        &mut (*s).nKKT,
        (*(*s).KKT).x,
        (*s).KKT_p,
        (*s).KKT_i,
        &mut (*s).idum,
        &mut (*s).nrhs,
        ((*s).iparm).as_mut_ptr(),
        &mut (*s).msglvl,
        b,
        (*s).sol,
        &mut (*s).error,
    );
    if (*s).error != 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 21],
                &[libc::c_char; 21],
            >(b"solve_linsys_pardiso\0"))
                .as_ptr(),
        );
        printf(
            b"Error during linear system solution: %d\0" as *const u8
                as *const libc::c_char,
            (*s).error as libc::c_int,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return 1 as libc::c_int as c_int;
    }
    if (*s).polish == 0 {
        j = 0 as libc::c_int as c_int;
        while j < (*s).n {
            *b.offset(j as isize) = *((*s).sol).offset(j as isize);
            j += 1;
        }
        j = 0 as libc::c_int as c_int;
        while j < (*s).m {
            let ref mut fresh15 = *b.offset((j + (*s).n) as isize);
            *fresh15
                += *((*s).rho_inv_vec).offset(j as isize)
                    * *((*s).sol).offset((j + (*s).n) as isize);
            j += 1;
        }
    }
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn update_linsys_solver_matrices_pardiso(
    mut s: *mut pardiso_solver,
    mut P: *const csc,
    mut A: *const csc,
) -> c_int {
    update_KKT_P((*s).KKT, P, (*s).PtoKKT, (*s).sigma, (*s).Pdiag_idx, (*s).Pdiag_n);
    update_KKT_A((*s).KKT, A, (*s).AtoKKT);
    (*s).phase = PARDISO_NUMERIC as c_int;
    pardiso(
        ((*s).pt).as_mut_ptr(),
        &mut (*s).maxfct,
        &mut (*s).mnum,
        &mut (*s).mtype,
        &mut (*s).phase,
        &mut (*s).nKKT,
        (*(*s).KKT).x,
        (*s).KKT_p,
        (*s).KKT_i,
        &mut (*s).idum,
        &mut (*s).nrhs,
        ((*s).iparm).as_mut_ptr(),
        &mut (*s).msglvl,
        &mut (*s).fdum,
        &mut (*s).fdum,
        &mut (*s).error,
    );
    return (*s).error;
}
#[no_mangle]
pub unsafe extern "C" fn update_linsys_solver_rho_vec_pardiso(
    mut s: *mut pardiso_solver,
    mut rho_vec: *const c_float,
) -> c_int {
    let mut i: c_int = 0;
    i = 0 as libc::c_int as c_int;
    while i < (*s).m {
        *((*s).rho_inv_vec).offset(i as isize) = 1.0f64 / *rho_vec.offset(i as isize);
        i += 1;
    }
    update_KKT_param2((*s).KKT, (*s).rho_inv_vec, (*s).rhotoKKT, (*s).m);
    (*s).phase = PARDISO_NUMERIC as c_int;
    pardiso(
        ((*s).pt).as_mut_ptr(),
        &mut (*s).maxfct,
        &mut (*s).mnum,
        &mut (*s).mtype,
        &mut (*s).phase,
        &mut (*s).nKKT,
        (*(*s).KKT).x,
        (*s).KKT_p,
        (*s).KKT_i,
        &mut (*s).idum,
        &mut (*s).nrhs,
        ((*s).iparm).as_mut_ptr(),
        &mut (*s).msglvl,
        &mut (*s).fdum,
        &mut (*s).fdum,
        &mut (*s).error,
    );
    return (*s).error;
}
