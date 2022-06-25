use ::libc;
extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
    fn QDLDL_etree(
        n: QDLDL_int,
        Ap: *const QDLDL_int,
        Ai: *const QDLDL_int,
        work: *mut QDLDL_int,
        Lnz: *mut QDLDL_int,
        etree: *mut QDLDL_int,
    ) -> QDLDL_int;
    fn QDLDL_factor(
        n: QDLDL_int,
        Ap: *const QDLDL_int,
        Ai: *const QDLDL_int,
        Ax: *const QDLDL_float,
        Lp: *mut QDLDL_int,
        Li: *mut QDLDL_int,
        Lx: *mut QDLDL_float,
        D: *mut QDLDL_float,
        Dinv: *mut QDLDL_float,
        Lnz: *const QDLDL_int,
        etree: *const QDLDL_int,
        bwork: *mut QDLDL_bool,
        iwork: *mut QDLDL_int,
        fwork: *mut QDLDL_float,
    ) -> QDLDL_int;
    fn QDLDL_solve(
        n: QDLDL_int,
        Lp: *const QDLDL_int,
        Li: *const QDLDL_int,
        Lx: *const QDLDL_float,
        Dinv: *const QDLDL_float,
        x: *mut QDLDL_float,
    );
    fn amd_l_order(
        n: libc::c_longlong,
        Ap: *const libc::c_longlong,
        Ai: *const libc::c_longlong,
        P: *mut libc::c_longlong,
        Control: *mut c_float,
        Info: *mut c_float,
    ) -> libc::c_longlong;
    fn csc_spfree(A: *mut csc);
    fn csc_pinv(p: *const c_int, n: c_int) -> *mut c_int;
    fn csc_symperm(
        A: *const csc,
        pinv: *const c_int,
        AtoC: *mut c_int,
        values: c_int,
    ) -> *mut csc;
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
}
pub type c_int = libc::c_longlong;
pub type c_float = libc::c_double;
pub type QDLDL_int = libc::c_longlong;
pub type QDLDL_float = libc::c_double;
pub type QDLDL_bool = libc::c_uchar;
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
pub struct qdldl {
    pub type_0: linsys_solver_type,
    pub solve: Option::<unsafe extern "C" fn(*mut qdldl, *mut c_float) -> c_int>,
    pub free: Option::<unsafe extern "C" fn(*mut qdldl) -> ()>,
    pub update_matrices: Option::<
        unsafe extern "C" fn(*mut qdldl, *const csc, *const csc) -> c_int,
    >,
    pub update_rho_vec: Option::<
        unsafe extern "C" fn(*mut qdldl, *const c_float) -> c_int,
    >,
    pub nthreads: c_int,
    pub L: *mut csc,
    pub Dinv: *mut c_float,
    pub P: *mut c_int,
    pub bp: *mut c_float,
    pub sol: *mut c_float,
    pub rho_inv_vec: *mut c_float,
    pub sigma: c_float,
    pub polish: c_int,
    pub n: c_int,
    pub m: c_int,
    pub Pdiag_idx: *mut c_int,
    pub Pdiag_n: c_int,
    pub KKT: *mut csc,
    pub PtoKKT: *mut c_int,
    pub AtoKKT: *mut c_int,
    pub rhotoKKT: *mut c_int,
    pub D: *mut QDLDL_float,
    pub etree: *mut QDLDL_int,
    pub Lnz: *mut QDLDL_int,
    pub iwork: *mut QDLDL_int,
    pub bwork: *mut QDLDL_bool,
    pub fwork: *mut QDLDL_float,
}
pub type qdldl_solver = qdldl;
pub const c_calloc: unsafe extern "C" fn(
    libc::c_ulong,
    libc::c_ulong,
) -> *mut libc::c_void = calloc;
pub const c_malloc: unsafe extern "C" fn(libc::c_ulong) -> *mut libc::c_void = malloc;
pub const c_print: unsafe extern "C" fn(*const libc::c_char, ...) -> libc::c_int = printf;
pub const c_free: unsafe extern "C" fn(*mut libc::c_void) -> () = free;
pub const OSQP_NULL: libc::c_int = 0 as libc::c_int;
pub const AMD_INFO: libc::c_int = 20 as libc::c_int;
#[no_mangle]
pub unsafe extern "C" fn free_linsys_solver_qdldl(mut s: *mut qdldl_solver) {
    if !s.is_null() {
        if !((*s).L).is_null() {
            csc_spfree((*s).L);
        }
        if !((*s).P).is_null() {
            free((*s).P as *mut libc::c_void);
        }
        if !((*s).Dinv).is_null() {
            free((*s).Dinv as *mut libc::c_void);
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
        if !((*s).KKT).is_null() {
            csc_spfree((*s).KKT);
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
        if !((*s).D).is_null() {
            free((*s).D as *mut libc::c_void);
        }
        if !((*s).etree).is_null() {
            free((*s).etree as *mut libc::c_void);
        }
        if !((*s).Lnz).is_null() {
            free((*s).Lnz as *mut libc::c_void);
        }
        if !((*s).iwork).is_null() {
            free((*s).iwork as *mut libc::c_void);
        }
        if !((*s).bwork).is_null() {
            free((*s).bwork as *mut libc::c_void);
        }
        if !((*s).fwork).is_null() {
            free((*s).fwork as *mut libc::c_void);
        }
        free(s as *mut libc::c_void);
    }
}
unsafe extern "C" fn LDL_factor(
    mut A: *mut csc,
    mut p: *mut qdldl_solver,
    mut nvar: c_int,
) -> c_int {
    let mut sum_Lnz: c_int = 0;
    let mut factor_status: c_int = 0;
    sum_Lnz = QDLDL_etree((*A).n, (*A).p, (*A).i, (*p).iwork, (*p).Lnz, (*p).etree);
    if sum_Lnz < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"LDL_factor\0"))
                .as_ptr(),
        );
        printf(
            b"Error in KKT matrix LDL factorization when computing the elimination tree.\0"
                as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        if sum_Lnz == -(1 as libc::c_int) as libc::c_longlong {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"LDL_factor\0"))
                    .as_ptr(),
            );
            printf(
                b"Matrix is not perfectly upper triangular.\0" as *const u8
                    as *const libc::c_char,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
        } else if sum_Lnz == -(2 as libc::c_int) as libc::c_longlong {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"LDL_factor\0"))
                    .as_ptr(),
            );
            printf(
                b"Integer overflow in L nonzero count.\0" as *const u8
                    as *const libc::c_char,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
        }
        return sum_Lnz;
    }
    let ref mut fresh0 = (*(*p).L).i;
    *fresh0 = malloc(
        (::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(sum_Lnz as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut c_int;
    let ref mut fresh1 = (*(*p).L).x;
    *fresh1 = malloc(
        (::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(sum_Lnz as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut c_float;
    (*(*p).L).nzmax = sum_Lnz;
    factor_status = QDLDL_factor(
        (*A).n,
        (*A).p,
        (*A).i,
        (*A).x,
        (*(*p).L).p,
        (*(*p).L).i,
        (*(*p).L).x,
        (*p).D,
        (*p).Dinv,
        (*p).Lnz,
        (*p).etree,
        (*p).bwork,
        (*p).iwork,
        (*p).fwork,
    );
    if factor_status < 0 as libc::c_int as libc::c_longlong {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<&[u8; 11], &[libc::c_char; 11]>(b"LDL_factor\0"))
                .as_ptr(),
        );
        printf(
            b"Error in KKT matrix LDL factorization when computing the nonzero elements. There are zeros in the diagonal matrix\0"
                as *const u8 as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        return factor_status;
    } else {
        if factor_status < nvar {
            printf(
                b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
                (*::std::mem::transmute::<
                    &[u8; 11],
                    &[libc::c_char; 11],
                >(b"LDL_factor\0"))
                    .as_ptr(),
            );
            printf(
                b"Error in KKT matrix LDL factorization when computing the nonzero elements. The problem seems to be non-convex\0"
                    as *const u8 as *const libc::c_char,
            );
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            return -(2 as libc::c_int) as c_int;
        }
    }
    return 0 as libc::c_int as c_int;
}
unsafe extern "C" fn permute_KKT(
    mut KKT: *mut *mut csc,
    mut p: *mut qdldl_solver,
    mut Pnz: c_int,
    mut Anz: c_int,
    mut m: c_int,
    mut PtoKKT: *mut c_int,
    mut AtoKKT: *mut c_int,
    mut rhotoKKT: *mut c_int,
) -> c_int {
    let mut info: *mut c_float = 0 as *mut c_float;
    let mut amd_status: c_int = 0;
    let mut Pinv: *mut c_int = 0 as *mut c_int;
    let mut KKT_temp: *mut csc = 0 as *mut csc;
    let mut KtoPKPt: *mut c_int = 0 as *mut c_int;
    let mut i: c_int = 0;
    info = malloc(
        (AMD_INFO as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<c_float>() as libc::c_ulong),
    ) as *mut c_float;
    amd_status = amd_l_order(
        (**KKT).n,
        (**KKT).p as *const libc::c_longlong,
        (**KKT).i as *const libc::c_longlong,
        (*p).P,
        0 as *mut c_float,
        info,
    );
    if amd_status < 0 as libc::c_int as libc::c_longlong {
        free(info as *mut libc::c_void);
        return amd_status;
    }
    Pinv = csc_pinv((*p).P, (**KKT).n);
    if PtoKKT.is_null() && AtoKKT.is_null() && rhotoKKT.is_null() {
        KKT_temp = csc_symperm(
            *KKT,
            Pinv,
            OSQP_NULL as *mut c_int,
            1 as libc::c_int as c_int,
        );
    } else {
        KtoPKPt = malloc(
            (*((**KKT).p).offset((**KKT).n as isize) as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        KKT_temp = csc_symperm(*KKT, Pinv, KtoPKPt, 1 as libc::c_int as c_int);
        if !PtoKKT.is_null() {
            i = 0 as libc::c_int as c_int;
            while i < Pnz {
                *PtoKKT
                    .offset(
                        i as isize,
                    ) = *KtoPKPt.offset(*PtoKKT.offset(i as isize) as isize);
                i += 1;
            }
        }
        if !AtoKKT.is_null() {
            i = 0 as libc::c_int as c_int;
            while i < Anz {
                *AtoKKT
                    .offset(
                        i as isize,
                    ) = *KtoPKPt.offset(*AtoKKT.offset(i as isize) as isize);
                i += 1;
            }
        }
        if !rhotoKKT.is_null() {
            i = 0 as libc::c_int as c_int;
            while i < m {
                *rhotoKKT
                    .offset(
                        i as isize,
                    ) = *KtoPKPt.offset(*rhotoKKT.offset(i as isize) as isize);
                i += 1;
            }
        }
        free(KtoPKPt as *mut libc::c_void);
    }
    csc_spfree(*KKT);
    *KKT = KKT_temp;
    free(Pinv as *mut libc::c_void);
    free(info as *mut libc::c_void);
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn init_linsys_solver_qdldl(
    mut sp: *mut *mut qdldl_solver,
    mut P: *const csc,
    mut A: *const csc,
    mut sigma: c_float,
    mut rho_vec: *const c_float,
    mut polish: c_int,
) -> c_int {
    let mut KKT_temp: *mut csc = 0 as *mut csc;
    let mut i: c_int = 0;
    let mut n_plus_m: c_int = 0;
    let mut s: *mut qdldl_solver = 0 as *mut qdldl_solver;
    s = calloc(
        1 as libc::c_int as libc::c_ulong,
        ::std::mem::size_of::<qdldl_solver>() as libc::c_ulong,
    ) as *mut qdldl_solver;
    *sp = s;
    (*s).n = (*P).n;
    (*s).m = (*A).m;
    n_plus_m = (*s).n + (*s).m;
    (*s).sigma = sigma;
    (*s).polish = polish;
    let ref mut fresh2 = (*s).solve;
    *fresh2 = Some(
        solve_linsys_qdldl
            as unsafe extern "C" fn(*mut qdldl_solver, *mut c_float) -> c_int,
    );
    let ref mut fresh3 = (*s).free;
    *fresh3 = Some(
        free_linsys_solver_qdldl as unsafe extern "C" fn(*mut qdldl_solver) -> (),
    );
    let ref mut fresh4 = (*s).update_matrices;
    *fresh4 = Some(
        update_linsys_solver_matrices_qdldl
            as unsafe extern "C" fn(*mut qdldl_solver, *const csc, *const csc) -> c_int,
    );
    let ref mut fresh5 = (*s).update_rho_vec;
    *fresh5 = Some(
        update_linsys_solver_rho_vec_qdldl
            as unsafe extern "C" fn(*mut qdldl_solver, *const c_float) -> c_int,
    );
    (*s).type_0 = QDLDL_SOLVER;
    (*s).nthreads = 1 as libc::c_int as c_int;
    let ref mut fresh6 = (*s).L;
    *fresh6 = malloc(::std::mem::size_of::<csc>() as libc::c_ulong) as *mut csc;
    (*(*s).L).m = n_plus_m;
    (*(*s).L).n = n_plus_m;
    (*(*s).L).nz = -(1 as libc::c_int) as c_int;
    let ref mut fresh7 = (*s).Dinv;
    *fresh7 = malloc(
        (::std::mem::size_of::<QDLDL_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_float;
    let ref mut fresh8 = (*s).D;
    *fresh8 = malloc(
        (::std::mem::size_of::<QDLDL_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_float;
    let ref mut fresh9 = (*s).P;
    *fresh9 = malloc(
        (::std::mem::size_of::<QDLDL_int>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_int;
    let ref mut fresh10 = (*s).bp;
    *fresh10 = malloc(
        (::std::mem::size_of::<QDLDL_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_float;
    let ref mut fresh11 = (*s).sol;
    *fresh11 = malloc(
        (::std::mem::size_of::<QDLDL_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_float;
    let ref mut fresh12 = (*s).rho_inv_vec;
    *fresh12 = malloc(
        (::std::mem::size_of::<c_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul((*s).m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut c_float;
    let ref mut fresh13 = (*s).etree;
    *fresh13 = malloc(
        (n_plus_m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<QDLDL_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut QDLDL_int;
    let ref mut fresh14 = (*s).Lnz;
    *fresh14 = malloc(
        (n_plus_m as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<QDLDL_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut QDLDL_int;
    let ref mut fresh15 = (*(*s).L).p;
    *fresh15 = malloc(
        ((n_plus_m + 1 as libc::c_int as libc::c_longlong) as libc::c_ulonglong)
            .wrapping_mul(
                ::std::mem::size_of::<QDLDL_int>() as libc::c_ulong as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut c_int;
    let ref mut fresh16 = (*(*s).L).i;
    *fresh16 = OSQP_NULL as *mut c_int;
    let ref mut fresh17 = (*(*s).L).x;
    *fresh17 = OSQP_NULL as *mut c_float;
    let ref mut fresh18 = (*s).iwork;
    *fresh18 = malloc(
        (::std::mem::size_of::<QDLDL_int>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(
                (3 as libc::c_int as libc::c_longlong * n_plus_m) as libc::c_ulonglong,
            ) as libc::c_ulong,
    ) as *mut QDLDL_int;
    let ref mut fresh19 = (*s).bwork;
    *fresh19 = malloc(
        (::std::mem::size_of::<QDLDL_bool>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_bool;
    let ref mut fresh20 = (*s).fwork;
    *fresh20 = malloc(
        (::std::mem::size_of::<QDLDL_float>() as libc::c_ulong as libc::c_ulonglong)
            .wrapping_mul(n_plus_m as libc::c_ulonglong) as libc::c_ulong,
    ) as *mut QDLDL_float;
    if polish != 0 {
        i = 0 as libc::c_int as c_int;
        while i < (*A).m {
            *((*s).rho_inv_vec).offset(i as isize) = sigma;
            i += 1;
        }
        KKT_temp = form_KKT(
            P,
            A,
            0 as libc::c_int as c_int,
            sigma,
            (*s).rho_inv_vec,
            OSQP_NULL as *mut c_int,
            OSQP_NULL as *mut c_int,
            OSQP_NULL as *mut *mut c_int,
            OSQP_NULL as *mut c_int,
            OSQP_NULL as *mut c_int,
        );
        if !KKT_temp.is_null() {
            permute_KKT(
                &mut KKT_temp,
                s,
                OSQP_NULL as c_int,
                OSQP_NULL as c_int,
                OSQP_NULL as c_int,
                OSQP_NULL as *mut c_int,
                OSQP_NULL as *mut c_int,
                OSQP_NULL as *mut c_int,
            );
        }
    } else {
        let ref mut fresh21 = (*s).PtoKKT;
        *fresh21 = malloc(
            (*((*P).p).offset((*P).n as isize) as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        let ref mut fresh22 = (*s).AtoKKT;
        *fresh22 = malloc(
            (*((*A).p).offset((*A).n as isize) as libc::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as libc::c_ulong as libc::c_ulonglong,
                ) as libc::c_ulong,
        ) as *mut c_int;
        let ref mut fresh23 = (*s).rhotoKKT;
        *fresh23 = malloc(
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
        KKT_temp = form_KKT(
            P,
            A,
            0 as libc::c_int as c_int,
            sigma,
            (*s).rho_inv_vec,
            (*s).PtoKKT,
            (*s).AtoKKT,
            &mut (*s).Pdiag_idx,
            &mut (*s).Pdiag_n,
            (*s).rhotoKKT,
        );
        if !KKT_temp.is_null() {
            permute_KKT(
                &mut KKT_temp,
                s,
                *((*P).p).offset((*P).n as isize),
                *((*A).p).offset((*A).n as isize),
                (*A).m,
                (*s).PtoKKT,
                (*s).AtoKKT,
                (*s).rhotoKKT,
            );
        }
    }
    if KKT_temp.is_null() {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const libc::c_char,
            (*::std::mem::transmute::<
                &[u8; 25],
                &[libc::c_char; 25],
            >(b"init_linsys_solver_qdldl\0"))
                .as_ptr(),
        );
        printf(
            b"Error forming and permuting KKT matrix\0" as *const u8
                as *const libc::c_char,
        );
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        free_linsys_solver_qdldl(s);
        *sp = OSQP_NULL as *mut qdldl_solver;
        return OSQP_LINSYS_SOLVER_INIT_ERROR as libc::c_int as c_int;
    }
    if LDL_factor(KKT_temp, s, (*P).n) < 0 as libc::c_int as libc::c_longlong {
        csc_spfree(KKT_temp);
        free_linsys_solver_qdldl(s);
        *sp = OSQP_NULL as *mut qdldl_solver;
        return OSQP_NONCVX_ERROR as libc::c_int as c_int;
    }
    if polish != 0 {
        csc_spfree(KKT_temp);
    } else {
        let ref mut fresh24 = (*s).KKT;
        *fresh24 = KKT_temp;
    }
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn permute_x(
    mut n: c_int,
    mut x: *mut c_float,
    mut b: *const c_float,
    mut P: *const c_int,
) {
    let mut j: c_int = 0;
    j = 0 as libc::c_int as c_int;
    while j < n {
        *x.offset(j as isize) = *b.offset(*P.offset(j as isize) as isize);
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn permutet_x(
    mut n: c_int,
    mut x: *mut c_float,
    mut b: *const c_float,
    mut P: *const c_int,
) {
    let mut j: c_int = 0;
    j = 0 as libc::c_int as c_int;
    while j < n {
        *x.offset(*P.offset(j as isize) as isize) = *b.offset(j as isize);
        j += 1;
    }
}
unsafe extern "C" fn LDLSolve(
    mut x: *mut c_float,
    mut b: *mut c_float,
    mut L: *const csc,
    mut Dinv: *const c_float,
    mut P: *const c_int,
    mut bp: *mut c_float,
) {
    permute_x((*L).n, bp, b, P);
    QDLDL_solve((*L).n, (*L).p, (*L).i, (*L).x, Dinv, bp);
    permutet_x((*L).n, x, bp, P);
}
#[no_mangle]
pub unsafe extern "C" fn solve_linsys_qdldl(
    mut s: *mut qdldl_solver,
    mut b: *mut c_float,
) -> c_int {
    let mut j: c_int = 0;
    if (*s).polish != 0 {
        LDLSolve(b, b, (*s).L, (*s).Dinv, (*s).P, (*s).bp);
    } else {
        LDLSolve((*s).sol, b, (*s).L, (*s).Dinv, (*s).P, (*s).bp);
        j = 0 as libc::c_int as c_int;
        while j < (*s).n {
            *b.offset(j as isize) = *((*s).sol).offset(j as isize);
            j += 1;
        }
        j = 0 as libc::c_int as c_int;
        while j < (*s).m {
            let ref mut fresh25 = *b.offset((j + (*s).n) as isize);
            *fresh25
                += *((*s).rho_inv_vec).offset(j as isize)
                    * *((*s).sol).offset((j + (*s).n) as isize);
            j += 1;
        }
    }
    return 0 as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn update_linsys_solver_matrices_qdldl(
    mut s: *mut qdldl_solver,
    mut P: *const csc,
    mut A: *const csc,
) -> c_int {
    update_KKT_P((*s).KKT, P, (*s).PtoKKT, (*s).sigma, (*s).Pdiag_idx, (*s).Pdiag_n);
    update_KKT_A((*s).KKT, A, (*s).AtoKKT);
    return (QDLDL_factor(
        (*(*s).KKT).n,
        (*(*s).KKT).p,
        (*(*s).KKT).i,
        (*(*s).KKT).x,
        (*(*s).L).p,
        (*(*s).L).i,
        (*(*s).L).x,
        (*s).D,
        (*s).Dinv,
        (*s).Lnz,
        (*s).etree,
        (*s).bwork,
        (*s).iwork,
        (*s).fwork,
    ) < 0 as libc::c_int as libc::c_longlong) as libc::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn update_linsys_solver_rho_vec_qdldl(
    mut s: *mut qdldl_solver,
    mut rho_vec: *const c_float,
) -> c_int {
    let mut i: c_int = 0;
    i = 0 as libc::c_int as c_int;
    while i < (*s).m {
        *((*s).rho_inv_vec).offset(i as isize) = 1.0f64 / *rho_vec.offset(i as isize);
        i += 1;
    }
    update_KKT_param2((*s).KKT, (*s).rho_inv_vec, (*s).rhotoKKT, (*s).m);
    return (QDLDL_factor(
        (*(*s).KKT).n,
        (*(*s).KKT).p,
        (*(*s).KKT).i,
        (*(*s).KKT).x,
        (*(*s).L).p,
        (*(*s).L).i,
        (*(*s).L).x,
        (*s).D,
        (*s).Dinv,
        (*s).Lnz,
        (*s).etree,
        (*s).bwork,
        (*s).iwork,
        (*s).fwork,
    ) < 0 as libc::c_int as libc::c_longlong) as libc::c_int as c_int;
}
