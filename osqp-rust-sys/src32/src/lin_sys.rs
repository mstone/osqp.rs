extern "C" {
    fn init_linsys_solver_qdldl(
        sp: *mut *mut qdldl_solver,
        P: *const csc,
        A: *const csc,
        sigma: c_float,
        rho_vec: *const c_float,
        polish: c_int,
    ) -> c_int;
    fn init_linsys_solver_pardiso(
        sp: *mut *mut pardiso_solver,
        P: *const csc,
        A: *const csc,
        sigma: c_float,
        rho_vec: *const c_float,
        polish: c_int,
    ) -> c_int;
    fn lh_load_pardiso(libname: *const ::std::os::raw::c_char) -> c_int;
    fn lh_unload_pardiso() -> c_int;
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
pub type qdldl_solver = qdldl;
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
pub type QDLDL_float = ::std::os::raw::c_double;
pub type QDLDL_bool = ::std::os::raw::c_uchar;
pub type QDLDL_int = ::std::os::raw::c_int;
pub type pardiso_solver = pardiso;
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
    pub pt: [*mut ::std::os::raw::c_void; 64],
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
pub const OSQP_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
#[no_mangle]
pub static mut LINSYS_SOLVER_NAME: [*const ::std::os::raw::c_char; 2] = [
    b"qdldl\0" as *const u8 as *const ::std::os::raw::c_char,
    b"mkl pardiso\0" as *const u8 as *const ::std::os::raw::c_char,
];
#[no_mangle]
pub unsafe extern "C" fn load_linsys_solver(
    mut linsys_solver: linsys_solver_type,
) -> c_int {
    match linsys_solver as ::std::os::raw::c_uint {
        0 => return 0 as ::std::os::raw::c_int,
        1 => return lh_load_pardiso(OSQP_NULL as *const ::std::os::raw::c_char),
        _ => return 0 as ::std::os::raw::c_int,
    };
}
#[no_mangle]
pub unsafe extern "C" fn unload_linsys_solver(
    mut linsys_solver: linsys_solver_type,
) -> c_int {
    match linsys_solver as ::std::os::raw::c_uint {
        0 => return 0 as ::std::os::raw::c_int,
        1 => return lh_unload_pardiso(),
        _ => return 0 as ::std::os::raw::c_int,
    };
}
#[no_mangle]
pub unsafe extern "C" fn init_linsys_solver(
    mut s: *mut *mut LinSysSolver,
    mut P: *const csc,
    mut A: *const csc,
    mut sigma: c_float,
    mut rho_vec: *const c_float,
    mut linsys_solver: linsys_solver_type,
    mut polish: c_int,
) -> c_int {
    match linsys_solver as ::std::os::raw::c_uint {
        0 => {
            return init_linsys_solver_qdldl(
                s as *mut *mut qdldl_solver,
                P,
                A,
                sigma,
                rho_vec,
                polish,
            );
        }
        1 => {
            return init_linsys_solver_pardiso(
                s as *mut *mut pardiso_solver,
                P,
                A,
                sigma,
                rho_vec,
                polish,
            );
        }
        _ => {
            return init_linsys_solver_qdldl(
                s as *mut *mut qdldl_solver,
                P,
                A,
                sigma,
                rho_vec,
                polish,
            );
        }
    };
}
