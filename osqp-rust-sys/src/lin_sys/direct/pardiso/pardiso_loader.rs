extern "C" {
    fn lh_unload_lib(libhandle: soHandle_t) -> c_int;
    fn lh_load_lib(libname: *const ::std::os::raw::c_char) -> soHandle_t;
    fn printf(_: *const ::std::os::raw::c_char, _: ...) -> ::std::os::raw::c_int;
    fn lh_load_sym(h: soHandle_t, symName: *const ::std::os::raw::c_char) -> voidfun;
}
pub type c_int = ::std::os::raw::c_longlong;
pub type c_float = ::std::os::raw::c_double;
pub type soHandle_t = *mut ::std::os::raw::c_void;
pub type mkl_get_mt_t = Option::<unsafe extern "C" fn() -> ::std::os::raw::c_int>;
pub type voidfun = Option::<unsafe extern "C" fn() -> ()>;
pub type mkl_set_ifl_t = Option::<unsafe extern "C" fn(::std::os::raw::c_int) -> ::std::os::raw::c_int>;
pub type pardiso_t = Option::<
    unsafe extern "C" fn(
        *mut *mut ::std::os::raw::c_void,
        *const c_int,
        *const c_int,
        *const c_int,
        *const c_int,
        *const c_int,
        *const c_float,
        *const c_int,
        *const c_int,
        *mut c_int,
        *const c_int,
        *mut c_int,
        *const c_int,
        *mut c_float,
        *mut c_float,
        *mut c_int,
    ) -> (),
>;
pub const c_print: unsafe extern "C" fn(*const ::std::os::raw::c_char, ...) -> ::std::os::raw::c_int = printf;
pub const OSQP_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
static mut Pardiso_handle: soHandle_t = OSQP_NULL as soHandle_t;
static mut func_pardiso: pardiso_t = unsafe {
    ::std::mem::transmute::<isize, pardiso_t>(OSQP_NULL as isize)
};
static mut func_mkl_set_interface_layer: mkl_set_ifl_t = unsafe {
    ::std::mem::transmute::<isize, mkl_set_ifl_t>(OSQP_NULL as isize)
};
static mut func_mkl_get_max_threads: mkl_get_mt_t = unsafe {
    ::std::mem::transmute::<isize, mkl_get_mt_t>(OSQP_NULL as isize)
};
#[no_mangle]
pub unsafe extern "C" fn pardiso(
    mut pt: *mut *mut ::std::os::raw::c_void,
    mut maxfct: *const c_int,
    mut mnum: *const c_int,
    mut mtype: *const c_int,
    mut phase: *const c_int,
    mut n: *const c_int,
    mut a: *const c_float,
    mut ia: *const c_int,
    mut ja: *const c_int,
    mut perm: *mut c_int,
    mut nrhs: *const c_int,
    mut iparm: *mut c_int,
    mut msglvl: *const c_int,
    mut b: *mut c_float,
    mut x: *mut c_float,
    mut error: *mut c_int,
) {
    if func_pardiso.is_some() {
        func_pardiso
            .expect(
                "non-null function pointer",
            )(
            pt,
            maxfct,
            mnum,
            mtype,
            phase,
            n,
            a,
            ia,
            ja,
            perm,
            nrhs,
            iparm,
            msglvl,
            b,
            x,
            error,
        );
    } else {
        printf(
            b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
            (*::std::mem::transmute::<&[u8; 8], &[::std::os::raw::c_char; 8]>(b"pardiso\0"))
                .as_ptr(),
        );
        printf(b"Pardiso not loaded correctly\0" as *const u8 as *const ::std::os::raw::c_char);
        printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
    };
}
#[no_mangle]
pub unsafe extern "C" fn mkl_set_interface_layer(mut code: c_int) -> c_int {
    return func_mkl_set_interface_layer
        .expect("non-null function pointer")(code as ::std::os::raw::c_int) as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn mkl_get_max_threads() -> c_int {
    return ::std::mem::transmute::<
        _,
        fn() -> ::std::os::raw::c_int,
    >(func_mkl_get_max_threads.expect("non-null function pointer"))() as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn lh_load_pardiso(mut libname: *const ::std::os::raw::c_char) -> c_int {
    if !libname.is_null() {
        Pardiso_handle = lh_load_lib(libname);
    } else {
        Pardiso_handle = lh_load_lib(
            b"libmkl_rt.dylib\0" as *const u8 as *const ::std::os::raw::c_char,
        );
    }
    if Pardiso_handle.is_null() {
        return 1 as ::std::os::raw::c_int as c_int;
    }
    func_pardiso = ::std::mem::transmute::<
        voidfun,
        pardiso_t,
    >(lh_load_sym(Pardiso_handle, b"pardiso\0" as *const u8 as *const ::std::os::raw::c_char));
    if func_pardiso.is_none() {
        return 1 as ::std::os::raw::c_int as c_int;
    }
    func_mkl_set_interface_layer = ::std::mem::transmute::<
        voidfun,
        mkl_set_ifl_t,
    >(
        lh_load_sym(
            Pardiso_handle,
            b"MKL_Set_Interface_Layer\0" as *const u8 as *const ::std::os::raw::c_char,
        ),
    );
    if func_mkl_set_interface_layer.is_none() {
        return 1 as ::std::os::raw::c_int as c_int;
    }
    func_mkl_get_max_threads = ::std::mem::transmute::<
        voidfun,
        mkl_get_mt_t,
    >(
        lh_load_sym(
            Pardiso_handle,
            b"MKL_Get_Max_Threads\0" as *const u8 as *const ::std::os::raw::c_char,
        ),
    );
    if func_mkl_get_max_threads.is_none() {
        return 1 as ::std::os::raw::c_int as c_int;
    }
    return 0 as ::std::os::raw::c_int as c_int;
}
#[no_mangle]
pub unsafe extern "C" fn lh_unload_pardiso() -> c_int {
    if Pardiso_handle.is_null() {
        return 0 as ::std::os::raw::c_int as c_int;
    }
    return lh_unload_lib(Pardiso_handle);
}
