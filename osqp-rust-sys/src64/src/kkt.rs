extern "C" {
    fn realloc(_: *mut ::std::os::raw::c_void, _: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn free(_: *mut ::std::os::raw::c_void);
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn triplet_to_csr(T: *const csc, TtoC: *mut c_int) -> *mut csc;
    fn triplet_to_csc(T: *const csc, TtoC: *mut c_int) -> *mut csc;
    fn csc_spfree(A: *mut csc);
    fn csc_spalloc(
        m: c_int,
        n: c_int,
        nzmax: c_int,
        values: c_int,
        triplet: c_int,
    ) -> *mut csc;
}
pub type c_int = ::std::os::raw::c_longlong;
pub type c_float = ::std::os::raw::c_double;
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
pub const c_realloc: unsafe extern "C" fn(
    *mut ::std::os::raw::c_void,
    ::std::os::raw::c_ulong,
) -> *mut ::std::os::raw::c_void = realloc;
pub const c_malloc: unsafe extern "C" fn(::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void = malloc;
pub const OSQP_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const c_free: unsafe extern "C" fn(*mut ::std::os::raw::c_void) -> () = free;
#[no_mangle]
pub unsafe extern "C" fn form_KKT(
    mut P: *const csc,
    mut A: *const csc,
    mut format: c_int,
    mut param1: c_float,
    mut param2: *mut c_float,
    mut PtoKKT: *mut c_int,
    mut AtoKKT: *mut c_int,
    mut Pdiag_idx: *mut *mut c_int,
    mut Pdiag_n: *mut c_int,
    mut param2toKKT: *mut c_int,
) -> *mut csc {
    let mut nKKT: c_int = 0;
    let mut nnzKKTmax: c_int = 0;
    let mut KKT_trip: *mut csc = 0 as *mut csc;
    let mut KKT: *mut csc = 0 as *mut csc;
    let mut ptr: c_int = 0;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut zKKT: c_int = 0 as ::std::os::raw::c_int as c_int;
    let mut KKT_TtoC: *mut c_int = 0 as *mut c_int;
    nKKT = (*P).m + (*A).m;
    nnzKKTmax = *((*P).p).offset((*P).n as isize) + (*P).m
        + *((*A).p).offset((*A).n as isize) + (*A).m;
    KKT_trip = csc_spalloc(
        nKKT,
        nKKT,
        nnzKKTmax,
        1 as ::std::os::raw::c_int as c_int,
        1 as ::std::os::raw::c_int as c_int,
    );
    if KKT_trip.is_null() {
        return OSQP_NULL as *mut csc;
    }
    if !Pdiag_idx.is_null() {
        *Pdiag_idx = malloc(
            ((*P).m as ::std::os::raw::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as ::std::os::raw::c_ulong as ::std::os::raw::c_ulonglong,
                ) as ::std::os::raw::c_ulong,
        ) as *mut c_int;
        *Pdiag_n = 0 as ::std::os::raw::c_int as c_int;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*P).n {
        if *((*P).p).offset(j as isize)
            == *((*P).p).offset((j + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong) as isize)
        {
            *((*KKT_trip).i).offset(zKKT as isize) = j;
            *((*KKT_trip).p).offset(zKKT as isize) = j;
            *((*KKT_trip).x).offset(zKKT as isize) = param1;
            zKKT += 1;
        }
        ptr = *((*P).p).offset(j as isize);
        while ptr < *((*P).p).offset((j + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong) as isize)
        {
            i = *((*P).i).offset(ptr as isize);
            *((*KKT_trip).i).offset(zKKT as isize) = i;
            *((*KKT_trip).p).offset(zKKT as isize) = j;
            *((*KKT_trip).x).offset(zKKT as isize) = *((*P).x).offset(ptr as isize);
            if !PtoKKT.is_null() {
                *PtoKKT.offset(ptr as isize) = zKKT;
            }
            if i == j {
                let ref mut fresh0 = *((*KKT_trip).x).offset(zKKT as isize);
                *fresh0 += param1;
                if !Pdiag_idx.is_null() {
                    *(*Pdiag_idx).offset(*Pdiag_n as isize) = ptr;
                    *Pdiag_n += 1;
                }
            }
            zKKT += 1;
            if i < j
                && ptr + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
                    == *((*P).p)
                        .offset((j + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong) as isize)
            {
                *((*KKT_trip).i).offset(zKKT as isize) = j;
                *((*KKT_trip).p).offset(zKKT as isize) = j;
                *((*KKT_trip).x).offset(zKKT as isize) = param1;
                zKKT += 1;
            }
            ptr += 1;
        }
        j += 1;
    }
    if !Pdiag_idx.is_null() {
        *Pdiag_idx = realloc(
            *Pdiag_idx as *mut ::std::os::raw::c_void,
            (*Pdiag_n as ::std::os::raw::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as ::std::os::raw::c_ulong as ::std::os::raw::c_ulonglong,
                ) as ::std::os::raw::c_ulong,
        ) as *mut c_int;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*A).n {
        ptr = *((*A).p).offset(j as isize);
        while ptr < *((*A).p).offset((j + 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong) as isize)
        {
            *((*KKT_trip).p)
                .offset(zKKT as isize) = (*P).m + *((*A).i).offset(ptr as isize);
            *((*KKT_trip).i).offset(zKKT as isize) = j;
            *((*KKT_trip).x).offset(zKKT as isize) = *((*A).x).offset(ptr as isize);
            if !AtoKKT.is_null() {
                *AtoKKT.offset(ptr as isize) = zKKT;
            }
            zKKT += 1;
            ptr += 1;
        }
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int as c_int;
    while j < (*A).m {
        *((*KKT_trip).i).offset(zKKT as isize) = j + (*P).n;
        *((*KKT_trip).p).offset(zKKT as isize) = j + (*P).n;
        *((*KKT_trip).x).offset(zKKT as isize) = -*param2.offset(j as isize);
        if !param2toKKT.is_null() {
            *param2toKKT.offset(j as isize) = zKKT;
        }
        zKKT += 1;
        j += 1;
    }
    (*KKT_trip).nz = zKKT;
    if PtoKKT.is_null() && AtoKKT.is_null() && param2toKKT.is_null() {
        if format == 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
            KKT = triplet_to_csc(KKT_trip, OSQP_NULL as *mut c_int);
        } else {
            KKT = triplet_to_csr(KKT_trip, OSQP_NULL as *mut c_int);
        }
    } else {
        KKT_TtoC = malloc(
            (zKKT as ::std::os::raw::c_ulonglong)
                .wrapping_mul(
                    ::std::mem::size_of::<c_int>() as ::std::os::raw::c_ulong as ::std::os::raw::c_ulonglong,
                ) as ::std::os::raw::c_ulong,
        ) as *mut c_int;
        if KKT_TtoC.is_null() {
            csc_spfree(KKT_trip);
            free(*Pdiag_idx as *mut ::std::os::raw::c_void);
            return OSQP_NULL as *mut csc;
        }
        if format == 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
            KKT = triplet_to_csc(KKT_trip, KKT_TtoC);
        } else {
            KKT = triplet_to_csr(KKT_trip, KKT_TtoC);
        }
        if !PtoKKT.is_null() {
            i = 0 as ::std::os::raw::c_int as c_int;
            while i < *((*P).p).offset((*P).n as isize) {
                *PtoKKT
                    .offset(
                        i as isize,
                    ) = *KKT_TtoC.offset(*PtoKKT.offset(i as isize) as isize);
                i += 1;
            }
        }
        if !AtoKKT.is_null() {
            i = 0 as ::std::os::raw::c_int as c_int;
            while i < *((*A).p).offset((*A).n as isize) {
                *AtoKKT
                    .offset(
                        i as isize,
                    ) = *KKT_TtoC.offset(*AtoKKT.offset(i as isize) as isize);
                i += 1;
            }
        }
        if !param2toKKT.is_null() {
            i = 0 as ::std::os::raw::c_int as c_int;
            while i < (*A).m {
                *param2toKKT
                    .offset(
                        i as isize,
                    ) = *KKT_TtoC.offset(*param2toKKT.offset(i as isize) as isize);
                i += 1;
            }
        }
        free(KKT_TtoC as *mut ::std::os::raw::c_void);
    }
    csc_spfree(KKT_trip);
    return KKT;
}
#[no_mangle]
pub unsafe extern "C" fn update_KKT_P(
    mut KKT: *mut csc,
    mut P: *const csc,
    mut PtoKKT: *const c_int,
    param1: c_float,
    mut Pdiag_idx: *const c_int,
    Pdiag_n: c_int,
) {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    i = 0 as ::std::os::raw::c_int as c_int;
    while i < *((*P).p).offset((*P).n as isize) {
        *((*KKT).x)
            .offset(*PtoKKT.offset(i as isize) as isize) = *((*P).x).offset(i as isize);
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int as c_int;
    while i < Pdiag_n {
        j = *Pdiag_idx.offset(i as isize);
        let ref mut fresh1 = *((*KKT).x).offset(*PtoKKT.offset(j as isize) as isize);
        *fresh1 += param1;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn update_KKT_A(
    mut KKT: *mut csc,
    mut A: *const csc,
    mut AtoKKT: *const c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int as c_int;
    while i < *((*A).p).offset((*A).n as isize) {
        *((*KKT).x)
            .offset(*AtoKKT.offset(i as isize) as isize) = *((*A).x).offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn update_KKT_param2(
    mut KKT: *mut csc,
    mut param2: *const c_float,
    mut param2toKKT: *const c_int,
    m: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int as c_int;
    while i < m {
        *((*KKT).x)
            .offset(
                *param2toKKT.offset(i as isize) as isize,
            ) = -*param2.offset(i as isize);
        i += 1;
    }
}
