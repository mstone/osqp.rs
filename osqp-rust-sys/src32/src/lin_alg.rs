extern "C" {
    fn sqrt(_: ::std::os::raw::c_double) -> ::std::os::raw::c_double;
    fn malloc(_: ::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void;
    fn printf(_: *const ::std::os::raw::c_char, _: ...) -> ::std::os::raw::c_int;
}
pub type c_int = ::std::os::raw::c_int;
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
pub const OSQP_NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const c_malloc: unsafe extern "C" fn(::std::os::raw::c_ulong) -> *mut ::std::os::raw::c_void = malloc;
pub const c_print: unsafe extern "C" fn(*const ::std::os::raw::c_char, ...) -> ::std::os::raw::c_int = printf;
pub const c_sqrt: unsafe extern "C" fn(::std::os::raw::c_double) -> ::std::os::raw::c_double = sqrt;
#[no_mangle]
pub unsafe extern "C" fn vec_add_scaled(
    mut c: *mut c_float,
    mut a: *const c_float,
    mut b: *const c_float,
    mut n: c_int,
    mut sc: c_float,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *c.offset(i as isize) = *a.offset(i as isize) + sc * *b.offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_scaled_norm_inf(
    mut S: *const c_float,
    mut v: *const c_float,
    mut l: c_int,
) -> c_float {
    let mut i: c_int = 0;
    let mut abs_Sv_i: c_float = 0.;
    let mut max: c_float = 0.0f64;
    i = 0 as ::std::os::raw::c_int;
    while i < l {
        abs_Sv_i = if *S.offset(i as isize) * *v.offset(i as isize)
            < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
        {
            -(*S.offset(i as isize) * *v.offset(i as isize))
        } else {
            *S.offset(i as isize) * *v.offset(i as isize)
        };
        if abs_Sv_i > max {
            max = abs_Sv_i;
        }
        i += 1;
    }
    return max;
}
#[no_mangle]
pub unsafe extern "C" fn vec_norm_inf(mut v: *const c_float, mut l: c_int) -> c_float {
    let mut i: c_int = 0;
    let mut abs_v_i: c_float = 0.;
    let mut max: c_float = 0.0f64;
    i = 0 as ::std::os::raw::c_int;
    while i < l {
        abs_v_i = if *v.offset(i as isize) < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double {
            -*v.offset(i as isize)
        } else {
            *v.offset(i as isize)
        };
        if abs_v_i > max {
            max = abs_v_i;
        }
        i += 1;
    }
    return max;
}
#[no_mangle]
pub unsafe extern "C" fn vec_norm_inf_diff(
    mut a: *const c_float,
    mut b: *const c_float,
    mut l: c_int,
) -> c_float {
    let mut nmDiff: c_float = 0.0f64;
    let mut tmp: c_float = 0.;
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < l {
        tmp = if *a.offset(i as isize) - *b.offset(i as isize)
            < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
        {
            -(*a.offset(i as isize) - *b.offset(i as isize))
        } else {
            *a.offset(i as isize) - *b.offset(i as isize)
        };
        if tmp > nmDiff {
            nmDiff = tmp;
        }
        i += 1;
    }
    return nmDiff;
}
#[no_mangle]
pub unsafe extern "C" fn vec_mean(mut a: *const c_float, mut n: c_int) -> c_float {
    let mut mean: c_float = 0.0f64;
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        mean += *a.offset(i as isize);
        i += 1;
    }
    mean /= n as c_float;
    return mean;
}
#[no_mangle]
pub unsafe extern "C" fn int_vec_set_scalar(
    mut a: *mut c_int,
    mut sc: c_int,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *a.offset(i as isize) = sc;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_set_scalar(
    mut a: *mut c_float,
    mut sc: c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *a.offset(i as isize) = sc;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_add_scalar(
    mut a: *mut c_float,
    mut sc: c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        let ref mut fresh0 = *a.offset(i as isize);
        *fresh0 += sc;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_mult_scalar(
    mut a: *mut c_float,
    mut sc: c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        let ref mut fresh1 = *a.offset(i as isize);
        *fresh1 *= sc;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_copy(mut a: *mut c_float, mut n: c_int) -> *mut c_float {
    let mut b: *mut c_float = 0 as *mut c_float;
    let mut i: c_int = 0;
    b = malloc(
        (n as ::std::os::raw::c_ulong)
            .wrapping_mul(::std::mem::size_of::<c_float>() as ::std::os::raw::c_ulong),
    ) as *mut c_float;
    if b.is_null() {
        return OSQP_NULL as *mut c_float;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *b.offset(i as isize) = *a.offset(i as isize);
        i += 1;
    }
    return b;
}
#[no_mangle]
pub unsafe extern "C" fn prea_int_vec_copy(
    mut a: *const c_int,
    mut b: *mut c_int,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *b.offset(i as isize) = *a.offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn prea_vec_copy(
    mut a: *const c_float,
    mut b: *mut c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *b.offset(i as isize) = *a.offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_recipr(
    mut a: *const c_float,
    mut b: *mut c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *b.offset(i as isize) = 1.0f64 / *a.offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_prod(
    mut a: *const c_float,
    mut b: *const c_float,
    mut n: c_int,
) -> c_float {
    let mut prod: c_float = 0.0f64;
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        prod += *a.offset(i as isize) * *b.offset(i as isize);
        i += 1;
    }
    return prod;
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_prod(
    mut a: *const c_float,
    mut b: *const c_float,
    mut c: *mut c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *c.offset(i as isize) = *b.offset(i as isize) * *a.offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_sqrt(mut a: *mut c_float, mut n: c_int) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *a.offset(i as isize) = sqrt(*a.offset(i as isize));
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_max(
    mut a: *mut c_float,
    mut n: c_int,
    mut max_val: c_float,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *a
            .offset(
                i as isize,
            ) = if *a.offset(i as isize) > max_val {
            *a.offset(i as isize)
        } else {
            max_val
        };
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_min(
    mut a: *mut c_float,
    mut n: c_int,
    mut min_val: c_float,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *a
            .offset(
                i as isize,
            ) = if *a.offset(i as isize) < min_val {
            *a.offset(i as isize)
        } else {
            min_val
        };
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_max_vec(
    mut a: *const c_float,
    mut b: *const c_float,
    mut c: *mut c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *c
            .offset(
                i as isize,
            ) = if *a.offset(i as isize) > *b.offset(i as isize) {
            *a.offset(i as isize)
        } else {
            *b.offset(i as isize)
        };
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn vec_ew_min_vec(
    mut a: *const c_float,
    mut b: *const c_float,
    mut c: *mut c_float,
    mut n: c_int,
) {
    let mut i: c_int = 0;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *c
            .offset(
                i as isize,
            ) = if *a.offset(i as isize) < *b.offset(i as isize) {
            *a.offset(i as isize)
        } else {
            *b.offset(i as isize)
        };
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn mat_mult_scalar(mut A: *mut csc, mut sc: c_float) {
    let mut i: c_int = 0;
    let mut nnzA: c_int = 0;
    nnzA = *((*A).p).offset((*A).n as isize);
    i = 0 as ::std::os::raw::c_int;
    while i < nnzA {
        let ref mut fresh2 = *((*A).x).offset(i as isize);
        *fresh2 *= sc;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn mat_premult_diag(mut A: *mut csc, mut d: *const c_float) {
    let mut j: c_int = 0;
    let mut i: c_int = 0;
    j = 0 as ::std::os::raw::c_int;
    while j < (*A).n {
        i = *((*A).p).offset(j as isize);
        while i < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            let ref mut fresh3 = *((*A).x).offset(i as isize);
            *fresh3 *= *d.offset(*((*A).i).offset(i as isize) as isize);
            i += 1;
        }
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn mat_postmult_diag(mut A: *mut csc, mut d: *const c_float) {
    let mut j: c_int = 0;
    let mut i: c_int = 0;
    j = 0 as ::std::os::raw::c_int;
    while j < (*A).n {
        i = *((*A).p).offset(j as isize);
        while i < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            let ref mut fresh4 = *((*A).x).offset(i as isize);
            *fresh4 *= *d.offset(j as isize);
            i += 1;
        }
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn mat_vec(
    mut A: *const csc,
    mut x: *const c_float,
    mut y: *mut c_float,
    mut plus_eq: c_int,
) {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    if plus_eq == 0 {
        i = 0 as ::std::os::raw::c_int;
        while i < (*A).m {
            *y.offset(i as isize) = 0 as ::std::os::raw::c_int as c_float;
            i += 1;
        }
    }
    if *((*A).p).offset((*A).n as isize) == 0 as ::std::os::raw::c_int {
        return;
    }
    if plus_eq == -(1 as ::std::os::raw::c_int) {
        j = 0 as ::std::os::raw::c_int;
        while j < (*A).n {
            i = *((*A).p).offset(j as isize);
            while i < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
                let ref mut fresh5 = *y.offset(*((*A).i).offset(i as isize) as isize);
                *fresh5 -= *((*A).x).offset(i as isize) * *x.offset(j as isize);
                i += 1;
            }
            j += 1;
        }
    } else {
        j = 0 as ::std::os::raw::c_int;
        while j < (*A).n {
            i = *((*A).p).offset(j as isize);
            while i < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
                let ref mut fresh6 = *y.offset(*((*A).i).offset(i as isize) as isize);
                *fresh6 += *((*A).x).offset(i as isize) * *x.offset(j as isize);
                i += 1;
            }
            j += 1;
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn mat_tpose_vec(
    mut A: *const csc,
    mut x: *const c_float,
    mut y: *mut c_float,
    mut plus_eq: c_int,
    mut skip_diag: c_int,
) {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut k: c_int = 0;
    if plus_eq == 0 {
        i = 0 as ::std::os::raw::c_int;
        while i < (*A).n {
            *y.offset(i as isize) = 0 as ::std::os::raw::c_int as c_float;
            i += 1;
        }
    }
    if *((*A).p).offset((*A).n as isize) == 0 as ::std::os::raw::c_int {
        return;
    }
    if plus_eq == -(1 as ::std::os::raw::c_int) {
        if skip_diag != 0 {
            j = 0 as ::std::os::raw::c_int;
            while j < (*A).n {
                k = *((*A).p).offset(j as isize);
                while k < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
                    i = *((*A).i).offset(k as isize);
                    let ref mut fresh7 = *y.offset(j as isize);
                    *fresh7
                        -= if i == j {
                            0 as ::std::os::raw::c_int as ::std::os::raw::c_double
                        } else {
                            *((*A).x).offset(k as isize) * *x.offset(i as isize)
                        };
                    k += 1;
                }
                j += 1;
            }
        } else {
            j = 0 as ::std::os::raw::c_int;
            while j < (*A).n {
                k = *((*A).p).offset(j as isize);
                while k < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
                    let ref mut fresh8 = *y.offset(j as isize);
                    *fresh8
                        -= *((*A).x).offset(k as isize)
                            * *x.offset(*((*A).i).offset(k as isize) as isize);
                    k += 1;
                }
                j += 1;
            }
        }
    } else if skip_diag != 0 {
        j = 0 as ::std::os::raw::c_int;
        while j < (*A).n {
            k = *((*A).p).offset(j as isize);
            while k < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
                i = *((*A).i).offset(k as isize);
                let ref mut fresh9 = *y.offset(j as isize);
                *fresh9
                    += if i == j {
                        0 as ::std::os::raw::c_int as ::std::os::raw::c_double
                    } else {
                        *((*A).x).offset(k as isize) * *x.offset(i as isize)
                    };
                k += 1;
            }
            j += 1;
        }
    } else {
        j = 0 as ::std::os::raw::c_int;
        while j < (*A).n {
            k = *((*A).p).offset(j as isize);
            while k < *((*A).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
                let ref mut fresh10 = *y.offset(j as isize);
                *fresh10
                    += *((*A).x).offset(k as isize)
                        * *x.offset(*((*A).i).offset(k as isize) as isize);
                k += 1;
            }
            j += 1;
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn mat_inf_norm_cols(mut M: *const csc, mut E: *mut c_float) {
    let mut j: c_int = 0;
    let mut ptr: c_int = 0;
    j = 0 as ::std::os::raw::c_int;
    while j < (*M).n {
        *E.offset(j as isize) = 0.0f64;
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < (*M).n {
        ptr = *((*M).p).offset(j as isize);
        while ptr < *((*M).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            *E
                .offset(
                    j as isize,
                ) = if (if *((*M).x).offset(ptr as isize)
                < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
            {
                -*((*M).x).offset(ptr as isize)
            } else {
                *((*M).x).offset(ptr as isize)
            }) > *E.offset(j as isize)
            {
                if *((*M).x).offset(ptr as isize) < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double {
                    -*((*M).x).offset(ptr as isize)
                } else {
                    *((*M).x).offset(ptr as isize)
                }
            } else {
                *E.offset(j as isize)
            };
            ptr += 1;
        }
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn mat_inf_norm_rows(mut M: *const csc, mut E: *mut c_float) {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut ptr: c_int = 0;
    j = 0 as ::std::os::raw::c_int;
    while j < (*M).m {
        *E.offset(j as isize) = 0.0f64;
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < (*M).n {
        ptr = *((*M).p).offset(j as isize);
        while ptr < *((*M).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            i = *((*M).i).offset(ptr as isize);
            *E
                .offset(
                    i as isize,
                ) = if (if *((*M).x).offset(ptr as isize)
                < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
            {
                -*((*M).x).offset(ptr as isize)
            } else {
                *((*M).x).offset(ptr as isize)
            }) > *E.offset(i as isize)
            {
                if *((*M).x).offset(ptr as isize) < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double {
                    -*((*M).x).offset(ptr as isize)
                } else {
                    *((*M).x).offset(ptr as isize)
                }
            } else {
                *E.offset(i as isize)
            };
            ptr += 1;
        }
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn mat_inf_norm_cols_sym_triu(
    mut M: *const csc,
    mut E: *mut c_float,
) {
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut ptr: c_int = 0;
    let mut abs_x: c_float = 0.;
    j = 0 as ::std::os::raw::c_int;
    while j < (*M).n {
        *E.offset(j as isize) = 0.0f64;
        j += 1;
    }
    j = 0 as ::std::os::raw::c_int;
    while j < (*M).n {
        ptr = *((*M).p).offset(j as isize);
        while ptr < *((*M).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            i = *((*M).i).offset(ptr as isize);
            abs_x = if *((*M).x).offset(ptr as isize)
                < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double
            {
                -*((*M).x).offset(ptr as isize)
            } else {
                *((*M).x).offset(ptr as isize)
            };
            *E
                .offset(
                    j as isize,
                ) = if abs_x > *E.offset(j as isize) {
                abs_x
            } else {
                *E.offset(j as isize)
            };
            if i != j {
                *E
                    .offset(
                        i as isize,
                    ) = if abs_x > *E.offset(i as isize) {
                    abs_x
                } else {
                    *E.offset(i as isize)
                };
            }
            ptr += 1;
        }
        j += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn quad_form(mut P: *const csc, mut x: *const c_float) -> c_float {
    let mut quad_form_0: c_float = 0.0f64;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut ptr: c_int = 0;
    j = 0 as ::std::os::raw::c_int;
    while j < (*P).n {
        ptr = *((*P).p).offset(j as isize);
        while ptr < *((*P).p).offset((j + 1 as ::std::os::raw::c_int) as isize) {
            i = *((*P).i).offset(ptr as isize);
            if i == j {
                quad_form_0
                    += 0.5f64 * *((*P).x).offset(ptr as isize) * *x.offset(i as isize)
                        * *x.offset(i as isize);
            } else if i < j {
                quad_form_0
                    += *((*P).x).offset(ptr as isize) * *x.offset(i as isize)
                        * *x.offset(j as isize);
            } else {
                printf(
                    b"ERROR in %s: \0" as *const u8 as *const ::std::os::raw::c_char,
                    (*::std::mem::transmute::<
                        &[u8; 10],
                        &[::std::os::raw::c_char; 10],
                    >(b"quad_form\0"))
                        .as_ptr(),
                );
                printf(
                    b"quad_form matrix is not upper triangular\0" as *const u8
                        as *const ::std::os::raw::c_char,
                );
                printf(b"\n\0" as *const u8 as *const ::std::os::raw::c_char);
                return OSQP_NULL as c_float;
            }
            ptr += 1;
        }
        j += 1;
    }
    return quad_form_0;
}
