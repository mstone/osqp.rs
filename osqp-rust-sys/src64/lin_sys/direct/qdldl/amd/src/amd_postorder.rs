extern "C" {
    fn amd_l_post_tree(
        root: ::std::os::raw::c_longlong,
        k: ::std::os::raw::c_longlong,
        Child: *mut ::std::os::raw::c_longlong,
        Sibling: *const ::std::os::raw::c_longlong,
        Order: *mut ::std::os::raw::c_longlong,
        Stack: *mut ::std::os::raw::c_longlong,
    ) -> ::std::os::raw::c_longlong;
}
pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const AMD_post_tree: unsafe extern "C" fn(
    ::std::os::raw::c_longlong,
    ::std::os::raw::c_longlong,
    *mut ::std::os::raw::c_longlong,
    *const ::std::os::raw::c_longlong,
    *mut ::std::os::raw::c_longlong,
    *mut ::std::os::raw::c_longlong,
) -> ::std::os::raw::c_longlong = amd_l_post_tree;
#[no_mangle]
pub unsafe extern "C" fn amd_l_postorder(
    mut nn: ::std::os::raw::c_longlong,
    mut Parent: *mut ::std::os::raw::c_longlong,
    mut Nv: *mut ::std::os::raw::c_longlong,
    mut Fsize: *mut ::std::os::raw::c_longlong,
    mut Order: *mut ::std::os::raw::c_longlong,
    mut Child: *mut ::std::os::raw::c_longlong,
    mut Sibling: *mut ::std::os::raw::c_longlong,
    mut Stack: *mut ::std::os::raw::c_longlong,
) {
    let mut i: ::std::os::raw::c_longlong = 0;
    let mut j: ::std::os::raw::c_longlong = 0;
    let mut k: ::std::os::raw::c_longlong = 0;
    let mut parent: ::std::os::raw::c_longlong = 0;
    let mut frsize: ::std::os::raw::c_longlong = 0;
    let mut f: ::std::os::raw::c_longlong = 0;
    let mut fprev: ::std::os::raw::c_longlong = 0;
    let mut maxfrsize: ::std::os::raw::c_longlong = 0;
    let mut bigfprev: ::std::os::raw::c_longlong = 0;
    let mut bigf: ::std::os::raw::c_longlong = 0;
    let mut fnext: ::std::os::raw::c_longlong = 0;
    j = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    while j < nn {
        *Child.offset(j as isize) = EMPTY as ::std::os::raw::c_longlong;
        *Sibling.offset(j as isize) = EMPTY as ::std::os::raw::c_longlong;
        j += 1;
    }
    j = nn - 1 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    while j >= 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
        if *Nv.offset(j as isize) > 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
            parent = *Parent.offset(j as isize);
            if parent != EMPTY as ::std::os::raw::c_longlong {
                *Sibling.offset(j as isize) = *Child.offset(parent as isize);
                *Child.offset(parent as isize) = j;
            }
        }
        j -= 1;
    }
    i = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    while i < nn {
        if *Nv.offset(i as isize) > 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
            && *Child.offset(i as isize) != EMPTY as ::std::os::raw::c_longlong
        {
            fprev = EMPTY as ::std::os::raw::c_longlong;
            maxfrsize = EMPTY as ::std::os::raw::c_longlong;
            bigfprev = EMPTY as ::std::os::raw::c_longlong;
            bigf = EMPTY as ::std::os::raw::c_longlong;
            f = *Child.offset(i as isize);
            while f != EMPTY as ::std::os::raw::c_longlong {
                frsize = *Fsize.offset(f as isize);
                if frsize >= maxfrsize {
                    maxfrsize = frsize;
                    bigfprev = fprev;
                    bigf = f;
                }
                fprev = f;
                f = *Sibling.offset(f as isize);
            }
            fnext = *Sibling.offset(bigf as isize);
            if fnext != EMPTY as ::std::os::raw::c_longlong {
                if bigfprev == EMPTY as ::std::os::raw::c_longlong {
                    *Child.offset(i as isize) = fnext;
                } else {
                    *Sibling.offset(bigfprev as isize) = fnext;
                }
                *Sibling.offset(bigf as isize) = EMPTY as ::std::os::raw::c_longlong;
                *Sibling.offset(fprev as isize) = bigf;
            }
        }
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    while i < nn {
        *Order.offset(i as isize) = EMPTY as ::std::os::raw::c_longlong;
        i += 1;
    }
    k = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    i = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    while i < nn {
        if *Parent.offset(i as isize) == EMPTY as ::std::os::raw::c_longlong
            && *Nv.offset(i as isize) > 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong
        {
            k = amd_l_post_tree(
                i,
                k,
                Child,
                Sibling as *const ::std::os::raw::c_longlong,
                Order,
                Stack,
            );
        }
        i += 1;
    }
}
