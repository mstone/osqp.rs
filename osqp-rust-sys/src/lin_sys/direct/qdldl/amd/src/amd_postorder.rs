use ::libc;
extern "C" {
    fn amd_l_post_tree(
        root: libc::c_longlong,
        k: libc::c_longlong,
        Child: *mut libc::c_longlong,
        Sibling: *const libc::c_longlong,
        Order: *mut libc::c_longlong,
        Stack: *mut libc::c_longlong,
    ) -> libc::c_longlong;
}
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
pub const AMD_post_tree: unsafe extern "C" fn(
    libc::c_longlong,
    libc::c_longlong,
    *mut libc::c_longlong,
    *const libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
) -> libc::c_longlong = amd_l_post_tree;
#[no_mangle]
pub unsafe extern "C" fn amd_l_postorder(
    mut nn: libc::c_longlong,
    mut Parent: *mut libc::c_longlong,
    mut Nv: *mut libc::c_longlong,
    mut Fsize: *mut libc::c_longlong,
    mut Order: *mut libc::c_longlong,
    mut Child: *mut libc::c_longlong,
    mut Sibling: *mut libc::c_longlong,
    mut Stack: *mut libc::c_longlong,
) {
    let mut i: libc::c_longlong = 0;
    let mut j: libc::c_longlong = 0;
    let mut k: libc::c_longlong = 0;
    let mut parent: libc::c_longlong = 0;
    let mut frsize: libc::c_longlong = 0;
    let mut f: libc::c_longlong = 0;
    let mut fprev: libc::c_longlong = 0;
    let mut maxfrsize: libc::c_longlong = 0;
    let mut bigfprev: libc::c_longlong = 0;
    let mut bigf: libc::c_longlong = 0;
    let mut fnext: libc::c_longlong = 0;
    j = 0 as libc::c_int as libc::c_longlong;
    while j < nn {
        *Child.offset(j as isize) = EMPTY as libc::c_longlong;
        *Sibling.offset(j as isize) = EMPTY as libc::c_longlong;
        j += 1;
    }
    j = nn - 1 as libc::c_int as libc::c_longlong;
    while j >= 0 as libc::c_int as libc::c_longlong {
        if *Nv.offset(j as isize) > 0 as libc::c_int as libc::c_longlong {
            parent = *Parent.offset(j as isize);
            if parent != EMPTY as libc::c_longlong {
                *Sibling.offset(j as isize) = *Child.offset(parent as isize);
                *Child.offset(parent as isize) = j;
            }
        }
        j -= 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < nn {
        if *Nv.offset(i as isize) > 0 as libc::c_int as libc::c_longlong
            && *Child.offset(i as isize) != EMPTY as libc::c_longlong
        {
            fprev = EMPTY as libc::c_longlong;
            maxfrsize = EMPTY as libc::c_longlong;
            bigfprev = EMPTY as libc::c_longlong;
            bigf = EMPTY as libc::c_longlong;
            f = *Child.offset(i as isize);
            while f != EMPTY as libc::c_longlong {
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
            if fnext != EMPTY as libc::c_longlong {
                if bigfprev == EMPTY as libc::c_longlong {
                    *Child.offset(i as isize) = fnext;
                } else {
                    *Sibling.offset(bigfprev as isize) = fnext;
                }
                *Sibling.offset(bigf as isize) = EMPTY as libc::c_longlong;
                *Sibling.offset(fprev as isize) = bigf;
            }
        }
        i += 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < nn {
        *Order.offset(i as isize) = EMPTY as libc::c_longlong;
        i += 1;
    }
    k = 0 as libc::c_int as libc::c_longlong;
    i = 0 as libc::c_int as libc::c_longlong;
    while i < nn {
        if *Parent.offset(i as isize) == EMPTY as libc::c_longlong
            && *Nv.offset(i as isize) > 0 as libc::c_int as libc::c_longlong
        {
            k = amd_l_post_tree(
                i,
                k,
                Child,
                Sibling as *const libc::c_longlong,
                Order,
                Stack,
            );
        }
        i += 1;
    }
}
