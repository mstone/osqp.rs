extern "C" {
    fn amd_post_tree(
        root: ::std::os::raw::c_int,
        k: ::std::os::raw::c_int,
        Child: *mut ::std::os::raw::c_int,
        Sibling: *const ::std::os::raw::c_int,
        Order: *mut ::std::os::raw::c_int,
        Stack: *mut ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}
pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const AMD_post_tree: unsafe extern "C" fn(
    ::std::os::raw::c_int,
    ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *const ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
) -> ::std::os::raw::c_int = amd_post_tree;
#[no_mangle]
pub unsafe extern "C" fn amd_postorder(
    mut nn: ::std::os::raw::c_int,
    mut Parent: *mut ::std::os::raw::c_int,
    mut Nv: *mut ::std::os::raw::c_int,
    mut Fsize: *mut ::std::os::raw::c_int,
    mut Order: *mut ::std::os::raw::c_int,
    mut Child: *mut ::std::os::raw::c_int,
    mut Sibling: *mut ::std::os::raw::c_int,
    mut Stack: *mut ::std::os::raw::c_int,
) {
    let mut i: ::std::os::raw::c_int = 0;
    let mut j: ::std::os::raw::c_int = 0;
    let mut k: ::std::os::raw::c_int = 0;
    let mut parent: ::std::os::raw::c_int = 0;
    let mut frsize: ::std::os::raw::c_int = 0;
    let mut f: ::std::os::raw::c_int = 0;
    let mut fprev: ::std::os::raw::c_int = 0;
    let mut maxfrsize: ::std::os::raw::c_int = 0;
    let mut bigfprev: ::std::os::raw::c_int = 0;
    let mut bigf: ::std::os::raw::c_int = 0;
    let mut fnext: ::std::os::raw::c_int = 0;
    j = 0 as ::std::os::raw::c_int;
    while j < nn {
        *Child.offset(j as isize) = EMPTY;
        *Sibling.offset(j as isize) = EMPTY;
        j += 1;
    }
    j = nn - 1 as ::std::os::raw::c_int;
    while j >= 0 as ::std::os::raw::c_int {
        if *Nv.offset(j as isize) > 0 as ::std::os::raw::c_int {
            parent = *Parent.offset(j as isize);
            if parent != EMPTY {
                *Sibling.offset(j as isize) = *Child.offset(parent as isize);
                *Child.offset(parent as isize) = j;
            }
        }
        j -= 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < nn {
        if *Nv.offset(i as isize) > 0 as ::std::os::raw::c_int
            && *Child.offset(i as isize) != EMPTY
        {
            fprev = EMPTY;
            maxfrsize = EMPTY;
            bigfprev = EMPTY;
            bigf = EMPTY;
            f = *Child.offset(i as isize);
            while f != EMPTY {
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
            if fnext != EMPTY {
                if bigfprev == EMPTY {
                    *Child.offset(i as isize) = fnext;
                } else {
                    *Sibling.offset(bigfprev as isize) = fnext;
                }
                *Sibling.offset(bigf as isize) = EMPTY;
                *Sibling.offset(fprev as isize) = bigf;
            }
        }
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < nn {
        *Order.offset(i as isize) = EMPTY;
        i += 1;
    }
    k = 0 as ::std::os::raw::c_int;
    i = 0 as ::std::os::raw::c_int;
    while i < nn {
        if *Parent.offset(i as isize) == EMPTY
            && *Nv.offset(i as isize) > 0 as ::std::os::raw::c_int
        {
            k = amd_post_tree(i, k, Child, Sibling as *const ::std::os::raw::c_int, Order, Stack);
        }
        i += 1;
    }
}
