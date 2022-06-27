pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
#[no_mangle]
pub unsafe extern "C" fn amd_l_post_tree(
    mut root: ::std::os::raw::c_longlong,
    mut k: ::std::os::raw::c_longlong,
    mut Child: *mut ::std::os::raw::c_longlong,
    mut Sibling: *const ::std::os::raw::c_longlong,
    mut Order: *mut ::std::os::raw::c_longlong,
    mut Stack: *mut ::std::os::raw::c_longlong,
) -> ::std::os::raw::c_longlong {
    let mut f: ::std::os::raw::c_longlong = 0;
    let mut head: ::std::os::raw::c_longlong = 0;
    let mut h: ::std::os::raw::c_longlong = 0;
    let mut i: ::std::os::raw::c_longlong = 0;
    head = 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong;
    *Stack.offset(0 as ::std::os::raw::c_int as isize) = root;
    while head >= 0 as ::std::os::raw::c_int as ::std::os::raw::c_longlong {
        i = *Stack.offset(head as isize);
        if *Child.offset(i as isize) != EMPTY as ::std::os::raw::c_longlong {
            f = *Child.offset(i as isize);
            while f != EMPTY as ::std::os::raw::c_longlong {
                head += 1;
                f = *Sibling.offset(f as isize);
            }
            h = head;
            f = *Child.offset(i as isize);
            while f != EMPTY as ::std::os::raw::c_longlong {
                let fresh0 = h;
                h = h - 1;
                *Stack.offset(fresh0 as isize) = f;
                f = *Sibling.offset(f as isize);
            }
            *Child.offset(i as isize) = EMPTY as ::std::os::raw::c_longlong;
        } else {
            head -= 1;
            let fresh1 = k;
            k = k + 1;
            *Order.offset(i as isize) = fresh1;
        }
    }
    return k;
}
