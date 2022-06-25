use ::libc;
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
#[no_mangle]
pub unsafe extern "C" fn amd_l_post_tree(
    mut root: libc::c_longlong,
    mut k: libc::c_longlong,
    mut Child: *mut libc::c_longlong,
    mut Sibling: *const libc::c_longlong,
    mut Order: *mut libc::c_longlong,
    mut Stack: *mut libc::c_longlong,
) -> libc::c_longlong {
    let mut f: libc::c_longlong = 0;
    let mut head: libc::c_longlong = 0;
    let mut h: libc::c_longlong = 0;
    let mut i: libc::c_longlong = 0;
    head = 0 as libc::c_int as libc::c_longlong;
    *Stack.offset(0 as libc::c_int as isize) = root;
    while head >= 0 as libc::c_int as libc::c_longlong {
        i = *Stack.offset(head as isize);
        if *Child.offset(i as isize) != EMPTY as libc::c_longlong {
            f = *Child.offset(i as isize);
            while f != EMPTY as libc::c_longlong {
                head += 1;
                f = *Sibling.offset(f as isize);
            }
            h = head;
            f = *Child.offset(i as isize);
            while f != EMPTY as libc::c_longlong {
                let fresh0 = h;
                h = h - 1;
                *Stack.offset(fresh0 as isize) = f;
                f = *Sibling.offset(f as isize);
            }
            *Child.offset(i as isize) = EMPTY as libc::c_longlong;
        } else {
            head -= 1;
            let fresh1 = k;
            k = k + 1;
            *Order.offset(i as isize) = fresh1;
        }
    }
    return k;
}
