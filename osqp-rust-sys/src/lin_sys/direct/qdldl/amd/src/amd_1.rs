use ::libc;
extern "C" {
    fn amd_l2(
        n: libc::c_longlong,
        Pe: *mut libc::c_longlong,
        Iw: *mut libc::c_longlong,
        Len: *mut libc::c_longlong,
        iwlen: libc::c_longlong,
        pfree: libc::c_longlong,
        Nv: *mut libc::c_longlong,
        Next: *mut libc::c_longlong,
        Last: *mut libc::c_longlong,
        Head: *mut libc::c_longlong,
        Elen: *mut libc::c_longlong,
        Degree: *mut libc::c_longlong,
        W: *mut libc::c_longlong,
        Control: *mut c_float,
        Info: *mut c_float,
    );
}
pub type c_float = libc::c_double;
pub const AMD_2: unsafe extern "C" fn(
    libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    libc::c_longlong,
    libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut c_float,
    *mut c_float,
) -> () = amd_l2;
#[no_mangle]
pub unsafe extern "C" fn amd_l1(
    mut n: libc::c_longlong,
    mut Ap: *const libc::c_longlong,
    mut Ai: *const libc::c_longlong,
    mut P: *mut libc::c_longlong,
    mut Pinv: *mut libc::c_longlong,
    mut Len: *mut libc::c_longlong,
    mut slen: libc::c_longlong,
    mut S: *mut libc::c_longlong,
    mut Control: *mut c_float,
    mut Info: *mut c_float,
) {
    let mut i: libc::c_longlong = 0;
    let mut j: libc::c_longlong = 0;
    let mut k: libc::c_longlong = 0;
    let mut p: libc::c_longlong = 0;
    let mut pfree: libc::c_longlong = 0;
    let mut iwlen: libc::c_longlong = 0;
    let mut pj: libc::c_longlong = 0;
    let mut p1: libc::c_longlong = 0;
    let mut p2: libc::c_longlong = 0;
    let mut pj2: libc::c_longlong = 0;
    let mut Iw: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Pe: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Nv: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Head: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Elen: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Degree: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut s: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut W: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Sp: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    let mut Tp: *mut libc::c_longlong = 0 as *mut libc::c_longlong;
    iwlen = slen - 6 as libc::c_int as libc::c_longlong * n;
    s = S;
    Pe = s;
    s = s.offset(n as isize);
    Nv = s;
    s = s.offset(n as isize);
    Head = s;
    s = s.offset(n as isize);
    Elen = s;
    s = s.offset(n as isize);
    Degree = s;
    s = s.offset(n as isize);
    W = s;
    s = s.offset(n as isize);
    Iw = s;
    s = s.offset(iwlen as isize);
    Sp = Nv;
    Tp = W;
    pfree = 0 as libc::c_int as libc::c_longlong;
    j = 0 as libc::c_int as libc::c_longlong;
    while j < n {
        *Pe.offset(j as isize) = pfree;
        *Sp.offset(j as isize) = pfree;
        pfree += *Len.offset(j as isize);
        j += 1;
    }
    k = 0 as libc::c_int as libc::c_longlong;
    while k < n {
        p1 = *Ap.offset(k as isize);
        p2 = *Ap.offset((k + 1 as libc::c_int as libc::c_longlong) as isize);
        p = p1;
        while p < p2 {
            j = *Ai.offset(p as isize);
            if j < k {
                let ref mut fresh0 = *Sp.offset(j as isize);
                let fresh1 = *fresh0;
                *fresh0 = *fresh0 + 1;
                *Iw.offset(fresh1 as isize) = k;
                let ref mut fresh2 = *Sp.offset(k as isize);
                let fresh3 = *fresh2;
                *fresh2 = *fresh2 + 1;
                *Iw.offset(fresh3 as isize) = j;
                p += 1;
                pj2 = *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize);
                pj = *Tp.offset(j as isize);
                while pj < pj2 {
                    i = *Ai.offset(pj as isize);
                    if i < k {
                        let ref mut fresh4 = *Sp.offset(i as isize);
                        let fresh5 = *fresh4;
                        *fresh4 = *fresh4 + 1;
                        *Iw.offset(fresh5 as isize) = j;
                        let ref mut fresh6 = *Sp.offset(j as isize);
                        let fresh7 = *fresh6;
                        *fresh6 = *fresh6 + 1;
                        *Iw.offset(fresh7 as isize) = i;
                        pj += 1;
                    } else if i == k {
                        pj += 1;
                        break;
                    } else {
                        break;
                    }
                }
                *Tp.offset(j as isize) = pj;
            } else if j == k {
                p += 1;
                break;
            } else {
                break;
            }
        }
        *Tp.offset(k as isize) = p;
        k += 1;
    }
    j = 0 as libc::c_int as libc::c_longlong;
    while j < n {
        pj = *Tp.offset(j as isize);
        while pj < *Ap.offset((j + 1 as libc::c_int as libc::c_longlong) as isize) {
            i = *Ai.offset(pj as isize);
            let ref mut fresh8 = *Sp.offset(i as isize);
            let fresh9 = *fresh8;
            *fresh8 = *fresh8 + 1;
            *Iw.offset(fresh9 as isize) = j;
            let ref mut fresh10 = *Sp.offset(j as isize);
            let fresh11 = *fresh10;
            *fresh10 = *fresh10 + 1;
            *Iw.offset(fresh11 as isize) = i;
            pj += 1;
        }
        j += 1;
    }
    amd_l2(
        n,
        Pe,
        Iw,
        Len,
        iwlen,
        pfree,
        Nv,
        Pinv,
        P,
        Head,
        Elen,
        Degree,
        W,
        Control,
        Info,
    );
}
