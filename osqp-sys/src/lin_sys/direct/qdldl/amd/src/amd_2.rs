use ::libc;
extern "C" {
    fn sqrt(_: libc::c_double) -> libc::c_double;
    fn amd_l_postorder(
        nn: libc::c_longlong,
        Parent: *mut libc::c_longlong,
        Npiv: *mut libc::c_longlong,
        Fsize: *mut libc::c_longlong,
        Order: *mut libc::c_longlong,
        Child: *mut libc::c_longlong,
        Sibling: *mut libc::c_longlong,
        Stack: *mut libc::c_longlong,
    );
}
pub type c_float = libc::c_double;
pub const EMPTY: libc::c_int = -(1 as libc::c_int);
pub const AMD_postorder: unsafe extern "C" fn(
    libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
    *mut libc::c_longlong,
) -> () = amd_l_postorder;
pub const AMD_OK: libc::c_int = 0 as libc::c_int;
pub const AMD_STATUS: libc::c_int = 0 as libc::c_int;
pub const AMD_NCMPA: libc::c_int = 8 as libc::c_int;
pub const AMD_DMAX: libc::c_int = 13 as libc::c_int;
pub const AMD_NDENSE: libc::c_int = 6 as libc::c_int;
pub const AMD_NMULTSUBS_LU: libc::c_int = 12 as libc::c_int;
pub const AMD_NMULTSUBS_LDL: libc::c_int = 11 as libc::c_int;
pub const AMD_NDIV: libc::c_int = 10 as libc::c_int;
pub const AMD_LNZ: libc::c_int = 9 as libc::c_int;
pub const NULL: libc::c_int = 0 as libc::c_int;
pub const LONG_MAX: libc::c_long = 0x7fffffffffffffff as libc::c_long;
pub const SuiteSparse_long_max: libc::c_long = 0x7fffffffffffffff as libc::c_long;
pub const Int_MAX: libc::c_long = 0x7fffffffffffffff as libc::c_long;
pub const AMD_DEFAULT_AGGRESSIVE: libc::c_int = 1 as libc::c_int;
pub const AMD_DEFAULT_DENSE: libc::c_double = 10.0f64;
pub const AMD_AGGRESSIVE: libc::c_int = 1 as libc::c_int;
pub const AMD_DENSE: libc::c_int = 0 as libc::c_int;
unsafe extern "C" fn clear_flag(
    mut wflg: libc::c_longlong,
    mut wbig: libc::c_longlong,
    mut W: *mut libc::c_longlong,
    mut n: libc::c_longlong,
) -> libc::c_longlong {
    let mut x: libc::c_longlong = 0;
    if wflg < 2 as libc::c_int as libc::c_longlong || wflg >= wbig {
        x = 0 as libc::c_int as libc::c_longlong;
        while x < n {
            if *W.offset(x as isize) != 0 as libc::c_int as libc::c_longlong {
                *W.offset(x as isize) = 1 as libc::c_int as libc::c_longlong;
            }
            x += 1;
        }
        wflg = 2 as libc::c_int as libc::c_longlong;
    }
    return wflg;
}
#[no_mangle]
pub unsafe extern "C" fn amd_l2(
    mut n: libc::c_longlong,
    mut Pe: *mut libc::c_longlong,
    mut Iw: *mut libc::c_longlong,
    mut Len: *mut libc::c_longlong,
    mut iwlen: libc::c_longlong,
    mut pfree: libc::c_longlong,
    mut Nv: *mut libc::c_longlong,
    mut Next: *mut libc::c_longlong,
    mut Last: *mut libc::c_longlong,
    mut Head: *mut libc::c_longlong,
    mut Elen: *mut libc::c_longlong,
    mut Degree: *mut libc::c_longlong,
    mut W: *mut libc::c_longlong,
    mut Control: *mut c_float,
    mut Info: *mut c_float,
) {
    let mut deg: libc::c_longlong = 0;
    let mut degme: libc::c_longlong = 0;
    let mut dext: libc::c_longlong = 0;
    let mut lemax: libc::c_longlong = 0;
    let mut e: libc::c_longlong = 0;
    let mut elenme: libc::c_longlong = 0;
    let mut eln: libc::c_longlong = 0;
    let mut i: libc::c_longlong = 0;
    let mut ilast: libc::c_longlong = 0;
    let mut inext: libc::c_longlong = 0;
    let mut j: libc::c_longlong = 0;
    let mut jlast: libc::c_longlong = 0;
    let mut jnext: libc::c_longlong = 0;
    let mut k: libc::c_longlong = 0;
    let mut knt1: libc::c_longlong = 0;
    let mut knt2: libc::c_longlong = 0;
    let mut knt3: libc::c_longlong = 0;
    let mut lenj: libc::c_longlong = 0;
    let mut ln: libc::c_longlong = 0;
    let mut me: libc::c_longlong = 0;
    let mut mindeg: libc::c_longlong = 0;
    let mut nel: libc::c_longlong = 0;
    let mut nleft: libc::c_longlong = 0;
    let mut nvi: libc::c_longlong = 0;
    let mut nvj: libc::c_longlong = 0;
    let mut nvpiv: libc::c_longlong = 0;
    let mut slenme: libc::c_longlong = 0;
    let mut wbig: libc::c_longlong = 0;
    let mut we: libc::c_longlong = 0;
    let mut wflg: libc::c_longlong = 0;
    let mut wnvi: libc::c_longlong = 0;
    let mut ok: libc::c_longlong = 0;
    let mut ndense: libc::c_longlong = 0;
    let mut ncmpa: libc::c_longlong = 0;
    let mut dense: libc::c_longlong = 0;
    let mut aggressive: libc::c_longlong = 0;
    let mut hash: libc::c_ulonglong = 0;
    let mut f: c_float = 0.;
    let mut r: c_float = 0.;
    let mut ndiv: c_float = 0.;
    let mut s: c_float = 0.;
    let mut nms_lu: c_float = 0.;
    let mut nms_ldl: c_float = 0.;
    let mut dmax: c_float = 0.;
    let mut alpha: c_float = 0.;
    let mut lnz: c_float = 0.;
    let mut lnzme: c_float = 0.;
    let mut p: libc::c_longlong = 0;
    let mut p1: libc::c_longlong = 0;
    let mut p2: libc::c_longlong = 0;
    let mut p3: libc::c_longlong = 0;
    let mut p4: libc::c_longlong = 0;
    let mut pdst: libc::c_longlong = 0;
    let mut pend: libc::c_longlong = 0;
    let mut pj: libc::c_longlong = 0;
    let mut pme: libc::c_longlong = 0;
    let mut pme1: libc::c_longlong = 0;
    let mut pme2: libc::c_longlong = 0;
    let mut pn: libc::c_longlong = 0;
    let mut psrc: libc::c_longlong = 0;
    lnz = 0 as libc::c_int as c_float;
    ndiv = 0 as libc::c_int as c_float;
    nms_lu = 0 as libc::c_int as c_float;
    nms_ldl = 0 as libc::c_int as c_float;
    dmax = 1 as libc::c_int as c_float;
    me = EMPTY as libc::c_longlong;
    mindeg = 0 as libc::c_int as libc::c_longlong;
    ncmpa = 0 as libc::c_int as libc::c_longlong;
    nel = 0 as libc::c_int as libc::c_longlong;
    lemax = 0 as libc::c_int as libc::c_longlong;
    if !Control.is_null() {
        alpha = *Control.offset(AMD_DENSE as isize);
        aggressive = (*Control.offset(AMD_AGGRESSIVE as isize)
            != 0 as libc::c_int as libc::c_double) as libc::c_int as libc::c_longlong;
    } else {
        alpha = AMD_DEFAULT_DENSE;
        aggressive = AMD_DEFAULT_AGGRESSIVE as libc::c_longlong;
    }
    if alpha < 0 as libc::c_int as libc::c_double {
        dense = n - 2 as libc::c_int as libc::c_longlong;
    } else {
        dense = (alpha * sqrt(n as c_float)) as libc::c_longlong;
    }
    dense = if 16 as libc::c_int as libc::c_longlong > dense {
        16 as libc::c_int as libc::c_longlong
    } else {
        dense
    };
    dense = if n < dense { n } else { dense };
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        *Last.offset(i as isize) = EMPTY as libc::c_longlong;
        *Head.offset(i as isize) = EMPTY as libc::c_longlong;
        *Next.offset(i as isize) = EMPTY as libc::c_longlong;
        *Nv.offset(i as isize) = 1 as libc::c_int as libc::c_longlong;
        *W.offset(i as isize) = 1 as libc::c_int as libc::c_longlong;
        *Elen.offset(i as isize) = 0 as libc::c_int as libc::c_longlong;
        *Degree.offset(i as isize) = *Len.offset(i as isize);
        i += 1;
    }
    wbig = Int_MAX as libc::c_longlong - n;
    wflg = clear_flag(0 as libc::c_int as libc::c_longlong, wbig, W, n);
    ndense = 0 as libc::c_int as libc::c_longlong;
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        deg = *Degree.offset(i as isize);
        if deg == 0 as libc::c_int as libc::c_longlong {
            *Elen
                .offset(
                    i as isize,
                ) = (-(1 as libc::c_int) - 2 as libc::c_int) as libc::c_longlong;
            nel += 1;
            *Pe.offset(i as isize) = EMPTY as libc::c_longlong;
            *W.offset(i as isize) = 0 as libc::c_int as libc::c_longlong;
        } else if deg > dense {
            ndense += 1;
            *Nv.offset(i as isize) = 0 as libc::c_int as libc::c_longlong;
            *Elen.offset(i as isize) = EMPTY as libc::c_longlong;
            nel += 1;
            *Pe.offset(i as isize) = EMPTY as libc::c_longlong;
        } else {
            inext = *Head.offset(deg as isize);
            if inext != EMPTY as libc::c_longlong {
                *Last.offset(inext as isize) = i;
            }
            *Next.offset(i as isize) = inext;
            *Head.offset(deg as isize) = i;
        }
        i += 1;
    }
    while nel < n {
        deg = mindeg;
        while deg < n {
            me = *Head.offset(deg as isize);
            if me != EMPTY as libc::c_longlong {
                break;
            }
            deg += 1;
        }
        mindeg = deg;
        inext = *Next.offset(me as isize);
        if inext != EMPTY as libc::c_longlong {
            *Last.offset(inext as isize) = EMPTY as libc::c_longlong;
        }
        *Head.offset(deg as isize) = inext;
        elenme = *Elen.offset(me as isize);
        nvpiv = *Nv.offset(me as isize);
        nel += nvpiv;
        *Nv.offset(me as isize) = -nvpiv;
        degme = 0 as libc::c_int as libc::c_longlong;
        if elenme == 0 as libc::c_int as libc::c_longlong {
            pme1 = *Pe.offset(me as isize);
            pme2 = pme1 - 1 as libc::c_int as libc::c_longlong;
            p = pme1;
            while p
                <= pme1 + *Len.offset(me as isize) - 1 as libc::c_int as libc::c_longlong
            {
                i = *Iw.offset(p as isize);
                nvi = *Nv.offset(i as isize);
                if nvi > 0 as libc::c_int as libc::c_longlong {
                    degme += nvi;
                    *Nv.offset(i as isize) = -nvi;
                    pme2 += 1;
                    *Iw.offset(pme2 as isize) = i;
                    ilast = *Last.offset(i as isize);
                    inext = *Next.offset(i as isize);
                    if inext != EMPTY as libc::c_longlong {
                        *Last.offset(inext as isize) = ilast;
                    }
                    if ilast != EMPTY as libc::c_longlong {
                        *Next.offset(ilast as isize) = inext;
                    } else {
                        *Head.offset(*Degree.offset(i as isize) as isize) = inext;
                    }
                }
                p += 1;
            }
        } else {
            p = *Pe.offset(me as isize);
            pme1 = pfree;
            slenme = *Len.offset(me as isize) - elenme;
            knt1 = 1 as libc::c_int as libc::c_longlong;
            while knt1 <= elenme + 1 as libc::c_int as libc::c_longlong {
                if knt1 > elenme {
                    e = me;
                    pj = p;
                    ln = slenme;
                } else {
                    let fresh0 = p;
                    p = p + 1;
                    e = *Iw.offset(fresh0 as isize);
                    pj = *Pe.offset(e as isize);
                    ln = *Len.offset(e as isize);
                }
                knt2 = 1 as libc::c_int as libc::c_longlong;
                while knt2 <= ln {
                    let fresh1 = pj;
                    pj = pj + 1;
                    i = *Iw.offset(fresh1 as isize);
                    nvi = *Nv.offset(i as isize);
                    if nvi > 0 as libc::c_int as libc::c_longlong {
                        if pfree >= iwlen {
                            *Pe.offset(me as isize) = p;
                            *Len.offset(me as isize) -= knt1;
                            if *Len.offset(me as isize)
                                == 0 as libc::c_int as libc::c_longlong
                            {
                                *Pe.offset(me as isize) = EMPTY as libc::c_longlong;
                            }
                            *Pe.offset(e as isize) = pj;
                            *Len.offset(e as isize) = ln - knt2;
                            if *Len.offset(e as isize)
                                == 0 as libc::c_int as libc::c_longlong
                            {
                                *Pe.offset(e as isize) = EMPTY as libc::c_longlong;
                            }
                            ncmpa += 1;
                            j = 0 as libc::c_int as libc::c_longlong;
                            while j < n {
                                pn = *Pe.offset(j as isize);
                                if pn >= 0 as libc::c_int as libc::c_longlong {
                                    *Pe.offset(j as isize) = *Iw.offset(pn as isize);
                                    *Iw
                                        .offset(
                                            pn as isize,
                                        ) = -j - 2 as libc::c_int as libc::c_longlong;
                                }
                                j += 1;
                            }
                            psrc = 0 as libc::c_int as libc::c_longlong;
                            pdst = 0 as libc::c_int as libc::c_longlong;
                            pend = pme1 - 1 as libc::c_int as libc::c_longlong;
                            while psrc <= pend {
                                let fresh2 = psrc;
                                psrc = psrc + 1;
                                j = -*Iw.offset(fresh2 as isize)
                                    - 2 as libc::c_int as libc::c_longlong;
                                if j >= 0 as libc::c_int as libc::c_longlong {
                                    *Iw.offset(pdst as isize) = *Pe.offset(j as isize);
                                    let fresh3 = pdst;
                                    pdst = pdst + 1;
                                    *Pe.offset(j as isize) = fresh3;
                                    lenj = *Len.offset(j as isize);
                                    knt3 = 0 as libc::c_int as libc::c_longlong;
                                    while knt3 <= lenj - 2 as libc::c_int as libc::c_longlong {
                                        let fresh4 = psrc;
                                        psrc = psrc + 1;
                                        let fresh5 = pdst;
                                        pdst = pdst + 1;
                                        *Iw.offset(fresh5 as isize) = *Iw.offset(fresh4 as isize);
                                        knt3 += 1;
                                    }
                                }
                            }
                            p1 = pdst;
                            psrc = pme1;
                            while psrc <= pfree - 1 as libc::c_int as libc::c_longlong {
                                let fresh6 = pdst;
                                pdst = pdst + 1;
                                *Iw.offset(fresh6 as isize) = *Iw.offset(psrc as isize);
                                psrc += 1;
                            }
                            pme1 = p1;
                            pfree = pdst;
                            pj = *Pe.offset(e as isize);
                            p = *Pe.offset(me as isize);
                        }
                        degme += nvi;
                        *Nv.offset(i as isize) = -nvi;
                        let fresh7 = pfree;
                        pfree = pfree + 1;
                        *Iw.offset(fresh7 as isize) = i;
                        ilast = *Last.offset(i as isize);
                        inext = *Next.offset(i as isize);
                        if inext != EMPTY as libc::c_longlong {
                            *Last.offset(inext as isize) = ilast;
                        }
                        if ilast != EMPTY as libc::c_longlong {
                            *Next.offset(ilast as isize) = inext;
                        } else {
                            *Head.offset(*Degree.offset(i as isize) as isize) = inext;
                        }
                    }
                    knt2 += 1;
                }
                if e != me {
                    *Pe.offset(e as isize) = -me - 2 as libc::c_int as libc::c_longlong;
                    *W.offset(e as isize) = 0 as libc::c_int as libc::c_longlong;
                }
                knt1 += 1;
            }
            pme2 = pfree - 1 as libc::c_int as libc::c_longlong;
        }
        *Degree.offset(me as isize) = degme;
        *Pe.offset(me as isize) = pme1;
        *Len.offset(me as isize) = pme2 - pme1 + 1 as libc::c_int as libc::c_longlong;
        *Elen
            .offset(
                me as isize,
            ) = -(nvpiv + degme) - 2 as libc::c_int as libc::c_longlong;
        wflg = clear_flag(wflg, wbig, W, n);
        pme = pme1;
        while pme <= pme2 {
            i = *Iw.offset(pme as isize);
            eln = *Elen.offset(i as isize);
            if eln > 0 as libc::c_int as libc::c_longlong {
                nvi = -*Nv.offset(i as isize);
                wnvi = wflg - nvi;
                p = *Pe.offset(i as isize);
                while p
                    <= *Pe.offset(i as isize) + eln
                        - 1 as libc::c_int as libc::c_longlong
                {
                    e = *Iw.offset(p as isize);
                    we = *W.offset(e as isize);
                    if we >= wflg {
                        we -= nvi;
                    } else if we != 0 as libc::c_int as libc::c_longlong {
                        we = *Degree.offset(e as isize) + wnvi;
                    }
                    *W.offset(e as isize) = we;
                    p += 1;
                }
            }
            pme += 1;
        }
        pme = pme1;
        while pme <= pme2 {
            i = *Iw.offset(pme as isize);
            p1 = *Pe.offset(i as isize);
            p2 = p1 + *Elen.offset(i as isize) - 1 as libc::c_int as libc::c_longlong;
            pn = p1;
            hash = 0 as libc::c_int as libc::c_ulonglong;
            deg = 0 as libc::c_int as libc::c_longlong;
            if aggressive != 0 {
                p = p1;
                while p <= p2 {
                    e = *Iw.offset(p as isize);
                    we = *W.offset(e as isize);
                    if we != 0 as libc::c_int as libc::c_longlong {
                        dext = we - wflg;
                        if dext > 0 as libc::c_int as libc::c_longlong {
                            deg += dext;
                            let fresh8 = pn;
                            pn = pn + 1;
                            *Iw.offset(fresh8 as isize) = e;
                            hash = hash.wrapping_add(e as libc::c_ulonglong);
                        } else {
                            *Pe
                                .offset(
                                    e as isize,
                                ) = -me - 2 as libc::c_int as libc::c_longlong;
                            *W.offset(e as isize) = 0 as libc::c_int as libc::c_longlong;
                        }
                    }
                    p += 1;
                }
            } else {
                p = p1;
                while p <= p2 {
                    e = *Iw.offset(p as isize);
                    we = *W.offset(e as isize);
                    if we != 0 as libc::c_int as libc::c_longlong {
                        dext = we - wflg;
                        deg += dext;
                        let fresh9 = pn;
                        pn = pn + 1;
                        *Iw.offset(fresh9 as isize) = e;
                        hash = hash.wrapping_add(e as libc::c_ulonglong);
                    }
                    p += 1;
                }
            }
            *Elen.offset(i as isize) = pn - p1 + 1 as libc::c_int as libc::c_longlong;
            p3 = pn;
            p4 = p1 + *Len.offset(i as isize);
            p = p2 + 1 as libc::c_int as libc::c_longlong;
            while p < p4 {
                j = *Iw.offset(p as isize);
                nvj = *Nv.offset(j as isize);
                if nvj > 0 as libc::c_int as libc::c_longlong {
                    deg += nvj;
                    let fresh10 = pn;
                    pn = pn + 1;
                    *Iw.offset(fresh10 as isize) = j;
                    hash = hash.wrapping_add(j as libc::c_ulonglong);
                }
                p += 1;
            }
            if *Elen.offset(i as isize) == 1 as libc::c_int as libc::c_longlong
                && p3 == pn
            {
                *Pe.offset(i as isize) = -me - 2 as libc::c_int as libc::c_longlong;
                nvi = -*Nv.offset(i as isize);
                degme -= nvi;
                nvpiv += nvi;
                nel += nvi;
                *Nv.offset(i as isize) = 0 as libc::c_int as libc::c_longlong;
                *Elen.offset(i as isize) = EMPTY as libc::c_longlong;
            } else {
                *Degree
                    .offset(
                        i as isize,
                    ) = if *Degree.offset(i as isize) < deg {
                    *Degree.offset(i as isize)
                } else {
                    deg
                };
                *Iw.offset(pn as isize) = *Iw.offset(p3 as isize);
                *Iw.offset(p3 as isize) = *Iw.offset(p1 as isize);
                *Iw.offset(p1 as isize) = me;
                *Len.offset(i as isize) = pn - p1 + 1 as libc::c_int as libc::c_longlong;
                hash = hash.wrapping_rem(n as libc::c_ulonglong);
                j = *Head.offset(hash as isize);
                if j <= EMPTY as libc::c_longlong {
                    *Next.offset(i as isize) = -j - 2 as libc::c_int as libc::c_longlong;
                    *Head
                        .offset(
                            hash as isize,
                        ) = -i - 2 as libc::c_int as libc::c_longlong;
                } else {
                    *Next.offset(i as isize) = *Last.offset(j as isize);
                    *Last.offset(j as isize) = i;
                }
                *Last.offset(i as isize) = hash as libc::c_longlong;
            }
            pme += 1;
        }
        *Degree.offset(me as isize) = degme;
        lemax = if lemax > degme { lemax } else { degme };
        wflg += lemax;
        wflg = clear_flag(wflg, wbig, W, n);
        pme = pme1;
        while pme <= pme2 {
            i = *Iw.offset(pme as isize);
            if *Nv.offset(i as isize) < 0 as libc::c_int as libc::c_longlong {
                hash = *Last.offset(i as isize) as libc::c_ulonglong;
                j = *Head.offset(hash as isize);
                if j == EMPTY as libc::c_longlong {
                    i = EMPTY as libc::c_longlong;
                } else if j < EMPTY as libc::c_longlong {
                    i = -j - 2 as libc::c_int as libc::c_longlong;
                    *Head.offset(hash as isize) = EMPTY as libc::c_longlong;
                } else {
                    i = *Last.offset(j as isize);
                    *Last.offset(j as isize) = EMPTY as libc::c_longlong;
                }
                while i != EMPTY as libc::c_longlong
                    && *Next.offset(i as isize) != EMPTY as libc::c_longlong
                {
                    ln = *Len.offset(i as isize);
                    eln = *Elen.offset(i as isize);
                    p = *Pe.offset(i as isize) + 1 as libc::c_int as libc::c_longlong;
                    while p
                        <= *Pe.offset(i as isize) + ln
                            - 1 as libc::c_int as libc::c_longlong
                    {
                        *W.offset(*Iw.offset(p as isize) as isize) = wflg;
                        p += 1;
                    }
                    jlast = i;
                    j = *Next.offset(i as isize);
                    while j != EMPTY as libc::c_longlong {
                        ok = (*Len.offset(j as isize) == ln
                            && *Elen.offset(j as isize) == eln) as libc::c_int
                            as libc::c_longlong;
                        p = *Pe.offset(j as isize)
                            + 1 as libc::c_int as libc::c_longlong;
                        while ok != 0
                            && p
                                <= *Pe.offset(j as isize) + ln
                                    - 1 as libc::c_int as libc::c_longlong
                        {
                            if *W.offset(*Iw.offset(p as isize) as isize) != wflg {
                                ok = 0 as libc::c_int as libc::c_longlong;
                            }
                            p += 1;
                        }
                        if ok != 0 {
                            *Pe
                                .offset(
                                    j as isize,
                                ) = -i - 2 as libc::c_int as libc::c_longlong;
                            *Nv.offset(i as isize) += *Nv.offset(j as isize);
                            *Nv
                                .offset(j as isize) = 0 as libc::c_int as libc::c_longlong;
                            *Elen.offset(j as isize) = EMPTY as libc::c_longlong;
                            j = *Next.offset(j as isize);
                            *Next.offset(jlast as isize) = j;
                        } else {
                            jlast = j;
                            j = *Next.offset(j as isize);
                        }
                    }
                    wflg += 1;
                    i = *Next.offset(i as isize);
                }
            }
            pme += 1;
        }
        p = pme1;
        nleft = n - nel;
        pme = pme1;
        while pme <= pme2 {
            i = *Iw.offset(pme as isize);
            nvi = -*Nv.offset(i as isize);
            if nvi > 0 as libc::c_int as libc::c_longlong {
                *Nv.offset(i as isize) = nvi;
                deg = *Degree.offset(i as isize) + degme - nvi;
                deg = if deg < nleft - nvi { deg } else { nleft - nvi };
                inext = *Head.offset(deg as isize);
                if inext != EMPTY as libc::c_longlong {
                    *Last.offset(inext as isize) = i;
                }
                *Next.offset(i as isize) = inext;
                *Last.offset(i as isize) = EMPTY as libc::c_longlong;
                *Head.offset(deg as isize) = i;
                mindeg = if mindeg < deg { mindeg } else { deg };
                *Degree.offset(i as isize) = deg;
                let fresh11 = p;
                p = p + 1;
                *Iw.offset(fresh11 as isize) = i;
            }
            pme += 1;
        }
        *Nv.offset(me as isize) = nvpiv;
        *Len.offset(me as isize) = p - pme1;
        if *Len.offset(me as isize) == 0 as libc::c_int as libc::c_longlong {
            *Pe.offset(me as isize) = EMPTY as libc::c_longlong;
            *W.offset(me as isize) = 0 as libc::c_int as libc::c_longlong;
        }
        if elenme != 0 as libc::c_int as libc::c_longlong {
            pfree = p;
        }
        if !Info.is_null() {
            f = nvpiv as c_float;
            r = (degme + ndense) as c_float;
            dmax = if dmax > f + r { dmax } else { f + r };
            lnzme = f * r
                + (f - 1 as libc::c_int as libc::c_double) * f
                    / 2 as libc::c_int as libc::c_double;
            lnz += lnzme;
            ndiv += lnzme;
            s = f * r * r + r * (f - 1 as libc::c_int as libc::c_double) * f
                + (f - 1 as libc::c_int as libc::c_double) * f
                    * (2 as libc::c_int as libc::c_double * f
                        - 1 as libc::c_int as libc::c_double)
                    / 6 as libc::c_int as libc::c_double;
            nms_lu += s;
            nms_ldl += (s + lnzme) / 2 as libc::c_int as libc::c_double;
        }
    }
    if !Info.is_null() {
        f = ndense as c_float;
        dmax = if dmax > ndense as c_float { dmax } else { ndense as c_float };
        lnzme = (f - 1 as libc::c_int as libc::c_double) * f
            / 2 as libc::c_int as libc::c_double;
        lnz += lnzme;
        ndiv += lnzme;
        s = (f - 1 as libc::c_int as libc::c_double) * f
            * (2 as libc::c_int as libc::c_double * f
                - 1 as libc::c_int as libc::c_double)
            / 6 as libc::c_int as libc::c_double;
        nms_lu += s;
        nms_ldl += (s + lnzme) / 2 as libc::c_int as libc::c_double;
        *Info.offset(AMD_LNZ as isize) = lnz;
        *Info.offset(AMD_NDIV as isize) = ndiv;
        *Info.offset(AMD_NMULTSUBS_LDL as isize) = nms_ldl;
        *Info.offset(AMD_NMULTSUBS_LU as isize) = nms_lu;
        *Info.offset(AMD_NDENSE as isize) = ndense as c_float;
        *Info.offset(AMD_DMAX as isize) = dmax;
        *Info.offset(AMD_NCMPA as isize) = ncmpa as c_float;
        *Info.offset(AMD_STATUS as isize) = AMD_OK as c_float;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        *Pe
            .offset(
                i as isize,
            ) = -*Pe.offset(i as isize) - 2 as libc::c_int as libc::c_longlong;
        i += 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        *Elen
            .offset(
                i as isize,
            ) = -*Elen.offset(i as isize) - 2 as libc::c_int as libc::c_longlong;
        i += 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        if *Nv.offset(i as isize) == 0 as libc::c_int as libc::c_longlong {
            j = *Pe.offset(i as isize);
            if !(j == EMPTY as libc::c_longlong) {
                while *Nv.offset(j as isize) == 0 as libc::c_int as libc::c_longlong {
                    j = *Pe.offset(j as isize);
                }
                e = j;
                j = i;
                while *Nv.offset(j as isize) == 0 as libc::c_int as libc::c_longlong {
                    jnext = *Pe.offset(j as isize);
                    *Pe.offset(j as isize) = e;
                    j = jnext;
                }
            }
        }
        i += 1;
    }
    amd_l_postorder(n, Pe, Nv, Elen, W, Head, Next, Last);
    k = 0 as libc::c_int as libc::c_longlong;
    while k < n {
        *Head.offset(k as isize) = EMPTY as libc::c_longlong;
        *Next.offset(k as isize) = EMPTY as libc::c_longlong;
        k += 1;
    }
    e = 0 as libc::c_int as libc::c_longlong;
    while e < n {
        k = *W.offset(e as isize);
        if k != EMPTY as libc::c_longlong {
            *Head.offset(k as isize) = e;
        }
        e += 1;
    }
    nel = 0 as libc::c_int as libc::c_longlong;
    k = 0 as libc::c_int as libc::c_longlong;
    while k < n {
        e = *Head.offset(k as isize);
        if e == EMPTY as libc::c_longlong {
            break;
        }
        *Next.offset(e as isize) = nel;
        nel += *Nv.offset(e as isize);
        k += 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        if *Nv.offset(i as isize) == 0 as libc::c_int as libc::c_longlong {
            e = *Pe.offset(i as isize);
            if e != EMPTY as libc::c_longlong {
                *Next.offset(i as isize) = *Next.offset(e as isize);
                let ref mut fresh12 = *Next.offset(e as isize);
                *fresh12 += 1;
            } else {
                let fresh13 = nel;
                nel = nel + 1;
                *Next.offset(i as isize) = fresh13;
            }
        }
        i += 1;
    }
    i = 0 as libc::c_int as libc::c_longlong;
    while i < n {
        k = *Next.offset(i as isize);
        *Last.offset(k as isize) = i;
        i += 1;
    }
}
