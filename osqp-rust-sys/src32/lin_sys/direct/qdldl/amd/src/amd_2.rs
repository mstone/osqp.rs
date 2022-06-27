extern "C" {
    fn sqrt(_: ::std::os::raw::c_double) -> ::std::os::raw::c_double;
    fn amd_postorder(
        nn: ::std::os::raw::c_int,
        Parent: *mut ::std::os::raw::c_int,
        Npiv: *mut ::std::os::raw::c_int,
        Fsize: *mut ::std::os::raw::c_int,
        Order: *mut ::std::os::raw::c_int,
        Child: *mut ::std::os::raw::c_int,
        Sibling: *mut ::std::os::raw::c_int,
        Stack: *mut ::std::os::raw::c_int,
    );
}
pub type c_float = ::std::os::raw::c_double;
pub const AMD_DENSE: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_AGGRESSIVE: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const AMD_DEFAULT_DENSE: ::std::os::raw::c_double = 10.0f64;
pub const AMD_DEFAULT_AGGRESSIVE: ::std::os::raw::c_int = 1 as ::std::os::raw::c_int;
pub const Int_MAX: ::std::os::raw::c_int = INT_MAX;
pub const INT_MAX: ::std::os::raw::c_int = 2147483647 as ::std::os::raw::c_int;
pub const NULL: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_LNZ: ::std::os::raw::c_int = 9 as ::std::os::raw::c_int;
pub const AMD_NDIV: ::std::os::raw::c_int = 10 as ::std::os::raw::c_int;
pub const AMD_NMULTSUBS_LDL: ::std::os::raw::c_int = 11 as ::std::os::raw::c_int;
pub const AMD_NMULTSUBS_LU: ::std::os::raw::c_int = 12 as ::std::os::raw::c_int;
pub const AMD_NDENSE: ::std::os::raw::c_int = 6 as ::std::os::raw::c_int;
pub const AMD_DMAX: ::std::os::raw::c_int = 13 as ::std::os::raw::c_int;
pub const AMD_NCMPA: ::std::os::raw::c_int = 8 as ::std::os::raw::c_int;
pub const AMD_STATUS: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const AMD_OK: ::std::os::raw::c_int = 0 as ::std::os::raw::c_int;
pub const EMPTY: ::std::os::raw::c_int = -(1 as ::std::os::raw::c_int);
pub const AMD_postorder: unsafe extern "C" fn(
    ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
    *mut ::std::os::raw::c_int,
) -> () = amd_postorder;
unsafe extern "C" fn clear_flag(
    mut wflg: ::std::os::raw::c_int,
    mut wbig: ::std::os::raw::c_int,
    mut W: *mut ::std::os::raw::c_int,
    mut n: ::std::os::raw::c_int,
) -> ::std::os::raw::c_int {
    let mut x: ::std::os::raw::c_int = 0;
    if wflg < 2 as ::std::os::raw::c_int || wflg >= wbig {
        x = 0 as ::std::os::raw::c_int;
        while x < n {
            if *W.offset(x as isize) != 0 as ::std::os::raw::c_int {
                *W.offset(x as isize) = 1 as ::std::os::raw::c_int;
            }
            x += 1;
        }
        wflg = 2 as ::std::os::raw::c_int;
    }
    return wflg;
}
#[no_mangle]
pub unsafe extern "C" fn amd_2(
    mut n: ::std::os::raw::c_int,
    mut Pe: *mut ::std::os::raw::c_int,
    mut Iw: *mut ::std::os::raw::c_int,
    mut Len: *mut ::std::os::raw::c_int,
    mut iwlen: ::std::os::raw::c_int,
    mut pfree: ::std::os::raw::c_int,
    mut Nv: *mut ::std::os::raw::c_int,
    mut Next: *mut ::std::os::raw::c_int,
    mut Last: *mut ::std::os::raw::c_int,
    mut Head: *mut ::std::os::raw::c_int,
    mut Elen: *mut ::std::os::raw::c_int,
    mut Degree: *mut ::std::os::raw::c_int,
    mut W: *mut ::std::os::raw::c_int,
    mut Control: *mut c_float,
    mut Info: *mut c_float,
) {
    let mut deg: ::std::os::raw::c_int = 0;
    let mut degme: ::std::os::raw::c_int = 0;
    let mut dext: ::std::os::raw::c_int = 0;
    let mut lemax: ::std::os::raw::c_int = 0;
    let mut e: ::std::os::raw::c_int = 0;
    let mut elenme: ::std::os::raw::c_int = 0;
    let mut eln: ::std::os::raw::c_int = 0;
    let mut i: ::std::os::raw::c_int = 0;
    let mut ilast: ::std::os::raw::c_int = 0;
    let mut inext: ::std::os::raw::c_int = 0;
    let mut j: ::std::os::raw::c_int = 0;
    let mut jlast: ::std::os::raw::c_int = 0;
    let mut jnext: ::std::os::raw::c_int = 0;
    let mut k: ::std::os::raw::c_int = 0;
    let mut knt1: ::std::os::raw::c_int = 0;
    let mut knt2: ::std::os::raw::c_int = 0;
    let mut knt3: ::std::os::raw::c_int = 0;
    let mut lenj: ::std::os::raw::c_int = 0;
    let mut ln: ::std::os::raw::c_int = 0;
    let mut me: ::std::os::raw::c_int = 0;
    let mut mindeg: ::std::os::raw::c_int = 0;
    let mut nel: ::std::os::raw::c_int = 0;
    let mut nleft: ::std::os::raw::c_int = 0;
    let mut nvi: ::std::os::raw::c_int = 0;
    let mut nvj: ::std::os::raw::c_int = 0;
    let mut nvpiv: ::std::os::raw::c_int = 0;
    let mut slenme: ::std::os::raw::c_int = 0;
    let mut wbig: ::std::os::raw::c_int = 0;
    let mut we: ::std::os::raw::c_int = 0;
    let mut wflg: ::std::os::raw::c_int = 0;
    let mut wnvi: ::std::os::raw::c_int = 0;
    let mut ok: ::std::os::raw::c_int = 0;
    let mut ndense: ::std::os::raw::c_int = 0;
    let mut ncmpa: ::std::os::raw::c_int = 0;
    let mut dense: ::std::os::raw::c_int = 0;
    let mut aggressive: ::std::os::raw::c_int = 0;
    let mut hash: ::std::os::raw::c_uint = 0;
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
    let mut p: ::std::os::raw::c_int = 0;
    let mut p1: ::std::os::raw::c_int = 0;
    let mut p2: ::std::os::raw::c_int = 0;
    let mut p3: ::std::os::raw::c_int = 0;
    let mut p4: ::std::os::raw::c_int = 0;
    let mut pdst: ::std::os::raw::c_int = 0;
    let mut pend: ::std::os::raw::c_int = 0;
    let mut pj: ::std::os::raw::c_int = 0;
    let mut pme: ::std::os::raw::c_int = 0;
    let mut pme1: ::std::os::raw::c_int = 0;
    let mut pme2: ::std::os::raw::c_int = 0;
    let mut pn: ::std::os::raw::c_int = 0;
    let mut psrc: ::std::os::raw::c_int = 0;
    lnz = 0 as ::std::os::raw::c_int as c_float;
    ndiv = 0 as ::std::os::raw::c_int as c_float;
    nms_lu = 0 as ::std::os::raw::c_int as c_float;
    nms_ldl = 0 as ::std::os::raw::c_int as c_float;
    dmax = 1 as ::std::os::raw::c_int as c_float;
    me = EMPTY;
    mindeg = 0 as ::std::os::raw::c_int;
    ncmpa = 0 as ::std::os::raw::c_int;
    nel = 0 as ::std::os::raw::c_int;
    lemax = 0 as ::std::os::raw::c_int;
    if !Control.is_null() {
        alpha = *Control.offset(AMD_DENSE as isize);
        aggressive = (*Control.offset(AMD_AGGRESSIVE as isize)
            != 0 as ::std::os::raw::c_int as ::std::os::raw::c_double) as ::std::os::raw::c_int;
    } else {
        alpha = AMD_DEFAULT_DENSE;
        aggressive = AMD_DEFAULT_AGGRESSIVE;
    }
    if alpha < 0 as ::std::os::raw::c_int as ::std::os::raw::c_double {
        dense = n - 2 as ::std::os::raw::c_int;
    } else {
        dense = (alpha * sqrt(n as c_float)) as ::std::os::raw::c_int;
    }
    dense = if 16 as ::std::os::raw::c_int > dense { 16 as ::std::os::raw::c_int } else { dense };
    dense = if n < dense { n } else { dense };
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *Last.offset(i as isize) = EMPTY;
        *Head.offset(i as isize) = EMPTY;
        *Next.offset(i as isize) = EMPTY;
        *Nv.offset(i as isize) = 1 as ::std::os::raw::c_int;
        *W.offset(i as isize) = 1 as ::std::os::raw::c_int;
        *Elen.offset(i as isize) = 0 as ::std::os::raw::c_int;
        *Degree.offset(i as isize) = *Len.offset(i as isize);
        i += 1;
    }
    wbig = Int_MAX - n;
    wflg = clear_flag(0 as ::std::os::raw::c_int, wbig, W, n);
    ndense = 0 as ::std::os::raw::c_int;
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        deg = *Degree.offset(i as isize);
        if deg == 0 as ::std::os::raw::c_int {
            *Elen.offset(i as isize) = -(1 as ::std::os::raw::c_int) - 2 as ::std::os::raw::c_int;
            nel += 1;
            *Pe.offset(i as isize) = EMPTY;
            *W.offset(i as isize) = 0 as ::std::os::raw::c_int;
        } else if deg > dense {
            ndense += 1;
            *Nv.offset(i as isize) = 0 as ::std::os::raw::c_int;
            *Elen.offset(i as isize) = EMPTY;
            nel += 1;
            *Pe.offset(i as isize) = EMPTY;
        } else {
            inext = *Head.offset(deg as isize);
            if inext != EMPTY {
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
            if me != EMPTY {
                break;
            }
            deg += 1;
        }
        mindeg = deg;
        inext = *Next.offset(me as isize);
        if inext != EMPTY {
            *Last.offset(inext as isize) = EMPTY;
        }
        *Head.offset(deg as isize) = inext;
        elenme = *Elen.offset(me as isize);
        nvpiv = *Nv.offset(me as isize);
        nel += nvpiv;
        *Nv.offset(me as isize) = -nvpiv;
        degme = 0 as ::std::os::raw::c_int;
        if elenme == 0 as ::std::os::raw::c_int {
            pme1 = *Pe.offset(me as isize);
            pme2 = pme1 - 1 as ::std::os::raw::c_int;
            p = pme1;
            while p <= pme1 + *Len.offset(me as isize) - 1 as ::std::os::raw::c_int {
                i = *Iw.offset(p as isize);
                nvi = *Nv.offset(i as isize);
                if nvi > 0 as ::std::os::raw::c_int {
                    degme += nvi;
                    *Nv.offset(i as isize) = -nvi;
                    pme2 += 1;
                    *Iw.offset(pme2 as isize) = i;
                    ilast = *Last.offset(i as isize);
                    inext = *Next.offset(i as isize);
                    if inext != EMPTY {
                        *Last.offset(inext as isize) = ilast;
                    }
                    if ilast != EMPTY {
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
            knt1 = 1 as ::std::os::raw::c_int;
            while knt1 <= elenme + 1 as ::std::os::raw::c_int {
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
                knt2 = 1 as ::std::os::raw::c_int;
                while knt2 <= ln {
                    let fresh1 = pj;
                    pj = pj + 1;
                    i = *Iw.offset(fresh1 as isize);
                    nvi = *Nv.offset(i as isize);
                    if nvi > 0 as ::std::os::raw::c_int {
                        if pfree >= iwlen {
                            *Pe.offset(me as isize) = p;
                            *Len.offset(me as isize) -= knt1;
                            if *Len.offset(me as isize) == 0 as ::std::os::raw::c_int {
                                *Pe.offset(me as isize) = EMPTY;
                            }
                            *Pe.offset(e as isize) = pj;
                            *Len.offset(e as isize) = ln - knt2;
                            if *Len.offset(e as isize) == 0 as ::std::os::raw::c_int {
                                *Pe.offset(e as isize) = EMPTY;
                            }
                            ncmpa += 1;
                            j = 0 as ::std::os::raw::c_int;
                            while j < n {
                                pn = *Pe.offset(j as isize);
                                if pn >= 0 as ::std::os::raw::c_int {
                                    *Pe.offset(j as isize) = *Iw.offset(pn as isize);
                                    *Iw.offset(pn as isize) = -j - 2 as ::std::os::raw::c_int;
                                }
                                j += 1;
                            }
                            psrc = 0 as ::std::os::raw::c_int;
                            pdst = 0 as ::std::os::raw::c_int;
                            pend = pme1 - 1 as ::std::os::raw::c_int;
                            while psrc <= pend {
                                let fresh2 = psrc;
                                psrc = psrc + 1;
                                j = -*Iw.offset(fresh2 as isize) - 2 as ::std::os::raw::c_int;
                                if j >= 0 as ::std::os::raw::c_int {
                                    *Iw.offset(pdst as isize) = *Pe.offset(j as isize);
                                    let fresh3 = pdst;
                                    pdst = pdst + 1;
                                    *Pe.offset(j as isize) = fresh3;
                                    lenj = *Len.offset(j as isize);
                                    knt3 = 0 as ::std::os::raw::c_int;
                                    while knt3 <= lenj - 2 as ::std::os::raw::c_int {
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
                            while psrc <= pfree - 1 as ::std::os::raw::c_int {
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
                        if inext != EMPTY {
                            *Last.offset(inext as isize) = ilast;
                        }
                        if ilast != EMPTY {
                            *Next.offset(ilast as isize) = inext;
                        } else {
                            *Head.offset(*Degree.offset(i as isize) as isize) = inext;
                        }
                    }
                    knt2 += 1;
                }
                if e != me {
                    *Pe.offset(e as isize) = -me - 2 as ::std::os::raw::c_int;
                    *W.offset(e as isize) = 0 as ::std::os::raw::c_int;
                }
                knt1 += 1;
            }
            pme2 = pfree - 1 as ::std::os::raw::c_int;
        }
        *Degree.offset(me as isize) = degme;
        *Pe.offset(me as isize) = pme1;
        *Len.offset(me as isize) = pme2 - pme1 + 1 as ::std::os::raw::c_int;
        *Elen.offset(me as isize) = -(nvpiv + degme) - 2 as ::std::os::raw::c_int;
        wflg = clear_flag(wflg, wbig, W, n);
        pme = pme1;
        while pme <= pme2 {
            i = *Iw.offset(pme as isize);
            eln = *Elen.offset(i as isize);
            if eln > 0 as ::std::os::raw::c_int {
                nvi = -*Nv.offset(i as isize);
                wnvi = wflg - nvi;
                p = *Pe.offset(i as isize);
                while p <= *Pe.offset(i as isize) + eln - 1 as ::std::os::raw::c_int {
                    e = *Iw.offset(p as isize);
                    we = *W.offset(e as isize);
                    if we >= wflg {
                        we -= nvi;
                    } else if we != 0 as ::std::os::raw::c_int {
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
            p2 = p1 + *Elen.offset(i as isize) - 1 as ::std::os::raw::c_int;
            pn = p1;
            hash = 0 as ::std::os::raw::c_int as ::std::os::raw::c_uint;
            deg = 0 as ::std::os::raw::c_int;
            if aggressive != 0 {
                p = p1;
                while p <= p2 {
                    e = *Iw.offset(p as isize);
                    we = *W.offset(e as isize);
                    if we != 0 as ::std::os::raw::c_int {
                        dext = we - wflg;
                        if dext > 0 as ::std::os::raw::c_int {
                            deg += dext;
                            let fresh8 = pn;
                            pn = pn + 1;
                            *Iw.offset(fresh8 as isize) = e;
                            hash = hash.wrapping_add(e as ::std::os::raw::c_uint);
                        } else {
                            *Pe.offset(e as isize) = -me - 2 as ::std::os::raw::c_int;
                            *W.offset(e as isize) = 0 as ::std::os::raw::c_int;
                        }
                    }
                    p += 1;
                }
            } else {
                p = p1;
                while p <= p2 {
                    e = *Iw.offset(p as isize);
                    we = *W.offset(e as isize);
                    if we != 0 as ::std::os::raw::c_int {
                        dext = we - wflg;
                        deg += dext;
                        let fresh9 = pn;
                        pn = pn + 1;
                        *Iw.offset(fresh9 as isize) = e;
                        hash = hash.wrapping_add(e as ::std::os::raw::c_uint);
                    }
                    p += 1;
                }
            }
            *Elen.offset(i as isize) = pn - p1 + 1 as ::std::os::raw::c_int;
            p3 = pn;
            p4 = p1 + *Len.offset(i as isize);
            p = p2 + 1 as ::std::os::raw::c_int;
            while p < p4 {
                j = *Iw.offset(p as isize);
                nvj = *Nv.offset(j as isize);
                if nvj > 0 as ::std::os::raw::c_int {
                    deg += nvj;
                    let fresh10 = pn;
                    pn = pn + 1;
                    *Iw.offset(fresh10 as isize) = j;
                    hash = hash.wrapping_add(j as ::std::os::raw::c_uint);
                }
                p += 1;
            }
            if *Elen.offset(i as isize) == 1 as ::std::os::raw::c_int && p3 == pn {
                *Pe.offset(i as isize) = -me - 2 as ::std::os::raw::c_int;
                nvi = -*Nv.offset(i as isize);
                degme -= nvi;
                nvpiv += nvi;
                nel += nvi;
                *Nv.offset(i as isize) = 0 as ::std::os::raw::c_int;
                *Elen.offset(i as isize) = EMPTY;
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
                *Len.offset(i as isize) = pn - p1 + 1 as ::std::os::raw::c_int;
                hash = hash.wrapping_rem(n as ::std::os::raw::c_uint);
                j = *Head.offset(hash as isize);
                if j <= EMPTY {
                    *Next.offset(i as isize) = -j - 2 as ::std::os::raw::c_int;
                    *Head.offset(hash as isize) = -i - 2 as ::std::os::raw::c_int;
                } else {
                    *Next.offset(i as isize) = *Last.offset(j as isize);
                    *Last.offset(j as isize) = i;
                }
                *Last.offset(i as isize) = hash as ::std::os::raw::c_int;
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
            if *Nv.offset(i as isize) < 0 as ::std::os::raw::c_int {
                hash = *Last.offset(i as isize) as ::std::os::raw::c_uint;
                j = *Head.offset(hash as isize);
                if j == EMPTY {
                    i = EMPTY;
                } else if j < EMPTY {
                    i = -j - 2 as ::std::os::raw::c_int;
                    *Head.offset(hash as isize) = EMPTY;
                } else {
                    i = *Last.offset(j as isize);
                    *Last.offset(j as isize) = EMPTY;
                }
                while i != EMPTY && *Next.offset(i as isize) != EMPTY {
                    ln = *Len.offset(i as isize);
                    eln = *Elen.offset(i as isize);
                    p = *Pe.offset(i as isize) + 1 as ::std::os::raw::c_int;
                    while p <= *Pe.offset(i as isize) + ln - 1 as ::std::os::raw::c_int {
                        *W.offset(*Iw.offset(p as isize) as isize) = wflg;
                        p += 1;
                    }
                    jlast = i;
                    j = *Next.offset(i as isize);
                    while j != EMPTY {
                        ok = (*Len.offset(j as isize) == ln
                            && *Elen.offset(j as isize) == eln) as ::std::os::raw::c_int;
                        p = *Pe.offset(j as isize) + 1 as ::std::os::raw::c_int;
                        while ok != 0
                            && p <= *Pe.offset(j as isize) + ln - 1 as ::std::os::raw::c_int
                        {
                            if *W.offset(*Iw.offset(p as isize) as isize) != wflg {
                                ok = 0 as ::std::os::raw::c_int;
                            }
                            p += 1;
                        }
                        if ok != 0 {
                            *Pe.offset(j as isize) = -i - 2 as ::std::os::raw::c_int;
                            *Nv.offset(i as isize) += *Nv.offset(j as isize);
                            *Nv.offset(j as isize) = 0 as ::std::os::raw::c_int;
                            *Elen.offset(j as isize) = EMPTY;
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
            if nvi > 0 as ::std::os::raw::c_int {
                *Nv.offset(i as isize) = nvi;
                deg = *Degree.offset(i as isize) + degme - nvi;
                deg = if deg < nleft - nvi { deg } else { nleft - nvi };
                inext = *Head.offset(deg as isize);
                if inext != EMPTY {
                    *Last.offset(inext as isize) = i;
                }
                *Next.offset(i as isize) = inext;
                *Last.offset(i as isize) = EMPTY;
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
        if *Len.offset(me as isize) == 0 as ::std::os::raw::c_int {
            *Pe.offset(me as isize) = EMPTY;
            *W.offset(me as isize) = 0 as ::std::os::raw::c_int;
        }
        if elenme != 0 as ::std::os::raw::c_int {
            pfree = p;
        }
        if !Info.is_null() {
            f = nvpiv as c_float;
            r = (degme + ndense) as c_float;
            dmax = if dmax > f + r { dmax } else { f + r };
            lnzme = f * r
                + (f - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double) * f
                    / 2 as ::std::os::raw::c_int as ::std::os::raw::c_double;
            lnz += lnzme;
            ndiv += lnzme;
            s = f * r * r + r * (f - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double) * f
                + (f - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double) * f
                    * (2 as ::std::os::raw::c_int as ::std::os::raw::c_double * f
                        - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double)
                    / 6 as ::std::os::raw::c_int as ::std::os::raw::c_double;
            nms_lu += s;
            nms_ldl += (s + lnzme) / 2 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        }
    }
    if !Info.is_null() {
        f = ndense as c_float;
        dmax = if dmax > ndense as c_float { dmax } else { ndense as c_float };
        lnzme = (f - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double) * f
            / 2 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        lnz += lnzme;
        ndiv += lnzme;
        s = (f - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double) * f
            * (2 as ::std::os::raw::c_int as ::std::os::raw::c_double * f
                - 1 as ::std::os::raw::c_int as ::std::os::raw::c_double)
            / 6 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        nms_lu += s;
        nms_ldl += (s + lnzme) / 2 as ::std::os::raw::c_int as ::std::os::raw::c_double;
        *Info.offset(AMD_LNZ as isize) = lnz;
        *Info.offset(AMD_NDIV as isize) = ndiv;
        *Info.offset(AMD_NMULTSUBS_LDL as isize) = nms_ldl;
        *Info.offset(AMD_NMULTSUBS_LU as isize) = nms_lu;
        *Info.offset(AMD_NDENSE as isize) = ndense as c_float;
        *Info.offset(AMD_DMAX as isize) = dmax;
        *Info.offset(AMD_NCMPA as isize) = ncmpa as c_float;
        *Info.offset(AMD_STATUS as isize) = AMD_OK as c_float;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *Pe.offset(i as isize) = -*Pe.offset(i as isize) - 2 as ::std::os::raw::c_int;
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        *Elen.offset(i as isize) = -*Elen.offset(i as isize) - 2 as ::std::os::raw::c_int;
        i += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        if *Nv.offset(i as isize) == 0 as ::std::os::raw::c_int {
            j = *Pe.offset(i as isize);
            if !(j == EMPTY) {
                while *Nv.offset(j as isize) == 0 as ::std::os::raw::c_int {
                    j = *Pe.offset(j as isize);
                }
                e = j;
                j = i;
                while *Nv.offset(j as isize) == 0 as ::std::os::raw::c_int {
                    jnext = *Pe.offset(j as isize);
                    *Pe.offset(j as isize) = e;
                    j = jnext;
                }
            }
        }
        i += 1;
    }
    amd_postorder(n, Pe, Nv, Elen, W, Head, Next, Last);
    k = 0 as ::std::os::raw::c_int;
    while k < n {
        *Head.offset(k as isize) = EMPTY;
        *Next.offset(k as isize) = EMPTY;
        k += 1;
    }
    e = 0 as ::std::os::raw::c_int;
    while e < n {
        k = *W.offset(e as isize);
        if k != EMPTY {
            *Head.offset(k as isize) = e;
        }
        e += 1;
    }
    nel = 0 as ::std::os::raw::c_int;
    k = 0 as ::std::os::raw::c_int;
    while k < n {
        e = *Head.offset(k as isize);
        if e == EMPTY {
            break;
        }
        *Next.offset(e as isize) = nel;
        nel += *Nv.offset(e as isize);
        k += 1;
    }
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        if *Nv.offset(i as isize) == 0 as ::std::os::raw::c_int {
            e = *Pe.offset(i as isize);
            if e != EMPTY {
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
    i = 0 as ::std::os::raw::c_int;
    while i < n {
        k = *Next.offset(i as isize);
        *Last.offset(k as isize) = i;
        i += 1;
    }
}
