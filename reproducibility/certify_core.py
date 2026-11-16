#!/usr/bin/env python3
"""Directed-interval verification of the FGUM appendix inequalities.

This is a reconstruction from the manuscript's equations, not recovered
original certificate code.  Floating point calculations only propose root
brackets and Chernoff tilts.  Every accepted inequality is checked with
mpmath.iv; failed or inconclusive checks never count as certificates.

The analytic peeling, conditional-uniformity, and first-moment arguments
remain hypotheses of this numerical verification.  In particular this
program is not a formal verification of the entire rank theorem.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

import mpmath
from mpmath import mp, iv


SCHEMA = "fgum-core-interval-certificate-v1"
PROVENANCE = "Reconstructed from the manuscript; not recovered original output"
PAIRS = tuple((k, D) for k in range(3, 8) for D in range(k + 1, 9))
ORIGINAL_TARGETS = {
    (3,4): ("-1.082e-4", "-1.287e-1", ".339728"),
    (3,5): ("-7.759e-5", "-1.358e-1", ".421490"),
    (3,6): ("-6.080e-5", "-1.091e-1", ".463277"),
    (3,7): ("-5.017e-5", "-3.863e-4", ".491235"),
    (3,8): ("-4.280e-5", "-3.686e-1", ".512585"),
    (4,5): ("-1.202e-3", "-1.373e-3", ".101722"),
    (4,6): ("-1.093e-3", "-5.389e-4", ".179980"),
    (4,7): ("-1.012e-3", "-2.107e-1", ".232882"),
    (4,8): ("-9.495e-4", "-6.869e-2", ".271464"),
    (5,6): ("-2.851e-3", "-1.497e-3", ".019509"),
    (5,7): ("-2.720e-3", "-2.500e-3", ".051637"),
    (5,8): ("-2.616e-3", "-2.195e-3", ".079882"),
    (6,7): ("-4.987e-3", "-3.154e-5", ".002578"),
    (6,8): ("-1.713e-3", "-4.959e-4", ".010985"),
    (7,8): ("-7.589e-3", "-7.686e-4", ".000255"),
}
# The original manuscript rounded these spectral radii down even though the
# column was labeled an upper bound.  These replacements are rounded up.
SPECTRAL_CORRECTIONS = {
    (3,5): ".421491", (3,6): ".463278", (3,7): ".491236", (3,8): ".512586",
    (4,5): ".101723", (4,6): ".179981", (4,8): ".271465",
    (5,7): ".051638", (7,8): ".000256",
}
TARGETS = {pair:(values[0],values[1],SPECTRAL_CORRECTIONS.get(pair,values[2]))
           for pair,values in ORIGINAL_TARGETS.items()}
DELTA = Fraction(1, 10000)
CENTRAL_LO, CENTRAL_HI = Fraction(2, 5), Fraction(3, 5)
ROOT_RADIUS = Fraction(1, 10**45)
ANALYTIC_COVERAGE = {
    "small_support":"0 < alpha <= 1/10000; beta < 1",
    "small_alpha_exponent":"Phi <= c alpha [k/2 (1+log(D-1)) + (k/2-1) log(alpha)]",
    "large_support":"alpha > alpha_max + o(1) excluded by odd rows",
    "central_zero":"Gamma(0)=Gamma'(0)=0 analytically at exact balance",
    "exterior_convexity":"U''=(k-1)c/[alpha(1-alpha)] > 0",
}


class CertificateFailure(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise CertificateFailure(message)


def exact_bound(raw):
    """Convert a finite internal mpf endpoint to its exact rational value."""
    sign, mantissa, exponent, bits = raw
    require(bits >= 0, "Nonfinite interval endpoint")
    value = Fraction(int(mantissa))
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << -exponent)
    return -value if sign else value


def lower(value):
    return exact_bound(value._mpi_[0])


def upper(value):
    return exact_bound(value._mpi_[1])


def I(value):
    """An outward interval containing an exact integer, decimal, or fraction."""
    if isinstance(value, Fraction):
        return iv.mpf(value.numerator) / value.denominator
    return iv.mpf(value)


def box(lo, hi):
    lo, hi = Fraction(lo), Fraction(hi)
    require(lo <= hi, "Reversed input interval")
    return iv.mpf([I(lo).a, I(hi).b])


def decimal_out(value, upper_bound, digits=50):
    """Directed decimal conversion, preserving the direction of every bound."""
    value = Fraction(value)
    with localcontext() as ctx:
        ctx.prec = digits
        ctx.rounding = ROUND_CEILING if upper_bound else ROUND_FLOOR
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def enclosure(value):
    return {"lower": decimal_out(lower(value), False),
            "upper": decimal_out(upper(value), True)}


def ceil_decimal(value, places):
    scale = 10**places
    scaled = Fraction(value) * scale
    return Fraction(-((-scaled.numerator) // scaled.denominator), scale)


def phi(z, D):
    return 1 - (1-z)**(D-1)


def pi(z, D):
    return 1-(1-z)**D-D*z*(1-z)**(D-1)


def parameters(y, k, D, context):
    p = 1-(1-y)**(context.mpf(1)/(D-1))
    q = y**(k-1)
    a = phi(q,D)
    e = p/(a*y**(k-2))
    c = D*p*y
    balance = e*pi(q,D)+(k-1)*pi(p,D)-c
    return {"y":y,"p":p,"q":q,"a":a,"e":e,"c":c,"balance":balance}


def propose_root(k, D, dps):
    """Proposal only.  Certification subsequently checks the bracket signs."""
    with mp.workdps(dps+25):
        lo,hi=mp.mpf('0.1'),mp.mpf(1)
        require(parameters(lo,k,D,mp)["balance"] > 0, "Root proposal needs a new lower seed")
        require(parameters(hi,k,D,mp)["balance"] < 0, "Root proposal needs a new upper seed")
        for _ in range(int((dps+15)*3.4)):
            mid=(lo+hi)/2
            if parameters(mid,k,D,mp)["balance"] > 0: lo=mid
            else: hi=mid
        center=Fraction(mp.nstr((lo+hi)/2,dps+10))
    return center-ROOT_RADIUS,center+ROOT_RADIUS


def root_data(k,D,lo,hi):
    require(0 < lo < hi < 1, "Root bracket outside (0,1)")
    require(hi-lo <= 2*ROOT_RADIUS, "Root bracket wider than stated radius")
    at_lo=parameters(I(lo),k,D,iv)["balance"]
    at_hi=parameters(I(hi),k,D,iv)["balance"]
    require(lower(at_lo)>0 and upper(at_hi)<0, "Balance endpoint signs not certified")
    data=parameters(box(lo,hi),k,D,iv)
    y,p,q,a,e,c=(data[key] for key in ("y","p","q","a","e","c"))
    require(lower(e)>0 and upper(e)<1, "Erasure probability outside (0,1)")
    pprime=(1-y)**(I(1)/(D-1)-1)/(D-1)
    aprime=(D-1)*(1-q)**(D-2)*(k-1)*y**(k-2)
    eprime=e*(pprime/p-aprime/a-(k-2)/y)
    balance_prime=eprime*(pi(q,D)-D*q*a)
    require(lower(eprime)>0, "Local increasing-branch condition failed")
    require(upper(balance_prime)<0, "Balance derivative not strictly negative")
    rhoA=[e*math.comb(D,d)*q**d*(1-q)**(D-d) for d in range(2,D+1)]
    rhoB=[math.comb(D,d)*p**d*(1-p)**(D-d) for d in range(2,D+1)]
    VA=sum(r*(I(1) if d==2 else I(d)/4) for d,r in enumerate(rhoA,2))
    VB=sum(r*(I(1) if d==2 else I(d)/4) for d,r in enumerate(rhoB,2))
    require(lower(VA)>0 and lower(VB)>0, "Nonpositive tilt variance")
    oddA=sum(rhoA[d-2] for d in range(3,D+1,2))
    oddB=sum(rhoB[d-2] for d in range(3,D+1,2))
    odd_max=box(max(lower(oddA),lower(oddB)),max(upper(oddA),upper(oddB)))
    alpha_max=1-odd_max/c
    require(lower(odd_max)>0, "Odd-degree cutoff not strict")
    require(lower(alpha_max)>CENTRAL_HI and upper(alpha_max)<1, "Unexpected cutoff range")
    exterior_end=ceil_decimal(upper(alpha_max),45)+Fraction(1,1000)
    require(exterior_end<1, "Exterior cover would cross alpha=1")
    beta=iv.exp(1)*(D-1)**(I(k)/2)*I(DELTA)**(I(k)/2-1)*(1-I(DELTA))**(-I(k)/2)
    require(upper(beta)<1, "Small-support beta is not below one")
    small_alpha_factor=I(k)/2*(1+iv.log(D-1))+(I(k)/2-1)*iv.log(I(DELTA))
    require(upper(small_alpha_factor)<0, "Small-alpha variational bound is not negative")
    J01=aprime
    J10=(D-1)*(1-p)**(D-2)*e*y**(k-2)
    J11=(D-1)*(1-p)**(D-2)*e*a*(k-2)*y**(k-3)
    spectral=(J11+iv.sqrt(J11**2+4*J01*J10))/2
    require(upper(spectral)<1, "Jacobian stability not certified")
    require(upper(spectral)<=Fraction(TARGETS[k,D][2]), "Corrected spectral-radius upper bound not verified")
    x_end=ceil_decimal(upper(1-p),45)
    data.update({"rhoA":rhoA,"rhoB":rhoB,"VA":VA,"VB":VB,
                 "oddA":oddA,"oddB":oddB,"alpha_max":alpha_max,
                 "exterior_end":exterior_end,"beta":beta,"spectral":spectral,
                 "small_alpha_factor":small_alpha_factor,
                 "x_end":x_end,"eprime":eprime,"balance_prime":balance_prime,
                 "balance_at_lower":at_lo,"balance_at_upper":at_hi})
    return data


def log_polynomial(rho,t):
    total=I(0)
    for d,r in enumerate(rho,2):
        z=sum(math.comb(d,j)*iv.exp(j*t) for j in range(0,d+1,2))
        total+=r*iv.log(z)
    return total


def second_derivative(rho,t):
    """Positive pair-sum variance identity avoids subtracting two moments."""
    result=I(0)
    for d,r in enumerate(rho,2):
        js=list(range(0,d+1,2))
        weights=[math.comb(d,j)*iv.exp(j*t) for j in js]
        z=sum(weights)
        numerator=sum(weights[i]*weights[j]*(js[i]-js[j])**2
                      for i in range(len(js)) for j in range(i+1,len(js)))
        result+=r*numerator/z**2
    return result


def gamma_second(data,k,alpha):
    c,VA,VB=data["c"],data["VA"],data["VB"]
    u=c*(alpha-I(Fraction(1,2)))
    return (second_derivative(data["rhoA"],u/VA)/VA**2-2/VA
            +(k-1)*(second_derivative(data["rhoB"],u/VB)/VB**2-2/VB)
            +(k-1)/(c*alpha*(1-alpha)))


def high_branch_polynomial(x,k,D):
    """T=Q+(1-x)Q_x; e_x=-T/Q^2, regular also at x=0."""
    y=1-x**(D-1)
    q=y**(k-1)
    a=phi(q,D)
    da_dy=(D-1)*(1-q)**(D-2)*(k-1)*y**(k-2)
    Q=a*y**(k-2)
    Q_y=da_dy*y**(k-2)+a*(k-2)*y**(k-3)
    return Q-(1-x)*(D-1)*x**(D-2)*Q_y


def numeric_profile(rho):
    return [float((lower(r)+upper(r))/2) for r in rho]


def propose_tilt(rho,target):
    """Binary64 proposal; no inequality is accepted on this calculation."""
    def mean(t):
        total=0.0
        for d,r in enumerate(rho,2):
            js=range(0,d+1,2)
            weights=[math.comb(d,j)*math.exp(j*t) for j in js]
            total+=r*sum(j*w for j,w in zip(js,weights))/sum(weights)
        return total
    lo,hi=-64.0,64.0
    if mean(hi)<=target:return repr(hi)
    require(mean(lo)<target, "Tilt proposal failed to bracket target")
    for _ in range(55):
        mid=(lo+hi)/2
        if mean(mid)<target:lo=mid
        else:hi=mid
    return repr((lo+hi)/2)


def exterior_bounds(data,k,lo,hi,tA,tB):
    tA,tB=I(tA),I(tB)
    const=log_polynomial(data["rhoA"],tA)+(k-1)*log_polynomial(data["rhoB"],tB)
    slope=tA+(k-1)*tB
    def evaluate(alpha):
        alpha=I(alpha)
        h=-alpha*iv.log(alpha)-(1-alpha)*iv.log(1-alpha)
        return const-data["c"]*alpha*slope-(k-1)*data["c"]*h
    return evaluate(lo),evaluate(hi)


def adaptive_cover(lo,hi,test,max_depth=40):
    cells=[]
    stack=[(Fraction(lo),Fraction(hi),0)]
    while stack:
        left,right,depth=stack.pop()
        record=test(left,right)
        if record is not None:
            cells.append({"lower":str(left),"upper":str(right),**record})
        else:
            require(depth<max_depth,"Adaptive cover did not certify its cell")
            mid=(left+right)/2
            stack.append((mid,right,depth+1))
            stack.append((left,mid,depth+1))
        require(len(cells)+len(stack)<100000,"Adaptive cover exceeded resource guard")
    return cells


def summarize(k,D,data,branch,central,exterior):
    gamma_max=max(Fraction(cell["gamma_second_derivative"]["upper"]) for cell in central)
    U_max=max(Fraction(cell[end]["upper"]) for part in exterior for cell in part["cells"]
              for end in ("U_at_lower","U_at_upper"))
    ihat_num=sum(Fraction(math.comb(D,w),2 if D%2==0 and w==D//2 else 1)*w
                 for w in range(1,D//2+1))
    ihat=ihat_num/(2**(D-1)-1)
    alpha=(1-data["e"])*I(1-Fraction(1,2**(D-1)))
    score=1-alpha*I(ihat)/D
    return {"e_core":enclosure(data["e"]),"y":enclosure(data["y"]),
            "core_columns_per_b":enclosure(data["c"]),
            "beta":enclosure(data["beta"]),"spectral_radius":enclosure(data["spectral"]),
            "small_alpha_phi_factor":enclosure(data["small_alpha_factor"]),
            "alpha_max":enclosure(data["alpha_max"]),
            "exterior_U_upper":decimal_out(U_max,True),
            "central_Gamma_second_upper":decimal_out(gamma_max,True),
            "branch_T_lower":decimal_out(min(Fraction(cell["T"]["lower"]) for cell in branch),False),
            "original_manuscript_spectral_upper_valid":upper(data["spectral"])<=Fraction(ORIGINAL_TARGETS[k,D][2]),
            "cells":{"branch":len(branch),"central":len(central),
                     "exterior":sum(len(part["cells"]) for part in exterior)},
            "ihat_star_exact":str(ihat),"alpha_min":enclosure(alpha),
            "fgum_score":enclosure(score)}


def certify_pair(k,D,dps=80):
    require((k,D) in PAIRS,"Unsupported parameter pair")
    require(dps>=70,"At least 70 decimal interval digits are required")
    require(mpmath.__version__=="1.3.0","This reconstruction is verified for mpmath 1.3.0")
    iv.dps=dps
    root_lo,root_hi=propose_root(k,D,dps)
    data=root_data(k,D,root_lo,root_hi)
    U_target,gamma_target,_=(Fraction(x) for x in TARGETS[k,D])
    def branch_test(lo,hi):
        T=high_branch_polynomial(box(lo,hi),k,D)
        return {"T":enclosure(T)} if lower(T)>0 else None
    branch=adaptive_cover(0,data["x_end"],branch_test)
    def central_test(lo,hi):
        result=gamma_second(data,k,box(lo,hi))
        return {"gamma_second_derivative":enclosure(result)} if upper(result)<=gamma_target else None
    central=adaptive_cover(CENTRAL_LO,CENTRAL_HI,central_test)
    c_mid=float((lower(data["c"])+upper(data["c"]))/2)
    A_num,B_num=numeric_profile(data["rhoA"]),numeric_profile(data["rhoB"])
    def exterior_test(lo,hi):
        target=c_mid*float((lo+hi)/2)
        tA,tB=propose_tilt(A_num,target),propose_tilt(B_num,target)
        u_lo,u_hi=exterior_bounds(data,k,lo,hi,tA,tB)
        if max(upper(u_lo),upper(u_hi))<=U_target:
            return {"tilt_A":tA,"tilt_B":tB,
                    "U_at_lower":enclosure(u_lo),"U_at_upper":enclosure(u_hi)}
        return None
    exterior=[]
    for lo,hi in ((DELTA,CENTRAL_LO),(CENTRAL_HI,data["exterior_end"])):
        exterior.append({"lower":str(lo),"upper":str(hi),
                         "cells":adaptive_cover(lo,hi,exterior_test)})
    summary=summarize(k,D,data,branch,central,exterior)
    return {"schema":SCHEMA,"status":"passed","provenance":PROVENANCE,
            "k":k,"D":D,"interval_decimal_precision":dps,
            "mpmath_version":mpmath.__version__,"published_upper_targets":{
                "exterior_U":TARGETS[k,D][0],"central_Gamma_second":TARGETS[k,D][1],
                "spectral_radius":TARGETS[k,D][2]},
            "original_manuscript_spectral_upper":ORIGINAL_TARGETS[k,D][2],
            "root":{"variable":"y","lower":str(root_lo),"upper":str(root_hi),
                    "radius":"1e-45","balance_at_lower":enclosure(data["balance_at_lower"]),
                    "balance_at_upper":enclosure(data["balance_at_upper"]),
                    "balance_derivative":enclosure(data["balance_prime"]),
                    "e_derivative":enclosure(data["eprime"])},
            "profiles":{"A":[enclosure(x) for x in data["rhoA"]],
                        "B":[enclosure(x) for x in data["rhoB"]],
                        "V_A":enclosure(data["VA"]),"V_B":enclosure(data["VB"]),
                        "odd_A":enclosure(data["oddA"]),"odd_B":enclosure(data["oddB"])},
            "high_branch":{"variable":"x=(1-y)^(1/(D-1))","lower":"0",
                           "upper":str(data["x_end"]),"cells":branch},
            "central":{"lower":str(CENTRAL_LO),"upper":str(CENTRAL_HI),"cells":central},
            "exterior":exterior,"summary":summary,"analytic_coverage":dict(ANALYTIC_COVERAGE)}


def checked_cells(record,expected_lo,expected_hi):
    require(Fraction(record["lower"])==expected_lo and Fraction(record["upper"])==expected_hi,
            "Cover has wrong endpoints")
    cells=record["cells"]
    require(bool(cells),"Empty interval cover")
    position=expected_lo
    for cell in cells:
        lo,hi=Fraction(cell["lower"]),Fraction(cell["upper"])
        require(lo==position and hi>lo and hi<=expected_hi,"Cover contains a gap, overlap, or reversed cell")
        position=hi
        yield cell,lo,hi
    require(position==expected_hi,"Incomplete interval cover")


def verify_certificate(document):
    """Replay the proof inputs, recomputing every inequality and full cover.

    Stored interval endpoints and summary claims are additionally checked to
    enclose the recomputed values, rather than trusted as proof inputs.
    """
    require(document["schema"]==SCHEMA,"Unknown certificate schema")
    require(document["status"]=="passed","Certificate does not declare passed status")
    require(document["provenance"]==PROVENANCE,"Incorrect reconstruction provenance")
    k,D=int(document["k"]),int(document["D"])
    require((k,D) in PAIRS,"Unsupported parameter pair")
    dps=int(document["interval_decimal_precision"])
    require(dps>=70,"Insufficient recorded precision")
    require(document["mpmath_version"]==mpmath.__version__=="1.3.0","mpmath version does not match verified runtime")
    expected_targets={"exterior_U":TARGETS[k,D][0],"central_Gamma_second":TARGETS[k,D][1],
                      "spectral_radius":TARGETS[k,D][2]}
    require(document["published_upper_targets"]==expected_targets,"Recorded numerical targets changed")
    require(document["original_manuscript_spectral_upper"]==ORIGINAL_TARGETS[k,D][2],"Original manuscript target changed")
    require(document["analytic_coverage"]==ANALYTIC_COVERAGE,"Recorded analytic coverage changed")
    iv.dps=dps
    root=document["root"]
    require(root["variable"]=="y" and Fraction(root["radius"])==ROOT_RADIUS,"Root convention changed")
    data=root_data(k,D,Fraction(root["lower"]),Fraction(root["upper"]))
    U_target,gamma_target,_=(Fraction(x) for x in TARGETS[k,D])
    def check_enclosure(stored,actual):
        require(Fraction(stored["lower"])<=lower(actual) and Fraction(stored["upper"])>=upper(actual),
                "Stored enclosure does not contain replayed enclosure")
    check_enclosure(root["balance_at_lower"],data["balance_at_lower"])
    check_enclosure(root["balance_at_upper"],data["balance_at_upper"])
    check_enclosure(root["balance_derivative"],data["balance_prime"])
    check_enclosure(root["e_derivative"],data["eprime"])
    require(Fraction(root["balance_at_lower"]["lower"])>0 and
            Fraction(root["balance_at_upper"]["upper"])<0 and
            Fraction(root["balance_derivative"]["upper"])<0 and
            Fraction(root["e_derivative"]["lower"])>0,"Stored root bounds fail to exhibit their strict signs")
    require(document["high_branch"]["variable"]=="x=(1-y)^(1/(D-1))","High-branch coordinate changed")
    for cell,lo,hi in checked_cells(document["high_branch"],Fraction(0),data["x_end"]):
        T=high_branch_polynomial(box(lo,hi),k,D)
        require(lower(T)>0,"High-branch monotonicity failed")
        check_enclosure(cell["T"],T)
        require(Fraction(cell["T"]["lower"])>0,"Stored branch bound is not positive")
    for cell,lo,hi in checked_cells(document["central"],CENTRAL_LO,CENTRAL_HI):
        result=gamma_second(data,k,box(lo,hi))
        require(upper(result)<=gamma_target,"Central curvature bound failed")
        check_enclosure(cell["gamma_second_derivative"],result)
        require(Fraction(cell["gamma_second_derivative"]["upper"])<=gamma_target,"Stored central bound exceeds target")
    require(len(document["exterior"])==2,"Expected exactly two exterior ranges")
    for part,expected in zip(document["exterior"],((DELTA,CENTRAL_LO),(CENTRAL_HI,data["exterior_end"]))):
        for cell,lo,hi in checked_cells(part,*expected):
            u_lo,u_hi=exterior_bounds(data,k,lo,hi,cell["tilt_A"],cell["tilt_B"])
            require(max(upper(u_lo),upper(u_hi))<=U_target,"Exterior endpoint bound failed")
            check_enclosure(cell["U_at_lower"],u_lo)
            check_enclosure(cell["U_at_upper"],u_hi)
            require(max(Fraction(cell["U_at_lower"]["upper"]),Fraction(cell["U_at_upper"]["upper"]))<=U_target,
                    "Stored exterior bound exceeds target")
    require(len(document["profiles"]["A"])==D-1 and len(document["profiles"]["B"])==D-1,
            "Wrong number of degree-profile entries")
    for stored,actual in zip(document["profiles"]["A"],data["rhoA"]):check_enclosure(stored,actual)
    for stored,actual in zip(document["profiles"]["B"],data["rhoB"]):check_enclosure(stored,actual)
    for saved,actual in (("V_A","VA"),("V_B","VB"),("odd_A","oddA"),("odd_B","oddB")):
        check_enclosure(document["profiles"][saved],data[actual])
    summary=summarize(k,D,data,document["high_branch"]["cells"],document["central"]["cells"],document["exterior"])
    require(document["summary"]==summary,"Stored summary does not match replayed certificate")
    return {"verified":True,"k":k,"D":D,"summary":summary}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k",type=int)
    parser.add_argument("--D",type=int)
    parser.add_argument("--dps",type=int,default=80)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--verify",type=Path)
    args=parser.parse_args()
    if args.verify:
        document=json.loads(args.verify.read_text())
        result=verify_certificate(document)
        print(json.dumps(result,sort_keys=True))
    else:
        require(args.k is not None and args.D is not None,"Specify --k and --D")
        document=certify_pair(args.k,args.D,args.dps)
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(document,indent=2,sort_keys=True)+"\n")
            print(json.dumps({"k":args.k,"D":args.D,"summary":document["summary"]},sort_keys=True))
        else:
            print(json.dumps(document,indent=2,sort_keys=True))


if __name__=="__main__":
    try:
        main()
    except (CertificateFailure,ValueError,KeyError,ArithmeticError) as exc:
        print("CERTIFICATE FAILED: "+str(exc),file=sys.stderr)
        raise SystemExit(1)
