"""Build validation tables and figures from the seeded study outputs."""
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, PercentFormatter


def read(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def read_many(root, pattern):
    rows = []
    for path in sorted(root.glob(pattern)):
        rows.extend(read(path))
    return rows


def read_calibration_rows(root):
    rows = []
    for path in sorted(root.glob("cal_*_rows.csv")):
        if re.fullmatch(r"cal_\d+_(?:null|weak|strong)_rows\.csv", path.name):
            rows.extend(read(path))
    return rows


def write(dest, name, rows):
    with (dest / (name + ".csv")).open("w", newline="") as file:
        w = csv.DictWriter(file, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def num(r, key):
    try:
        return float(r[key])
    except (ValueError, TypeError):
        return math.nan


def integers(r, key):
    value = num(r, key)
    return int(value) if math.isfinite(value) else 0


def finite(vals):
    return [x for x in vals if math.isfinite(x)]


def med(vals):
    vals = finite(vals)
    return median(vals) if vals else math.nan


def avg(vals):
    vals = finite(vals)
    return mean(vals) if vals else math.nan


def group(rows, *keys):
    groups = defaultdict(list)
    for r in rows:
        name = tuple(k(r) if callable(k) else r[k] for k in keys)
        groups[name].append(r)
    return sorted(groups.items())


def term_family(r):
    return r["term"].split(".", 1)[0]


def outcome(rows):
    estimable = [r for r in rows if all(math.isfinite(num(r, k)) for k in ("low", "high", "truth", "estimate"))]
    null = [r for r in estimable if num(r, "truth") == 0]
    signal = [r for r in estimable if num(r, "truth") != 0]
    errors = finite(num(r, "estimate") - num(r, "truth") for r in rows)
    return dict(attempted=len(rows), available=len(estimable),
                coverage=avg(num(r, "low") <= num(r, "truth") <= num(r, "high") for r in estimable),
                false_positive_rate=avg(num(r, "low") > 0 or num(r, "high") < 0 for r in null),
                null_intervals=len(null), directional_power=avg(
                    num(r, "low") > 0 if num(r, "truth") > 0 else num(r, "high") < 0 for r in signal),
                nonzero_intervals=len(signal), bias=avg(errors),
                rmse=math.sqrt(avg(x*x for x in errors)), median_absolute_error=med(map(abs, errors)),
                median_width=med(num(r, "high") - num(r, "low") for r in estimable),
                mean_interval_midpoint_error=avg((num(r, "high") + num(r, "low")) / 2 - num(r, "truth")
                                                  for r in estimable))


def calibration(raw, dest, complete=True):
    estimates = read_calibration_rows(raw)
    changes = read_many(raw, "cal_*_changes.csv")
    diag = read_many(raw, "cal_*_diagnostics.csv")
    expected = {(str(n), signal) for n in (70, 100, 150, 250) for signal in ("null", "weak", "strong")}
    if complete and {(x["n"], x["signal"]) for x in estimates} != expected:
        raise ValueError("The 12 calibration cells must be complete.")
    write(dest, "calibration", [dict(n=n, signal=s, order=k, term=t, method=m, **outcome(v))
        for (n, s, k, t, m), v in group(estimates, lambda r: integers(r, "n"), "signal",
                                       lambda r: integers(r, "order"), "term", "method")])
    selected = [r for r in estimates if .02 <= num(r, "ratio") <= .20]
    main = [dict(n=n, family=f, method=m, **outcome(v))
        for (n, f, m), v in group(selected, lambda r: integers(r, "n"), term_family, "method")]
    write(dest, "calibration_main_range", main)
    write(dest, "calibration_changes", [dict(n=n, signal=s, from_order=low, to_order=k,
        term=t, method=m, **outcome(v))
        for (n, s, low, k, t, m), v in group(changes, lambda r: integers(r, "n"), "signal",
            lambda r: integers(r, "from"), lambda r: integers(r, "order"), "term", "method")])
    write(dest, "calibration_information", [dict(n=n, signal=s, order=k, budget=budget, n_fit=len(v),
        full_rank=sum(r["rank"] == r["parameters"] for r in v),
        converged=sum(r["converged"] == "TRUE" for r in v),
        exact=sum(r["exact"] == "TRUE" for r in v),
        median_condition=med(num(r, "condition") for r in v),
        median_min_eigenvalue=med(num(r, "lambda_min") for r in v),
        median_min_exclusion_ess=med(num(r, "omission_ess") for r in v),
        jackknife_cpu=med(num(r, "jack_cpu") for r in v),
        deletion_score_cpu=med(num(r, "multiplier_cpu") for r in v),
        bayesian_cpu=med(num(r, "bb_cpu") for r in v))
        for (n,s,k,budget),v in group(diag,lambda r: integers(r,"n"),"signal",lambda r: integers(r,"order"),"budget")])
    raw_mc = read_many(raw, "cal_*_mc_rows.csv")
    if raw_mc:
        write(dest, "calibration_mc_precision", [dict(n=n, order=k, budget=b, term=term,
            n_samples=len(v), sd=stdev(finite(num(r, "estimate") for r in v))
                 if len(finite(num(r, "estimate") for r in v)) > 1 else math.nan,
            exact=v[0]["exact"])
            for (n,k,b,term),v in group(raw_mc,lambda r: integers(r,"n"),lambda r: integers(r,"order"),
                                      lambda r: integers(r,"budget"),"term")])
    spatial = read_many(raw, "spatial_*_rows.csv")
    if spatial:
        write(dest, "calibration_spatial_stress", [dict(n=n, order=k, family=f, method=m, **outcome(v))
            for (n,k,f,m),v in group(spatial,lambda r: integers(r,"n"),lambda r: integers(r,"order"),
                                     term_family,"method")])
    return main, selected


def mc_empirical(raw_mc, diagnostics, uncertainty, meta):
    exact = {(r["case"], r["sampler"], r["order"], r["budget"]): r["exact"] == "TRUE"
             for r in diagnostics}
    jk = {(r["case"], r["sampler"], r["order"], r["from"], r["response"], r["term"]): num(r, "se")
          for r in uncertainty if r["method"] == "jackknife"}
    fields = defaultdict(dict)
    for r in raw_mc:
        k = (r["case"], r["sampler"], r["order"], r["budget"], r["response"], r["term"])
        fields[k][r["replicate"]] = num(r, "estimate")

    def sd(steps, all_exact):
        z = finite(steps)
        if all_exact and z:
            return (0.0, len(z))
        return (stdev(z) if len(z) > 1 else math.nan, len(z))

    def result(case, sampler, k, from_order, b, r, t, d):
        se = jk.get((case, sampler, k, from_order, r, t), math.nan)
        return dict(case=case, sampler=sampler, order=k, from_order=from_order,
                    budget=b, response=r, term=t, mc_sd=d[0], n_replicates=d[1],
                    jk_se_at_10000=se, mc_to_jk=d[0] / se if se > 0 else math.nan)

    out = []
    for (case,a,k,b,r,t), values in sorted(fields.items()):
        is_exact = exact[case, a, k, b]
        out.append(result(case,a,k,"NA",b,r,t,sd(values.values(), is_exact)))
        orders = [str(meta[case][name]) for name in ("low", "mid", "high")]
        if k == orders[0]:
            continue
        for lo in dict.fromkeys((orders[0], orders[orders.index(k) - 1])):
            prior = fields[case, a, lo, b, r, t]
            other_exact = exact[case, a, lo, b]
            if is_exact and values:
                vals = [next(iter(values.values())) - x for x in prior.values()]
            elif other_exact and prior:
                vals = [x - next(iter(prior.values())) for x in values.values()]
            else:
                vals = [values[j] - prior[j] for j in values.keys() & prior.keys()]
            out.append(result(case,a,k,lo,b,r,t,sd(vals, is_exact and other_exact)))
    return out


def empirical(raw, dest, meta):
    uncertainty = read_many(raw, "*_all_uncertainty.csv")
    profiles = read_many(raw, "*_profiles.csv")
    diagnostics = read_many(raw, "*_all_diagnostics.csv")
    assert {r["case"] for r in uncertainty} == set(meta)
    assert {r["case"] for r in profiles} == set(meta)
    support = read(raw / "contrast_support.csv")
    sup = {(r["case"],r["sampler"],r["order"],r["term"]): r["supported"] == "TRUE"
           for r in support if r["budget"] == "10000"}
    for r in uncertainty:
        keys = (r["case"], r["sampler"], r["order"], r["term"])
        lo = (r["case"], r["sampler"], r["from"], r["term"])
        r["marginal_support"] = (sup[keys] and sup.get(lo, True)) if keys in sup else ""
    write(dest, "empirical_intervals", uncertainty)
    write(dest, "empirical_profiles", profiles)
    write(dest, "empirical_support", support)
    write(dest, "empirical_information", [dict(case=c, sampler=a, order=k, budget=b, response=r,
        n_runs=len(v), sets=integers(v[0],"sets"), exact=v[0]["exact"],
        base_sharing=med(num(x,"target_share") for x in v),
        proposal_sharing=med(num(x,"proposal_share") for x in v),
        expected_shared_species=med(num(x,"expected_shared") for x in v),
        median_ess=med(num(x,"ess") for x in v), min_exclusion_ess=min(num(x,"min_exclusion_ess") for x in v),
        converged=sum(x["converged"]=="TRUE" for x in v),
        full_rank=sum(x["rank"]==x["parameters"] for x in v),
        median_condition=med(num(x,"condition") for x in v),
        median_min_eigenvalue=med(num(x,"lambda_min") for x in v),
        zero_active=sum(num(x,"active")==0 for x in v),
        point_active=integers(next(x for x in v if x["replicate"]=="1"),"active"))
        for (c,a,k,b,r),v in group(diagnostics,"case","sampler","order","budget","response")])
    precision = mc_empirical(read_many(raw,"*_all_mc_rows.csv"), diagnostics, uncertainty, meta)
    write(dest,"empirical_mc_precision",precision)
    write(dest,"empirical_timing",read_many(raw,"*_all_timings.csv"))
    scoring = read(raw / "holdout.csv")
    write(dest,"empirical_held_out",[dict(case=c,sampler=a,order=k,response=r,splits=len(v),
        converged=sum(x["converged"]=="TRUE" for x in v),
        median_full_rmse=med(num(x,"rmse") for x in v), median_spread_rmse=med(num(x,"rmse_spread") for x in v),
        median_rmse_gain=med(num(x,"gain") for x in v),
        gain_q10=quantile(num(x,"gain") for x in v), gain_q90=quantile((num(x,"gain") for x in v),.9))
        for (c,a,k,r),v in group(scoring,"case","sampler","order","response")])
    primary = [r for r in uncertainty if r["term"] == "share.position" and
               r["order"] == str(meta[r["case"]]["high"]) and r["from"] == "2"]
    write(dest, "position_share_changes", primary)
    refit = read_many(raw, "*_refit50000.csv")
    expected_refits = {"N78TPLANT", "S7MFI", "N68PDIAT"}
    assert {r["case"] for r in refit} == expected_refits
    write(dest, "position_share_changes_50000", refit)
    return primary


def quantile(vals, p=.1):
    z = sorted(finite(vals))
    if not z:
        return math.nan
    i = (len(z) - 1) * p
    j = int(i)
    return z[j] + (z[min(j+1,len(z)-1)] - z[j]) * (i-j)


METHODS = ("bayesian", "jackknife", "multiplier")
COLORS = ("#265889", "#b04a25", "#187463")
MARKERS = ("o", "s", "D")
LINES = ("-", "--", "-.")
LABELS = ("Conditional Poisson Bayesian", "Delete-one jackknife", "Deletion-score multiplier")


def diagrams(cal, empirical_rows, meta, dest):
    plt.rcParams.update({"axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,
                         "axes.axisbelow":True,"grid.alpha":.22,"font.size":11,"figure.facecolor":"white"})
    fig, axs = plt.subplots(1, 2, figsize=(10.6,4.6), sharey=True)
    fig.subplots_adjust(left=.09,right=.985,top=.81,bottom=.28,wspace=.13)
    for ax,f in zip(axs,("spread","position")):
        ax.axhline(.95,color=".42",linewidth=1)
        for j,m in enumerate(METHODS):
            data=sorted((r["n"],r["coverage"]) for r in cal if r["family"]==f and r["method"]==m)
            ax.plot(*zip(*data),color=COLORS[j],marker=MARKERS[j],linestyle=LINES[j],label=LABELS[j])
        ax.set(title=f"{f.capitalize()} effects",xticks=(70,100,150,250),ylim=(0,1.025),xlabel="Number of sites")
        ax.yaxis.set_major_formatter(PercentFormatter(1))
    axs[0].set_ylabel("Generating-parameter coverage")
    h,l=axs[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.5,.035),ncol=3,frameon=False,fontsize=9.5)
    fig.suptitle("Calibration at manuscript-scale sample sizes",fontsize=15,weight="bold",y=.96)
    fig.savefig(dest/"calibration_coverage.png",dpi=150)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11.6,5.8),sharex=True,sharey=True)
    fig.subplots_adjust(left=.21,right=.985,top=.84,bottom=.23,wspace=.07)
    codes=list(meta)
    label=[meta[c]["community"] for c in codes]
    for ax,response,title in zip(axs,("total","relative"),("Incidence","Relative abundance")):
        offsets={("conventional","bayesian"):-.24,("conventional","jackknife"):-.08,
                 ("species","bayesian"):.08,("species","jackknife"):.24}
        ax.axvline(0,color=".43",linewidth=1)
        for j,r in enumerate(empirical_rows):
            if r["response"]!=response:
                continue
            y=codes.index(r["case"])+offsets[r["sampler"],r["method"]]
            color=COLORS[0 if r["method"]=="bayesian" else 1]
            marker="o" if r["sampler"]=="conventional" else "s"
            lo,hi=num(r,"low"),num(r,"high")
            if math.isfinite(lo+hi):
                ax.hlines(y,lo,hi,color=color,alpha=.9,lw=1.7)
            ax.scatter(num(r,"estimate"),y,color=color,marker=marker,s=23)
        ax.set(title=title,yticks=list(range(len(codes))),yticklabels=label,ylim=(len(codes)-.45,-.52),
               xlabel="Change in position share (percentage points)")
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{100 * x:.0f}"))
    for a,m in (("conventional","o"),("species","s")):
        for method,j in (("Bayesian",0),("Jackknife",1)):
            axs[0].scatter([],[],color=COLORS[j],marker=m,s=26,label=f"{a.capitalize()} · {method}")
    h,l=axs[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.56,.03),ncol=2,frameon=False,fontsize=9)
    fig.suptitle("Position allocation at the breakpoint versus pairs",fontsize=15,weight="bold",y=.96)
    fig.savefig(dest/"position_share_changes.png",dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw",type=Path)
    parser.add_argument("output",type=Path)
    opt=parser.parse_args()
    opt.output.mkdir(parents=True,exist_ok=True)
    meta={r["id"]:r for r in read(Path(__file__).with_name("communities.csv"))}
    c,_=calibration(opt.raw,opt.output)
    e=empirical(opt.raw,opt.output,meta)
    diagrams(c,e,meta,opt.output)
