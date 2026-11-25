"""Create review tables and figures from the output of simulate.R."""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, PercentFormatter


def load(root, name):
    with (root / f"{name}.csv").open() as file:
        return list(csv.DictReader(file))


def export(root, name, rows):
    with (root / f"{name}.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def number(row, key):
    return float(row[key]) if row[key] not in ("NA", "NaN", "") else math.nan


def finite(values):
    return [value for value in values if math.isfinite(value)]


def med(values):
    values = finite(values)
    return median(values) if values else math.nan


def average(values):
    values = finite(values)
    return mean(values) if values else math.nan


def grouped(rows, *keys):
    result = defaultdict(list)
    for row in rows:
        values = tuple(key(row) if callable(key) else row[key] for key in keys)
        result[values].append(row)
    return sorted(result.items())


def family(row):
    return row["effect"].split(".")[0]


def order_label(row):
    ratio = int(row.get("order", row.get("to"))) / int(row["n"])
    return min(range(5), key=lambda i: abs(ratio - (.09, .3, .6, .85, .97)[i]))


def evaluate(data, change=False):
    keys = ("spatial", "signal", order_label, family, "method")
    if change:
        data = [row for row in data if row["effect"] == "spread.gradient"]
    result = []
    for (spatial, signal, j, effect, method), rows in grouped(data, *keys):
        intervals = [row for row in rows if all(math.isfinite(number(row, key)) for key in ("low", "high", "estimate"))]
        signed = [number(row, "estimate") - number(row, "truth") for row in rows]
        exact = sum(row.get("exact") == "TRUE" for row in rows)
        result.append(dict(spatial=spatial, signal=signal, stage=j, effect=effect, method=method,
            attempted=len(rows), valid=len(intervals), exact=exact, bias=average(signed),
            abs_error=med(map(abs, signed)), coverage=average(
                number(row, "low") <= number(row, "truth") <= number(row, "high") for row in intervals),
            excludes_zero=average(number(row, "low") > 0 or number(row, "high") < 0 for row in intervals),
            median_width=med(number(row, "high") - number(row, "low") for row in intervals)))
    return result


COLORS = {"direct_bb": "#264b80", "product_bb": "#8c6641", "jackknife": "#b94d29",
          "multiplier": "#187464", "multiplier_raw": "#725290"}
LABELS = {"direct_bb": "Direct Bayesian", "product_bb": "Reweighted Bayesian",
          "jackknife": "Delete-one jackknife", "multiplier": "Deletion-score multiplier",
          "multiplier_raw": "Raw inclusion multiplier"}
MARKERS = ["o", "s", "^", "D", "v"]
LINES = ["-", "--", "-.", ":", (0, (6, 2, 1, 2))]
ORDER = ["Pairs", "≈ .30", "≈ .60", "≈ .85", "n − 1"]
FAMILIES = ("spread", "position", "importance", "shapley")


def configure():
    plt.rcParams.update({"figure.dpi": 135, "savefig.dpi": 150, "font.size": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": .22, "axes.axisbelow": True,
        "figure.facecolor": "white", "axes.facecolor": "white"})


def calibration_plot(data, dest):
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 8), sharex=True, sharey=True)
    fig.subplots_adjust(left=.115, right=.98, bottom=.19, top=.88, hspace=.33, wspace=.13)
    for i, spatial in enumerate(("FALSE", "TRUE")):
        for j, effect in enumerate(("spread", "position")):
            ax = axs[i, j]
            ax.axhline(.95, color=".4", lw=1)
            for m, method in enumerate(COLORS):
                pairs = grouped([row for row in data if row["spatial"] == spatial and family(row) == effect and row["method"] == method], order_label)
                y = []
                for _, rows in pairs:
                    valid = [row for row in rows if math.isfinite(number(row, "low")) and math.isfinite(number(row, "high"))]
                    y.append(average(number(row, "low") <= number(row, "truth") <= number(row, "high") for row in valid))
                ax.plot(range(len(y)), y, color=COLORS[method], linestyle=LINES[m],
                        marker=MARKERS[m], ms=5, lw=1.7, label=LABELS[method])
            ax.set(title=f"{'Independent' if i == 0 else 'Spatially correlated'} sites · {effect}",
                   ylim=(0, 1.025), xticks=range(5), xticklabels=ORDER)
            ax.yaxis.set_major_formatter(PercentFormatter(1))
            if j == 0:
                ax.set_ylabel("Generating-parameter coverage")
            if i == 1:
                ax.set_xlabel("Order as a share of the sites")
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .025), ncol=3,
               frameon=False, fontsize=10)
    fig.suptitle("Coverage of nominal 95% intervals by order", fontsize=16, weight="bold", y=.974)
    fig.savefig(dest / "coverage.png")
    plt.close(fig)


def precision_plot(data, dest):
    palette = ["#264b80", "#b94d29", "#187464", "#725290"]
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.7), sharex=True, sharey=True, layout="constrained")
    for ax, changes in zip(axs, (False, True)):
        for m, name in enumerate(FAMILIES):
            pairs = grouped([row for row in data if row["confirmation"] == "FALSE" and
                             family(row) == name and (" - " in row["contrast"]) == changes], lambda row: int(row["budget"]))
            ax.plot([key[0] for key, _ in pairs], [med(number(row, "ratio") for row in rows) for _, rows in pairs],
                    label=name.capitalize(), color=palette[m], marker=MARKERS[m], linestyle=LINES[m], lw=1.8)
        ax.axhline(.2, color=".35", lw=1.1, label="20% precision target")
        ax.set(xscale="log", yscale="log", title="Adjacent-order changes" if changes else "Single-order effects",
               xlabel="Combinations per independent fit")
        ax.set_xticks([500, 2000, 10000, 50000], ["500", "2,000", "10,000", "50,000"])
        ax.yaxis.set_major_locator(FixedLocator([.025, .05, .10, .20, .40, .80]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{100 * value:g}%"))
        ax.yaxis.set_minor_formatter(FuncFormatter(lambda value, _: ""))
    axs[0].set_ylabel("Median Monte Carlo SD / delete-one SD")
    axs[1].legend(loc="best", frameon=False, fontsize=9)
    fig.suptitle("Larger samples reduce numerical variation across orders", fontsize=15, weight="bold")
    fig.savefig(dest / "precision.png")
    plt.close(fig)


def ecology_plot(data, dest):
    fig, axs = plt.subplots(1, 3, figsize=(12.9, 4.1), layout="constrained", sharex=True)
    for m, signal in enumerate(("null", "weak", "strong")):
        bins = grouped([r for r in data if r["signal"] == signal], order_label)
        position = [key[0] for key, _ in bins]
        sharing = [med(number(r, "probability_sharing") for r in rows) for _, rows in bins]
        usable = [average(r["converged"] == "TRUE" and r["rank"] == r["dimension"] and
                          math.isfinite(number(r, "D2")) for r in rows) for _, rows in bins]
        def ranking(rows):
            ranks = [number(r, "importance.gradient") > number(r, "importance.control")
                     for r in rows if math.isfinite(number(r, "importance.gradient"))]
            return mean(ranks) if len(ranks) >= 5 else math.nan
        rank = [ranking(rows) for _, rows in bins]
        for ax, y in zip(axs, (sharing, usable, rank)):
            ax.plot(position, y, color=list(COLORS.values())[m], marker=MARKERS[m],
                    linestyle=LINES[m], lw=1.8, label=signal.capitalize())
            ax.set(ylim=(-.01, 1.03), xticks=range(5), xticklabels=ORDER)
            ax.yaxis.set_major_formatter(PercentFormatter(1))
    for ax, title in zip(axs, ("Base probability of sharing", "Identified, converged fits", "Gradient ranked first\n(at least five estimable fits)")):
        ax.set_title(title, fontsize=11, weight="bold")
    axs[2].legend(loc="best", frameon=False, fontsize=9)
    fig.supxlabel("Order as a share of the sites")
    fig.suptitle("Sharing and inferential diagnostics across orders", fontsize=15, weight="bold")
    fig.savefig(dest / "ecological.png")
    plt.close(fig)


def pilot_plot(data, dest):
    fig, ax = plt.subplots(figsize=(8.8, 4.9), layout="constrained")
    groups = grouped(data, "scenario", lambda row: int(row["order"]))
    for j, ((scenario, order), rows) in enumerate(groups):
        baseline = next(number(r, "cost") for r in rows if number(r, "alpha") == 0)
        if not math.isfinite(baseline) or baseline <= 0:
            continue
        xs = [number(r, "alpha") for r in rows]
        ys = [number(r, "cost") / baseline for r in rows]
        selected = number(rows[0], "selected")
        m = (int(scenario) - 1) % 3
        color = list(COLORS.values())[m]
        ax.plot(xs, ys, color=color, linestyle=LINES[j % 3], marker=MARKERS[j % 3], lw=1.25,
                label=f"Scenario {scenario} · order {order}", alpha=.78)
        for x, y in zip(xs, ys):
            if x == selected:
                ax.scatter([x], [y], color=color, edgecolor="white", s=90, zorder=5, linewidth=.8)
    ax.axhline(1, color=".37", lw=1)
    ax.set(yscale="log", xlabel="Targeted mixture fraction",
           ylabel="Validation variance × CPU time / conventional",
           xticks=(0, .1, .25, .5, .75, .9))
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9)
    fig.suptitle("Pilot choices and fresh-fit efficiency", fontsize=15, weight="bold")
    fig.savefig(dest / "pilot.png")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    configure()
    calibration = load(args.raw, "calibration")
    changes = load(args.raw, "changes")
    precision = load(args.raw, "precision")
    ecology = load(args.raw, "ecological")
    export(args.output, "coverage", evaluate(calibration))
    export(args.output, "changes", evaluate(changes, change=True))
    diagnostics = load(args.raw, "precision_diagnostics")
    rows = []
    for (contrast, effect, budget), records in grouped(
        [r for r in precision if r["confirmation"] == "FALSE"],
        lambda r: "adjacent" if " - " in r["contrast"] else "single", family, lambda r: int(r["budget"])):
        ratios = finite(number(r, "ratio") for r in records)
        rows.append(dict(contrast=contrast, effect=effect, budget=budget, evaluated=len(records),
             valid=len(ratios), median_ratio=med(ratios), max_ratio=max(ratios, default=math.nan),
             passed=sum(r["pass"] == "TRUE" for r in records)))
    export(args.output, "precision", rows)
    export(args.output, "diagnostics", [dict(
        scenario=s, order=k, budget=m, alpha=number(rows[0], "alpha"),
        sharing=med(number(r, "sharing_estimate") for r in rows),
        proposal_sharing=med(number(r, "shared") for r in rows),
        expected_shared_species=number(rows[0], "expected_shared_species"),
        sharing_lower=number(rows[0], "sharing_lower"), sharing_upper=number(rows[0], "sharing_upper"),
        ess=med(number(r, "ess") for r in rows), min_delete_ess=number(rows[0], "min_delete_ess"),
        condition=med(number(r, "condition") for r in rows),
        convergence=average(r["converged"] == "TRUE" for r in rows),
        active=average(number(r, "active") > 0 for r in rows),
        full_rank=average(r["rank"] == r["parameters"] for r in rows))
        for (s, k, m), rows in grouped(diagnostics, "scenario", lambda r: int(r["order"]), lambda r: int(r["budget"]))])
    export(args.output, "ecology", [dict(signal=signal, stage=j, fits=len(rows),
        median_sharing=med(number(r, "probability_sharing") for r in rows),
        median_ess=med(number(r, "ess") for r in rows),
        converged_full_rank=average(r["converged"] == "TRUE" and r["rank"] == r["dimension"] and
                                    math.isfinite(number(r, "D2")) for r in rows),
        rankable=sum(math.isfinite(number(r, "importance.gradient")) for r in rows),
        gradient_top=average(number(r, "importance.gradient") > number(r, "importance.control")
                             for r in rows if math.isfinite(number(r, "importance.gradient"))),
        median_D2=med(number(r, "D2") for r in rows),
        median_test_D2=med(number(r, "D2_test") for r in rows),
        median_test_RMSE=med(number(r, "RMSE_test") for r in rows))
        for (signal, j), rows in grouped(ecology, "signal", order_label)])
    export(args.output, "permutation", [dict(signal=signal, stage=j, attempted=len(rows),
        valid=sum(math.isfinite(number(r, "p")) for r in rows),
        rate=average(number(r, "p") <= .05 for r in rows if math.isfinite(number(r, "p"))))
        for (signal, j), rows in grouped(load(args.raw, "permutation"), "signal", order_label)])
    raw_finite = load(args.raw, "finite_reference")
    export(args.output, "finite_reference", [dict(stage=j, effect=effect, attempted=len(rows),
        valid=len([v for v in errors if math.isfinite(v)]), mean_bias=average(errors),
        median_absolute_error=med(map(abs, errors)))
        for (j, effect), rows in grouped(raw_finite, order_label, family)
        for errors in [[number(r, "sampled") - number(r, "reference") for r in rows]]])
    identification = load(args.raw, "identification")
    export(args.output, "identification", [dict(n=n, stage=j, fits=len(rows),
        convergence=average(r["converged"] == "TRUE" for r in rows),
        full_rank=average(r["rank"] == r["dimension"] for r in rows),
        condition=med(number(r, "condition") for r in rows),
        min_delete_ess=med(number(r, "min_delete_ess") for r in rows),
        constrained_active=average(number(r, "constrained_active") > 0 for r in rows))
        for (n, j), rows in grouped(identification, lambda r: int(r["n"]), order_label)])
    export(args.output, "additive", [dict(order=k, experiments=len(rows),
        bb_median_se=med(number(r, "bb_se") for r in rows),
        bb_coverage=average(r["bb_covered"] == "TRUE" for r in rows),
        jackknife_median_se=med(number(r, "jack_se") for r in rows),
        jackknife_coverage=average(r["jack_covered"] == "TRUE" for r in rows if r["jack_covered"] != "NA"),
        true_se=number(rows[0], "exact_sampling_se"))
        for (k,), rows in grouped(load(args.raw, "additive"), lambda r: int(r["order"]))])
    sensitivity = load(args.raw, "sensitivity")
    export(args.output, "sensitivity", [dict(n=n, occupancy=occ, breadth=b, spatial=spatial,
        order=k, fits=len(rows), sharing=med(number(r, "sharing") for r in rows),
        identified=average(r["identified"] == "TRUE" for r in rows),
        importance_rank=average(r["gradient_first"] == "TRUE" for r in rows if r["gradient_first"] != "NA"))
        for (n, occ, b, spatial, k), rows in grouped(sensitivity, "n", "occupancy", "breadth", "spatial", "order")])
    for name in ("pilot", "pilot_validation", "thresholds", "focused"):
        export(args.output, name, load(args.raw, name))
    calibration_plot(calibration, args.output)
    precision_plot(precision, args.output)
    ecology_plot(ecology, args.output)
    pilot_plot(load(args.raw, "pilot_validation"), args.output)


if __name__ == "__main__":
    main()
