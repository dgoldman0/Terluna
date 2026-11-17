"""The figures of the vacuum economy's scoping study and of the proposal for the fund before the build, and both built
as PDFs.

Reads the period model's scoping product (research/studies/vacuum_economy_periods/results/scoping.json) and draws its
figures; then builds the study's text (scoping_study.md) and the fund proposal's (rd_fund.md), both in
research/studies/vacuum_economy_periods/, with them into PDFs through pandoc and xelatex. Calculates nothing: every
plotted value is read from the product and listed in the manifest with the source, renderer, header and text hashes.

usage (repository root):
  python visualization/vacuum-economy-periods/report.py [--out DIR] [--figures-only]

Writes DIR/figures/*.pdf, DIR/manifest.json, DIR/scoping_study.pdf and DIR/rd_fund.pdf; DIR defaults to the ignored
research/runs/vacuum_economy_periods/report/.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'research/studies/vacuum_economy_periods'
SOURCE = STUDY / 'results/scoping.json'
TEXTS = (STUDY / 'scoping_study.md', STUDY / 'rd_fund.md')
HEADER = Path(__file__).resolve().parent / 'header.tex'
OUT = ROOT / 'research/runs/vacuum_economy_periods/report'

# The validated categorical order, ink and chrome (light print surface).
S = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
INK, INK2, MUTED, GRID, AXIS, SURFACE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#ffffff'
YEARS = list(range(2041, 2342))
W = 6.3          # inches, the text width of the page


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def style():
    plt.rcParams.update({
        'font.family': 'Inter', 'font.size': 8.5, 'axes.titlesize': 9, 'axes.labelsize': 8.5, 'xtick.labelsize': 8,
        'ytick.labelsize': 8, 'legend.fontsize': 8, 'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': MUTED,
        'ytick.color': MUTED, 'axes.edgecolor': AXIS, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.5, 'axes.axisbelow': True, 'lines.linewidth': 1.6,
        'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round', 'legend.frameon': False, 'figure.facecolor': SURFACE,
        'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE, 'pdf.fonttype': 42})


def millions(x, _):
    return ('%.0f million' % (x / 1e6)) if x >= 1e8 else ('%.1f million' % (x / 1e6)) if x >= 1e6 else ('{:,.0f}'.format(x))


def save(fig, out, name, drawn, values):
    fig.savefig(out / 'figures' / ('%s.pdf' % name), bbox_inches='tight')
    plt.close(fig)
    drawn[name] = values


def end_label(ax, x, y, text, color, dy=0.0):
    ax.annotate(text, (x, y), xytext=(4, dy), textcoords='offset points', color=INK2, fontsize=8, va='center')
    ax.plot([x], [y], 'o', ms=4, color=color, mec=SURFACE, mew=1.2, zorder=5)


def figures(p, out):
    drawn = {}
    base, dollar, ramp = p['paths']['base']['series'], p['paths']['dollar']['series'], p['paths']['ramp']['series']

    # 1. Residents.
    fig, ax = plt.subplots(figsize=(W, 2.9))
    for s, label, color, dash in ((ramp, "The programme's ramp, carried on", S[2], (4, 2)), (dollar, 'Growth as a result, dollars alone', S[1], ()),
                                  (base, 'Growth as a result, U (base case)', S[0], ())):
        n = [x if x > 0 else None for x in s['residents']]
        ax.plot(YEARS, n, color=color, dashes=dash, label=label)
    ax.set_yscale('log')
    ax.set_ylim(1e4, 3e8)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: {1e4: '10,000', 1e5: '100,000', 1e6: '1 million', 1e7: '10 million', 1e8: '100 million'}.get(x, '')))
    ax.set_xlim(2041, 2341)
    ax.axhline(6.47e6, color=AXIS, lw=0.6)
    ax.annotate('The plan, 6.47 million', (2045, 6.47e6), xytext=(0, 3), textcoords='offset points', color=MUTED, fontsize=7.5)
    for s, color, dy in ((base, S[0], 5), (dollar, S[1], -5), (ramp, S[2], 0)):
        end_label(ax, 2341, s['residents'][-1], millions(s['residents'][-1], None), color, dy)
    ax.legend(loc='upper left')
    ax.set_ylabel('Residents')
    save(fig, out, 'residents', drawn, dict(base=base['residents'][::10], dollar=dollar['residents'][::10], ramp=ramp['residents'][::10]))

    # 2. The funders' position.
    fig, ax = plt.subplots(figsize=(W, 2.6))
    for s, label, color in ((dollar, 'Dollars alone', S[1]), (base, 'U (base case)', S[0])):
        ax.plot(YEARS, [x / 1e12 for x in s['funders_cumulative']], color=color, label=label)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xlim(2041, 2341)
    ax.set_ylabel("Funders' cumulative position,\ntrillions of 2041 dollars")
    ax.legend(loc='upper left')
    save(fig, out, 'funders', drawn, dict(base=[x / 1e12 for x in base['funders_cumulative'][::10]], dollar=[x / 1e12 for x in dollar['funders_cumulative'][::10]]))

    # 3. U's rate and the reserve.
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(W, 3.6), sharex=True, gridspec_kw=dict(height_ratios=(1, 1.3)))
    a1.plot(YEARS, base['u_rate'], color=S[0])
    a1.axhline(0.75, color=S[7], lw=0.8, dashes=(3, 2))
    a1.annotate('Floor, 75% of par', (2045, 0.75), xytext=(0, 3), textcoords='offset points', color=INK2, fontsize=7.5)
    a1.set_ylim(0.6, 1.05)
    a1.set_ylabel("U's rate,\ndollars a U")
    a2.plot(YEARS, [x / 1e12 for x in base['u_reserve']], color=S[0], label='Reserve')
    a2.plot(YEARS, [x / 1e12 for x in base['u_run_need']], color=S[1], label='What the run test asks')
    a2.set_ylabel('Trillions of\n2041 dollars')
    a2.set_xlim(2041, 2341)
    a2.legend(loc='upper left')
    save(fig, out, 'rate_reserve', drawn, dict(rate=base['u_rate'][::10], reserve=[x / 1e12 for x in base['u_reserve'][::10]],
                                               run_need=[x / 1e12 for x in base['u_run_need'][::10]]))

    # 4. Payments to Earth and receipts from Earth, by kind.
    pay_labels = [('residents_purchases', "Residents' purchases from Earth"), ('partners_and_services', "Partners' and services' shares of sales"),
                  ('supplies', 'Supplies'), ('equipment_plant_power_tugs', 'Equipment, plant, power and tugs'),
                  ('fares_teleoperation_and_other', 'Fares, teleoperation and other')]
    rec_groups = [('computing', ('computing',), 'Computing'), ('visits', ('visits',), 'Visits'), ('culture', ('culture',), 'Sport and entertainment'),
                  ('institutions', ('institutions',), 'Hosted research and professional work'),
                  ('depot', ('propellant', 'materials', 'metals'), 'Propellant, material and metals')]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, 3.7), sharey=True)
    ys = [[x / 1e12 for x in base['payments'][k]] for k, _ in pay_labels]
    a1.stackplot(YEARS, *ys, colors=S[:5], labels=[l for _, l in pay_labels], edgecolor=SURFACE, linewidth=0.4)
    a1.set_title('Paid to Earth', loc='left', color=INK)
    rs = [[sum(base['receipts_by_line'][k][i] for k in keys) / 1e12 for i in range(len(YEARS))] for _, keys, _ in rec_groups]
    a2.stackplot(YEARS, *rs, colors=S[:5], labels=[l for _, _, l in rec_groups], edgecolor=SURFACE, linewidth=0.4)
    a2.set_title("Received from Earth's customers", loc='left', color=INK)
    a1.set_ylabel('Trillions of 2041 dollars a year')
    for a in (a1, a2):
        a.set_xlim(2041, 2341)
        handles, labels = a.get_legend_handles_labels()     # listed from the top band down, as the stack reads
        a.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0, -0.11), fontsize=7, handlelength=1.2, borderaxespad=0)
    save(fig, out, 'flows', drawn, dict(paid={k: [x / 1e12 for x in base['payments'][k][::10]] for k, _ in pay_labels},
                                        received={k: [r[i] for i in range(0, len(YEARS), 10)] for (k, _, _), r in zip(rec_groups, rs)}))

    # 5. Sensitivities, as the change in residents by 2341.
    base_n = p['paths']['base']['measures']['residents_2341']

    def tornado(rows, name_key, name):
        fig, ax = plt.subplots(figsize=(W, 0.32 * len(rows) + 0.7))
        labels, values = [], {}
        for i, row in enumerate(reversed(rows)):
            ends = row['ends']
            same_side = len(ends) == 2 and (ends[0]['change'] >= 0) == (ends[1]['change'] >= 0)
            for j, e in enumerate(ends):
                change = e['change'] / 1e6
                color = S[0] if change >= 0 else S[7]
                y, h = (i + (0.16 if j == 0 else -0.16), 0.28) if same_side else (i, 0.55)
                ax.barh(y, change, height=h, color=color, edgecolor=SURFACE, linewidth=0.6)
                ax.annotate(e['setting'], (change, y), xytext=(3 if change >= 0 else -3, 0), textcoords='offset points',
                            ha='left' if change >= 0 else 'right', va='center', fontsize=7, color=INK2)
                values['%s: %s' % (row[name_key], e['setting'])] = change
            labels.append(row[name_key])
        ax.set_yticks(range(len(rows)), labels, color=INK2)
        ax.axvline(0, color=AXIS, lw=0.8)
        ax.grid(axis='y', visible=False)
        ax.set_xlabel('Change in residents by 2341 from the base case\'s %.1f million, millions' % (base_n / 1e6))
        lo = min(min(values.values()), 0)
        hi = max(max(values.values()), 0)
        ax.set_xlim(lo - 3.5, hi + 3.5)
        save(fig, out, name, drawn, values)

    tornado(p['sensitivities'], 'input', 'sensitivities')

    # 6. Residents by 2341 across demand and technology.
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.4), sharey=True)
    grid_values = {}
    for ax, scenario in zip(axes, ('low', 'central', 'high')):
        for k, (tcase, x) in enumerate((('slow', 0), ('central', 1), ('fast', 2))):
            for money, dx, color in ((':U', -0.17, S[0]), (':dollar', 0.17, S[1])):
                g = p['grid']['%s:%s%s' % (scenario, tcase, money)]
                v = (g['residents_2341'] or 0) / 1e6
                ax.bar(x + dx, v, width=0.32, color=color, edgecolor=SURFACE, linewidth=0.6, label=('U' if money == ':U' else 'Dollars alone') if k == 0 else None)
                if v < 0.1:
                    ax.annotate('none', (x + dx, 0), xytext=(0, 2), textcoords='offset points', ha='center', va='bottom', fontsize=6.5, color=MUTED, rotation=90)
                grid_values['%s:%s%s' % (scenario, tcase, money)] = v
        ax.set_xticks(range(3), ('slow', 'central', 'fast'))
        ax.set_title('%s demand' % scenario.capitalize(), loc='left', color=INK)
        ax.grid(axis='x', visible=False)
    axes[0].set_ylabel('Residents by 2341, millions')
    axes[0].legend(loc='upper left')
    save(fig, out, 'grid', drawn, grid_values)

    # 7. Earth's holdings of U by group.
    hold = [('working_balances', 'Working balances'), ('holders_for_use', 'Holders for use'), ('funders', 'The funders'),
            ('bonds', "Investors' bonds")]
    fig, ax = plt.subplots(figsize=(W, 2.5))
    ax.stackplot(YEARS, *[[x / 1e12 for x in base['holdings_USD'][k]] for k, _ in hold], colors=S[:5], labels=[l for _, l in hold],
                 edgecolor=SURFACE, linewidth=0.4)
    ax.set_xlim(2041, 2341)
    ax.set_ylabel('Trillions of 2041 dollars')
    ax.legend(loc='upper left')
    save(fig, out, 'holdings', drawn, {k: [x / 1e12 for x in base['holdings_USD'][k][::10]] for k, _ in hold})
    return drawn


def front(text):
    """The simple `key: "value"` lines of a document's front matter, which set its running head (runninghead), the
    heading of its title-page note (abouttitle) and whether it has a page of contents (contents: false)."""
    lines = text.read_text().split('\n')
    end = lines.index('---', 1)
    out = {}
    for line in lines[1:end]:
        if ':' in line and not line.startswith(' '):
            key, value = line.split(':', 1)
            out[key.strip()] = value.strip().strip('"')
    return out


def build(out, text):
    """One document as a PDF: the scoping study, or the shorter proposal for the fund that precedes the build. Its own
    running head and note heading are set in a small header file after the shared one, since pandoc's command-line
    header replaces any header-includes a document carries."""
    name = text.stem
    shutil.copy(text, out / text.name)
    meta = front(text)
    own = out / ('%s_header.tex' % name)
    own.write_text(''.join('\\def\\%s{%s}\n' % (key, meta[key]) for key in ('runninghead', 'abouttitle') if key in meta))
    contents = ['--toc', '--toc-depth=2'] if meta.get('contents', 'true') != 'false' else []
    cmd = ['pandoc', str(out / text.name), '-o', str(out / ('%s.pdf' % name)), '--pdf-engine=xelatex', '-H', str(HEADER), '-H', str(own),
           '--resource-path', str(out), '--number-sections'] + contents + [
           '-V', 'papersize=a4', '-V', 'geometry:margin=24mm', '-V', 'fontsize=10.5pt', '-V', 'linkcolor=[HTML]{184F95}',
           '-V', 'urlcolor=[HTML]{184F95}', '-V', 'toccolor=[HTML]{0B0B0B}']
    subprocess.run(cmd, check=True, cwd=out)
    return out / ('%s.pdf' % name)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--out', default=str(OUT))
    parser.add_argument('--figures-only', action='store_true')
    args = parser.parse_args()
    out = Path(args.out)
    (out / 'figures').mkdir(parents=True, exist_ok=True)
    p = json.loads(SOURCE.read_text())
    if p['schema'] != 'terluna.research.vacuum-economy-periods.scoping/1':
        raise ValueError('unexpected product %s' % p['schema'])
    style()
    drawn = figures(p, out)
    texts = [t for t in TEXTS if t.exists()]
    manifest = dict(source=str(SOURCE.relative_to(ROOT)), source_sha256=digest(SOURCE), renderer_sha256=digest(__file__),
                    header_sha256=digest(HEADER), texts={str(t.relative_to(ROOT)): digest(t) for t in texts}, figures=drawn)
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    if not args.figures_only:
        for t in texts:
            print(build(out, t))
    else:
        print(out / 'figures')


if __name__ == '__main__':
    main()
