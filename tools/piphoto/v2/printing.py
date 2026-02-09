import math

from FFSlib import NoResult, nuclei


def _format_parens(val):
    """Format np.array([mean, unc]) as parenthetical notation, e.g. '1.551(78)'.
    Returns '---' for NoResult."""
    if isinstance(val, NoResult):
        return "---"
    mean, unc = val
    if unc <= 0:
        return f"{mean:.3f}"
    mag = math.floor(math.log10(abs(unc)))
    n_dec = -(mag - 1)  # decimal places for 2 sig figs of unc
    unc_r = round(unc, n_dec)
    unc_int = int(round(unc_r * 10**n_dec))
    # Trim trailing zero (e.g. 20 -> use 1 sig fig instead)
    if unc_int % 10 == 0 and unc_int >= 10:
        n_dec -= 1
        unc_r = round(unc, n_dec)
        unc_int = int(round(unc_r * 10**n_dec))
    n_dec = max(0, n_dec)
    mean_r = round(mean, n_dec)
    mean_str = f"{mean_r:.{n_dec}f}"
    return f"{mean_str}({unc_int})"


def _print_nucleus_table(d, nuc, title, rows, authors):
    """Print a table for one nucleus.

    rows: list of quantity key strings (one per row)
    authors: list of (label, author_key) — columns
    """
    w = 18
    row_w = max(16, max(len(q) for q in rows) + 2)

    col_labels = [label for label, _ in authors]
    total_width = row_w + len(authors) * w

    print(title)
    print("=" * total_width)

    # Column headers
    header = f"{'':<{row_w}}"
    for label in col_labels:
        header += f"{label:>{w}}"
    print(header)
    print("-" * total_width)

    # Data rows
    for q in rows:
        row = f"{q:<{row_w}}"
        for _, author in authors:
            row += f"{_format_parens(d[author][nuc][q]):>{w}}"
        print(row)

    print("-" * total_width)


def printLit(d=None):
    """Print literature results per nucleus, with rows=quantities and columns=authors."""
    if d is None:
        from FFSlib import get_literature_results

        d = get_literature_results()

    authors = [("TDA", "TDA"), ("Lenkewitz", "lenke"), ("Braun", "braun")]

    onebody_rows = ["F_T^S+V", "F_T^S-V", "F_L^S+V", "F_L^S-V", "E1N", "L1N"]

    twobody_rows = [
        "F^{(a)}_T-F^{(b)}_T O(q^3)",
        "F^{(a)}_L-F^{(b)}_L O(q^3)",
        "E2N O(q^3)",
        "L2N O(q^3)",
        # "Static F^{(a)}_T-F^{(b)}_T",
        # "Static F^{(a)}_L-F^{(b)}_L",
        # "Static E2N",
        # "Static L2N",
        "F^{(a)}_T-F^{(b)}_T O(q^4)",
        "F^{(a)}_L-F^{(b)}_L O(q^4)",
        "E2N O(q^4)",
        "L2N O(q^4)",
    ]

    for nuc in nuclei:
        print(f"{nuc}")
        _print_nucleus_table(d, nuc, "One body", onebody_rows, authors)
        print()
        _print_nucleus_table(d, nuc, "Two body", twobody_rows, authors)
        print("\n")


if __name__ == "__main__":
    printLit()
