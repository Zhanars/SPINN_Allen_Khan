"""Assemble the revised manuscript from manuscript/r1_src/*.tex and compile it.
usage: python build_r1.py [hl|clean|both]
 * inserts the highlight switch, flattens \\input of the number/abstract files,
 * highlighted version: everything new or changed with respect to the first version (v1_submitted/manuscript.tex) is
   blue; paragraphs, headings, title-page data and references carried over verbatim are black,
 * re-orders the bibliography by order of first citation (MDPI requirement),
 * checks for undefined labels/citations and unresolved numbers (??),
 * writes manuscript/manuscript_r1_highlighted.tex / manuscript_r1_clean.tex and compiles them (pdflatex x2)."""
import os
import re
import sys
import subprocess

BS = chr(92)
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "r1_src")


def read(name):
    return open(os.path.join(SRC, name), encoding="utf-8").read()


def flatten_inputs(text):
    def rep(m):
        fn = m.group(1)
        p = os.path.join(SRC, fn)
        return open(p, encoding="utf-8").read() if os.path.exists(p) else m.group(0)
    for _ in range(5):                      # nested inputs (appendix -> 73_app_orig -> tab_orig)
        new = re.sub(re.escape(BS) + r"input\{([^}]+)\}", rep, text)
        if new == text:
            break
        text = new
    return text


def reorder_bib(body, refs_text):
    head, rest = refs_text.split(BS + "begin{thebibliography}", 1)
    rest_head, rest_body = rest.split("\n", 1)
    rest_body, tail = rest_body.split(BS + "end{thebibliography}", 1)
    items = re.split(r"(?m)^(?=" + re.escape(BS) + r"bibitem\{)", rest_body)
    items = [i for i in items if i.strip().startswith(BS + "bibitem")]
    bykey = {}
    for it in items:
        k = re.match(re.escape(BS) + r"bibitem\{([^}]+)\}", it.strip()).group(1)
        bykey[k] = it.strip() + "\n"
    cited = []
    for m in re.finditer(re.escape(BS) + r"cite[a-z]*\{([^}]+)\}", body):
        for k in m.group(1).split(","):
            k = k.strip()
            if k not in cited:
                cited.append(k)
    missing = [k for k in cited if k not in bykey]
    unused = [k for k in bykey if k not in cited]
    ordered = [bykey[k] for k in cited if k in bykey]
    new = head + BS + "begin{thebibliography}" + rest_head + "\n" + "".join(ordered) + BS + "end{thebibliography}" + tail
    return new, missing, unused


V1 = os.path.join(HERE, "v1_submitted", "manuscript.tex")
HL_PREAMBLE = "\n".join([
    "% highlighted version: every passage that is new or changed with respect to the first version is blue;",
    "% paragraphs, headings, title-page data and references carried over verbatim stay black (build_r1.py)",
    BS + "renewcommand{" + BS + "rev}[1]{#1}   % blue is the default; a colour group at a paragraph start could allow a break after a heading",
    BS + "DeclareCaptionFont{hlblue}{" + BS + "color{blue}}",
    BS + "AtBeginEnvironment{tabular}{" + BS + "color{blue}}",
    BS + "AtBeginEnvironment{tabularx}{" + BS + "color{blue}}",
    BS + "AtBeginDocument{",
    "  " + BS + "captionsetup[table]{labelfont+=hlblue,textfont+=hlblue}",
    "  " + BS + "captionsetup[figure]{labelfont+=hlblue,textfont+=hlblue}}",
    ""])
HEADING = re.compile(re.escape(BS) + r"(?:part|section|subsection|subsubsection|paragraph|reftitle)\*?[\[{]")
VERBATIM = ("lstlisting", "verbatim")


FRONT = ("Title", "Author", "address", "corres", "conference", "abstract", "keyword", "MSC")
SIMPLE_HEADING = re.compile(re.escape(BS) + r"((?:sub)*section)\{([^{}]*)\}")


def _norm(t):
    t = re.sub(r"(?<!" + re.escape(BS) + r")%.*", "", t)
    t = t.replace(BS + " ", " ")          # a control space ('Y.Zh.\ and') changes the spacing, not the text
    return " ".join(t.split())


def _arg(s, cmd):
    """(start, end) of the braced argument of the first \\cmd{...} in s (balanced braces), or None."""
    m = re.search(re.escape(BS + cmd) + r"\{", s)
    if not m:
        return None
    i, depth = m.end(), 1
    while depth:
        c = s[i]
        if c == BS:
            i += 2
            continue
        depth += (c == "{") - (c == "}")
        i += 1
    return m.end(), i - 1


def highlight_front(pre):
    """Title-page data: blue only where the content differs from the first version (the template stays black).
    \\texorpdfstring keeps the colour out of the PDF metadata (title, subject, keywords)."""
    v1 = open(V1, encoding="utf-8").read()
    changed = []
    for cmd in FRONT:
        a, b = _arg(pre, cmd), _arg(v1, cmd)
        if a is None:
            continue
        if b is None or _norm(pre[a[0]:a[1]]) != _norm(v1[b[0]:b[1]]):
            # a group, so that the colour does not leak into the template text printed after the item
            pre = (pre[:a[0]] + BS + "texorpdfstring{" + BS + "begingroup" + BS + "color{blue}}{}" + pre[a[0]:a[1]]
                   + BS + "texorpdfstring{" + BS + "endgroup}{}" + pre[a[1]:])
            changed.append(cmd)
    return pre, changed


def _v1_body():
    s = open(V1, encoding="utf-8").read()
    return s[s.index(BS + "begin{document}") + len(BS + "begin{document}"):s.index(BS + "end{document}")]


def _env_delta(block):
    b = len(re.findall(re.escape(BS) + r"begin\{(?!document)", block))
    e = len(re.findall(re.escape(BS) + r"end\{(?!document)", block))
    return b - e


def highlight_body(body):
    """Blue is the default colour of the highlighted version.  Blocks (paragraphs separated by blank lines) whose text
    occurs verbatim in the first version are switched to black, provided they are outside every environment, close
    every environment they open, contain no verbatim material and are not headings (a colour switch right after a
    heading could allow a page break there).  Everything else stays blue, so a change is never shown in black."""
    v1 = _norm(_v1_body())
    # every heading starts its own block (a blank line before \section is harmless: the heading ends the paragraph)
    body = re.sub(r"(?m)^(?=" + re.escape(BS) + r"(?:sub)*section\{)", "\n", body)
    parts = re.split(r"(\n\s*\n)", body)
    depth, verb, cur, out, n_black, prev_heading = 0, False, "blue", [BS + "color{blue}%\n"], 0, False
    for i in range(0, len(parts), 2):
        blk = parts[i]
        sep = parts[i + 1] if i + 1 < len(parts) else ""
        nb = _norm(blk)
        is_heading = bool(HEADING.match(nb))
        if is_heading:
            # a heading carried over from the first version: colour its text black inside the heading command
            # (a colour switch between the heading and the next paragraph could allow a page break there)
            first_line = next((ln for ln in blk.split("\n") if ln.strip() and not ln.lstrip().startswith("%")), "")
            m = SIMPLE_HEADING.match(first_line.lstrip())
            if m and _norm(m.group(0)) in v1:
                t = m.group(2)
                new = BS + m.group(1) + "{" + BS + "texorpdfstring{" + BS + "protect" + BS + "textcolor{black}{" + t + "}}{" + t + "}}"
                blk = blk.replace(m.group(0), new, 1)
                n_black += 1
        eligible = (nb and depth == 0 and not verb and not is_heading and _env_delta(blk) == 0
                    and not any(BS + "begin{" + v in blk for v in VERBATIM) and nb in v1)
        want = "black" if eligible else "blue"
        if nb and want != cur and depth == 0 and not verb:
            lead = blk[:len(blk) - len(blk.lstrip())]
            starts_text = not nb.startswith(BS) or nb.startswith(BS + "emph") or nb.startswith(BS + "textbf")
            if starts_text:
                blk = lead + BS + "leavevmode" + BS + "color{" + want + "}" + blk.lstrip()
                cur = want
            elif not prev_heading:
                blk = lead + BS + "color{" + want + "}" + blk.lstrip()
                cur = want
        if want == "black" and cur == "black":
            n_black += 1
        for v in VERBATIM:
            if BS + "begin{" + v in blk:
                verb = True
            if BS + "end{" + v in blk:
                verb = False
        depth += _env_delta(blk)
        if nb:
            prev_heading = is_heading
        out.append(blk + sep)
    if cur != "blue":
        out.append(BS + "color{blue}\n")
    return "".join(out), n_black


def highlight_bib(refs2):
    v1 = _v1_body()
    old = {}
    for it in re.split(r"(?m)^(?=" + re.escape(BS) + r"bibitem\{)", v1):
        m = re.match(re.escape(BS) + r"bibitem\{([^}]+)\}", it.strip())
        if m:
            old[m.group(1)] = _norm(it.split(BS + "end{thebibliography}")[0])
    n, first = 0, []

    def rep(m):
        # the colour is switched before \bibitem (at the end of the previous entry's paragraph), so that the
        # reference number has the colour of its own entry; the first entry's colour is set before the list
        nonlocal n
        key, text = m.group(1), m.group(2)
        same = old.get(key) == _norm(BS + "bibitem{" + key + "}" + text)
        n += same
        col = "black" if same else "blue"
        if not first:
            first.append(col)
            return BS + "bibitem{" + key + "}" + text
        return BS + "ifhmode" + BS + "unskip" + BS + "fi" + BS + "color{" + col + "}" + BS + "bibitem{" + key + "}" + text
    out = re.sub(r"(?ms)^" + re.escape(BS) + r"bibitem\{([^}]+)\}(.*?)(?=^" + re.escape(BS) + r"bibitem\{|^" + re.escape(BS) + r"end\{thebibliography\})", rep, refs2)
    if first:     # the heading 'References' (unchanged) and the first entry
        out = out.replace(BS + "reftitle{", BS + "color{" + first[0] + "}" + BS + "reftitle{", 1)
    return out, n


def build(mode):
    pre = read("00_preamble.tex").replace(BS + "HLSWITCH", BS + ("hltrue" if mode == "hl" else "hlfalse"))
    body = "".join(read(f) for f in ["10_intro.tex", "20_methods.tex", "30_results.tex", "40_discussion.tex", "50_conclusions.tex", "60_backmatter.tex", "70_appendix.tex"] if os.path.exists(os.path.join(SRC, f)))
    refs = read("80_refs.tex")
    full_body = flatten_inputs(body)
    refs2, missing, unused = reorder_bib(full_body, refs)
    pre = flatten_inputs(pre)
    if mode == "hl":
        # blue is the default now; a manual switch back to black would hide the changes that follow it
        full_body = full_body.replace(BS + "ifhl" + BS + "color{black}" + BS + "fi", "")
        full_body, n_black = highlight_body(full_body)
        refs2, n_bib = highlight_bib(refs2)
        pre, front = highlight_front(pre)
        pre = pre + "\n" + HL_PREAMBLE
        print(f"[hl] unchanged blocks/headings shown in black: {n_black}; unchanged references shown in black: {n_bib}; changed title-page items: {front}")
    doc = pre + "\n" + BS + "begin{document}\n\n" + full_body + "\n" + refs2
    name = "manuscript_r1_highlighted" if mode == "hl" else "manuscript_r1_clean"
    open(os.path.join(HERE, name + ".tex"), "w", encoding="utf-8").write(doc)
    log = ""
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", name + ".tex"], cwd=HERE, capture_output=True, text=True, encoding="utf-8", errors="replace")
        log = r.stdout
    errs = re.findall(r"(?m)^! .*", log)
    undef = re.findall(r"(?m)^.*(?:undefined|Undefined).*$", log)
    used = set(re.findall(re.escape(BS) + r"N\{([^}]+)\}", doc))
    defined = set(re.findall(re.escape(BS) + r"csname num@(.+?)" + re.escape(BS) + "endcsname", doc))
    qq = len(used - defined)                # \N{key} without a generated value would print ??[key]
    pages = re.findall(r"Output written on .*\((\d+) pages", log)
    print(f"[{mode}] pages={pages} errors={len(errs)} undefined-warnings={len(undef)} unresolved-numbers={qq} missing-bib={missing} unused-bib={len(unused)}")
    for e in errs[:8]:
        print("   ", e)
    for u in undef[:8]:
        print("   ", u.strip()[:160])
    if unused:
        print("   unused bib keys:", unused)
    return errs


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    for m in (["hl", "clean"] if which == "both" else [which]):
        build(m)
