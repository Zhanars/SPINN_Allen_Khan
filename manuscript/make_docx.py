"""Build a Word version of the revised manuscript (clean) from manuscript/r1_src.
usage: python make_docx.py          -> manuscript_r1_clean.docx
Needs: pandoc, python-docx. Run build_r1.py first (the .aux file gives the numbers of tables, figures, theorems, equations).
Cross-references, citations, \\N{...} numbers, theorem/figure/table numbers are resolved here; pandoc does the rest (equations become Word equations)."""
import os
import re
import subprocess
import sys

BS = chr(92)
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "r1_src")
sys.path.insert(0, HERE)
from build_r1 import flatten_inputs  # noqa: E402


def read(fn):
    return open(os.path.join(SRC, fn), encoding="utf-8").read()


# ------------------------------------------------------------------ helpers
def grab(s, i):
    """s[i] == '{': return (inner, index after the matching '}')."""
    assert s[i] == "{", s[i:i + 40]
    d = 0
    j = i
    while j < len(s):
        c = s[j]
        if c == BS:
            j += 2
            continue
        if c == "{":
            d += 1
        elif c == "}":
            d -= 1
            if d == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ValueError("unbalanced: " + s[i:i + 60])


def unwrap(s, name):
    """\\name{X} -> X (all occurrences)."""
    key = BS + name + "{"
    out = []
    i = 0
    while True:
        k = s.find(key, i)
        if k < 0:
            out.append(s[i:])
            break
        out.append(s[i:k])
        inner, j = grab(s, k + len(key) - 1)
        out.append(unwrap(inner, name))
        i = j
    return "".join(out)


def remove_macro(s, name, nargs=1, keep_arg=None):
    """remove \\name{..}{..} (nargs groups); keep_arg: index of the argument to keep."""
    key = BS + name
    out = []
    i = 0
    while True:
        k = s.find(key, i)
        if k < 0:
            out.append(s[i:])
            break
        e = k + len(key)
        if e < len(s) and (s[e].isalpha()):
            out.append(s[i:e])
            i = e
            continue
        out.append(s[i:k])
        j = e
        args = []
        for _ in range(nargs):
            while j < len(s) and s[j] in " \n":
                j += 1
            if j < len(s) and s[j] == "{":
                a, j = grab(s, j)
                args.append(a)
            else:
                break
        if keep_arg is not None and keep_arg < len(args):
            out.append(args[keep_arg])
        i = j
    return "".join(out)


# ------------------------------------------------------------------ data from the build
def read_aux():
    aux = open(os.path.join(HERE, "manuscript_r1_clean.aux"), encoding="utf-8").read()
    lab = {}
    for m in re.finditer(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}", aux):
        lab[m.group(1)] = m.group(2)
    return lab


def read_numbers():
    d = {}
    txt = read("numbers_r1.tex")
    for line in txt.split("\n"):
        m = re.match(r"\\expandafter\\def\\csname num@(.*?)\\endcsname\{(.*)\}\s*$", line)
        if m:
            d[m.group(1)] = m.group(2)
    return d


def read_bib():
    clean = open(os.path.join(HERE, "manuscript_r1_clean.tex"), encoding="utf-8").read()
    k = clean.index(BS + "begin{thebibliography}")
    body = clean[k:]
    items = re.split(r"(?m)^" + re.escape(BS) + r"bibitem\{", body)[1:]
    keys, texts = [], []
    for it in items:
        key, rest = it.split("}", 1)
        rest = rest.split(BS + "end{thebibliography}")[0].strip()
        keys.append(key)
        texts.append(rest)
    return keys, texts


# ------------------------------------------------------------------ text transformations
def expand_numbers(s, nums):
    out = []
    i = 0
    key = BS + "N{"
    while True:
        k = s.find(key, i)
        if k < 0:
            out.append(s[i:])
            break
        out.append(s[i:k])
        inner, j = grab(s, k + 2)
        out.append(nums.get(inner, "??"))
        i = j
    return "".join(out)


MATH_ON = [BS + "[", BS + "begin{equation}", BS + "begin{align}"]
MATH_OFF = [BS + "]", BS + "end{equation}", BS + "end{align}"]


def ensuremath_modes(s):
    """\\ensuremath{X} -> $X$ in text mode, X in math mode."""
    out = []
    i = 0
    in_math = False
    key = BS + "ensuremath{"
    n = len(s)
    while i < n:
        if s.startswith(key, i):
            inner, j = grab(s, i + len(key) - 1)
            inner = ensuremath_modes(inner)
            out.append(inner if in_math else "$" + inner + "$")
            i = j
            continue
        if s[i] == BS and i + 1 < n and s[i + 1] == BS:
            out.append(BS + BS)
            i += 2
            continue
        if s[i] == BS and i + 1 < n and s[i + 1] == "$":
            out.append(BS + "$")
            i += 2
            continue
        handled = False
        for t in MATH_ON:
            if s.startswith(t, i):
                in_math = True
                out.append(t)
                i += len(t)
                handled = True
                break
        if handled:
            continue
        for t in MATH_OFF:
            if s.startswith(t, i):
                in_math = False
                out.append(t)
                i += len(t)
                handled = True
                break
        if handled:
            continue
        if s[i] == "$":
            in_math = not in_math
        out.append(s[i])
        i += 1
    return "".join(out)


def fmt_cites(nums_list):
    nums_list = sorted(set(nums_list))
    parts = []
    i = 0
    while i < len(nums_list):
        j = i
        while j + 1 < len(nums_list) and nums_list[j + 1] == nums_list[j] + 1:
            j += 1
        if j - i >= 2:
            parts.append(f"{nums_list[i]}\u2013{nums_list[j]}")
        else:
            parts.extend(str(x) for x in nums_list[i:j + 1])
        i = j + 1
    return "[" + ",".join(parts) + "]"


def resolve_cites(s, bibkeys):
    def rep(m):
        ks = [k.strip() for k in m.group(1).split(",")]
        return fmt_cites([bibkeys.index(k) + 1 for k in ks if k in bibkeys])
    return re.sub(r"~?" + re.escape(BS) + r"cite[a-z]*\{([^}]+)\}", lambda m: " " + rep(m), s)


def resolve_refs(s, lab):
    def eq(m):
        return "(" + lab.get(m.group(1), "??") + ")"

    def rf(m):
        return lab.get(m.group(1), "??")
    s = re.sub(re.escape(BS) + r"eqref\{([^}]+)\}", eq, s)
    s = re.sub(re.escape(BS) + r"ref\{([^}]+)\}", rf, s)
    return s


def convert_envs(s, lab):
    # ---- theorem-like
    def thm(m):
        title = m.group(2) or ""
        body = m.group(3)
        lm = re.search(re.escape(BS) + r"label\{([^}]+)\}", body)
        num = lab.get(lm.group(1), "") if lm else ""
        body = re.sub(re.escape(BS) + r"label\{[^}]+\}", "", body)
        head = BS + "textbf{Theorem " + num + "}" + (" (" + title + ")." if title else ".")
        return "\n\n" + head + " " + body.strip() + "\n\n"
    s = re.sub(re.escape(BS) + r"begin\{Theorem\}(\[(.*?)\])?(.*?)" + re.escape(BS) + r"end\{Theorem\}", thm, s, flags=re.S)
    s = re.sub(re.escape(BS) + r"begin\{proof\}(.*?)" + re.escape(BS) + r"end\{proof\}",
               lambda m: "\n\n" + BS + "textit{Proof.} " + m.group(1).strip() + " $\\square$\n\n", s, flags=re.S)

    # ---- figures
    def fig(m):
        body = m.group(1)
        im = re.search(re.escape(BS) + r"includegraphics(\[[^\]]*\])?\{([^}]+)\}", body)
        path = im.group(2)
        opt = im.group(1) or ""
        ck = body.index(BS + "caption{")
        cap, _ = grab(body, ck + len(BS + "caption"))
        lm = re.search(re.escape(BS) + r"label\{([^}]+)\}", cap)
        num = lab.get(lm.group(1), "") if lm else ""
        cap = re.sub(re.escape(BS) + r"label\{[^}]+\}", "", cap)
        w = "width=15cm"
        mm = re.search(r"width=([0-9.]+)cm", opt)
        if mm:
            w = "width=" + mm.group(1) + "cm"
        return ("\n\n" + BS + "begin{center}" + BS + "includegraphics[" + w + "]{" + path + "}" + BS + "end{center}\n\n"
                + BS + "textbf{Figure " + num + ".} " + cap.strip() + "\n\n")
    s = re.sub(re.escape(BS) + r"begin\{figure\}(?:\[[^\]]*\])?(.*?)" + re.escape(BS) + r"end\{figure\}", fig, s, flags=re.S)

    # ---- tables
    def tab(m):
        body = m.group(1)
        ck = body.index(BS + "caption{")
        cap, e = grab(body, ck + len(BS + "caption"))
        lm = re.search(re.escape(BS) + r"label\{([^}]+)\}", cap)
        num = lab.get(lm.group(1), "") if lm else ""
        cap = re.sub(re.escape(BS) + r"label\{[^}]+\}", "", cap)
        rest = body[:ck] + body[e:]
        rest = re.sub(re.escape(BS) + r"begin\{adjustwidth\}\{[^}]*\}\{[^}]*\}", "", rest)
        rest = rest.replace(BS + "end{adjustwidth}", "")
        rest = rest.replace(BS + "tablesize{" + BS + "footnotesize}", "")
        rest = remove_macro(rest, "newcolumntype", 2)
        rest = rest.replace(BS + "small", "").replace(BS + "footnotesize{}", "")

        def tabx(mm):
            spec_raw = mm.group(2)
            spec = tabspec(spec_raw)
            return BS + "begin{tabular}{" + spec + "}" + mm.group(3) + BS + "end{tabular}"
        rest = re.sub(re.escape(BS) + r"begin\{tabularx\}\{[^}]*\}(\{((?:[^{}]|\{[^{}]*\})*)\})(.*?)" + re.escape(BS) + r"end\{tabularx\}",
                      lambda mm: BS + "begin{tabular}{" + tabspec(mm.group(2)) + "}" + mm.group(3) + BS + "end{tabular}", rest, flags=re.S)
        rest = rest.replace(BS + "noindent{", "{")
        block = "\n\n" + BS + "textbf{Table " + num + ".} " + cap.strip() + "\n\n" + rest.strip() + "\n\n"
        km = re.search(re.escape(BS) + r"begin\{tabular\}\{([^}]*)\}", rest)
        if km and len(km.group(1)) >= 6:
            block = "\n\nLANDSCAPESTART\n\n" + block + "LANDSCAPEEND\n\n"
        return block
    s = re.sub(re.escape(BS) + r"begin\{table\}(?:\[[^\]]*\])?(.*?)" + re.escape(BS) + r"end\{table\}", tab, s, flags=re.S)

    # ---- equations: number at the end of the (first) row
    def eqn(m):
        env, body = m.group(1), m.group(2)
        lm = re.search(re.escape(BS) + r"label\{([^}]+)\}", body)
        body = re.sub(re.escape(BS) + r"label\{[^}]+\}", "", body)
        body = body.replace(BS + "notag", "")
        if lm and lm.group(1) in lab:
            tag = BS + "qquad (" + lab[lm.group(1)] + ")"
            if env == "align" and (BS + BS) in body:
                k = body.index(BS + BS)
                body = body[:k] + tag + body[k:]
            else:
                body = body.rstrip() + tag
        if env == "align":
            return BS + "begin{aligned}" + body + BS + "end{aligned}" if False else BS + "begin{align*}" + body + BS + "end{align*}"
        return BS + "begin{equation*}" + body + BS + "end{equation*}"
    s = re.sub(re.escape(BS) + r"begin\{(equation|align)\}(.*?)" + re.escape(BS) + r"end\{\1\}", eqn, s, flags=re.S)
    return s


def tabspec(raw):
    out = []
    i = 0
    while i < len(raw):
        c = raw[i]
        if c == "p" and i + 1 < len(raw) and raw[i + 1] == "{":
            _, j = grab(raw, i + 1)
            out.append("l")
            i = j
        elif c == ">" and i + 1 < len(raw) and raw[i + 1] == "{":
            _, j = grab(raw, i + 1)
            i = j
        elif c in "lcr":
            out.append(c)
            i += 1
        elif c in "XC":
            out.append("c" if c == "C" else "l")
            i += 1
        else:
            i += 1
    return "".join(out)


def clean_misc(s):
    s = re.sub(r"\{" + re.escape(BS) + r"rm\s+([^{}]*)\}", lambda m: BS + "mathrm{" + m.group(1) + "}", s)
    s = s.replace(BS + "mbox{", BS + "text{")
    s = s.replace(BS + "tfrac", BS + "frac").replace(BS + "dfrac", BS + "frac")
    for name in ["label"]:
        s = re.sub(re.escape(BS) + name + r"\{[^}]*\}", "", s)
    s = re.sub(re.escape(BS) + r"(needspace|vspace|hspace)\{[^}]*\}", "", s)
    s = s.replace(BS + "noindent", "").replace(BS + "medskip", "").replace(BS + "newpage", "")
    s = re.sub(re.escape(BS) + r"texorpdfstring\{([^}]*)\}\{[^}]*\}", r"\1", s)
    s = s.replace(BS + "sloppy", "")
    s = re.sub(re.escape(BS) + r"def" + re.escape(BS) + r"theH[^\n]*", "", s)
    s = s.replace(BS + "appendixtitles{yes}", "").replace(BS + "appendixstart", "")
    s = s.replace("~", " ")
    s = s.replace("---", "\u2014").replace("--", "\u2013")
    return s


def number_headings(s):
    """\\section etc. -> unnumbered pandoc headings with manual numbers."""
    sec = sub = subsub = 0
    appendix = False
    letters = "ABCDEFGH"
    out = []
    i = 0
    pat = re.compile(re.escape(BS) + r"(section|subsection|subsubsection|appendix)(\*?)(\{)?")
    pos = 0
    res = []
    while True:
        m = pat.search(s, pos)
        if not m:
            res.append(s[pos:])
            break
        res.append(s[pos:m.start()])
        kind = m.group(1)
        if kind == "appendix":
            appendix = True
            sec = 0
            pos = m.end() - (1 if m.group(3) else 0)
            continue
        title, j = grab(s, m.end() - 1)
        if m.group(2) == "*":
            res.append(BS + kind + "*{" + title + "}")
            pos = j
            continue
        if kind == "section":
            sec += 1
            sub = subsub = 0
            num = (f"Appendix {letters[sec - 1]}." if appendix else f"{sec}.")
            res.append(BS + "section*{" + num + " " + title + "}")
        elif kind == "subsection":
            sub += 1
            subsub = 0
            pre = letters[sec - 1] if appendix else str(sec)
            res.append(BS + "subsection*{" + f"{pre}.{sub}. " + title + "}")
        else:
            subsub += 1
            pre = letters[sec - 1] if appendix else str(sec)
            res.append(BS + "subsubsection*{" + f"{pre}.{sub}.{subsub}. " + title + "}")
        pos = j
    return "".join(res)


def backmatter(s):
    heads = [("supplementary", "Supplementary Materials"), ("authorcontributions", "Author Contributions"), ("funding", "Funding"),
             ("institutionalreview", "Institutional Review Board Statement"), ("informedconsent", "Informed Consent Statement"),
             ("dataavailability", "Data Availability Statement"), ("acknowledgments", "Acknowledgments"),
             ("conflictsofinterest", "Conflicts of Interest")]
    for name, title in heads:
        key = BS + name + "{"
        k = s.find(key)
        if k < 0:
            continue
        inner, j = grab(s, k + len(key) - 1)
        inner = re.sub(re.escape(BS) + r"linksupplementary\{[^}]*\}", "the journal website", inner)
        s = s[:k] + "\n\n" + BS + "section*{" + title + "}\n\n" + inner + "\n\n" + s[j:]
    k = s.find(BS + "abbreviations{")
    if k >= 0:
        a, j = grab(s, k + len(BS + "abbreviations"))
        b, j2 = grab(s, j)
        kb = b.find(BS + "begin{tabular}")
        if kb >= 0:
            _, jb = grab(b, kb + len(BS + "begin{tabular}"))
            ke = b.index(BS + "end{tabular}", jb)
            rows = [r for r in b[jb:ke].split(BS + BS) if "&" in r]
            lines = [BS + "textbf{" + r.split("&")[0].strip() + "} \u2014 " + r.split("&")[1].strip() for r in rows]
            b = b[:kb] + "\n\n".join(lines) + b[ke + len(BS + "end{tabular}"):]
        b = b.replace(BS + BS, "\n\n")
        s = s[:k] + "\n\n" + BS + "section*{" + a + "}\n\n" + b + "\n\n" + s[j2:]
    return s


def front_matter(nums, lab):
    pre = read("00_preamble.tex")

    def arg(name):
        k = pre.find(BS + name + "{")
        a, _ = grab(pre, k + len(BS + name))
        return a
    title = arg("Title")
    title = re.sub(re.escape(BS) + r"texorpdfstring\{[^}]*\}\{[^}]*\}", "", title).replace("--", "\u2013")
    title = title.replace("$^" + BS + "dagger$", "")
    title = title.strip()
    sup = BS + "textsuperscript"
    authors = ("Zhuldyz Issa" + sup + "{1}, Adil Nauryzbay" + sup + "{1}, Zhanars Abdiramanov" + sup + "{1,2,*}, Yessimkhan Zhumazhanov" + sup
               + "{1} and Dastan Urnek" + sup + "{1}")
    abstract = unwrap(read("abstract_r1.tex"), "rev")
    abstract = expand_numbers(abstract, nums)
    abstract = ensuremath_modes(abstract)
    abstract = abstract.replace("--", "\u2013")
    kw = arg("keyword").replace("--", "\u2013")
    return title, authors, abstract, kw


def main():
    lab = read_aux()
    nums = read_numbers()
    bibkeys, bibtexts = read_bib()
    sup = BS + "textsuperscript"

    body = "".join(read(f) for f in ["10_intro.tex", "20_methods.tex", "30_results.tex", "40_discussion.tex", "50_conclusions.tex",
                                     "60_backmatter.tex", "70_appendix.tex"])
    body = flatten_inputs(body)
    body = re.sub(r"(?m)^\s*%.*$", "", body)
    body = re.sub(r"(?<!\\)%[^\n]*", "", body)
    body = unwrap(body, "rev")
    body = backmatter(body)
    body = expand_numbers(body, nums)
    body = ensuremath_modes(body)
    body = convert_envs(body, lab)
    body = resolve_cites(body, bibkeys)
    body = resolve_refs(body, lab)
    body = clean_misc(body)
    body = number_headings(body)

    title, authors, abstract, kw = front_matter(nums, lab)
    refs = []
    for i, t in enumerate(bibtexts, 1):
        t = re.sub(r"\s+", " ", t)
        t = t.replace("~", " ").replace("--", "\u2013")
        refs.append(f"{i}. " + t)
    refs_tex = "\n\n".join(refs)

    doc = []
    doc.append(BS + "documentclass{article}" + BS + "usepackage{graphicx}" + BS + "begin{document}")
    doc.append(BS + "title{" + title + "}")
    doc.append(BS + "maketitle")
    doc.append("\n\n" + authors + "\n\n")
    doc.append(sup + "{1} Department of Mathematics and Mathematical Modeling, Abai Kazakh National Pedagogical University, Almaty 050010, Kazakhstan\n\n"
               + sup + "{2} Institute of Information and Computational Technologies, Almaty 050010, Kazakhstan\n\n"
               + sup + "{*} Correspondence: zh.abdiramanov@abaiuniversity.edu.kz\n\n")
    doc.append(BS + "section*{Abstract}\n\n" + abstract + "\n\n" + BS + "textbf{Keywords:} " + kw + "\n\n" + BS + "textbf{MSC:} 65M99; 35K57; 68T07\n\n")
    doc.append(body)
    doc.append("\n\n" + BS + "section*{References}\n\n" + refs_tex + "\n\n")
    doc.append(BS + "end{document}")
    tex = "".join(doc)
    # sanity: leftovers
    left = sorted(set(re.findall(r"\\(?:N|ensuremath|rev|ref|eqref|cite)\{", tex)))
    if left:
        print("WARNING leftovers:", left)
    open(os.path.join(HERE, "_docx_input.tex"), "w", encoding="utf-8").write(tex)
    out = os.path.join(HERE, "manuscript_r1_clean.docx")
    r = subprocess.run(["pandoc", "_docx_input.tex", "-f", "latex", "-t", "docx", "-o", out, "--resource-path=" + HERE], cwd=HERE,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stderr[:1500])
    postprocess(out)
    print("written", out, os.path.getsize(out) // 1024, "KB")
    try:
        os.remove(os.path.join(HERE, "_docx_input.tex"))
    except OSError:
        pass


def postprocess(path):
    from docx import Document
    from docx.shared import Pt, Cm
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml.ns import qn
    from docx.oxml import OxmlElement

    d = Document(path)
    for sec in d.sections:
        sec.page_width = Cm(21.0)
        sec.page_height = Cm(29.7)
        sec.left_margin = sec.right_margin = Cm(2.2)
        sec.top_margin = sec.bottom_margin = Cm(2.2)

    def setfont(style, name, size=None, bold=None):
        style.font.name = name
        rpr = style.element.get_or_add_rPr()
        rf = rpr.find(qn("w:rFonts"))
        if rf is None:
            rf = OxmlElement("w:rFonts")
            rpr.append(rf)
        for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
            rf.set(qn(a), name)
        for a in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
            if rf.get(qn(a)) is not None:
                del rf.attrib[qn(a)]
        if size:
            style.font.size = Pt(size)
        if bold is not None:
            style.font.bold = bold
    for st in d.styles:
        try:
            if st.type == 1:
                nm = st.name
                if nm in ("Normal", "Body Text", "First Paragraph", "Compact", "Block Text", "Footnote Text", "Table Caption", "Image Caption"):
                    setfont(st, "Times New Roman", 11)
                elif nm.startswith("Heading") or nm in ("Title", "Subtitle", "Author", "Date", "Abstract"):
                    size = {"Heading 1": 14, "Heading 2": 12, "Heading 3": 11, "Title": 17}.get(nm, 11)
                    setfont(st, "Times New Roman", size, True)
                    st.font.color.rgb = None
                    try:
                        from docx.shared import RGBColor
                        st.font.color.rgb = RGBColor(0, 0, 0)
                    except Exception:
                        pass
                    st.font.italic = (nm == "Heading 3")
        except Exception:
            pass
    for p in d.paragraphs:
        pf = p.paragraph_format
        pf.space_after = Pt(6)
        pf.line_spacing = 1.15
        if p.style.name in ("Body Text", "First Paragraph", "Normal"):
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        if p._p.xpath(".//w:drawing") and not p.text.strip():
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    # ---- section breaks for wide tables (landscape pages)
    from copy import deepcopy
    body = d.element.body
    base_sect = body.find(qn("w:sectPr"))

    def make_sect(landscape):
        sp = deepcopy(base_sect)
        for el in sp.findall(qn("w:type")):
            sp.remove(el)
        sz = sp.find(qn("w:pgSz"))
        if sz is None:
            sz = OxmlElement("w:pgSz")
            sp.insert(0, sz)
        if landscape:
            sz.set(qn("w:w"), "16838")
            sz.set(qn("w:h"), "11906")
            sz.set(qn("w:orient"), "landscape")
        else:
            sz.set(qn("w:w"), "11906")
            sz.set(qn("w:h"), "16838")
            if sz.get(qn("w:orient")) is not None:
                del sz.attrib[qn("w:orient")]
        return sp
    for p in list(d.paragraphs):
        txt = p.text.strip()
        if txt in ("LANDSCAPESTART", "LANDSCAPEEND"):
            for r in list(p.runs):
                r._r.getparent().remove(r._r)
            pPr = p._p.get_or_add_pPr()
            pPr.append(make_sect(txt == "LANDSCAPEEND"))
            p.paragraph_format.space_after = Pt(0)

    # ---- tables: small font, booktabs-like rules, fixed column widths
    mm = 567.0  # twips per cm
    for t in d.tables:
        ncols = len(t.columns)
        landscape = ncols >= 6
        total = int((25.3 if landscape else 16.6) * mm)

        def cell_len(c):
            txt = "".join(x.text or "" for x in c._tc.iter(qn("m:t"), qn("w:t")))
            return len(txt)
        lens = [6] * ncols
        for row in t.rows:
            ci = 0
            for c in row.cells:
                if ci < ncols:
                    gs = c._tc.get_or_add_tcPr().find(qn("w:gridSpan"))
                    if gs is None:
                        lens[ci] = max(lens[ci], cell_len(c))
                ci += 1
        lens = [min(l, 46) if i == 0 else l for i, l in enumerate(lens)]
        tot_len = float(sum(lens))
        widths = [max(int(total * l / tot_len), int(1.3 * mm)) for l in lens]
        scale = total / float(sum(widths))
        widths = [int(w * scale) for w in widths]
        tbl = t._tbl
        tblPr = tbl.tblPr
        for el in tblPr.findall(qn("w:tblBorders")) + tblPr.findall(qn("w:tblLayout")) + tblPr.findall(qn("w:tblW")):
            tblPr.remove(el)
        tw = OxmlElement("w:tblW")
        tw.set(qn("w:w"), str(sum(widths)))
        tw.set(qn("w:type"), "dxa")
        tblPr.append(tw)
        b = OxmlElement("w:tblBorders")
        for edge in ("top", "bottom"):
            e = OxmlElement("w:" + edge)
            e.set(qn("w:val"), "single")
            e.set(qn("w:sz"), "8")
            e.set(qn("w:space"), "0")
            e.set(qn("w:color"), "000000")
            b.append(e)
        tblPr.append(b)
        lay = OxmlElement("w:tblLayout")
        lay.set(qn("w:type"), "fixed")
        tblPr.append(lay)
        grid = tbl.find(qn("w:tblGrid"))
        if grid is not None:
            for gc, w in zip(grid.findall(qn("w:gridCol")), widths):
                gc.set(qn("w:w"), str(w))
        for ri, row in enumerate(t.rows):
            for ci, c in enumerate(row.cells):
                tcPr = c._tc.get_or_add_tcPr()
                for el in tcPr.findall(qn("w:tcW")):
                    tcPr.remove(el)
                span = 1
                gs = tcPr.find(qn("w:gridSpan"))
                if gs is not None:
                    span = int(gs.get(qn("w:val")))
                tcw = OxmlElement("w:tcW")
                tcw.set(qn("w:w"), str(sum(widths[ci:ci + span])))
                tcw.set(qn("w:type"), "dxa")
                tcPr.insert(0, tcw)
                for p in c.paragraphs:
                    p.paragraph_format.space_after = Pt(0)
                    p.paragraph_format.line_spacing = 1.0
                    p.alignment = WD_ALIGN_PARAGRAPH.LEFT if ci == 0 or span > 1 else WD_ALIGN_PARAGRAPH.CENTER
                    for r in p.runs:
                        r.font.size = Pt(8)
                    for mr in p._p.iter(qn("m:r")):
                        rp = mr.find(qn("w:rPr"))
                        if rp is None:
                            rp = OxmlElement("w:rPr")
                            mr.insert(1 if mr.find(qn("m:rPr")) is not None else 0, rp)
                        for tag in ("w:sz", "w:szCs"):
                            for el in rp.findall(qn(tag)):
                                rp.remove(el)
                            el = OxmlElement(tag)
                            el.set(qn("w:val"), "16")
                            rp.append(el)
                if ri == 0:
                    bd = OxmlElement("w:tcBorders")
                    e = OxmlElement("w:bottom")
                    e.set(qn("w:val"), "single")
                    e.set(qn("w:sz"), "4")
                    e.set(qn("w:space"), "0")
                    e.set(qn("w:color"), "000000")
                    bd.append(e)
                    tcPr.append(bd)
    d.save(path)


if __name__ == "__main__":
    main()
