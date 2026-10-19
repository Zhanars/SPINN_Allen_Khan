"""Word version of the revised manuscript in the MDPI Symmetry Word template.
usage: python make_docx_mdpi.py [template.dot]   -> manuscript_r1_mdpi.docx
Pipeline: LaTeX sources -> (make_docx.py helpers) -> pandoc -> body paragraphs/tables/equations/images are re-styled and
inserted into the template (front matter filled in, MDPI paragraph styles, equation/number tables, numbered references).
Run build_r1.py first (the .aux file gives the numbers of tables, figures, theorems, equations)."""
import copy
import io
import os
import re
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_docx as md  # noqa: E402
from docx import Document  # noqa: E402
from docx.opc.constants import RELATIONSHIP_TYPE as RT  # noqa: E402
from docx.oxml import OxmlElement, parse_xml  # noqa: E402
from docx.oxml.ns import qn, nsdecls  # noqa: E402

BS = chr(92)
TEMPLATE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "symmetry-template.dot")
OUT = os.path.join(HERE, "manuscript_r1_mdpi.docx")
BACK = ["Supplementary Materials", "Author Contributions", "Funding", "Institutional Review Board Statement", "Informed Consent Statement",
        "Data Availability Statement", "Acknowledgments", "Conflicts of Interest"]
TEXT_W = 7857      # text column of the template (twips)
FULL_W = 10466     # full page width between the margins
IND = 2608         # left indent of the text column
MAX_IMG = int(13.4 * 360000)  # EMU


# ------------------------------------------------------------------ pandoc input
def build_body_docx():
    md.EQ_SENTINEL = True
    md.LANDSCAPE = False
    md.ABBR_TABLE = True
    lab = md.read_aux()
    nums = md.read_numbers()
    bibkeys, bibtexts = md.read_bib()
    body = "".join(md.read(f) for f in ["10_intro.tex", "20_methods.tex", "30_results.tex", "40_discussion.tex", "50_conclusions.tex",
                                        "60_backmatter.tex", "70_appendix.tex"])
    body = md.flatten_inputs(body)
    body = re.sub(r"(?m)^\s*%.*$", "", body)
    body = re.sub(r"(?<!\\)%[^\n]*", "", body)
    body = md.unwrap(body, "rev")
    body = md.backmatter(body)
    body = md.expand_numbers(body, nums)
    body = md.ensuremath_modes(body)
    body = md.convert_envs(body, lab)
    body = md.resolve_cites(body, bibkeys)
    body = md.resolve_refs(body, lab)
    body = md.clean_misc(body)
    body = md.number_headings(body)
    refs = []
    for t in bibtexts:
        t = re.sub(r"\s+", " ", t).replace("~", " ").replace("--", "–")
        refs.append(t)
    tex = (BS + "documentclass{article}" + BS + "usepackage{graphicx}" + BS + "begin{document}\n" + body + "\n\n" + BS + "section*{References}\n\n"
           + "\n\n".join(refs) + "\n\n" + BS + "end{document}")
    path = os.path.join(HERE, "_mdpi_input.tex")
    open(path, "w", encoding="utf-8").write(tex)
    out = os.path.join(HERE, "_mdpi_body.docx")
    r = subprocess.run(["pandoc", "_mdpi_input.tex", "-f", "latex", "-t", "docx", "-o", out, "--resource-path=" + HERE], cwd=HERE,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.stderr.strip():
        print(r.stderr[:1200])
    return out, nums, lab


# ------------------------------------------------------------------ template as a docx package
def template_as_docx():
    out = os.path.join(HERE, "_template_as.docx")
    zin = zipfile.ZipFile(TEMPLATE)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
        for n in zin.namelist():
            data = zin.read(n)
            if n == "[Content_Types].xml":
                data = data.replace(b"wordprocessingml.template.main+xml", b"wordprocessingml.document.main+xml")
            zout.writestr(n, data)
    return out


# ------------------------------------------------------------------ small XML helpers
def ptext(p):
    return "".join(t.text or "" for t in p.iter(qn("w:t")))


def pstyle(p):
    ps = p.find(qn("w:pPr") + "/" + qn("w:pStyle"))
    return ps.get(qn("w:val")) if ps is not None else None


def mk_p(style, ind=None, jc=None, before=None, after=None):
    p = OxmlElement("w:p")
    pPr = OxmlElement("w:pPr")
    ps = OxmlElement("w:pStyle")
    ps.set(qn("w:val"), style)
    pPr.append(ps)
    if before is not None or after is not None:
        sp = OxmlElement("w:spacing")
        if before is not None:
            sp.set(qn("w:before"), str(before))
        if after is not None:
            sp.set(qn("w:after"), str(after))
        pPr.append(sp)
    if ind is not None:
        i = OxmlElement("w:ind")
        i.set(qn("w:left"), str(ind))
        pPr.append(i)
    if jc is not None:
        j = OxmlElement("w:jc")
        j.set(qn("w:val"), jc)
        pPr.append(j)
    p.append(pPr)
    return p


def mk_run(text, bold=False, italic=False, sup=False):
    r = OxmlElement("w:r")
    rPr = OxmlElement("w:rPr")
    if bold:
        rPr.append(OxmlElement("w:b"))
    if italic:
        rPr.append(OxmlElement("w:i"))
    if sup:
        v = OxmlElement("w:vertAlign")
        v.set(qn("w:val"), "superscript")
        rPr.append(v)
    if len(rPr):
        r.append(rPr)
    t = OxmlElement("w:t")
    t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    t.text = text
    r.append(t)
    return r


def set_runs(p, runs):
    for ch in list(p):
        if ch.tag != qn("w:pPr"):
            p.remove(ch)
    for text, fmt in runs:
        p.append(mk_run(text, **fmt))


class Rels:
    """copy images / hyperlinks from the pandoc document to the template document."""

    def __init__(self, bdoc, tdoc):
        self.bdoc, self.tdoc = bdoc, tdoc
        self.cache = {}

    def remap(self, el):
        for blip in el.iter(qn("a:blip")):
            rid = blip.get(qn("r:embed"))
            if rid and rid in self.bdoc.part.related_parts:
                if rid not in self.cache:
                    blob = self.bdoc.part.related_parts[rid].blob
                    new, _ = self.tdoc.part.get_or_add_image(io.BytesIO(blob))
                    self.cache[rid] = new
                blip.set(qn("r:embed"), self.cache[rid])
        for h in el.iter(qn("w:hyperlink")):
            rid = h.get(qn("r:id"))
            if rid and rid in self.bdoc.part.rels:
                key = "h" + rid
                if key not in self.cache:
                    rel = self.bdoc.part.rels[rid]
                    self.cache[key] = self.tdoc.part.relate_to(rel.target_ref, RT.HYPERLINK, is_external=True)
                h.set(qn("r:id"), self.cache[key])


def fix_runs(el):
    """pandoc character styles -> direct formatting that exists in the template."""
    for r in el.iter(qn("w:r")):
        rPr = r.find(qn("w:rPr"))
        if rPr is None:
            continue
        rs = rPr.find(qn("w:rStyle"))
        if rs is not None:
            name = rs.get(qn("w:val"))
            rPr.remove(rs)
            if name == "VerbatimChar":
                f = OxmlElement("w:rFonts")
                for a in ("w:ascii", "w:hAnsi", "w:cs"):
                    f.set(qn(a), "Consolas")
                rPr.insert(0, f)
            elif name == "Hyperlink":
                c = OxmlElement("w:color")
                c.set(qn("w:val"), "0563C1")
                rPr.append(c)


# ------------------------------------------------------------------ equations
def equation_table(omathpara, num):
    tbl = parse_xml(
        '<w:tbl %s><w:tblPr><w:tblW w:w="7859" w:type="dxa"/><w:tblInd w:w="%d" w:type="dxa"/>'
        '<w:tblCellMar><w:left w:w="0" w:type="dxa"/><w:right w:w="0" w:type="dxa"/></w:tblCellMar>'
        '<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/></w:tblPr>'
        '<w:tblGrid><w:gridCol w:w="7428"/><w:gridCol w:w="431"/></w:tblGrid><w:tr>'
        '<w:tc><w:tcPr><w:tcW w:w="7428" w:type="dxa"/></w:tcPr><w:p><w:pPr><w:pStyle w:val="MDPI39equation"/></w:pPr></w:p></w:tc>'
        '<w:tc><w:tcPr><w:tcW w:w="431" w:type="dxa"/><w:vAlign w:val="center"/></w:tcPr><w:p><w:pPr><w:pStyle w:val="MDPI3aequationnumber"/></w:pPr></w:p></w:tc>'
        '</w:tr></w:tbl>' % (nsdecls("w"), IND))
    cells = tbl.findall(".//" + qn("w:tc"))
    cells[0].find(qn("w:p")).append(omathpara)
    if num:
        cells[1].find(qn("w:p")).append(mk_run("(" + num + ")"))
    return tbl


def extract_sentinel(omathpara):
    """find '@@n@@' in the (possibly fragmented) math text, remove the runs that carry it, return n."""
    runs = []
    pos = 0
    text = ""
    for t in omathpara.iter(qn("m:t")):
        s = t.text or ""
        runs.append((t.getparent(), pos, pos + len(s)))
        text += s
        pos += len(s)
    m = re.search(r"@@(\w+)@@", text)
    if not m:
        return None
    a, b = m.span()
    seen = set()
    for r, s, e in runs:
        if e > a and s < b and id(r) not in seen:
            seen.add(id(r))
            if r.getparent() is not None:
                r.getparent().remove(r)
    return m.group(1)


# ------------------------------------------------------------------ tables
class _Unsupported(Exception):
    pass


def _flat(el, sup=False):
    out = []
    for ch in el:
        tag = ch.tag
        if tag == qn("m:r"):
            t = "".join(x.text or "" for x in ch.iter(qn("m:t")))
            if t:
                out.append((t, sup))
        elif tag in (qn("m:e"), qn("m:oMath")):
            out += _flat(ch, sup)
        elif tag == qn("m:sSup"):
            out += _flat(ch.find(qn("m:e")), sup)
            out += _flat(ch.find(qn("m:sup")), True)
        elif tag == qn("m:d"):
            beg, end = "(", ")"
            pr = ch.find(qn("m:dPr"))
            if pr is not None:
                b = pr.find(qn("m:begChr"))
                e = pr.find(qn("m:endChr"))
                if b is not None:
                    beg = b.get(qn("m:val"), "")
                if e is not None:
                    end = e.get(qn("m:val"), "")
            out.append((beg, sup))
            for k, e in enumerate(ch.findall(qn("m:e"))):
                if k:
                    out.append((",", sup))
                out += _flat(e, sup)
            out.append((end, sup))
        elif tag.endswith("Pr") or tag == qn("m:ctrlPr"):
            continue
        else:
            raise _Unsupported(tag)
    return out


NUMERIC = re.compile(r"^[0-9().,\s±×−\-+<>%≤≥]*$")


def numeric_math_to_runs(par, size):
    """inline math that is just a number like (6.76 +- 1.39) x 10^-2 becomes ordinary text runs (sized like the table)."""
    for om in [c for c in par if c.tag == qn("m:oMath")]:
        try:
            parts = _flat(om)
        except _Unsupported:
            continue
        if not parts or not all(NUMERIC.match(t) for t, _ in parts):
            continue
        new = []
        for t, sup in parts:
            t = t.replace("±", " ± ").replace("×", " × ")
            r = OxmlElement("w:r")
            rPr = OxmlElement("w:rPr")
            for tag in ("w:sz", "w:szCs"):
                e = OxmlElement(tag)
                e.set(qn("w:val"), size)
                rPr.append(e)
            if sup:
                v = OxmlElement("w:vertAlign")
                v.set(qn("w:val"), "superscript")
                rPr.append(v)
            r.append(rPr)
            tt = OxmlElement("w:t")
            tt.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
            tt.text = t
            r.append(tt)
            new.append(r)
        for r in new:
            om.addprevious(r)
        par.remove(om)


def restyle_table(tbl, abbreviation=False):
    grid = tbl.find(qn("w:tblGrid"))
    ncols = len(grid.findall(qn("w:gridCol")))
    wide = ncols >= 6
    total = FULL_W if wide else TEXT_W
    rows = tbl.findall(qn("w:tr"))

    def clen(tc):
        return len("".join(x.text or "" for x in tc.iter(qn("m:t"), qn("w:t"))))
    lens = [6] * ncols
    for tr in rows:
        ci = 0
        for tc in tr.findall(qn("w:tc")):
            tcPr = tc.find(qn("w:tcPr"))
            gs = tcPr.find(qn("w:gridSpan")) if tcPr is not None else None
            span = int(gs.get(qn("w:val"))) if gs is not None else 1
            if span == 1 and ci < ncols:
                lens[ci] = max(lens[ci], clen(tc))
            ci += span
    if abbreviation:
        widths = [936, total - 936]
    else:
        lens = [min(x, 40) if i == 0 else x for i, x in enumerate(lens)]
        s = float(sum(lens))
        widths = [max(int(total * x / s), 620) for x in lens]
        k = total / float(sum(widths))
        widths = [int(w * k) for w in widths]
    for ch in list(tbl.find(qn("w:tblPr"))):
        tbl.find(qn("w:tblPr")).remove(ch)
    tblPr = tbl.find(qn("w:tblPr"))
    tblPr.append(parse_xml('<w:tblW %s w:w="%d" w:type="dxa"/>' % (nsdecls("w"), sum(widths))))
    if wide:
        tblPr.append(parse_xml('<w:jc %s w:val="center"/>' % nsdecls("w")))
    else:
        tblPr.append(parse_xml('<w:tblInd %s w:w="%d" w:type="dxa"/>' % (nsdecls("w"), IND)))
    if not abbreviation:
        tblPr.append(parse_xml('<w:tblBorders %s><w:top w:val="single" w:sz="8" w:space="0" w:color="auto"/>'
                               '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="auto"/></w:tblBorders>' % nsdecls("w")))
    tblPr.append(parse_xml('<w:tblLayout %s w:type="fixed"/>' % nsdecls("w")))
    tblPr.append(parse_xml('<w:tblCellMar %s><w:left w:w="40" w:type="dxa"/><w:right w:w="40" w:type="dxa"/></w:tblCellMar>' % nsdecls("w")))
    for gc, w in zip(grid.findall(qn("w:gridCol")), widths):
        gc.set(qn("w:w"), str(w))
    sz = "18" if abbreviation else "16"
    for ri, tr in enumerate(rows):
        for pr in tr.findall(qn("w:trPr")):
            tr.remove(pr)
        ci = 0
        for tc in tr.findall(qn("w:tc")):
            tcPr = tc.find(qn("w:tcPr"))
            if tcPr is None:
                tcPr = OxmlElement("w:tcPr")
                tc.insert(0, tcPr)
            span = 1
            gs = tcPr.find(qn("w:gridSpan"))
            if gs is not None:
                span = int(gs.get(qn("w:val")))
            for ch in list(tcPr):
                tcPr.remove(ch)
            tcPr.append(parse_xml('<w:tcW %s w:w="%d" w:type="dxa"/>' % (nsdecls("w"), sum(widths[ci:ci + span]))))
            if span > 1:
                tcPr.append(parse_xml('<w:gridSpan %s w:val="%d"/>' % (nsdecls("w"), span)))
            if ri == 0 and not abbreviation:
                tcPr.append(parse_xml('<w:tcBorders %s><w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/></w:tcBorders>' % nsdecls("w")))
            tcPr.append(parse_xml('<w:vAlign %s w:val="center"/>' % nsdecls("w")))
            left = (ci == 0) or span > 1 or abbreviation
            for p in tc.findall(qn("w:p")):
                old = p.find(qn("w:pPr"))
                if old is not None:
                    p.remove(old)
                pPr = parse_xml('<w:pPr %s><w:pStyle w:val="MDPI42tablebody"/><w:spacing w:line="240" w:lineRule="auto"/>%s'
                                '<w:rPr><w:sz w:val="%s"/><w:szCs w:val="%s"/></w:rPr></w:pPr>'
                                % (nsdecls("w"), '<w:jc w:val="left"/>' if left else "", sz, sz))
                p.insert(0, pPr)
                numeric_math_to_runs(p, sz)
                for r in p.iter(qn("w:r")):
                    rPr = r.find(qn("w:rPr"))
                    if rPr is None:
                        rPr = OxmlElement("w:rPr")
                        r.insert(1 if r.find(qn("m:rPr")) is not None else 0, rPr)
                    for tag in ("w:sz", "w:szCs"):
                        for e in rPr.findall(qn(tag)):
                            rPr.remove(e)
                        e = OxmlElement(tag)
                        e.set(qn("w:val"), sz)
                        rPr.append(e)
            ci += span
    fix_runs(tbl)
    return tbl


# ------------------------------------------------------------------ main assembly
def main():
    body_docx, nums, lab = build_body_docx()
    tdoc = Document(template_as_docx())
    bdoc = Document(body_docx)
    rels = Rels(bdoc, tdoc)
    tbody = tdoc.element.body
    kids = list(tbody)
    assert len(kids) == 100, len(kids)
    assert pstyle(kids[98]) == 'MDPI63notes'

    # ---------------- front matter
    title, authors_tex, abstract, kw = md.front_matter(nums, lab)
    pre = md.read("00_preamble.tex")
    k = pre.find(BS + "conference{")
    conf, _ = md.grab(pre, k + len(BS + "conference"))
    conf = conf.replace(BS + '"u', "ü").replace("p.~4", "p. 4").replace("--", "–").replace("~", " ")
    set_runs(kids[0], [("Article", {})])
    set_runs(kids[1], [(title, {}), ("†", {"sup": True})])
    set_runs(kids[2], [("Zhuldyz Issa ", {}), ("1", {"sup": True}), (", Adil Nauryzbay ", {}), ("1", {"sup": True}), (", Zhanars Abdiramanov ", {}),
                       ("1,2,", {"sup": True}), ("*", {}), (", Yessimkhan Zhumazhanov ", {}), ("1", {"sup": True}), (" and Dastan Urnek ", {}),
                       ("1", {"sup": True})])

    def aff(p, mark, text, bold=False, sup=True):
        for ch in list(p):
            if ch.tag != qn("w:pPr"):
                p.remove(ch)
        p.append(mk_run(mark, bold=bold, sup=sup))
        r = OxmlElement("w:r")
        r.append(OxmlElement("w:tab"))
        t = OxmlElement("w:t")
        t.text = text
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        r.append(t)
        p.append(r)
    aff(kids[4], "1", "Department of Mathematics and Mathematical Modeling, Abai Kazakh National Pedagogical University, Almaty 050010, Kazakhstan")
    aff(kids[5], "2", "Institute of Information and Computational Technologies, Almaty 050010, Kazakhstan")
    aff(kids[6], "*", "Correspondence: zh.abdiramanov@abaiuniversity.edu.kz", bold=True, sup=False)
    dag = copy.deepcopy(kids[6])
    aff(dag, "†", "This paper is an extended version of our paper published in " + conf, bold=True, sup=False)
    kids[6].addnext(dag)
    abs_text = re.sub(r"\$([^$]*)\$", r"\1", abstract).replace(BS + "%", "%").replace("~", " ")
    set_runs(kids[8], [(abs_text, {})])
    for ch in list(kids[8]):
        pass
    kw_p = kids[9]
    set_runs(kw_p, [("Keywords: ", {"bold": True}), (kw, {})])
    msc = copy.deepcopy(kw_p)
    set_runs(msc, [("MSC: ", {"bold": True}), ("65M99; 35K57; 68T07", {})])
    sp = msc.find(qn("w:pPr")).find(qn("w:spacing"))
    pPr = msc.find(qn("w:pPr"))
    if sp is None:
        sp = OxmlElement("w:spacing")
        pPr.append(sp)
    sp.set(qn("w:before"), "0")
    kw_p.addnext(msc)

    # ---------------- remove the template body (items 11..97), keep the Publisher's note (98) and sectPr (99)
    for el in kids[11:98]:
        tbody.remove(el)
    anchor = kids[98]

    # ---------------- convert the pandoc body
    out = []
    state = {"section": None, "prev_eq": False}
    src = [c for c in bdoc.element.body if c.tag != qn("w:sectPr")]
    i = 0

    def style_for(p):
        sid = pstyle(p)
        txt = ptext(p).strip()
        return sid, txt

    def emit(el):
        out.append(el)

    def copy_children(src_p, dst_p, children=None):
        for ch in (children if children is not None else list(src_p)):
            if ch.tag == qn("w:pPr"):
                continue
            c = copy.deepcopy(ch)
            rels.remap(c)
            fix_runs(c)
            dst_p.append(c)

    refs_mode = False
    ref_n = 0
    while i < len(src):
        el = src[i]
        i += 1
        if el.tag == qn("w:tbl"):
            grid = el.find(qn("w:tblGrid"))
            ncols = len(grid.findall(qn("w:gridCol")))
            abbr = state["section"] == "Abbreviations"
            t = copy.deepcopy(el)
            rels.remap(t)
            emit(restyle_table(t, abbreviation=abbr))
            state["prev_eq"] = False
            continue
        if el.tag != qn("w:p"):
            continue
        sid, txt = style_for(el)
        has_draw = bool(list(el.iter(qn("w:drawing"))))
        if sid in ("Heading1", "Heading2", "Heading3"):
            level = int(sid[-1])
            if level == 1:
                state["section"] = txt
                if txt in BACK:
                    # merge heading + first paragraph into one back-matter paragraph
                    j = i
                    while j < len(src) and src[j].tag == qn("w:p") and not (pstyle(src[j]) or "").startswith("Heading"):
                        j += 1
                    first = src[i]
                    np_ = mk_p("MDPI62backmatter", before=240 if txt == BACK[0] else None)
                    np_.append(mk_run(txt + ": ", bold=True))
                    copy_children(first, np_)
                    emit(np_)
                    i = j
                    continue
                if txt == "References":
                    refs_mode = True
                np_ = mk_p("MDPI21heading1", ind=0 if txt in ("References", "Abbreviations") or txt.startswith("Appendix") else None)
                np_.append(mk_run(txt))
                emit(np_)
            else:
                np_ = mk_p("MDPI22heading2" if level == 2 else "MDPI23heading3")
                np_.append(mk_run(txt))
                emit(np_)
            state["prev_eq"] = False
            continue
        if refs_mode:
            ref_n += 1
            text = re.sub(r"^\s*\d+\.\s*", "", "".join(t.text or "" for t in el.iter(qn("w:t"))))
            np_ = mk_p("MDPI81references")
            pPr = np_.find(qn("w:pPr"))
            numPr = parse_xml('<w:numPr %s><w:ilvl w:val="0"/><w:numId w:val="4"/></w:numPr>' % nsdecls("w"))
            pPr.append(numPr)
            # keep italic / bold runs of the reference: copy runs, strip the leading number from the first text run
            copy_children(el, np_)
            first_t = next((t for t in np_.iter(qn("w:t"))), None)
            if first_t is not None:
                first_t.text = re.sub(r"^\s*\d+\.\s*", "", first_t.text or "")
            emit(np_)
            continue
        if has_draw and not txt:
            np_ = mk_p("MDPI52figure", ind=IND, jc="center")
            copy_children(el, np_)
            for ext in np_.iter("{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}extent"):
                cx, cy = int(ext.get("cx")), int(ext.get("cy"))
                if cx > MAX_IMG:
                    cy = int(cy * MAX_IMG / cx)
                    cx = MAX_IMG
                    ext.set("cx", str(cx))
                    ext.set("cy", str(cy))
                    for a_ext in np_.iter(qn("a:ext")):
                        if a_ext.get("cx") is not None:
                            a_ext.set("cx", str(cx))
                            a_ext.set("cy", str(cy))
            emit(np_)
            state["prev_eq"] = False
            continue
        # classify body paragraph
        first_bold = None
        r0 = el.find(qn("w:r"))
        if r0 is not None:
            rp = r0.find(qn("w:rPr"))
            first_bold = rp is not None and rp.find(qn("w:b")) is not None
        if first_bold and re.match(r"Figure \d+\.", txt):
            style = "MDPI51figurecaption"
        elif first_bold and re.match(r"Table [A-Z]?\d+\.", txt):
            style = "MDPI41tablecaption"
        elif first_bold and re.match(r"Theorem \d+", txt):
            style = "MDPI82theorem"
        elif txt.startswith("Proof."):
            style = "MDPI83proof"
        elif state["section"] == "Abbreviations":
            style = "MDPI32textnoindent"
        else:
            style = "MDPI31text"
        after_style = {"MDPI31text": "MDPI32textnoindent", "MDPI32textnoindent": "MDPI32textnoindent", "MDPI82theorem": "MDPI82theorem",
                       "MDPI83proof": "MDPI83proof"}.get(style, "MDPI32textnoindent")
        omaths = [c for c in el if c.tag == qn("m:oMathPara")]
        if not omaths:
            np_ = mk_p(style if not state["prev_eq"] or style != "MDPI31text" else style)
            copy_children(el, np_)
            emit(np_)
            state["prev_eq"] = False
            continue
        # paragraph with display equations: split into text / equation table / text
        seg = []
        cur_style = style
        for ch in list(el):
            if ch.tag == qn("w:pPr"):
                continue
            if ch.tag == qn("m:oMathPara"):
                if seg and any((x.tag != qn("w:r")) or (ptext(x).strip() if x.tag == qn("w:r") else True) for x in seg):
                    np_ = mk_p(cur_style)
                    for x in seg:
                        c = copy.deepcopy(x)
                        rels.remap(c)
                        fix_runs(c)
                        np_.append(c)
                    emit(np_)
                    cur_style = after_style
                seg = []
                om = copy.deepcopy(ch)
                num = extract_sentinel(om)
                emit(equation_table(om, num))
                cur_style = after_style
                continue
            seg.append(ch)
        if seg and "".join(ptext(x) if x.tag == qn("w:r") else "x" for x in seg).strip():
            np_ = mk_p(cur_style)
            for x in seg:
                c = copy.deepcopy(x)
                rels.remap(c)
                fix_runs(c)
                np_.append(c)
            emit(np_)
        state["prev_eq"] = True

    for el in out:
        anchor.addprevious(el)

    tdoc.save(OUT)
    for f in ("_mdpi_input.tex", "_mdpi_body.docx", "_template_as.docx"):
        try:
            os.remove(os.path.join(HERE, f))
        except OSError:
            pass
    print("written", OUT, os.path.getsize(OUT) // 1024, "KB;", len(out), "elements;", ref_n, "references")


if __name__ == "__main__":
    main()
