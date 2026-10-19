"""Word version of the response to the reviewers.
usage: python make_docx_response.py   -> response_r1.docx
Needs build_r1.py to have run (numbers of sections/tables/figures come from build/manuscript_r1_clean.aux), pandoc, python-docx."""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_docx as md  # noqa: E402

BS = chr(92)
OUT = os.path.join(HERE, "response_r1.docx")


def main():
    lab = md.read_aux()
    nums = md.read_numbers()
    s = open(os.path.join(HERE, "response_r1.tex"), encoding="utf-8").read()
    s = s[s.index(BS + "begin{document}") + len(BS + "begin{document}"):s.index(BS + "end{document}")]
    s = re.sub(r"(?<!\\)%[^\n]*", "", s)
    # \rc{comment} -> Comment / Response paragraphs
    out = []
    key = BS + "rc{"
    i = 0
    while True:
        k = s.find(key, i)
        if k < 0:
            out.append(s[i:])
            break
        out.append(s[i:k])
        inner, j = md.grab(s, k + len(key) - 1)
        out.append("\n\n" + BS + "textbf{Comment.} " + BS + "emph{" + inner + "}\n\n" + BS + "textbf{Response.} ")
        i = j
    s = "".join(out)
    s = s.replace("{" + BS + "large" + BS + "bfseries Response to the Reviewers}", BS + "textbf{Response to the Reviewers}")
    s = re.sub(re.escape(BS + BS) + r"\[\d+pt\]", lambda m: BS + BS, s)
    s = re.sub(r"(\\begin\{(?:enumerate|itemize)\})\[[^\]]*\]", r"\1", s)
    s = md.expand_numbers(s, nums)
    s = md.ensuremath_modes(s)
    s = md.resolve_refs(s, lab)
    s = md.clean_misc(s)
    tex = BS + "documentclass{article}" + BS + "begin{document}\n" + s + "\n" + BS + "end{document}"
    path = os.path.join(HERE, "_resp_input.tex")
    open(path, "w", encoding="utf-8").write(tex)
    left = re.findall(r"\\(?:N|ensuremath|ref|eqref|rc)\{", tex)
    if left:
        print("WARNING leftovers:", sorted(set(left)))
    if "??" in tex:
        print("WARNING unresolved '??' in the text")
    r = subprocess.run(["pandoc", "_resp_input.tex", "-f", "latex", "-t", "docx", "-o", OUT], cwd=HERE, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    if r.stderr.strip():
        print(r.stderr[:1000])
    postprocess(OUT)
    try:
        os.remove(path)
    except OSError:
        pass
    print("written", OUT, os.path.getsize(OUT) // 1024, "KB")


def postprocess(path):
    from docx import Document
    from docx.shared import Pt, Cm, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml.ns import qn
    from docx.oxml import OxmlElement

    d = Document(path)
    for sec in d.sections:
        sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
        sec.left_margin = sec.right_margin = Cm(2.3)
        sec.top_margin = sec.bottom_margin = Cm(2.2)

    def setfont(st, size, bold=None):
        st.font.name = "Times New Roman"
        rpr = st.element.get_or_add_rPr()
        rf = rpr.find(qn("w:rFonts"))
        if rf is None:
            rf = OxmlElement("w:rFonts")
            rpr.append(rf)
        for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
            rf.set(qn(a), "Times New Roman")
        for a in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
            if rf.get(qn(a)) is not None:
                del rf.attrib[qn(a)]
        st.font.size = Pt(size)
        if bold is not None:
            st.font.bold = bold
    for st in d.styles:
        try:
            if st.type != 1:
                continue
            if st.name in ("Normal", "Body Text", "First Paragraph", "Compact", "Block Text"):
                setfont(st, 11)
            elif st.name.startswith("Heading"):
                setfont(st, 13 if st.name == "Heading 1" else 12, True)
                st.font.color.rgb = RGBColor(0, 0, 0)
                st.font.italic = False
        except Exception:
            pass
    header_done = False
    for p in d.paragraphs:
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 1.15
        if p.text.startswith("Dear Editor"):
            header_done = True
        if not header_done and p.text.strip():
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER     # title block
            if p.runs:
                p.runs[0].font.size = Pt(14)
                p.runs[0].font.bold = True
        elif p.style.name in ("Body Text", "First Paragraph", "Compact"):
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        if p.text.startswith("Comment."):
            p.paragraph_format.left_indent = Cm(0.6)
            p.paragraph_format.space_before = Pt(10)
    d.save(path)


if __name__ == "__main__":
    main()
