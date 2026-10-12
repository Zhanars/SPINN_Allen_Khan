"""Assemble the revised manuscript from manuscript/r1_src/*.tex and compile it.
usage: python build_r1.py [hl|clean|both]
 * inserts the highlight switch, flattens \\input of the number/abstract files,
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


def build(mode):
    pre = read("00_preamble.tex").replace(BS + "HLSWITCH", BS + ("hltrue" if mode == "hl" else "hlfalse"))
    body = "".join(read(f) for f in ["10_intro.tex", "20_methods.tex", "30_results.tex", "40_discussion.tex", "50_conclusions.tex", "60_backmatter.tex", "70_appendix.tex"] if os.path.exists(os.path.join(SRC, f)))
    refs = read("80_refs.tex")
    full_body = flatten_inputs(body)
    refs2, missing, unused = reorder_bib(full_body, refs)
    doc = flatten_inputs(pre) + "\n" + BS + "begin{document}\n\n" + full_body + "\n" + refs2
    name = "manuscript_r1_highlighted" if mode == "hl" else "manuscript_r1_clean"
    open(os.path.join(HERE, name + ".tex"), "w", encoding="utf-8").write(doc)
    log = ""
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", name + ".tex"], cwd=HERE, capture_output=True, text=True, encoding="utf-8", errors="replace")
        log = r.stdout
    errs = re.findall(r"(?m)^! .*", log)
    undef = re.findall(r"(?m)^.*(?:undefined|Undefined).*$", log)
    qq = doc.count("??[")
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
