"""
Parse Aldy .aldy output and normalise the major diplotype so that it can be
compared with the deCYPher gold standard using CyriPanel's own
normalize_genotype()/compare_genotype() logic.

Aldy-specific cleaning applied before handing off to CyriPanel's normaliser:
  1. Strip '.ALDY' suffixes      : *36.ALDY  -> *36
  2. Strip '+rsNNNNN' modifiers  : *1+rs28371706 -> *1
     ('+*NN' tandem components are kept)
"""
import re, os, sys

sys.path.insert(0, os.environ['SRC'] + '/benchmark')
from run_single_iteration import normalize_genotype, compare_genotype   # noqa


def clean_aldy_allele(hap):
    """
    Reduce one Aldy haplotype string to core star alleles.
      *36.ALDY -> *36     (Aldy approximate-allele marker)
      *4.028   -> *4      (sub-allele -> major allele, matching gold level)
      rs769258 -> dropped (attached variant, absent from gold standard)
    """
    kept = []
    for p in hap.split('+'):
        p = p.strip()
        if not p:
            continue
        m = re.match(r'^\*(\d+)', p)     # keep only *<digits>
        if m:
            kept.append('*' + m.group(1))
        # anything else (rsNNNN, stray tokens) is dropped
    return '+'.join(kept) if kept else None


def clean_aldy_diplotype(gt):
    if gt is None or gt in ('', 'None', 'NO_CALL'):
        return None
    gt = gt.strip()
    if '/' not in gt:
        return clean_aldy_allele(gt)
    return '/'.join(clean_aldy_allele(h) for h in gt.split('/'))


def parse_aldy_file(path):
    """Return (major_diplotype, n_solutions) from a .aldy file."""
    if not os.path.exists(path):
        return None, 0
    majors, nsol = [], 0
    with open(path) as f:
        for line in f:
            if line.startswith('#Solution'):
                nsol += 1
                continue
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) > 3 and fields[3]:
                majors.append(fields[3])
    if not majors:
        return None, nsol
    uniq = list(dict.fromkeys(majors))
    return uniq[0], nsol


if __name__ == '__main__':
    tests = [
        ('*36.ALDY+*10/*49',        '*49/*10+*36'),   # NA24631
        ('*1+rs28371706/*67',       '*1/*67'),
        ('*5/*13+rs59421388',       '*5/*13'),
        ('*1+rs769258/*1+rs769258', '*1/*1'),
        ('*1+*67/*36.ALDY',         '*36/*1+*67'),
        ('*1/*1+*1',                '*1/*1+*1'),
        ('*43/*43+*43',             '*43/*43+*43'),
        ('*1/*5',                   '*1/*5'),
        ('*4.028/*35',              '*35/*4'),
        ('*83.ALDY/*80',            '*80/*83'),
        ('*5/*74+rs766391487+rs28371704+rs769258', '*5/*74'),
        ('*39+rs59421388/*39+rs59421388',          '*39/*39'),
        ('*1+*67/*61',              '*61/*1+*67'),
        ('*10+*36/*10+*36+*36',     '*10+*36/*10+*36+*36'),
    ]
    print(f'{"Aldy raw":<28}{"cleaned":<22}{"normalised":<22}{"gold-norm":<22}match')
    print('-' * 100)
    for raw, gold in tests:
        cleaned = clean_aldy_diplotype(raw)
        na = normalize_genotype(cleaned)
        ng = normalize_genotype(gold)
        print(f'{raw:<28}{str(cleaned):<22}{str(na):<22}{str(ng):<22}{"OK" if na == ng else "**MISMATCH**"}')

    print()
    print('=== 實際檔案解析測試 ===')
    p = os.environ['NEW'] + '/aldy_bench/test/NA24631_bwa.cyp2d6.aldy'
    major, n = parse_aldy_file(p)
    cleaned = clean_aldy_diplotype(major)
    print(f'  raw major : {major}')
    print(f'  cleaned   : {cleaned}')
    print(f'  normalised: {normalize_genotype(cleaned)}')
    print(f'  gold      : {normalize_genotype("*49/*10+*36")}')
    print(f'  concordant: {compare_genotype(cleaned, "*49/*10+*36")}')
