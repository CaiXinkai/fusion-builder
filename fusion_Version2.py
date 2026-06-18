from pyfaidx import Fasta
from Bio.Seq import Seq

def parse_gtf(gtf_path):
    """
    解析 GTF，返回 gene -> transcript -> info 字典:
    { gene_name: { transcript_id: {'chr':..., 'strand':..., 'exons': [(start,end), ...]} } }
    注: 期望 GTF 的 attributes 包含 gene_name 和 transcript_id。
    """
    import gzip, io
    gene2transcripts = {}
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            chrom, src, feature, start, end, score, strand, frame, attr = cols
            if feature != 'exon':
                continue
            start = int(start); end = int(end)
            attrs = {}
            for kv in attr.strip().split(';'):
                kv = kv.strip()
                if not kv:
                    continue
                if ' ' in kv:
                    k, v = kv.split(' ', 1)
                    attrs[k] = v.strip().strip('"')
            gene_name = attrs.get('gene_name') or attrs.get('gene_id')
            tx_id = attrs.get('transcript_id')
            if gene_name is None or tx_id is None:
                continue
            gene2transcripts.setdefault(gene_name, {})
            tdict = gene2transcripts[gene_name].setdefault(tx_id, {'chr': chrom, 'strand': strand, 'exons': []})
            tdict['exons'].append((start, end))
    # 排序外显子
    for g, txs in gene2transcripts.items():
        for tx, info in txs.items():
            info['exons'].sort(key=lambda x: x[0])
    return gene2transcripts

def pick_canonical_transcript(transcripts_dict):
    """
    从 transcripts_dict 中挑选总外显子长度最长的 transcript_id
    """
    best = None; best_len = -1
    for tx, info in transcripts_dict.items():
        total = sum(e - s + 1 for s, e in info['exons'])
        if total > best_len:
            best = tx; best_len = total
    return best

def transcript_sequence(fasta: Fasta, chrom, exons, strand):
    """
    exons: list of (start,end) 排序为基因组升序
    返回转录本方向（5'->3'）的核苷酸序列（字符串）
    pyfaidx 使用 1-based 包含端点索引
    """
    seq_pieces = []
    for s, e in exons:
        piece = fasta[chrom][s-1:e].seq
        seq_pieces.append(piece)
    seq = ''.join(seq_pieces)
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    return seq

def split_transcript_at_genomic_pos(exons, strand, breakpoint_pos):
    """
    将 transcript 的外显子列表按基因组坐标 breakpoint_pos 划分为左右两部分（基因组升序）
    返回 (left_parts, right_parts)，每项为 exons 子列表（升序）
    解释：left_parts 表示从 transcript 起始到 breakpoint（包含 breakpoint 所在碱基），
           right_parts 表示 breakpoint+1 到 transcript 结束。
    """
    left_parts = []
    right_parts = []
    # 查找 breakpoint 落在哪个 exon
    found = False
    for idx, (s, e) in enumerate(exons):
        if s <= breakpoint_pos <= e:
            found = True
            # breakpoint 在该 exon 内，分割该 exon
            if s <= breakpoint_pos:
                left_parts.extend(exons[:idx])  # 之前全部属于左
                left_parts.append((s, breakpoint_pos))
            if breakpoint_pos + 1 <= e:
                right_parts.append((breakpoint_pos + 1, e))
                right_parts.extend(exons[idx+1:])  # 之后全部属于右
            else:
                right_parts.extend(exons[idx+1:])
            break
    if not found:
        # breakpoint 不在任何 exon（可能在内含子或转录本外）
        if breakpoint_pos < exons[0][0]:
            left_parts = []
            right_parts = exons.copy()
        elif breakpoint_pos > exons[-1][1]:
            left_parts = exons.copy()
            right_parts = []
        else:
            # breakpoint 在内含子：找到上游 exon 与下游 exon 的索引
            up = None
            down = None
            for i, (s, e) in enumerate(exons):
                if e < breakpoint_pos:
                    up = i
                elif s > breakpoint_pos and down is None:
                    down = i
            left_parts = exons[: (up + 1) if up is not None else 0]
            right_parts = exons[down:] if down is not None else []
    # 保证升序
    left_parts.sort(key=lambda x: x[0])
    right_parts.sort(key=lambda x: x[0])
    return left_parts, right_parts

def build_fusion(fasta_path, gtf_dict, left_gene, right_gene, left_breakpoint, right_breakpoint):
    """
    构造融合转录本与蛋白序列。
    left_breakpoint, right_breakpoint: 字符串 'chr:pos' 或 (chr,pos) 元组
    返回 dict 包含: left_gene,right_gene,left_tx,right_tx,fusion_nt,fusion_protein,notes
    """
    fasta = Fasta(fasta_path)
    def parse_bp(bp):
        if isinstance(bp, tuple):
            return bp
        if ':' in bp:
            chrom, pos = bp.split(':')
            return chrom, int(pos)
        raise ValueError("breakpoint format must be chr:pos")

    lchr, lpos = parse_bp(left_breakpoint)
    rchr, rpos = parse_bp(right_breakpoint)
    notes = []

    if left_gene not in gtf_dict:
        raise KeyError(f"left gene {left_gene} not found in GTF")
    if right_gene not in gtf_dict:
        raise KeyError(f"right gene {right_gene} not found in GTF")

    ltx = pick_canonical_transcript(gtf_dict[left_gene])
    rtx = pick_canonical_transcript(gtf_dict[right_gene])
    linfo = gtf_dict[left_gene][ltx]; rinfo = gtf_dict[right_gene][rtx]

    if linfo['chr'] != lchr:
        notes.append(f"Warning: left gene {left_gene} annotated on {linfo['chr']} but breakpoint uses {lchr}")
    if rinfo['chr'] != rchr:
        notes.append(f"Warning: right gene {right_gene} annotated on {rinfo['chr']} but breakpoint uses {rchr}")

    l_left_exons, l_right_exons = split_transcript_at_genomic_pos(linfo['exons'], linfo['strand'], lpos)
    r_left_exons, r_right_exons = split_transcript_at_genomic_pos(rinfo['exons'], rinfo['strand'], rpos)

    # 左基因取从转录本起点到断点（l_left_exons），右基因取断点之后到转录本末端（r_right_exons）
    left_seq_nt = transcript_sequence(fasta, linfo['chr'], l_left_exons, linfo['strand']) if l_left_exons else ''
    right_seq_nt = transcript_sequence(fasta, rinfo['chr'], r_right_exons, rinfo['strand']) if r_right_exons else ''

    if not left_seq_nt:
        notes.append("Left transcript portion is empty (breakpoint upstream of transcript start).")
    if not right_seq_nt:
        notes.append("Right transcript portion is empty (breakpoint downstream of transcript end).")

    fusion_nt = left_seq_nt + '/' + right_seq_nt

    concat = left_seq_nt + right_seq_nt
    seq_obj = Seq(concat)
    protein = str(seq_obj.translate(to_stop=False))
    left_nt_len = len(left_seq_nt)
    aa_left_full = left_nt_len // 3
    rem = left_nt_len % 3
    if concat == '':
        protein_with_slash = '/'
    else:
        # 在蛋白序列中插入 '/'
        protein_with_slash = protein[:aa_left_full] + '/' + protein[aa_left_full:]
        if rem != 0:
            notes.append(f"Frameshift: left part length {left_nt_len} not divisible by 3; codon split across junction.")

    return {
        'left_gene': left_gene,
        'right_gene': right_gene,
        'left_tx': ltx,
        'right_tx': rtx,
        'fusion_nt': fusion_nt,
        'fusion_protein': protein_with_slash,
        'notes': notes
    }