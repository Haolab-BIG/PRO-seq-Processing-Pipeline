__author__ = ["Liu Qian"]
__email__ = "18229048546@163.com"
__version__ = "0.0.1"

import os
import sys
import re
import time
import json, pickle
import gzip
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
import fisher
import traceback
import bisect 
import gc
from types import SimpleNamespace
from statsmodels.stats.multitest import multipletests
import statsmodels.stats.contingency_tables as contingency_tables
from scipy.stats import chi2
import warnings
import urllib.request
import hashlib
import tempfile
import logging
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.preprocessing import deseq2_norm_fit
import importlib.util
import getpass

def get_lb(fn):
    lb = os.path.basename(fn).replace('.gz', '').rsplit('.', 1)[0]
    return re.sub(r'[\._\-]*sort(ed)?\b', '', lb)

def check_is_sorted(fn_bed):
    chr_prev = None
    chr_set = set()
    chr_order = []
    with open(fn_bed) as f:
        last = 0
        for i in f:
            chr_, start = i.split('\t', 2)[:2]
            start = int(start)
            if chr_ != chr_prev:
                chr_order.append(chr_)
                if chr_ in chr_set:
                    return False
                chr_set.add(chr_)
                chr_prev = chr_
            elif start < last:
                return False
            last = start
    tmp = sorted(chr_order)
    if tmp == chr_order:
        return True
    return False

def process_input(pwout_raw, fls, respect_sorted=False):
    fls_new = []
    err = 0
    fn = os.path.realpath(fls)
    fls_new.append(fn)
    fls = fls_new
    fls = fls_new
    
    err = 0
    res = []
    lb_map = {}

    fn_lb = get_lb(fn)
    gz_suffix = ''
    fn_out_bed = f'{pwout_raw}/bed/{fn_lb}.sorted.bed'
    fn_out_bed_gz = f'{pwout_raw}/bed/{fn_lb}.sorted.bed{gz_suffix}'

    if os.path.exists(fn_out_bed) and os.path.getsize(fn_out_bed) > 10:
        res.append([fn_lb, fn_out_bed])
    
    fn_for_check = fn
    fn_converted_bed = f'{pwout_raw}/bed/{fn_lb}.converted.bed'
    if os.path.exists(fn_converted_bed) and os.path.getsize(fn_converted_bed) > 1000:
        fn_for_check = fn_converted_bed

    if fn_for_check.endswith('.bed') or fn_for_check.endswith('.bed.gz'):
        is_sorted = check_is_sorted(fn_for_check)

        ires = [fn_lb, fn_out_bed]

    return res

def download_ref_files(home, organism, need_download, ref_files):
    builtin_organisms = ['ce10', 'dm3', 'dm6', 'hg19', 'hg38', 'mm10', 'mm39', 'danrer10']
    url_base = 'https://bioinfo.vanderbilt.edu/NRSA/download/NRSA'
    err = 0
    
    url_md5sum = f'{url_base}/files_md5sum.json'
    
    for key, fn_dest in need_download.items():
        slug = {'gtf': organism, 'fa': 'fa', 'fantom': 'annotation', 'association': 'annotation', '4d': 'annotation'}[key]
        base_name = os.path.basename(fn_dest)
        url_full = f'{url_base}/{slug}/{base_name}'
    if err:
        sys.exit(1)
    return ref_files

def get_ref(home, organism, fn_gtf=None, fn_fa=None):
    ref_files = {}
    pw_code = os.path.dirname(os.path.realpath(__file__))
    pw_code_write = pw_code
    
    if fn_gtf:
        fn_gtf = os.path.abspath(fn_gtf)
    if fn_fa:
        fn_fa = os.path.abspath(fn_fa)

    builtin_organisms = ['ce10', 'dm3', 'dm6', 'hg19', 'hg38', 'mm10', 'mm39', 'danrer10']
    fn_extra_org = f'{pw_code_write}/ref/extra_organisms.txt'

    organism_lower = organism.lower()
    fn_gtf_exp = f'{pw_code}/ref/{organism}/RefSeq-{organism}.gtf'
    fn_fa_exp = f'{pw_code}/fa/{organism}/{organism}.fa'
    fn_gtf_exp_write = f'{pw_code_write}/ref/{organism}/RefSeq-{organism}.gtf'
    fn_fa_exp_write = f'{pw_code_write}/fa/{organism}/{organism}.fa'
    
    fn_protein_coding = f'{pw_code}/ref/{organism}/{organism}.protein_coding.txt'
    fn_protein_coding_write = f'{pw_code_write}/ref/{organism}/{organism}.protein_coding.txt'

    if fn_fa and not os.path.exists(fn_fa_exp_write):
        fn_fa_exp_write = fn_fa
    if fn_gtf and not os.path.exists(fn_gtf_exp_write):
        fn_gtf_exp_write = fn_gtf
    
    ref_files = {}
    need_download = {}
    for k, flist in [['gtf', [fn_gtf, fn_gtf_exp_write]], ['fa', [fn_fa, fn_fa_exp_write]]]:
        found = 0
        for fn in flist:
            if fn is None:
                continue
            if os.path.exists(fn):
                ref_files[k] = fn
                found = 1
                break
        if not found:
            tmp = [(fn, os.path.exists(str(fn))) for fn in flist]
            need_download[k] = flist[1]
        
    ref_files["fantom"], ref_files['association'] = {
        'hg19': ["human_permissive_enhancers_phase_1_and_2.bed", "human_enhancer_tss_associations.bed"],
        'mm10': ['mouse_permissive_enhancers_phase_1_and_2.bed', None],
        'hg38': ['human_permissive_enhancers_phase_1_and_2-hg38.bed', 'human_enhancer_tss_associations-hg38.bed']
        
    }.get(organism, [None, None])

    if organism in ['dm3', 'dm6', 'hg19', 'hg38', 'mm10']:
        ref_files['4d'] = f'data/4DGenome-{organism}.txt'
    else:
        ref_files['4d'] = None
    
    if ref_files['fantom']:
        ref_files["fantom"] = f'data/{ref_files["fantom"]}'
        if ref_files['association']:
            ref_files['association'] = f'data/{ref_files["association"]}'

    for key in ref_files:
        value = ref_files[key]
        if not value:
            continue
        if not os.path.exists(value):
            if os.path.exists(f'{value}.gz'):
                ref_files[key] = f'{value}.gz'
                continue
            fn_new = value.replace(f'{pw_code}/', f'{pw_code_write}/')
            if pw_code != pw_code_write and os.path.exists(fn_new):
                ref_files[key] = fn_new
                continue
            tmp = value.split('/')
            if len(tmp) > 1 and tmp[-2] == 'annotation':
                if organism in {'mm10', 'hg19', 'hg38'}:
                    need_download[key] = value
                ref_files[key] = None
    need_download = {k: v.replace(f'{pw_code}/', f'{pw_code_write}/') for k, v in need_download.items()}
    if need_download:
        folders = {os.path.dirname(fn) for fn in need_download.values()}
        ref_files = download_ref_files(home, organism, need_download, ref_files)

    tmp = json.dumps(ref_files, indent=3)
    return ref_files

class Analysis:
    def __init__(self, home, sample_id, fn_gtf, fa_in,is_long_eRNA=False, skip_get_mapped_reads=False, raw_input=True):
        
        self.bin_size = 200
        
        self.bin_dir = os.path.dirname(os.path.realpath(__file__))
        self.pwout = home or os.getcwd()
        self.pwout_raw = home
        self.pw_bed = f'{home}/bed'
        self.inter_dir = os.path.join(home, 'data')
        self.known_gene_dir = os.path.join(home, 'known_gene')
        self.longerna = is_long_eRNA
        self.skip_get_mapped_reads = skip_get_mapped_reads
        respect_sorted = False

        self.home = home

        self.status = 0

        fa_in = fa_in

        in1 = process_input(self.pwout_raw, home+'/bed/'+sample_id+'.sorted.bed', respect_sorted=respect_sorted)
        if in1 is None:
            self.status = 1
        self.control_bed = in1
        

        self.n_gene_cols = 2
        self.gene_cols = ['Transcript', 'Gene']
        self.longerna_flag_str = ''
        self.longerna_prefix = ''

        ref_fls = get_ref(home, 'hg38', fn_gtf=fn_gtf, fn_fa=fa_in)
        if ref_fls is None:
            self.status = 1
            return
        self.fa = ref_fls['fa']
        self.gtf = ref_fls['gtf']
        if not os.path.exists(self.gtf):
            self.status = 1
        
        if self.status:
            return
        
        self.ref = ref_fls
        self.input_fls = self.control_bed
        
        fn_count_pp_gb = f'{self.inter_dir}/{self.longerna_prefix}count_pp_gb.txt'
        self.out_fls = {
            'bed_peaks': {},
            'fdr': {},
            'count_pp_gb': fn_count_pp_gb
        }
        
        fn_lb_uniq = set()
        fn_lb_dup = {}

        fn_lb, fn_bed= self.input_fls[0][0],self.input_fls[0][1]
        if fn_lb not in fn_lb_uniq:
            fn_lb_uniq.add(fn_lb)
        else:
            fn_lb_dup.setdefault(fn_lb, []).append(fn_bed)
        header_bed_peaks = self.gene_cols + ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
        fn_bed_peaks = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_raw.txt')

        if not (self.skip_get_mapped_reads and os.path.exists(fn_count_pp_gb)):
            fh_bed_peaks = open(fn_bed_peaks, 'w')
            fh_bed_peaks.write('\t'.join(header_bed_peaks) + '\n')
        else:
            fh_bed_peaks = open(fn_bed_peaks, 'r')
        self.out_fls['bed_peaks'][fn_lb] = {'fn': fn_bed_peaks, 'fh': fh_bed_peaks, 'header': header_bed_peaks}
        
        header_fdr = header_bed_peaks + ['pvalue', 'FDR']
        fn_fdr = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_FDR.txt')
        self.out_fls['fdr'][fn_lb] = {'fn': fn_fdr, 'header': header_fdr}
        
        
        if fn_lb_dup:
            tmp = json.dumps(fn_lb_dup, indent=4)
            sys.exit(1)
        
        self.config = {
            'mapped_sites_only': 0,
            # 'pro_up': 500, 
            'pro_up': 1000, 
            # 'pro_down': 500, 
            'pro_down': 1000, 
            'gb_start': 1000, 
            'min_gene_len': 1000, 
            'window_size': 50, 
            'step': 5, 
            'tts_padding': 2000
        }

def build_tss(gtf_info, fn_tss, fn_tss_tts):
    with open(fn_tss, 'w') as f, open(fn_tss_tts, 'w') as o2:
        for k, v in gtf_info.items():
            itss, itts = [v['start'], v['end']] if v['strand'] == '+' else [v['end'], v['start']]
            
            f.write(f'{v["chr"]}\t{itss}\t{itss}\t{k}\t{v["strand"]}\n')
            o2.write(f'{v["chr"]}\t{itss}\t{itts}\t{k}\t{v["strand"]}\n')

def process_gtf(home, fn_gtf, pwout, force_rebuild=False, fake_gtf_path=None):
    custom_temp = f"{pwout}/temp"
    os.makedirs(custom_temp, exist_ok=True)
    tempfile.tempdir = custom_temp
    os.environ['TMPDIR'] = custom_temp

    err = {'no_transcript_id': 0, 'no_gene_name': 0, 'invalid_line_format': 0, 'invalid_strand': 0}
    fn_gtf = os.path.realpath(fn_gtf)
    fn_gtf_lb = os.path.basename(fn_gtf).replace('.gtf', '')
    
    fn_gtf_size = os.path.getsize(fn_gtf)
    fn_gtf_lb = f'{fn_gtf_lb}_size_{fn_gtf_size}'
    
    fn_tss = f'{pwout}/intermediate/{fn_gtf_lb}.tss.txt'
    fn_tss_tts = f'{pwout}/intermediate/{fn_gtf_lb}.tss_tts.txt'
    
    
    gtf_pwout = os.path.dirname(fn_gtf)
    
    gtf_col_idx = {
        'chr': 0,
        'source': 1,
        'feature': 2,
        'start': 3,
        'end': 4,
        'score': 5,
        'strand': 6,
        'attribute': 8,
    }
    res_raw = {}  
    
    last_exon = {}
    with open(fn_gtf) as f:
        skipped_transcript_on_random_chr = 0
        for i in f:
            line = i.strip().split('\t')
            line_err = None
            try:
                chr_, region_type, start, end, strand, attribute_raw = [line[_] for _ in [gtf_col_idx[_] for _ in ['chr', 'feature', 'start', 'end', 'strand', 'attribute']]]
                start = int(start)
                end = int(end)
            except:
                line_err = 'invalid_line_format'
            if region_type != 'exon':
                continue
            if strand not in {'+', '-'}:
                line_err = 'invalid_strand'

            m = re.search(r'transcript_id\s+"?([\w._\-]+)', attribute_raw)
            if m:
                transcript_id = m.group(1)
            else:
                line_err = 'no_transcript_id'
                
            m = re.search(r'gene_name\s+"?([\w._\-]+)', attribute_raw)
            if m:
                gene_name = m.group(1)
            else:
                line_err = 'gene_name_not_found'
            
            if line_err:
                tmp = err.setdefault(line_err, 0)
                err[line_err] += 1
                
                continue

            ires_exon = {'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end}
            ires = res_raw.setdefault(gene_name, {}).setdefault(transcript_id, ires_exon)
            if strand == '+':
                last_exon[transcript_id] = (start, end)
            elif transcript_id not in last_exon:
                last_exon[transcript_id] = (start, end)

            if chr_ != ires['chr'] or strand != ires['strand']:
                transcript_id_new = f'{transcript_id}@{chr_}@{strand}'
                ires = res_raw.setdefault(gene_name, {}).setdefault(transcript_id_new, ires_exon)
                if strand == '+':
                    last_exon[transcript_id_new] = (start, end)
                elif transcript_id_new not in last_exon:
                    last_exon[transcript_id_new] = (start, end)

            if start < ires['start']:
                ires['start'] = start
            if end > ires['end']:
                ires['end'] = end

    res = {}
    n_merged = 0
    meta = {'initial_n_transcripts': 0, 'n_merged': 0, 'final_n_transcripts': 0, 'skipped_transcript_on_random_chr': skipped_transcript_on_random_chr}

    for gn, v1 in res_raw.items():
        tmp = {}
        for transcript_id, v2 in v1.items():
            start, end, strand = v2['start'], v2['end'], v2['strand']
            unique_id = f'{v2["chr"]}_{start}_{end}_{strand}'
            v2['tss'] = start if strand == '+' else end
            v2['tts'] = end if strand == '+' else start
            tmp.setdefault(unique_id, []).append(transcript_id)
            meta['initial_n_transcripts'] += 1
        for transcript_list in tmp.values():
            if len(transcript_list) > 1:
                n_merged += len(transcript_list) - 1
                
            transcript_id_new = ';'.join(transcript_list)
            transcript_id_demo = transcript_list[0]
            res[transcript_id_new] = v1[transcript_id_demo]
            res[transcript_id_new]['last_exon'] = last_exon[transcript_id_demo]
    meta['n_merged'] = n_merged
    meta['final_n_transcripts'] = len(res)

    if pwout is not None:
        build_tss(res, fn_tss, fn_tss_tts)
    
    return res, fn_tss, fn_tss_tts, err

def add_value_to_gtf(gene_info, pro_up, pro_down, gb_down_distance, tts_padding, tts_down_length=5000, islongerna=False):
    strand = gene_info['strand']
    gene_raw_s, gene_raw_e = gene_info['start'], gene_info['end']
    if not islongerna:
        last_exon_s, last_exon_e = gene_info['last_exon']
    else:
        last_exon_s, last_exon_e = None, None
    if strand == '+':
        pp_start = gene_raw_s - pro_up
        pp_end = gene_raw_s + pro_down - 1
        gb_start = gene_raw_s + gb_down_distance
        gb_end = gene_raw_e
        tts = gene_raw_e
        strand_idx = 0
        tts_start = gene_raw_e
        tts_end = gene_raw_e + tts_padding
        tts_down_region_s, tts_down_region_e = gene_raw_e, gene_raw_e + tts_down_length
    else:
        pp_start = gene_raw_e - (pro_down - 1)
        pp_end = gene_raw_e + pro_up
        gb_start = gene_raw_s
        tts = gene_raw_s
        gb_end = gene_raw_e - gb_down_distance
        strand_idx = 1
        tts_start = gene_raw_s - tts_padding
        tts_end = gene_raw_s
        tts_down_region_s, tts_down_region_e = gene_raw_s - tts_down_length, gene_raw_s
    
    gene_seq = None
    gb_seq_N = 0
    
    gb_len_mappable = gb_end - gb_start + 1 - gb_seq_N

    new_info = dict(zip(
        ['pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq', 'tts_start', 'tts_end', 'last_exon_s', 'last_exon_e', 'tts_down_region_s', 'tts_down_region_e'], 
        [ pp_start,   pp_end,   gb_start ,  gb_end ,  strand_idx ,  gb_len_mappable,   gene_seq, tts_start, tts_end, last_exon_s, last_exon_e, tts_down_region_s, tts_down_region_e]
        ))
    gene_info.update(new_info)
    return gene_info

def pre_count_for_bed(fn_lb, fn_out_bed, pw_bed, bin_size=200, reuse=True):
    fn_n_lines = f'{pw_bed}/{fn_lb}.line_count.txt'
    
    fn_pre_count_list_json = f'{pw_bed}/{fn_lb}.pre_count_list.json'

    file_list = {}
    chr_map = {}
    n_lines = 0
    n_error = 0
    single_col = 0
    fatal = 0
    count_bin = {'bin_size': bin_size}
    count_per_base = {}
    chr_prev = None
    
    
    def dump_res(chr_prev, count_per_base, count_bin):
        chr_prev_new = chr_prev
        chr_map[chr_prev_new] = chr_prev
        fn_count_bin = f'{pw_bed}/{fn_lb}.count.bin_of_{bin_size}_chr{chr_prev_new}.pkl'
        fn_count_per_base = f'{pw_bed}/{fn_lb}.count.per_base_chr{chr_prev_new}.pkl'
        with open(fn_count_per_base, 'wb') as f:
            pickle.dump(count_per_base, f)
        with open(fn_count_bin, 'wb') as f:
            pickle.dump(count_bin, f)
        file_list[chr_prev_new] = (fn_count_per_base, fn_count_bin)
        count_bin = {'bin_size': bin_size}
        count_per_base = {}
        return count_per_base, count_bin
    
    
    with open(fn_out_bed, 'r') as fh_bed:
        for i in fh_bed:
            try:
                chr_, s, e, _, _, strand, *_ = i[:-1].split('\t', 6)
            except:
                nchr, ncol = len(i), len(i.split('\t'))
                if ncol == 1:
                    single_col += 1
                else:
                    n_error += 1
                    if n_error > 100:
                        fatal = 1
                        break
                continue
            strand_idx, read_end = (0, int(e)) if strand == '+' else (1, int(s) + 1)
            if chr_prev and chr_ != chr_prev:
                count_per_base, count_bin = dump_res(chr_prev, count_per_base, count_bin)


            chunk = read_end // bin_size
            count_bin.setdefault(chunk, [0, 0])
            count_bin[chunk][strand_idx] += 1

            count_per_base.setdefault(read_end, [0, 0])
            count_per_base[read_end][strand_idx] += 1
            n_lines += 1
            chr_prev = chr_
    if chr_prev:
        count_per_base, count_bin = dump_res(chr_prev, count_per_base, count_bin)

    if n_error:
        prefix = 'more than ' if fatal else ''
        sys.exit(1)

    with open(fn_n_lines, 'w') as f:
        f.write(f'{n_lines}\n')
    with open(fn_pre_count_list_json, 'w') as o:
        json.dump(file_list, o, indent=3)
    return file_list

def get_peak_method2(count_per_base, count_bin, chr_, strand_idx, s, e, bin_size):
    bin_start = s // bin_size
    bin_end = e // bin_size
    if bin_end - bin_start < 2:
        return sum([count_per_base.get(i, [0, 0])[strand_idx] for i in range(s, e + 1)])
    points = []
    left_mod = s % bin_size
    right_mod = e % bin_size
    if left_mod:
        points += range(s, (bin_start + 1) * bin_size)
    if right_mod != bin_size - 1:
        bin_list = range(bin_start + 1, bin_end)
        points += range(bin_end * bin_size, e + 1)
    else:
        bin_list = range(bin_start + 1, bin_end + 1)
    return sum([count_per_base.get(i, [0, 0])[strand_idx] for i in points]) + sum([count_bin.get(i, [0, 0])[strand_idx] for i in bin_list])


def get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, tts_start, tts_end, gb_len_mappable, gene_seq, window_size, step_size, bin_size, prev_pp_peak):
    gbc = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, gb_start, gb_end, bin_size=bin_size)
    if gbc == 0:
        gbd = 0
    else:
        gbd = gbc / gb_len_mappable

    tts_ct = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, tts_start, tts_end, bin_size=bin_size)
    
    pp_region_count = [count_per_base.get(i, [0, 0])[strand_idx] for i in range(pp_start, pp_end + 1)]
    window_ct = [(i, sum(pp_region_count[i: i + window_size])) for i in range(0, pp_end - pp_start - window_size + 2, step_size)]
    peak_window_start, ppc  = max(window_ct, key=lambda _: _[1]) 
    peak_window_ct = [(i, pp_region_count[i]) for i in range(peak_window_start, peak_window_start + window_size)]
    summit_pos, summit_count = max(peak_window_ct, key=lambda _: _[1])
    summit_pos += pp_start
    summit_pos_str = f'{chr_}{strand}:{summit_pos}'
    
    if ppc == 0:
        ppd = 0
    elif gene_seq:
        peak_windows_start_abs = peak_window_start + pp_start - gene_raw_s
        mappable_sites = window_size - gene_seq[peak_windows_start_abs: peak_windows_start_abs + window_size].count('N')
        ppd = ppc / mappable_sites
    else:
        ppd = ppc / window_size

    pp_res = {'ppc': ppc, 'ppd': ppd, 'mappable_sites': window_size, 'summit_pos': summit_pos_str, 'summit_count': summit_count}
    return gbc, gbd,  pp_res, tts_ct

def process_bed_files(analysis, fls, gtf_info, gtf_info_raw, fh_fa, reuse_pre_count=False, save_tts_count=True, islongerna=False, downstream_no_overlap_length=50000):
    invalid_chr_transcript = 0
    bin_size = analysis.bin_size
    fail_to_retrieve_seq = 0
    pwout = analysis.pwout
    window_size, step_size = analysis.config['window_size'], analysis.config['step']
    pw_bed = analysis.pw_bed
    n = 0
    prev = [time.time(), 0]
    section_size = 1000

    s_before_loop = time.time()
    s = s_before_loop
    ts_list = list(gtf_info)

    pp_str = {ts: [] for ts in ts_list}
    gb_str = {ts: [] for ts in ts_list}
    tts_str = {ts: [str(gtf_info[ts][_]) for _ in ['chr', 'strand', 'tts']] for ts in ts_list}
    tts_down_str = {}
    downstream_no_overlap_length_str = f'{downstream_no_overlap_length/1000:.0f}k'
    
    gtf_info_new = {}
    
    for ts, v in gtf_info.items():
        chr_ = v['chr']
        gtf_info_new.setdefault(chr_, {})[ts] = v

    invlude_all_ts = True
    ts_without_overlap = set(gtf_info_raw)
    n_ts_init = len(ts_list)
    err = 0
    fn_lb, fn_bed = fls[0][0],fls[0][1]
    fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']

    pre_count_flist = pre_count_for_bed(fn_lb, fn_bed, pw_bed, bin_size, reuse=reuse_pre_count)
    
    prev_peak = {}
    valid_chr = set(pre_count_flist)
    chr_excluded = {}
    ts_excluded = 0
    
    for ts, v in gtf_info.items():
        ts_chr = v['chr']
        if ts_chr not in valid_chr:
            ts_excluded += 1
            gtf_info_new[ts_chr].pop(ts)
            if ts in pp_str:
                del pp_str[ts]
            if ts in gb_str:
                del gb_str[ts]
            if ts in tts_str:
                del tts_str[ts]
            chr_excluded.setdefault(ts_chr, 0)
            chr_excluded[ts_chr] += 1

    n_current = n_ts_init - ts_excluded
    if n_current == 0:
        size = os.path.getsize(fn_bed)
        err = 1

    for chr_, ts_list_chr in gtf_info_new.items():
        if chr_ in pre_count_flist:
            fn_count_per_base, fn_count_bin = pre_count_flist[chr_]
            with open(fn_count_per_base, 'rb') as f:
                count_per_base = pickle.load(f)
            with open(fn_count_bin, 'rb') as f:
                count_bin = pickle.load(f)
            bin_size = count_bin['bin_size']
        
            for transcript_id in ts_list_chr:
                gene_info = gtf_info[transcript_id]
                chr_, strand, gene_raw_s, pp_start, pp_end, gb_start, gb_end, strand_idx, gb_len_mappable, gene_seq, tts_start, tts_end, tts_down_region_s, tts_down_region_e = [gene_info[_] for _ in ['chr', 'strand', 'start', 'pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq', 'tts_start', 'tts_end', 'tts_down_region_s', 'tts_down_region_e']]

                n += 1
                if n % section_size == 0:
                    now = time.time()
                    prev_time, prev_count = prev
                    time_gap = now - prev_time
                    speed = time_gap * 1000  / (n - prev_count)
                    prev = [now, n]

                s = time.time()
                gbc, gbd, pp_res, tts_ct = get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, tts_start, tts_end, gb_len_mappable, gene_seq,  window_size, step_size, bin_size, prev_peak)
                if islongerna == False and (invlude_all_ts or transcript_id in ts_without_overlap):
                    last_exon_s, last_exon_e = gene_info['last_exon']
                    last_exon_ct = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, last_exon_s, last_exon_e, bin_size)
                    tts_down_ct = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, tts_down_region_s, tts_down_region_e, bin_size)
                    last_exon_len = last_exon_e - last_exon_s + 1
                    if transcript_id not in tts_down_str:
                        tts_down_str[transcript_id] = [transcript_id, gene_info['gene_name'], chr_, strand, str(last_exon_len), str(last_exon_s), str(last_exon_e)]
                    tts_down_str[transcript_id] += [str(last_exon_ct), str(tts_down_ct), f'{tts_down_ct / last_exon_ct:.4f}' if last_exon_ct > 0 else 'NA']

                pp_str[transcript_id].append(str(pp_res['ppc']))
                gb_str[transcript_id] += [str(gbc), str(gbd)]
                tts_str[transcript_id].append(str(tts_ct))
                
                row = [transcript_id]
                if not analysis.longerna:
                    row.append(gene_info['gene_name'])
                row += [pp_res['ppc'], pp_res['mappable_sites'], pp_res['ppd'], pp_res['summit_pos'], gbc, gb_len_mappable, gbd]
                pro_vs_pb = round(pp_res['ppd'] / gbd, 4) if gbd else 'NA'
                row.append(pro_vs_pb)
                
                print('\t'.join(map(str, row)), file=fh_bed_peaks)
            del count_per_base
            del count_bin 
            gc.collect()
    
    if err:
        sys.exit(1)

    if islongerna == False:
        fn_tts_down = f'{pwout}/intermediate/tts_down_{downstream_no_overlap_length_str}.txt'
        with open(fn_tts_down, 'w') as o:
            header = ['Transcript', 'Gene', 'chr', 'strand', 'last_exon_len', 'last_exon_s', 'last_exon_e']
            fn_lb, _ = fls[0][0],fls[0][1]
            header += [f'last_exon_{fn_lb}', f'tts_down_{downstream_no_overlap_length_str}_{fn_lb}', f'ratio_{fn_lb}']
            print('\t'.join(header), file=o)
            for ts, v in tts_down_str.items():
                print('\t'.join(v), file=o)

    if save_tts_count:
        fn_tts_count = f'{pwout}/intermediate/count_tts.raw.txt'
        fn_tts_count_filtered = f'{pwout}/intermediate/count_tts.filtered.txt'
        tts_file_header = ['Transcript', 'Gene', 'chr', 'strand', 'TTS'] + [fls[0][0]]

        tss_list = {}
        gene_regions = {}
        for gene_info in gtf_info.values():
            chr_ = gene_info['chr']
            tss_list.setdefault(chr_, set()).add(gene_info['tss'])
            gene_regions.setdefault(chr_, []).append([gene_info['start'], gene_info['end'], gene_info['gene_name']])
        tss_list = {k: sorted(v) for k, v in tss_list.items()}
        gene_regions = {k: sorted(v, key=lambda x: x[0]) for k, v in gene_regions.items()}
        gene_regions_start = {k: [_[0] for _ in v] for k, v in gene_regions.items()}
        gene_regions_max_len = {k: max([_[1] - _[0] + 1 for _ in v]) for k, v in gene_regions.items()}

        tts_dowstream_no_tss_region_len = 2000
        n_excluded = 0
        
        with open(fn_tts_count, 'w') as o, open(fn_tts_count_filtered, 'w') as o1:
            print('\t'.join(tts_file_header), file=o)
            print('\t'.join(tts_file_header), file=o1)
            for ts in ts_list:
                gene_info = gtf_info[ts]
                gn = gene_info['gene_name']
                if ts in tts_str:
                    line = '\t'.join([ts, gn] + tts_str[ts])
                    print(line, file=o)
                    
                    
                    tts = gene_info['tts']
                    strand = gene_info['strand']
                    chr_ = gene_info['chr']
                    tss_list_chr = tss_list[chr_]
                    tts_region_s, tts_region_e = [tts, tts + tts_dowstream_no_tss_region_len] if strand == '+' else [tts - tts_dowstream_no_tss_region_len, tts]
                    left_idx =  bisect.bisect_left(tss_list_chr, tts_region_s)
                    right_idx = bisect.bisect_right(tss_list_chr, tts_region_e)
                        
                    if left_idx < right_idx:
                        n_excluded += 1
                        continue
                    
                    idx_tts_in_start_pos = bisect.bisect_right(gene_regions_start[chr_], tts)
                    excluded = False
                    gene_regions_max_len_chr = gene_regions_max_len[chr_]
                    
                    for iregion in gene_regions[chr_][:idx_tts_in_start_pos][::-1]:
                        if iregion[1] >= tts:
                            if gn == iregion[2]:
                                continue
                            n_excluded += 1
                            excluded = 1
                            break
                        if tts - iregion[0] > gene_regions_max_len_chr:
                            break
                    if excluded:
                        continue
                    print(line, file=o1)
        
    return pp_str, gb_str

if __name__=="__main__":
    home = sys.argv[1]
    sample_id = sys.argv[2]
    fn_gtf = sys.argv[3]
    fa_in = sys.argv[4]
    #工作目录
    analysis = Analysis(home, sample_id, fn_gtf, fa_in)
    if analysis.status:
        sys.exit(1)
    rep1 = len(analysis.control_bed)
    window_size = analysis.config['window_size']
    pro_up, pro_down, gb_down_distance, min_gene_len, tts_padding = [analysis.config[_] for _ in ['pro_up', 'pro_down', 'gb_start', 'min_gene_len', 'tts_padding']]
    fn_fa = analysis.ref['fa']
    fh_fa = open(fn_fa, 'r')
    gtf_info, fn_tss, fn_tss_tts, err = process_gtf(home, analysis.ref['gtf'], pwout=analysis.pwout)
    gtf_info_raw = gtf_info
    err_total = sum(err.values())
    len1 = len(gtf_info)
    gtf_info_new = {}
    short_genes = 0
    merged_transcripts = 0
    tmp = sorted(gtf_info)
    skipped_ts = {'invalid_chr': 0,  'short': 0}
    for k in tmp:
        v = gtf_info[k]
        if v['end'] - v['start'] + 1 > min_gene_len:
            merged_transcripts += k.count(';')
            gtf_info_new[k] = add_value_to_gtf(v, pro_up, pro_down, gb_down_distance, tts_padding)
        else:
            skipped_ts['short'] += 1
            short_genes += 1
    total_skipped = sum(skipped_ts.values())
    len1 = len(gtf_info_new)
    gtf_info = gtf_info_new
    fls = analysis.input_fls
    sam_order = [sample_id]
    n_gene_cols = analysis.n_gene_cols
    fn_count_pp_gb = analysis.out_fls['count_pp_gb']
    line_count = 0
    tts_down = 50000
    pp_str, gb_str = process_bed_files(analysis, fls, gtf_info, gtf_info_raw, fh_fa, reuse_pre_count=False, downstream_no_overlap_length=tts_down)
    fn_lb, fn_bed=fls[0][0],fls[0][1]
    analysis.out_fls['bed_peaks'][fn_lb]['fh'].close()
    fh_fa.close()