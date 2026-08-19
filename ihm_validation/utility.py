# -*- coding: utf-8 -*-
#
# utility.py - Various helper functions
#
# Copyright (C) 2019-2025 Arthur Zalevsky, Sai Ganesan, Benjamin M. Webb, Brinda Vallat
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Various helper functions
"""

import os
import glob
from pathlib import Path
from collections import Counter, defaultdict
from multiprocessing import Process
import numpy as np
import logging
import ihm, ihm.reader, ihm.model
import itertools
import time
import signal
import re
import requests
import base64
import json
import pypdf
import bokeh
from datetime import datetime, timezone
from bokeh.models import HoverTool

NA = 'Not available'

# Version of IHMValidation itself. It lives here rather than in report.py so
# that the assessment modules can stamp it into their caches without importing
# report, which imports them.
IHMV_VERSION = '3.2'

# Bump when the shape of a cache entry changes in a way older readers cannot
# understand.
CACHE_FORMAT = 1


# Directory caches cannot be wrapped the way a pickle can, so they carry their
# provenance in a sidecar of this name instead.
CACHE_METADATA_FILE = 'ihmv_cache.json'


def cache_metadata(software: dict = None) -> dict:
    """
    Provenance for a cache entry.

    Records when it was written, which IHMValidation wrote it, and the versions
    of any external tool whose output is baked into it - a MolProbity result is
    only meaningful next to the MolProbity that produced it.
    """
    return {
        'cache_format': CACHE_FORMAT,
        'created': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'ihmv_version': IHMV_VERSION,
        'software': dict(software or {}),
    }


def wrap_cache(data, software: dict = None) -> dict:
    """Bundle cached data together with its provenance, ready to pickle."""
    return {**cache_metadata(software), 'data': data}


def write_cache_metadata(dirname, software: dict = None) -> None:
    """Record provenance for a cache that is a directory rather than a pickle."""
    try:
        with open(Path(dirname, CACHE_METADATA_FILE), 'w') as f:
            json.dump(cache_metadata(software), f, indent=2)
    except OSError as e:
        logging.error(f'Could not write cache metadata to {dirname}: {e}')


def read_cache_metadata(dirname) -> dict:
    """Provenance for a directory cache, or None if it predates this."""
    path = Path(dirname, CACHE_METADATA_FILE)

    if not path.is_file():
        return None

    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError) as e:
        logging.error(f'Could not read cache metadata from {path}: {e}')
        return None


def unwrap_cache(raw):
    """
    Split a loaded cache entry into (data, metadata).

    Entries written before this metadata existed are bare payloads. They are
    still good data - some cost a 280 MB download or a long MolProbity run -
    so they come back with metadata None instead of being thrown away.
    """
    if isinstance(raw, dict) and 'cache_format' in raw and 'data' in raw:
        return raw['data'], {k: v for k, v in raw.items() if k != 'data'}

    return raw, None


def describe_cache(meta: dict) -> str:
    """One line of provenance for logging a cache hit."""
    if not meta:
        return 'no metadata, written before it was recorded'

    bits = [f"IHMValidation {meta.get('ihmv_version', NA)}",
            meta.get('created', NA)]
    bits += [f'{k} {v}' for k, v in sorted((meta.get('software') or {}).items())]

    return ', '.join(bits)

# The HTML report embeds plots as bokeh JSON and hands them to BokehJS loaded
# from the CDN, so the two have to be the same version - a 3.x document is not
# readable by 2.x BokehJS, and the plots silently come out blank. Take the
# version from the installed package rather than pinning it in the template,
# which is how they drifted apart in the first place.
BOKEH_VERSION = bokeh.__version__


def add_hover(plot, renderers, tooltips, mode: str = 'mouse') -> None:
    """
    Attach a HoverTool to `plot` showing `tooltips` for `renderers` only.

    `tooltips` is the usual bokeh list of (label, field) pairs, e.g.
    ``[('q', '@Q{0.0000}'), ('Log I(q)', '@logI{0.000}')]``.

    The tool is always tied to explicit renderers, because most plots in the
    report carry helper glyphs - error bars, fit lines, threshold markers -
    that have no values worth showing and would otherwise steal the hover.

    Use ``mode='vline'`` for line plots, where the pointer is rarely close
    enough to a vertex for the default 'mouse' mode to trigger.

    Tooltips only ever reach the HTML report; the PDF embeds SVG exports, and
    a HoverTool leaves those byte-for-byte identical.
    """
    if not isinstance(renderers, (list, tuple)):
        renderers = [renderers]

    plot.add_tools(HoverTool(renderers=list(renderers),
                             tooltips=list(tooltips),
                             mode=mode))


# Bokeh labels plots with font-family="helvetica", which is an alias rather than
# an installed font. The browser that measures the text to lay an axis out and
# the wkhtmltopdf that finally renders it resolve that alias independently, so
# text ends up narrower than the space reserved for it and right-aligned labels
# come out ragged. Name the font helvetica resolves to anyway: same typeface,
# but both sides now agree.
PLOT_FONT = 'Arial'


def set_plot_font(root, font: str = PLOT_FONT) -> None:
    """
    Name the font explicitly on every piece of text in a plot or layout.

    Walks the model tree rather than the caller listing properties, so it
    catches titles, axis and tick labels, legends, categorical group labels and
    anything iqplot built for us. Properties ending in `text_font_size` or
    `text_font_style` do not match and are left alone.
    """
    for obj in root.references():
        for prop in obj.properties():
            if prop.endswith('text_font'):
                try:
                    setattr(obj, prop, font)
                except (ValueError, AttributeError) as e:
                    logging.error(f'Could not set {prop} on {obj}: {e}')


def clear_legend_background(plot) -> None:
    """
    Stop a plot's legend painting itself onto an opaque panel.

    Every figure in the report is given a transparent background so the PDF
    watermark shows through, but a bokeh legend is not covered by that: it
    defaults to a white fill at 0.95 alpha and punches a white rectangle back
    through the page.

    Takes the plot rather than the legend so it also covers the plots that
    build their legend implicitly from `legend_label`. Safe on a plot with no
    legend at all - bokeh's splattable accessor makes it a no-op.
    """
    plot.legend.background_fill_color = None

# Bokeh 3.x emits a duplicate <text> element per label with stroke-opacity="0"
# (intended as an invisible outline). wkhtmltopdf ignores stroke-opacity and
# renders these as solid black outlines around every plot label.
_BOKEH_INVISIBLE_STROKE_TEXT_RE = re.compile(
    r'<text [^>]*stroke-opacity="0"[^>]*>[^<]*</text>'
)

# Fully transparent glyphs - we use them as invisible hover targets - are still
# exported, one empty <path> each. They have no geometry to draw, so they only
# add bulk to the SVG that ends up in the PDF.
_BOKEH_EMPTY_PATH_RE = re.compile(
    r'<path fill="none" stroke="none"/>'
)


def strip_bokeh_svg_noise(svg_path) -> None:
    """
    Drop elements bokeh's SVG export emits that cannot draw anything.

    Two kinds: the ghost text bokeh writes behind every label (which
    wkhtmltopdf turns into a black outline because it ignores stroke-opacity),
    and the empty paths left behind by fully transparent glyphs.
    """
    path = Path(svg_path)
    if not path.is_file():
        return
    text = path.read_text(encoding='utf-8')
    cleaned = _BOKEH_INVISIBLE_STROKE_TEXT_RE.sub('', text)
    cleaned = _BOKEH_EMPTY_PATH_RE.sub('', cleaned)
    if cleaned != text:
        path.write_text(cleaned, encoding='utf-8')


def dict_to_JSlist(d: dict) -> list:
    '''
    convert dictionary to list of lists
    '''
    output_list = []

    if bool(d) and len(list(d.keys())) > 0:
        # add headers for table, which are the keys of the dict
        header = list(d.keys())
        # Get number of columns
        N = len(header)
        # Get number of rows. +1 for header.
        M = len(d[header[0]]) + 1
        output_list = np.empty((M, N), dtype=object)
        output_list[0, :] = header
        # iterate over dict keys - columns
        for j, v in enumerate(d.values()):
            # iterate over values of every key - fill rows
            for i, el in enumerate(v, start=1):
                # Otherwise cast as str
                if isinstance(el, (type(None), type(ihm.unknown))):
                    el_ = NA

                # Check if int or float
                elif isinstance(el, int) or isinstance(el, float):
                    el_ = el

                # If string, try casting as int or float
                elif isinstance(el, str):
                    try:
                        el_ = int(el)
                    except (TypeError, ValueError):
                        el_ = str(el)

                    if isinstance(el_, str):
                        try:
                            el_ = float(el)
                        except (TypeError, ValueError):
                            el_ = str(el)

                else:
                    el_ = str(el)

                try:
                    output_list[i, j] = el_
                except IndexError:
                    logging.error(
                        'Dict has excessive elements. Ignoring them.')

        output_list = output_list.tolist()

    return output_list


def format_RB_text(tex: list) -> str:
    '''
    convert RB information to text for supp table
    '''
    val = []
    for el in tex:
        for subel in el:
            if subel == el[-1] and el == tex[-1]:
                val.append(str(subel))
            elif subel == el[-1] and el != tex[-1]:
                val.append(str(subel)+', ')
            else:
                val.append(str(subel)+':')

    if len(val) == 0 or len(tex) == 0 or val == '':
        val = ['-']
    return ''.join(val)


def format_flex_text(tex: list) -> str:
    '''
    convert flex information to text for supp table
    '''
    val = []
    for el in tex:
        for subel in el:
            if subel == el[-1] and el == tex[-1]:
                val.append(str(subel))
            else:
                val.append(str(subel)+', ')

    if len(val) == 0 or len(tex) == 0 or val == '':
        val = ['-']

    return ''.join(val)


def format_tuple(tex: list) -> str:
    return str(tex[0])+'-'+str(tex[1])


def dict_to_JSlist_rows(dict1: dict, dict2: dict) -> list:
    '''
    format rigid and flexible segments
    '''
    output_list = []
    output_list.append(['Chain ID', 'Rigid segments', 'Flexible segments'])
    for ind, el in dict1.items():
        output_list.append(
            [ind, format_RB_text(el), format_flex_text(dict2[ind])])
    return output_list


def islistempty(inlist: list) -> bool:
    '''
    minor func
    '''
    if isinstance(inlist, list):
        return all(map(islistempty, inlist))
    return False


def cat_list_string(listn: list) -> str:
    '''
    minor func
    '''
    result = ' '
    for ind in range(len(listn)):
        if ind == 0:
            result += str(listn[ind])
        else:
            result += ','
            result += str(listn[ind])
    return result


def get_key_from_val(dict1: dict, val1: str) -> list:
    '''
    minor func
    '''
    return dict1.keys()[dict1.values().index(val1)]


def get_val_from_key(dict1: dict, key1: str) -> list:
    '''
    minor func
    '''
    return dict1[key1]


def get_name(name) -> str:
    '''
    minor func
    '''
    return str(name)


def get_copy(name):
    '''
    minor func
    '''
    if str(name) == '?':
        copy = 'None listed'
    elif '.' in name:
        copy = (name.split('.')[1]).split('@')[0]
    else:
        copy = 1
    return copy


def get_unique_datasets(name: dict) -> list:
    '''
    get all datatypes that are yet to be validated
    the ones that can't or the ones that have already been validated
    are in the sub_data set
    '''
    all_data = set(name['Dataset type'])
    sub_data = {'Integrative model', 'Other', 'Comparative model',
                'Experimental model', 'De Novo model', 'SAS data', 'Crosslinking-MS data', '3DEM volume'}
    fin_data = list(all_data.difference(sub_data))
    output = list()
    for i in fin_data:
        min_list = [j for j in i.split() if j not in ['data']]
        output.append(' '.join(min_list))
    return output


def get_all_files(path_dir):
    '''
    minor func
    '''
    return glob.glob(path_dir)


def runInParallel(*fns):
    '''
    minor func
    '''
    proc = []
    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()


def runInParallel_noargs(*fns):
    '''
    minor func
    '''
    proc = []
    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)
    for p in proc:
        p.join()


def get_output_file_html(prefix: str) -> str:
    '''
    minor func
    '''
    return f'ValidationReport_{prefix}.html'


def get_supp_file_html(prefix: str) -> str:
    '''
    minor func
    '''
    return f'{prefix}_summary_validation.html'


def get_output_file_temp_html(prefix: str) -> str:
    '''
    minor func
    '''
    return f'{prefix}_full_validation.html'


def get_output_file_pdf(prefix: str) -> str:
    '''
    minor func
    '''
    return f'{prefix}_full_validation.pdf'


def get_output_file_json(prefix: str) -> str:
    '''
    minor func
    '''
    return f'ValidationReport_{prefix}.json'


def get_supp_file_pdf(prefix: str) -> str:
    '''
    minor func
    '''
    return f'{prefix}_summary_validation.pdf'


def get_subunits(sub_dict: dict) -> list:
    '''
    format chains for supplementary/summary table
    '''
    model_number = len(sub_dict['Model ID'])
    sublist = ['%s: Chain %s (%s residues)' % (sub_dict['Subunit name'][i], sub_dict['Chain ID']
                                               [i], str(sub_dict['Total residues'][i])) for i in range(model_number)]
    return list(set(sublist))


def get_datasets(data_dict: dict) -> list:
    '''
    format datasets for supplementary/summary table
    '''

    dataset_number = len(data_dict['ID'])
    datalist = ['%s, %s' % (data_dict['Dataset type'][i], data_dict['Details'][i])
                for i in range(dataset_number)]
    return datalist


def get_software(data_dict: dict) -> list:
    '''
    format software for supplementary/summary table
    '''

    if len(data_dict) > 0:
        dataset_number = len(data_dict['ID'])
        datalist = ['%s (version %s)' % (data_dict['Software name'][i],
                                         data_dict['Software version'][i]) for i in range(dataset_number)]
        return datalist
    else:
        return ['Software details not provided']


def get_RB(data_list: list) -> list:
    '''
    format RB for supplementary/summary table
    '''

    data_num = len(data_list)
    # datalist = ['%s: %s ' % (data_list[i][0], data_list[i][1],)
    #           for i in range(1, data_num)]
    datalist = []
    for i in range(1, data_num):
        if len(data_list[i][1]) < 1:
            data_list[i][1] = 'None'
        datalist.append('%s: %s ' % (data_list[i][0], data_list[i][1]))
    return datalist


def get_flex(data_list: list) -> list:
    '''
    format flexible regions for supplementary/summary table
    '''

    data_num = len(data_list)
    datalist = ['%s: %s ' % (data_list[i][0], data_list[i][2],)
                for i in range(1, data_num)]
    return datalist


def get_method_name(sample_dict: dict) -> str:
    '''
    format method name for supplementary/summary table
    '''

    datastr = '%s ' % (sample_dict['Method name'][0])
    return datastr.replace('monte carlo', 'Monte Carlo')


def get_method_type(sample_dict: dict) -> str:
    '''
    format method type  for supplementary/summary table
    '''

    datastr = '%s ' % (sample_dict['Method type'][0])
    return datastr.replace('monte carlo', 'Monte Carlo')


def get_restraints_info(restraints: dict) -> list:
    '''
    format restraints info for supplementary/summary table
    '''

    restraints_num = len(restraints['Restraint type'])
    datalist = []
    try:
        dataset = [(restraints['Restraint info'][i], restraints['Restraint type'][i])
                   for i in range(restraints_num)]
    except (ValueError, TypeError, IndexError):
        new_restraints = {key: list(set(val))
                          for key, val in restraints.items()}
        restraints_num = min(len(new_restraints['Restraint info']), len(
            new_restraints['Restraint type']))
        dataset = [(new_restraints['Restraint info'][i], new_restraints['Restraint type'][i])
                   for i in range(restraints_num)]
    for i, j in Counter(dataset).items():
        datalist.append(['%s unique %s: %s' % (j, i[1], i[0])])
    return datalist


def format_list_text(sublist: list) -> str:
    '''
    minor func
    '''
    val = ''
    for el in sublist:
        if el == sublist[-1]:
            val += str(el)+'. '
        else:
            val += str(el)+', '
    if val == '':
        val = '-'
    return val


def all_same(items: list):
    '''
    minor func
    '''
    return all(x == items[0] for x in items)


def mp_readable_format(mp: dict) -> list:
    '''
    Format MolProbity results for supplementary/summary table
    '''
    fin_string = []
    for ind, el in enumerate(mp['Models']):
        fin_string.append('Model-'+str(el)+': '+'Clashscore = ' +
                          str(mp['Clashscore'][ind]) + ', ' + 'Number of Ramachandran outliers = ' +
                          str(mp['Ramachandran outliers'][ind]) + ', '+'Number of sidechain outliers = ' +
                          str(mp['Sidechain outliers'][ind]))
    return fin_string


def get_rg_data(rg_dict: dict) -> list:
    '''
    format rg data for supplementary/summary table
    '''

    fin_rg = []
    for key, val in rg_dict.items():
        fin_rg.append(key+': Rg from Gunier is ' +
                      str(val[0])+' nm and Rg from p(r) is ' + str(val[1])+' nm')
    return fin_rg


def get_rg_data_fits(rg_dict: dict) -> list:
    '''
    format sas model fits for supplementary/summary table
    '''

    fin_rg = []
    for key, val in rg_dict.items():
        for ind, el in enumerate(val):
            count = ind+1
            fin_rg.append(key+': Fit ' + str(count) +
                          ' with &#x3A7;&#xb2; value ' + str(el))
    return fin_rg


def get_cx_data_fits(cx_dict: dict) -> list:
    '''
    format crosslinking-MS data for supplementary/summary table
    '''

    fin_cx = []
    count = 0
    for key, val in cx_dict.items():
        count += 1
        fin_cx.append('Crosslinking-MS Fit of medioid: model # ' + str(count) +
                      ', percentage satisfied ' + str(round(val, 2))+'%')
    return fin_cx


def clean_all(report=None):
    '''
    delete all generated files
    '''

    # dirname_ed = os.getcwd()
    # os.listdir('.')
    # for item in os.listdir('.'):
    #     if item.endswith('.txt'):
    #        os.remove(item)
    #    if item.endswith('.csv'):
    #        os.remove(item)
    #    if item.endswith('.json'):
    #        os.remove(item)
    #    if item.endswith('.sascif'):
    #        os.remove(item)

    if report:
        report.clean()


# Adapted from
# https://stackoverflow.com/questions/52859751/
# most-efficient-way-to-find-order-of-magnitude-of-float-in-python
def order_of_magnitude(value: float) -> float:
    '''
    calculate the order of magnitude for a given number

    >>> order_of_magnitude(135)
    2.0

    '''
    if value <= 0:
        raise(f'Wrong value: {value}. '
              'This function works only for positive values')
    return np.floor(np.log10(value))


def calc_optimal_range(counts: list) -> tuple:
    '''
    heuristics to find optimal range for plots

    >>> calc_optimal_range((10, 1567))
    (9.0, 1568.5669999999998)

    '''

    # Find min/max values
    upper = max(counts)
    lower = min(counts)

    # In peculiar cases add an arbitary offset to the range
    if upper == 0:
        upper = 10
        lower = 0

    # Find the data's order of magnitude and make offset.
    # Typically it would be .001%
    if upper > 0:
        oom = order_of_magnitude(upper)
        upper = upper * (1 + 10 ** (-oom))

    # if lower > 0:
    #     oom = order_of_magnitude(lower)
    #     # Do not allow the range to go below zero
    #     lower = max(0, lower * (1 - 10 ** (-oom)))

    lower = 0

    assert lower >= 0 and upper > 0

    return(lower, upper)

def compress_cx_stats(cx_stats: dict) -> list:
    '''Extract per-model satisfactions stats as a flat list'''
    out_stats = []
    for sg, sgv in cx_stats.items():
        for st, stv in sgv.items():
            for mg, mgv in stv.items():
                out_stats.append(mgv['cx_stats']['All']['Satisfied'])

    return out_stats

def get_python_ihm_version() -> str:
    """returns Python-IHM version"""
    import ihm
    return ihm.__version__

def get_hierarchy_from_atoms(atoms) -> dict:
    """Construct polymer hierarchy from a list of atoms"""
    def infinite_defaultdict(): return defaultdict(infinite_defaultdict)
    root = infinite_defaultdict()

    for a in atoms:
        root[a.asym_unit.id][a.seq_id][a.atom_id] = a

    return root

def get_hierarchy_from_model(model) -> dict:
    """Construct polymer hierarchy from atoms and beads in the model"""
    def infinite_defaultdict(): return defaultdict(infinite_defaultdict)
    root = infinite_defaultdict()

    for a in model.get_atoms():
        root[a.asym_unit.id][a.seq_id][a.atom_id] = a

    for r in model.representation:
        if r.granularity == 'by-residue':
            for i in range(r.asym_unit.seq_id_range[0],
                           r.asym_unit.seq_id_range[1] + 1):
                root[r.asym_unit.asym.id][i]['CA'] = None
                root[r.asym_unit.asym.id][i]['coarse-grained'] = None

        elif r.granularity == 'by-feature':
            for i in range(r.asym_unit.seq_id_range[0],
                           r.asym_unit.seq_id_range[1] + 1):
                root[r.asym_unit.asym.id][i]['coarse-grained'] = None

    for s in model.get_spheres():
        # Consider only by-residue spheres
        bs = get_bead_size(s)
        if bs == 1:

            seq_id = s.seq_id_range[0]

            if root[s.asym_unit.id][seq_id]['CA'] is None:
                root[s.asym_unit.id][seq_id]['CA'] = s
            if root[s.asym_unit.id][seq_id]['coarse-grained'] is None:
                root[s.asym_unit.id][seq_id]['coarse-grained'] = s

        else:

            for seq_id in range(s.seq_id_range[0], s.seq_id_range[1] + 1):

                if root[s.asym_unit.id][seq_id]['coarse-grained'] is None:
                    root[s.asym_unit.id][seq_id]['coarse-grained'] = s
                else:
                    s_ = root[s.asym_unit.id][seq_id]['coarse-grained']
                    # Select best possible resolution
                    if get_bead_size(s) < get_bead_size(s_):
                        root[s.asym_unit.id][seq_id]['coarse-grained'] = s


    return root

def get_bead_size(sphere: ihm.model.Sphere) -> int:
    """Number of residues per bead"""
    return sphere.seq_id_range[1] - sphere.seq_id_range[0] + 1

def pretty_print_representations(reprs: dict) -> list:
    """Pretty print information about representation scales"""
    pretty_reprs = []
    for reprs_ in reprs:
        out = ''
        if (reprs_['atomic'] and reprs_['coarse-grained']) or \
        (reprs_['coarse-grained'] and len(reprs_['coarse-grain_levels']) > 1):
            out += 'Multiscale: '

        if reprs_['atomic']:
            out += 'Atomic'

        if reprs_['coarse-grained']:
            if out != '':
                if out[-1] != ' ':
                    out += '; '

            out += 'Coarse-grained: '
            levels = reprs_['coarse-grain_levels']
            if len(levels) == 1:
                out += f'{levels[0]:d}'
            else:
                min_level = min(levels)
                max_level = max(levels)
                out += f'{min_level:d} - {max_level:d}'
            out += ' residue(s) per bead'

        pretty_reprs.append(out)

    return pretty_reprs

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]

def check_for_dataset_type(dataset_list: list=None, dataset_type=None) -> bool:
    """check if the specific dataset type is present in the dataset list"""
    flag = False

    for dataset in dataset_list:
        if isinstance(dataset, dataset_type):
            flag = True

    return flag


def summarize_entities(rep_info: dict) -> list:
    sum_entities = []

    for rep_ in rep_info:
        for k, v in rep_['Chains'].items():
            if isinstance(v['Total residues'], int):
                tr = f'{v["Total residues"]} residues'
            else:
                tr = None
            data_ = (v['Molecule name'], ', '.join(v['Chains']), tr)
            sum_entities.append(data_)

    sum_entities = sorted(set(sum_entities), key=lambda x: x[1])

    output = []

    for e in sum_entities:
        if e[2] is None:
            l = f"{e[0]}: chain(s) {e[1]}"
        else:
            l = f"{e[0]}: chain(s) {e[1]} ({e[2]})"

        output.append(l)

    return output

def summarize_segments(rep_info: dict) -> list:
    output = []

    for rep_ in rep_info:
        rigid, flexible = 0, 0
        for k, v in rep_['Chains'].items():
            rigid_ = len(v['Rigid segments']) * len(v['Chains'])
            flexible_ = len(v['Flexible segments']) * len(v['Chains'])

            rigid += rigid_
            flexible += flexible_

        output.append((f'{rigid}, {flexible}'))

    return output

def parse_ihm_cif(fname, encoding='utf8') -> tuple:
    try:
        with open(fname, encoding=encoding) as fh:
            system, = ihm.reader.read(fh)
    except UnicodeDecodeError:
        encoding = 'ascii'
        with open(fname, encoding=encoding, errors='ignore') as fh:
            system, = ihm.reader.read(fh)

    return(system, encoding)

def is_atomic(data: ihm.System|ihm.model.Model):
    flag = False

    if isinstance(data, ihm.System):
        for group, model in data._all_models():
            flag = is_atomic(model)
            if flag:
                break

    elif isinstance(data, ihm.model.Model):
        if len(data._atoms) > 0:
            flag = True

    return flag

def is_cg(data: ihm.System|ihm.model.Model):
    flag = False

    if isinstance(data, ihm.System):
        for group, model in data._all_models():
            flag = is_cg(model)
            if flag:
                break

    elif isinstance(data, ihm.model.Model):
        if len(data._spheres) > 0:
            flag = True

    return flag

# https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish
class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def format_range(data, format=".2f"):
    if len(data) > 1:
        min_ = np.nanmin(data)
        max_ = np.nanmax(data)
        r_ = f'{min_:{format}}-{max_:{format}}'
    else:
        min_ = np.nanmin(data)
        r_ = f'{min_:{format}}'

    return r_

def get_alphafolddb_link(acc: str) -> str|None:
    """Format link for AlphaFold DB"""
    url = None
    regexp = '^AF-(?P<uniprot>[0-9A-Za-z]+)-F1$'
    m = re.match(regexp, acc)
    if m:
        uid = m.groupdict()['uniprot']
        url_ = f'https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v6.cif'
        r = requests.head(url_)
        if r.status_code == 200:
            url = f"https://alphafold.ebi.ac.uk/entry/{uid}"

    return url

def encode_img(buffer):
    img_base64 = base64.b64encode(buffer)  # encode to base64 (bytes)
    img_base64 = img_base64.decode()    # convert bytes to string
    return img_base64

def load_img(path):
    with open(path, 'rb') as f:
        img_ = f.read() # read bytes from file

    img_base64 = encode_img(img_)
    return img_base64

def format_wwpdb_url(pdb_id: str) -> str:
    """Generate a url to the wwPDB entry page"""

    url = ''

    if is_pdb_id(pdb_id):
        url = f"https://dx.doi.org/10.2210/pdb{pdb_id.lower()}/pdb"
    elif is_pdbx_id(pdb_id):
        # TODO: Update in the future, when wwPDB will start issuing new DOIs
        url = f"https://www.wwpdb.org/pdb?id={pdb_id.lower()}"
    else:
        logging.error(f"Wrong PDB ID: {pdb_id}")

    return url

def is_pdb_id(pdb_id: str) -> bool:
    """Check if the PDB ID is in PDB format"""
    return pdb_id and len(pdb_id) == 4 and pdb_id[0].isdigit()

def is_pdbx_id(pdbx_id: str) -> bool:
    """Check if the PDB ID is in extended format"""
    return pdbx_id and len(pdbx_id) == 12 and pdbx_id[:4].lower() == 'pdb_'

def is_pdb_dev_id(pdb_dev_id: str) -> bool:
    """Check if the PDB-DEV ID"""
    return pdb_dev_id and len(pdb_dev_id) == 15 and pdb_dev_id[:7].upper() == 'PDBDEV_'

def format_pdb_id(pdb_id: str) -> str:
    """Format short PDB ID"""
    if is_pdb_id(pdb_id):
        pdb_id = pdb_id.upper()
    else:
        raise ValueError(f"Wrong PDB ID: {pdb_id}")

    return pdb_id

def format_pdbx_id(pdbx_id: str) -> str:
    """Format all PDB IDs to extended format"""
    if is_pdbx_id(pdbx_id):
        pdbx_id = pdbx_id.lower()
    elif is_pdb_id(pdbx_id):
        pdbx_id = f"pdb_0000{pdbx_id.lower()}"
    else:
        raise ValueError(f"Wrong PDBx ID: {pdbx_id}")

    return pdbx_id

def format_pdb_dev_id(pdb_dev_id: str) -> str:
    """Format PDB-DEV ID"""
    if is_pdb_dev_id(pdb_dev_id):
        pdb_dev_id = pdb_dev_id.upper()
    else:
        raise ValueError(f"Wrong PDB-DEV ID: {pdb_dev_id}")

    return pdb_dev_id

def get_datasets_summary(system: ihm.System) -> list:
    """Get counts for all data types used for modeling"""
    datasets = []
    datatypes = [d.data_type for d in set(system.orphan_datasets)]
    exps, models = [], []
    for k, v in Counter(datatypes).items():
        if re.search('model', k):
            models.append((k, 'Starting model', v))
        else:
            exps.append((k, 'Experimental data', v))
    exps = sorted(exps, key=lambda x: x[0])
    models = sorted(models, key=lambda x: x[1])
    datasets = exps + models
    if len(datasets) > 0:
        datasets.insert(0, ('Name', 'Type', 'Count'))

    return datasets

def add_pdf_watermark(input_pdf: str, watermark_pdf: str) -> str:
    """Add watermark to every page of the input file"""
    reader = pypdf.PdfReader(input_pdf)
    watermark = pypdf.PdfReader(watermark_pdf).pages[0]
    writer = pypdf.PdfWriter()

    for page in reader.pages:
        # Merge the original page ON TOP of the watermark
        page.merge_page(watermark, over=False)
        writer.add_page(page)

    with open(input_pdf, "wb") as f:
        writer.write(f)

    return input_pdf
