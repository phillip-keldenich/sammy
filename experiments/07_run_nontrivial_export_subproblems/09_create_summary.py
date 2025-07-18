if __name__ != "__main__":
    raise ImportError("Not meant to be imported!")

from pathlib import Path
import lzma
import json
import re
import tqdm

output_dir = Path(__file__).absolute().parent / "02_output"
summary_file = output_dir.parent / "10_summary_data.json.xz"
total_files = []
filename_re = re.compile(r'^solved-([^-]+)-(.+)-([0-9]+)\.json\.xz$')

solver_name_mapping = {
    'fixed_sat_kissat': 'non-incremental (kissat)',
    'fixed_sat_cadical': 'non-incremental (CaDiCaL)',
    'fixed_sat_cryptominisat': 'non-incremental (CryptoMS)',
    'fixed_sat_lingeling': 'non-incremental (Lingeling)',
    'incremental_sat_cadical': 'simple incremental (CaDiCaL)',
    'incremental_sat_cryptominisat': 'simple incremental (CryptoMS)',
    'incremental_sat_lingeling': 'simple incremental (Lingeling)',
    'satdsatur_cadical': 'alternating LB-UB (CaDiCaL)',
    'satdsatur_cryptominisat': 'alternating LB-UB (CryptoMS)',
    'satdsatur_lingeling': 'alternating LB-UB (Lingeling)',
    'newsatdsatur_cadical': 'greedy incremental (CaDiCaL)',
    'newsatdsatur_cryptominisat': 'greedy incremental (CryptoMS)',
    'newsatdsatur_lingeling': 'greedy incremental (Lingeling)'
}


def summary_from_data(file, data):
    instance_name = file.parent.name
    match = filename_re.match(file.name)
    if not match:
        raise 
    solver_name_raw = match[1]
    solver_name_mapped = solver_name_mapping[solver_name_raw]
    input_file = file.parent / f"amended-{match[2]}-{match[3]}.json.xz"
    with lzma.open(input_file, "rt") as inf:
        inf_data = json.load(inf)
    amend_events = inf_data['amendment']['events']
    if amend_events[1]['type'] != 'DONE_IMPLIED_VERTEX_ELIMINATION':
        raise RuntimeError("Unexpected events in amended file!")
    amend_start_time = amend_events[1]['time']
    return {
         'instance': instance_name,
         'solver': solver_name_mapped,
         'filename': str(file),
         'outcome': data['outcome'],
         'build_time': data['build_time'],
         'solve_time': data['solve_time'],
         'mes_time': inf_data['amendment']['events'][-1]['time'] - amend_start_time,
         'mes_size': inf_data['amendment']['mes_size'],
         'removal_size': inf_data['num_removed_configs'],
         'replace_size': len(data['assignments']) if data['assignments'] is not None\
                         else inf_data['num_removed_configs']
    }


summary = []
already_scanned = set()
if summary_file.is_file():
    with lzma.open(summary_file, "rt") as sfr:
        summary = json.load(sfr)
    for summary_entry in summary:
        already_scanned.add(summary_entry['filename'])
    print("Loaded existing summary file with", len(already_scanned), "entries")

try:
    for file in tqdm.tqdm(list(output_dir.glob("**/solved-*.json.xz"))):
        if str(file) in already_scanned:
            continue
        with lzma.open(file, "rt") as data_file:
            json_data = json.load(data_file)
        summary.append(summary_from_data(file, json_data))
finally:
    with lzma.open(summary_file, "wt") as sfw:
        json.dump(summary, sfw)

