import csv
from pathlib import Path
import uuid

# Path to the unharmonized cases
CASE_LIST_PTH = '../201904_locate_adhoc_data/tracked_results/unharmonized_cases.list'
# Path to the sequencing file mapping from Matt's catalog
MGI_MAP_PTH = '../matt_catalog/MGI.BamMap.dat'

# Read the list of all the unharmonized cases
unharmonized_cases = set(open(CASE_LIST_PTH).read().splitlines())

# Build the sample (case, sample_type) BAM map on MGI
# Data structure of the map
# (case, sample_type) -> bam_pth
wxs_bam_map = {}
wgs_bam_map = {}
with open(MGI_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        case = row['case']

        # Keep only the unharmonized cases of interest
        if not (case in unharmonized_cases \
                and row['reference'] == 'hg19' \
                and row['data_format'] == 'BAM'):
            continue

        # Get the corresponding bam map
        if row['experimental_strategy'] == 'WGS':
            bam_map = wgs_bam_map
        elif row['experimental_strategy'] == 'WXS':
            bam_map = wxs_bam_map

        sample = (case, row['sample_type'])
        bam_pth = Path(row['data_path'])
        bam_map[sample] = bam_pth

# At the moment we only test for the first 20 cases
all_wxs_samples = [f'{case}_{st}' for (case, st) in sorted(wxs_bam_map.keys())]


current_abs_pth = Path().resolve()

columns = [
    '# sample_name', 'case', 'disease', 'experimental_strategy', 'sample_type',
    'data_path', 'filesize', 'data_format', 'reference', 'UUID', 'system'
]
with open('processed_data/MGI.GBM_custom.BamMap.dat', 'w') as f:
    writer = csv.writer(f, dialect='excel-tab')
    writer.writerow(columns)  # Write header
    for sample in all_wxs_samples:
        case, sample_type = sample.split('_', 1)
        if sample_type == 'blood_normal':
            abbrv = 'N'
        elif sample_type == 'tumor':
            abbrv = 'T'
        else:
            raise ValueError(f'Unknown sample type {sample_type}')
        sample_name = f'{case}.WXS.{abbrv}.hg38-bobo'
        data_pth = Path(current_abs_pth,
                        'processed_data/wxs_mark_dup_bam',
                        f'{sample}.bam')
        data_id = str(uuid.uuid4())
        writer.writerow([
            sample_name, case, 'GBM', 'WXS', sample_type,
            str(data_pth), '0', 'BAM', 'hg38-bobo', data_id, 'MGI',
        ])
