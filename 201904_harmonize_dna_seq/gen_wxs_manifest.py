import csv
from pathlib import Path
import uuid

# Path to the unharmonized cases
CASE_LIST_PTH = '../201907_locate_data/tracked_results/wxs_or_wgs_unharmonized_cases.list'

CUSTOM_BAM_ROOT = 'processed_data/wxs_mark_dup_bam'

# Path to the output manifest file
PREV_OUT_MANIFEST_PTH = '../201904_locate_adhoc_data/tracked_results/MGI.GBM_custom_wxs.BamMap.dat'
OUT_MANIFEST_PTH = '../201907_locate_data/tracked_results/MGI.GBM_custom_wxs.BamMap.dat'

# Read the list of all the unharmonized cases
unharmonized_cases = set(open(CASE_LIST_PTH).read().splitlines())

# Build the sample (case, sample_type) BAM map on MGI
# Data structure of the map
# (case, sample_type) -> bam_pth
all_wxs_custom_bams = sorted(Path(CUSTOM_BAM_ROOT).resolve().glob('*.bam'))
wxs_bam_map = {}
for pth in all_wxs_custom_bams:
    case, sample_type = pth.stem.split('_', 1)
    wxs_bam_map[(case, sample_type)] = pth

# Read the previous manifest to maintain the UUID
previous_manifest_uuids = {}
with open(PREV_OUT_MANIFEST_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        previous_manifest_uuids[row['# sample_name']] = row['UUID']


# Generate the manifest
columns = [
    '# sample_name', 'case', 'disease', 'experimental_strategy', 'sample_type',
    'data_path', 'filesize', 'data_format', 'reference', 'UUID', 'system'
]
with open(OUT_MANIFEST_PTH, 'w') as f:
    writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
    writer.writerow(columns)  # Write header
    for (case, sample_type), data_pth in wxs_bam_map.items():
        if sample_type == 'blood_normal':
            abbrv = 'N'
        elif sample_type == 'tumor':
            abbrv = 'T'
        else:
            raise ValueError(f'Unknown sample type {sample_type}')
        sample_name = f'{case}.WXS.{abbrv}.hg38-bobo'
        file_size = data_pth.stat().st_size

        # Reuse the UUID
        try:
            data_id = previous_manifest_uuids[sample_name]
        except KeyError:
            # Generate a new UUID
            data_id = str(uuid.uuid4())

        writer.writerow([
            sample_name, case, 'GBM', 'WXS', sample_type,
            str(data_pth), str(file_size), 'BAM', 'hg38', data_id, 'MGI',
        ])
