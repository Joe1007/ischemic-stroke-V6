#!/usr/bin/env python3
"""Shared image/radiomics helpers used by V6 extraction (checkpoint IO, normalisation, ROI radiomics).

Merged 2026-09-30 from the earlier iterations so V6 runs on its own. Only the definitions V6
actually calls are kept; their logic is unchanged. Origin of each block is marked below.
"""
from __future__ import annotations
import gzip
import hashlib
import json
import logging
import math
import re
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import radiomics  # re-exported: extract_v6 records radiomics.__version__



# ---- from ischemic stroke analyseV3/extract_radiomics_v3.py -----------

REDUNDANT_GLCM_MESSAGE = (
    "GLCM is symmetrical, therefore Sum Average = 2 * Joint Average"
)

class KnownRedundantMessageFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        return REDUNDANT_GLCM_MESSAGE not in record.getMessage()

def configure_radiomics_logging():
    """Suppress one expected per-image redundancy notice, retaining other warnings."""
    logger = logging.getLogger("radiomics")
    known_filter = KnownRedundantMessageFilter()
    logger.addFilter(known_filter)
    for handler in logger.handlers:
        handler.addFilter(known_filter)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def canonical_hash(value: dict) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()

def same_geometry(left: sitk.Image, right: sitk.Image, tolerance: float = 1e-5) -> bool:
    return (
        left.GetSize() == right.GetSize()
        and np.allclose(left.GetSpacing(), right.GetSpacing(), atol=tolerance, rtol=0)
        and np.allclose(left.GetOrigin(), right.GetOrigin(), atol=tolerance, rtol=0)
        and np.allclose(left.GetDirection(), right.GetDirection(), atol=tolerance, rtol=0)
    )

def labels_from_file(path: Path) -> dict[int, str]:
    labels = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*(\d+)\s+(.+?)\s*$", line)
        if match and int(match.group(1)) != 0:
            labels[int(match.group(1))] = match.group(2)
    return labels

def clean_features(result: dict) -> dict[str, float]:
    cleaned = {}
    for key, value in result.items():
        key = str(key)
        if key.startswith("diagnostics_"):
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        cleaned[key] = numeric if math.isfinite(numeric) else None
    return cleaned


# ---- from ischemic stroke analyseV5/extract_v5.py ---------------------

def read(path):
    with gzip.open(path,'rt') as f:return json.load(f)

def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    with gzip.open(tmp,'wt') as f:json.dump(value,f,allow_nan=False)
    tmp.replace(path)

def normalize(values,brain,lesion,protocol):
    valid=np.isfinite(values)&(values>0)
    covered=brain&valid
    if covered.sum()<protocol['minimum_reference_voxels']:raise ValueError('Insufficient reference voxels')
    dil=sitk.GetArrayFromImage(sitk.BinaryDilate(sitk.GetImageFromArray(lesion.astype('uint8')),
                [protocol['normal_reference_exclusion_dilation_voxels']]*3))>0
    normal=covered&~dil
    if normal.sum()<protocol['minimum_reference_voxels']:raise ValueError('Insufficient nonlesional reference')
    p1,p99=np.percentile(values[covered],protocol['percentiles'])
    n1,n99=np.percentile(values[normal],protocol['percentiles'])
    if min(p99-p1,n99-n1)<=0:raise ValueError('Degenerate intensity reference')
    images={};qc={}
    for variant in protocol['variants']:
        lo,hi=(n1,n99) if variant=='normal_tissue_unclipped' else (p1,p99)
        arr=np.zeros_like(values,dtype='float32')
        scaled=np.maximum((values[valid]-lo)/(hi-lo),0)*protocol['scale']
        if variant=='clip_p99':scaled=np.minimum(scaled,protocol['scale'])
        arr[valid]=scaled
        images[variant]=arr
        lesion_values=arr[lesion&valid]
        qc[variant]={'reference_low':float(lo),'reference_high':float(hi),
            'reference_voxels':int(normal.sum() if variant=='normal_tissue_unclipped' else covered.sum()),
            'lesion_fraction_ge100':float(np.mean(lesion_values>=100)),
            'lesion_fraction_exactly100':float(np.mean(lesion_values==100)),
            'lesion_unique_intensities':int(len(np.unique(lesion_values))),
            'lesion_max':float(lesion_values.max()),'brain_max':float(arr[covered].max())}
    return images,qc

def crop_roi(mask,reference):
    pts=np.where(mask)
    if not len(pts[0]):raise ValueError('Empty ROI')
    # Keep a background border: PyRadiomics requires both zero and ROI labels.
    start=np.maximum(0,np.array([int(p.min()) for p in pts])-1)
    stop=np.minimum(mask.shape,np.array([int(p.max())+2 for p in pts]))
    slices=tuple(slice(a,b) for a,b in zip(start,stop))
    size=(stop-start)[::-1].tolist();index=start[::-1].tolist()
    cropped_ref=sitk.RegionOfInterest(reference,size,index)
    roi=sitk.GetImageFromArray(mask[slices].astype('uint8'));roi.CopyInformation(cropped_ref)
    return slices,cropped_ref,roi

def radiomics_roi(arr,mask,reference,extractor,protocol):
    # Original-only features are invariant to cropping that retains the whole ROI.
    inside=arr[mask]
    if (float(inside.max())-float(inside.min()))/protocol['radiomics']['binWidth']>protocol['maximum_gray_level_range_bins']:
        raise ValueError('Gray-level range exceeds predeclared resource/QC limit')
    slices,ref,roi=crop_roi(mask,reference)
    image=sitk.GetImageFromArray(arr[slices]);image.CopyInformation(ref)
    result=clean_features(extractor.execute(image,roi,label=1))
    if len(result)!=93 or any(v is None for v in result.values()):
        raise ValueError(f'Expected 93 finite features, got {len(result)}')
    return result
