#!/usr/bin/env python3

"""
Script to download both boltz1 and boltz2 caches in parallel
Adapted from https://github.com/jwohlwend/boltz/blob/main/src/boltz/main.py
Adapted by Berent Aldikacti
"""

import argparse
from pathlib import Path
import urllib.request
import tarfile

CCD_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl"
MOL_URL = "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar"

BOLTZ1_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz1_conf.ckpt",
    "https://huggingface.co/boltz-community/boltz-1/resolve/main/boltz1_conf.ckpt",
]

BOLTZ2_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz2_conf.ckpt",
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt",
]

BOLTZ2_AFFINITY_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz2_aff.ckpt",
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt",
]

def download_boltz1(cache: Path) -> None:
    """Download all the required data for Boltz1.

    Parameters
    ----------
    cache : Path
        The cache directory.

    """
    print("Starting Boltz1 downloads...")
    
    # Download CCD
    ccd = cache / "ccd.pkl"
    if not ccd.exists():
        print(f"Downloading the CCD dictionary to {ccd}")
        urllib.request.urlretrieve(CCD_URL, str(ccd))  # noqa: S310
        print("CCD dictionary download completed")
    else:
        print("CCD dictionary already exists, skipping")

    # Download model
    model = cache / "boltz1_conf.ckpt"
    if not model.exists():
        print(f"Downloading the Boltz1 model weights to {model}")
        for i, url in enumerate(BOLTZ1_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(model))  # noqa: S310
                print("Boltz1 model weights download completed")
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ1_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download Boltz1 model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                print(f"Failed to download from {url}, trying next URL...")
                continue
    else:
        print("Boltz1 model weights already exist, skipping")
    
    print("Boltz1 downloads completed!")


def download_boltz2(cache: Path) -> None:
    """Download all the required data for Boltz2.

    Parameters
    ----------
    cache : Path
        The cache directory.

    """
    print("Starting Boltz2 downloads...")
    
    # Download CCD
    mols = cache / "mols"
    tar_mols = cache / "mols.tar"
    if not tar_mols.exists():
        print(f"Downloading the molecular data to {tar_mols} (this may take a while)")
        urllib.request.urlretrieve(MOL_URL, str(tar_mols))  # noqa: S310
        print("Molecular data download completed")
    else:
        print("Molecular data tar already exists, skipping download")
        
    if not mols.exists():
        print(f"Extracting the molecular data to {mols} (this may take a while)")
        with tarfile.open(str(tar_mols), "r") as tar:
            tar.extractall(cache)  # noqa: S202
        print("Molecular data extraction completed")
    else:
        print("Molecular data already extracted, skipping")

    # Download model
    model = cache / "boltz2_conf.ckpt"
    if not model.exists():
        print(f"Downloading the Boltz2 model weights to {model}")
        for i, url in enumerate(BOLTZ2_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(model))  # noqa: S310
                print("Boltz2 model weights download completed")
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download Boltz2 model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                print(f"Failed to download from {url}, trying next URL...")
                continue
    else:
        print("Boltz2 model weights already exist, skipping")

    # Download affinity model
    affinity_model = cache / "boltz2_aff.ckpt"
    if not affinity_model.exists():
        print(f"Downloading the Boltz2 affinity weights to {affinity_model}")
        for i, url in enumerate(BOLTZ2_AFFINITY_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(affinity_model))  # noqa: S310
                print("Boltz2 affinity weights download completed")
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_AFFINITY_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download Boltz2 affinity model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                print(f"Failed to download from {url}, trying next URL...")
                continue
    else:
        print("Boltz2 affinity weights already exist, skipping")
    
    print("Boltz2 downloads completed!")


def main():
    parser = argparse.ArgumentParser(
        description="Download Boltz1 and Boltz2 model caches in parallel"
    )
    parser.add_argument(
        "cache",
        type=str,
        help="Path to the cache directory where models will be downloaded"
    )
    parser.add_argument(
        "mode",
        type=str,
        help="Boltz mode to use. Options: boltz1 or boltz2"
    )
    
    args = parser.parse_args()
    
    # Convert to Path object and create directory if it doesn't exist
    cache_path = Path(args.cache).expanduser()
    cache_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Using cache directory: {cache_path}")
    print("Starting parallel downloads...")
    
    if args.mode == 'boltz1':
        download_boltz1(cache_path)
    elif args.mode == 'boltz2':
        download_boltz2(cache_path)
    else:
        raise ValueError("Incorrect boltz mode argument. Options: boltz1 or boltz2")


if __name__ == "__main__":
    main()