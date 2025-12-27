#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# format_checker.py - Detect file format and check residue and atom names in IHMCIF file.
#
# Copyright (C) 2025 Arthur Zalevsky <aozalevsky@gmail.com>
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
Detect file format and check residue and atom names in IHMCIF file.

This module provides functionality to:
- Detect file format (PDB, mmCIF, or IHMCIF)
- Validate CIF files against PDBx dictionary
- Validate IHMCIF files against combined PDBx+IHM dictionary (following the approach
  from https://github.com/ihmwg/python-ihm/blob/main/examples/validate_pdb_ihm.py)
- Check residue and atom names in IHMCIF files
"""

import sys
import re
import logging
import argparse as ag
from pathlib import Path
from enum import Enum
from typing import Tuple

# Import ihm modules only when needed for validation
try:
    import ihm, ihm.reader, ihm.util.make_mmcif, ihm.dictionary
    import io
    import urllib.request
    IHM_AVAILABLE = True
except ImportError:
    IHM_AVAILABLE = False
    io = None
    urllib = None

# Non-standard histidine names (protonation states)
HISTIDINES = frozenset(('HIP', 'HID', 'HIE'))


class FileFormat(Enum):
    """Enumeration of supported file formats"""
    PDB = "PDB"
    MMCIF = "PDBx/mmCIF"
    IHMCIF = "IHMCIF"
    UNKNOWN = "UNKNOWN"


def detect_format(file_path: str, max_lines: int = 1000) -> Tuple[FileFormat, str]:
    """
    Detect the format of a structural biology file.
    
    Args:
        file_path: Path to the file to analyze
        max_lines: Maximum number of lines to read for detection (default: 1000)
    
    Returns:
        Tuple of (FileFormat enum, reason string)
    
    Raises:
        FileNotFoundError: If the file does not exist
        IOError: If the file cannot be read
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    if not file_path.is_file():
        raise IOError(f"Path is not a file: {file_path}")
    
    # Try different encodings
    encodings = ['utf-8', 'latin-1', 'ascii']
    content_lines = None
    
    for encoding in encodings:
        try:
            with open(file_path, 'r', encoding=encoding) as f:
                content_lines = [line for i, line in enumerate(f) if i < max_lines]
            break
        except UnicodeDecodeError:
            continue
    
    if content_lines is None:
        raise IOError(f"Could not read file with any supported encoding: {file_path}")
    
    if not content_lines:
        return FileFormat.UNKNOWN, "File is empty"
    
    # IMPORTANT: Check for CIF format FIRST before checking for PDB format
    # PDBx/mmCIF and IHMCIF files can have _atom_site tables that contain
    # ATOM/HETATM records, which would cause false positives for PDB format
    cif_pattern = re.compile(r'^data_')
    has_data_block = False
    for line in content_lines[:50]:  # Check first 50 lines
        if cif_pattern.match(line.strip()):
            has_data_block = True
            break
    
    # If it's a CIF file, process it as such (don't check for PDB format)
    if has_data_block:
        # Will be processed below as CIF format
        pass
    else:
        # Not a CIF file, check for PDB format (fixed-width with ATOM/HETATM records)
        # PDB files can have ATOM records later in the file (after HEADER, REMARK, etc.)
        pdb_pattern = re.compile(r'^(ATOM|HETATM)\s+\d+')
        
        # First, check the first 200 lines (most PDB files have ATOM records early)
        for i, line in enumerate(content_lines[:200]):
            if pdb_pattern.match(line.strip()):
                return FileFormat.PDB, f"Found ATOM/HETATM record at line {i+1}"
        
        # If not found, check the entire file (some PDB files have long headers)
        # But limit to avoid reading huge files unnecessarily
        if len(content_lines) > 200:
            for i, line in enumerate(content_lines[200:min(1000, len(content_lines))], start=200):
                if pdb_pattern.match(line.strip()):
                    return FileFormat.PDB, f"Found ATOM/HETATM record at line {i+1}"
        
        return FileFormat.UNKNOWN, "File does not appear to be PDB or PDBx/mmCIF format"
    
    # Now we know it's PDBx/mmCIF format, check if it's IHMCIF
    content = '\n'.join(content_lines)
    
    # First, check for IHM dictionary reference in _audit_conform section
    # This is the most reliable indicator (as per python-ihm library example)
    # IHMCIF files should reference both PDBx and IHM dictionaries
    ihm_dict_patterns = [
        r'mmcif_ihm_ext\.dic',  # Extended IHM dictionary
        r'mmcif_ihm\.dic',       # Standard IHM dictionary
    ]
    
    # Check in _audit_conform section (most reliable)
    # This follows the python-ihm library approach for detecting IHMCIF files
    # IHMCIF files reference the IHM dictionary in the audit_conform loop
    # Pattern 1: Check if _audit_conform.dict_name appears before the dictionary name
    audit_conform_field_pattern = re.compile(
        r'_audit_conform\.dict_name.*?(?:^|\n).*?(mmcif_ihm_ext\.dic|mmcif_ihm\.dic)',
        re.IGNORECASE | re.MULTILINE | re.DOTALL
    )
    
    # Pattern 2: Check for dictionary name in loop data (after _audit_conform.dict_name header)
    # The loop structure has dict_name as a column, then data rows with the dictionary name
    audit_conform_loop_pattern = re.compile(
        r'_audit_conform\.dict_name.*?\n(?:[^\n]*\n)*?[^\n]*(mmcif_ihm_ext\.dic|mmcif_ihm\.dic)',
        re.IGNORECASE | re.MULTILINE
    )
    
    # Pattern 3: Simple check for IHM dictionary name anywhere near audit_conform
    audit_conform_simple_pattern = re.compile(
        r'(?:^|\n).*?(?:mmcif_ihm_ext\.dic|mmcif_ihm\.dic).*?(?:^|\n).*?_audit_conform|'
        r'_audit_conform.*?(?:^|\n).*?(?:mmcif_ihm_ext\.dic|mmcif_ihm\.dic)',
        re.IGNORECASE | re.MULTILINE | re.DOTALL
    )
    
    if (audit_conform_field_pattern.search(content) or 
        audit_conform_loop_pattern.search(content) or
        audit_conform_simple_pattern.search(content)):
        return FileFormat.IHMCIF, "Found IHM dictionary reference in _audit_conform section"
    
    # Also check for IHM dictionary reference anywhere in the file
    for pattern in ihm_dict_patterns:
        if re.search(pattern, content, re.IGNORECASE):
            clean_pattern = pattern.replace('\\', '')
            return FileFormat.IHMCIF, f"Found IHM dictionary reference: {clean_pattern}"
    
    # Check for IHMCIF-specific categories (secondary check)
    ihm_category_patterns = [
        r'_ihm_entity_poly_segment',
        r'_ihm_struct_assembly',
        r'_ihm_model_representation',
        r'_ihm_modeling_protocol',
        r'_ihm_model_list',
        r'_ihm_model_group',
        r'_ihm_sphere_obj_site',
        r'_ihm_gaussian_obj_site',
        r'_ihm_cross_link_restraint',
        r'_ihm_dataset_group',
        r'_ihm_ensemble_info',
        r'_ihm_starting_model_details',
    ]
    
    ihm_categories_found = []
    for pattern in ihm_category_patterns:
        if re.search(pattern, content, re.IGNORECASE):
            clean_pattern = pattern.replace('\\', '').replace('^', '').replace('$', '')
            ihm_categories_found.append(clean_pattern)
    
    if ihm_categories_found:
        return FileFormat.IHMCIF, f"Found IHM-specific categories: {', '.join(ihm_categories_found[:3])}"
    
    # Check for mmCIF indicators
    mmcif_indicators = [
        r'_atom_site\.',
        r'_pdbx_',
        r'_entity\.',
        r'_struct\.',
        r'_cell\.',
        r'_symmetry\.',
    ]
    
    mmcif_found = False
    for pattern in mmcif_indicators:
        if re.search(pattern, content, re.IGNORECASE):
            mmcif_found = True
            break
    
    if mmcif_found:
        return FileFormat.MMCIF, "Found PDBx/mmCIF categories but no IHM-specific categories"
    
    # If it has data_ block but no clear indicators, it might be a minimal CIF
    # Check if it has any CIF-like structure
    if re.search(r'^_[\w.]+', content, re.MULTILINE):
        return FileFormat.MMCIF, "PDBx/mmCIF format detected but format type uncertain (assuming PDBx/mmCIF)"
    
    return FileFormat.UNKNOWN, "PDBx/mmCIF format detected but cannot determine type"

def parse_ihm_cif(fname, encoding='utf8') -> tuple:
    """Parse an IHMCIF file using the ihm library"""
    if not IHM_AVAILABLE:
        raise ImportError("ihm library is required for parsing IHMCIF files")
    
    try:
        with open(fname, encoding=encoding) as fh:
            system, = ihm.reader.read(fh)
    except UnicodeDecodeError:
        encoding = 'ascii'
        with open(fname, encoding=encoding, errors='ignore') as fh:
            system, = ihm.reader.read(fh)

    return(system, encoding)


# Dictionary cache to avoid reloading dictionaries multiple times
_DICT_CACHE = {}


def _load_pdbx_dictionary():
    """Load the PDBx/mmCIF dictionary from wwPDB"""
    if 'pdbx' in _DICT_CACHE:
        return _DICT_CACHE['pdbx']
    
    if not IHM_AVAILABLE:
        raise ImportError("ihm library is required for dictionary validation")
    
    try:
        fh = urllib.request.urlopen(
            'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic')
        d_pdbx = ihm.dictionary.read(fh)
        fh.close()
        _DICT_CACHE['pdbx'] = d_pdbx
        return d_pdbx
    except Exception as e:
        raise IOError(f"Failed to load PDBx dictionary: {e}")


def _load_ihm_dictionary():
    """Load the IHM dictionary from wwPDB"""
    if 'ihm' in _DICT_CACHE:
        return _DICT_CACHE['ihm']
    
    if not IHM_AVAILABLE:
        raise ImportError("ihm library is required for dictionary validation")
    
    try:
        fh = urllib.request.urlopen(
            'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ihm.dic')
        d_ihm = ihm.dictionary.read(fh)
        fh.close()
        _DICT_CACHE['ihm'] = d_ihm
        return d_ihm
    except Exception as e:
        raise IOError(f"Failed to load IHM dictionary: {e}")


def _load_flr_dictionary():
    """Load the FLRCIF dictionary from wwPDB (for FRET/fluorescence data)"""
    if 'flr' in _DICT_CACHE:
        return _DICT_CACHE['flr']
    
    if not IHM_AVAILABLE:
        raise ImportError("ihm library is required for dictionary validation")
    
    try:
        fh = urllib.request.urlopen(
            'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ihm_flr_ext.dic')
        d_flr = ihm.dictionary.read(fh)
        fh.close()
        _DICT_CACHE['flr'] = d_flr
        return d_flr
    except Exception as e:
        raise IOError(f"Failed to load FLRCIF dictionary: {e}")


def validate_cif_against_dictionary(file_path: str, dictionary) -> None:
    """
    Validate a CIF file against a dictionary.
    
    Args:
        file_path: Path to the CIF file to validate
        dictionary: Dictionary object to validate against
    
    Raises:
        ihm.dictionary.ValidatorError: If validation fails
        IOError: If file cannot be read
    """
    if not IHM_AVAILABLE:
        raise ImportError("ihm library is required for dictionary validation")
    
    # Read the file with proper encoding handling
    # The encoding for mmCIF files isn't strictly defined, so first try UTF-8
    # and if that fails, strip out any non-ASCII characters
    try:
        with open(file_path, 'rb') as f:
            cif_bytes = f.read()
        try:
            cif_text = cif_bytes.decode('utf-8')
        except UnicodeDecodeError:
            cif_text = cif_bytes.decode('ascii', errors='ignore')
    except Exception as e:
        raise IOError(f"Failed to read file {file_path}: {e}")
    
    # Validate against the dictionary
    fh = io.StringIO(cif_text)
    dictionary.validate(fh)


def validate_mmcif(file_path: str) -> None:
    """
    Validate a PDBx/mmCIF file against the PDBx dictionary.
    
    Args:
        file_path: Path to the PDBx/mmCIF file to validate
    
    Raises:
        ihm.dictionary.ValidatorError: If validation fails
    """
    d_pdbx = _load_pdbx_dictionary()
    validate_cif_against_dictionary(file_path, d_pdbx)


def validate_ihmcif(file_path: str) -> None:
    """
    Validate an IHMCIF file against the combined PDBx/mmCIF+IHMCIF dictionary.
    
    Deposited integrative models should conform to both the PDBx dictionary
    (used to define basic structural information such as residues and chains)
    and the IHM dictionary (used for information specific to integrative modeling).
    Some entries also use the FLRCIF dictionary for FRET/fluorescence data.
    
    Args:
        file_path: Path to the IHMCIF file to validate
    
    Raises:
        ihm.dictionary.ValidatorError: If validation fails
    """
    d_pdbx = _load_pdbx_dictionary()
    d_ihm = _load_ihm_dictionary()
    
    # Combine PDBx and IHM dictionaries using the + operator
    pdbx_ihm = d_pdbx + d_ihm
    
    # Check if the file references FLRCIF dictionary or contains FLRCIF categories
    # If so, also include it in validation
    try:
        with open(file_path, 'rb') as f:
            cif_bytes = f.read()
        try:
            cif_text = cif_bytes.decode('utf-8')
        except UnicodeDecodeError:
            cif_text = cif_bytes.decode('ascii', errors='ignore')
    except Exception as e:
        # If we can't read the file for FLRCIF detection, just do standard validation
        logging.warning(f"Error reading file for FLRCIF detection, using standard validation: {e}")
        validate_cif_against_dictionary(file_path, pdbx_ihm)
        return
    
    # Check for FLRCIF dictionary reference in _audit_conform section
    # First check _audit_conform section systematically (similar to IHM detection)
    audit_conform_field_pattern = re.compile(
        r'_audit_conform\.dict_name\s+(?:mmcif_ihm_flr_ext\.dic|mmcif_ihm_flr)',
        re.IGNORECASE | re.MULTILINE
    )
    
    audit_conform_loop_pattern = re.compile(
        r'_audit_conform\.dict_name.*?\n(?:[^\n]*\n)*?[^\n]*(mmcif_ihm_flr_ext\.dic|mmcif_ihm_flr)',
        re.IGNORECASE | re.MULTILINE
    )
    
    # Also check for FLRCIF dictionary reference anywhere or any _flr_ category
    flr_dict_pattern = re.compile(
        r'mmcif_ihm_flr_ext\.dic|mmcif_ihm_flr|_flr_',
        re.IGNORECASE
    )
    
    if (audit_conform_field_pattern.search(cif_text) or 
        audit_conform_loop_pattern.search(cif_text) or
        flr_dict_pattern.search(cif_text)):
        # File uses FLRCIF dictionary, include it in validation
        logging.info("FLRCIF dictionary detected, including in validation")
        try:
            d_flr = _load_flr_dictionary()
            pdbx_ihm_flr = pdbx_ihm + d_flr
            # Use validate_cif_against_dictionary which handles file reading properly
            validate_cif_against_dictionary(file_path, pdbx_ihm_flr)
            logging.info("IHMCIF file validated against PDBx/mmCIF+IHMCIF+FLRCIF dictionary")
        except (IOError, ImportError) as e:
            # If we can't load FLRCIF dictionary, fall back to standard validation
            logging.warning(f"Error loading FLRCIF dictionary, falling back to standard validation: {e}")
            validate_cif_against_dictionary(file_path, pdbx_ihm)
        # Don't catch validation errors - let them propagate
    else:
        # Standard IHMCIF validation without FLRCIF
        validate_cif_against_dictionary(file_path, pdbx_ihm)
        logging.info("IHMCIF file validated against PDBx/mmCIF+IHMCIF dictionary")


def check_entities_histidines(system: ihm.System, histidines=HISTIDINES):
    """Find any non-standard histidine chemical components"""
    out = []
    his = ihm.LPeptideAlphabet()['H']
    for e in system.entities:
        for c in e.sequence:
            if c.id in histidines:
                out.append(c.id)

    if len(out) > 0:
        raise(ValueError(f"Non-canonical histidine variant found: {', '.join(set(out))}"))

def check_models(system: ihm.System):
    """Find any non-standard histidine chemical components"""
    for state_group in system.state_groups:
        for state in state_group:
            for model_group in state:
                for model in model_group:
                    ihm.util.make_mmcif._check_atom_names(model, check_all=True)

def check_all_exception(system: ihm.System):
    """Perform all checks. Throw an exception if a check fails."""
    # Disable atom check until python-ihm fixes
    # checks = [check_entities_histidines, check_models]
    checks = [check_entities_histidines]

    for check in checks:
        check(system)

def check_all_log(system: ihm.System) -> int:
    """Perform all checks. Throw a message in the log if a check fails and return a non-zero exit code"""
    # Disable atom check until python-ihm fixes
    # checks = [check_entities_histidines, check_models]
    checks = [check_entities_histidines]
    exit_code = 0
    for check in checks:
        try:
            check(system)
        except ValueError as e:
            logging.error(e)
            exit_code = 127

    return exit_code

def check_file_format(fname: str, validate_dictionary: bool = True, raise_on_error: bool = True):
    """
    Check file format and validate that it is IHMCIF.
    
    Args:
        fname: Path to the file to check
        validate_dictionary: If True, validate IHMCIF files against dictionaries
        raise_on_error: If True, raise ValueError on format errors; if False, return error message
    
    Returns:
        If raise_on_error is False: (success: bool, error_msg: str or None)
        If raise_on_error is True: None (raises ValueError on error)
    
    Raises:
        ValueError: If file format is not IHMCIF (when raise_on_error is True)
    """
    format_type, reason = detect_format(fname)
    
    # This script only works with IHMCIF files
    if format_type == FileFormat.PDB:
        error_msg = f"File format is PDB, expected IHMCIF. {reason}"
        if raise_on_error:
            raise ValueError(error_msg)
        return False, error_msg
    elif format_type == FileFormat.MMCIF:
        error_msg = f"File format is PDBx/mmCIF, expected IHMCIF. {reason}"
        if raise_on_error:
            raise ValueError(error_msg)
        return False, error_msg
    elif format_type == FileFormat.UNKNOWN:
        error_msg = f"File format could not be determined. {reason}"
        if raise_on_error:
            raise ValueError(error_msg)
        return False, error_msg
    
    # Validate IHMCIF file against dictionaries if requested
    if validate_dictionary:
        if not IHM_AVAILABLE:
            logging.warning("Dictionary validation requested but ihm library is not available. Skipping dictionary validation.")
        elif format_type == FileFormat.IHMCIF:
            try:
                validate_ihmcif(fname)
                # Logging is done inside validate_ihmcif
            except Exception as e:
                error_msg = f"Dictionary validation failed: {e}"
                if raise_on_error:
                    raise ValueError(error_msg)
                return False, error_msg
    else:
        logging.info("Dictionary validation is disabled (use --no-dictionary-validation to explicitly disable)")
    
    if raise_on_error:
        return None
    return True, None

def check_file_exception(fname: str, check_format: bool = True, validate_dictionary: bool = True):
    """
    Parse a file, do all checks, throw an exception if a check fails.
    
    Args:
        fname: Path to the file to check
        check_format: If True, verify the file format before checking
        validate_dictionary: If True, validate CIF/IHMCIF files against dictionaries
    """
    if check_format:
        check_file_format(fname, validate_dictionary=validate_dictionary, raise_on_error=True)
    
    system, encoding = parse_ihm_cif(fname)
    check_all_exception(system)

def check_file_log(fname: str, check_format: bool = True, validate_dictionary: bool = True) -> int:
    """
    Parse a file, do all checks, throw a log message if a check fails and return a non-zero exit code
    
    Args:
        fname: Path to the file to check
        check_format: If True, verify the file format before checking
        validate_dictionary: If True, validate PDBx/mmCIF or IHMCIF files against dictionaries
    
    Returns:
        Exit code: 0 for success, non-zero for failure
    """
    if check_format:
        try:
            success, error_msg = check_file_format(fname, validate_dictionary=validate_dictionary, raise_on_error=False)
            if not success:
                logging.error(error_msg)
                return 1
        except Exception as e:
            logging.error(f"Format detection failed: {e}")
            return 1
    
    system, encoding = parse_ihm_cif(fname)
    exit_code = check_all_log(system)
    return exit_code

if __name__ == "__main__":
    parser = ag.ArgumentParser(
        description="Check residue and atom names in IHMCIF file and detect file format",
        formatter_class=ag.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i structure.cif          # Check IHMCIF file
  %(prog)s --detect-only structure.cif  # Only detect format
  %(prog)s --no-format-check -i file.cif  # Skip format check
        """
    )
    parser.add_argument("-i", "--input_file", help="Path to the input file")
    parser.add_argument(
        "file",
        nargs="?",
        help="Path to the input file (alternative to -i/--input_file)"
    )
    parser.add_argument(
        "--detect-only",
        action="store_true",
        help="Only detect file format, do not perform validation checks"
    )
    parser.add_argument(
        "--no-format-check",
        action="store_true",
        help="Skip format detection and validation (assume file is IHMCIF)"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only output the format name (for scripting, use with --detect-only)"
    )
    parser.add_argument(
        "--no-dictionary-validation",
        action="store_true",
        help="Skip dictionary validation (only check residue/atom names)"
    )

    args = parser.parse_args()
    
    # Get input file from either -i/--input_file or positional argument
    input_file = args.input_file or args.file
    
    if not input_file:
        parser.error("Input file is required (use -i/--input_file or provide as positional argument)")
    
    # If only detecting format
    if args.detect_only:
        try:
            format_type, reason = detect_format(input_file)
            if args.quiet:
                print(format_type.value)
            else:
                print(f"Format: {format_type.value}")
                print(f"Reason: {reason}")
            
            # Exit codes
            if format_type == FileFormat.UNKNOWN:
                sys.exit(2)
            else:
                sys.exit(0)
        except FileNotFoundError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        except IOError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Unexpected error: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Otherwise, perform validation checks
    check_format = not args.no_format_check
    validate_dictionary = not args.no_dictionary_validation
    try:
        check_file_exception(input_file, check_format=check_format, validate_dictionary=validate_dictionary)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)
