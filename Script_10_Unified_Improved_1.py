#!/usr/bin/env python3
"""
🌱 PLANT HYPOTHETICAL PROTEIN DRUG DISCOVERY PIPELINE v4.0
============================================================
BLAST + Foldseek + Secreted Localization + Druggable Pockets → PRIORITIZED TARGETS
Production-ready for Linux/VSCode bioinformatics workflows
"""

import os
import sys
import json
import time
import random
import requests
import subprocess
import logging
import argparse
import urllib.parse
import urllib.request
import re
from pathlib import Path
import re
import numpy as np
import math
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any
import pickle
import hashlib
import shutil
from dataclasses import dataclass, asdict
import threading
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# CONFIGURATION MANAGEMENT
# ============================================================================

@dataclass
class PipelineConfig:
    # BLAST options
    enable_blast: bool = True          # master on/off for BLAST
    blast_mode: str = "remote"         # "local" or "remote" (NCBI web API)

    # Core paths
    alphafold_output_dir: str = "Script_9_AlphaFold_Structure_Downloader_v2_output"
    fasta_dir: str = "hypothetical_proteins_fasta"
    output_dir: str = "./pipeline_output"
    
    # Analysis parameters
    plant_dbs: List[str] = None
    request_delay: float = 10.0
    poll_interval: float = 30.0
    timeout: int = 120
    max_retries: int = 5
    max_poll_attempts: int = 120
    blast_evalue: str = "1e-5"
    blast_max_targets: int = 5
    signalp_timeout: int = 45
    fpocket_timeout: int = 120
    min_tm_score: float = 0.15
    min_plddt: float = 70.0
    min_druggability_score: float = 0.5
    
    # Performance settings
    max_workers: int = 4
    enable_caching: bool = True
    enable_parallel: bool = True
    checkpoint_interval: int = 10
    
    # Visualization settings
    dpi: int = 300
    figsize_large: Tuple[int, int] = (12, 8)
    figsize_medium: Tuple[int, int] = (10, 6)
    figsize_small: Tuple[int, int] = (8, 8)
    
    def __post_init__(self):
        if self.plant_dbs is None:
            self.plant_dbs = ["afdb50", "afdb-swissprot", "pdb100", "cath50"]

class PipelineLogger:
    """Enhanced logging system for the pipeline"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.log_dir = Path(config.output_dir) / "logs"
        self.log_dir.mkdir(exist_ok=True, parents=True)
        
        # Setup logging
        log_file = self.log_dir / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def info(self, message: str):
        self.logger.info(message)
        
    def warning(self, message: str):
        self.logger.warning(message)
        
    def error(self, message: str):
        self.logger.error(message)
        
    def debug(self, message: str):
        self.logger.debug(message)

class CacheManager:
    """Caching system for expensive operations"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.cache_dir = Path(config.output_dir) / "cache"
        self.cache_dir.mkdir(exist_ok=True, parents=True)
        self.cache_file = self.cache_dir / "analysis_cache.pkl"
        self.cache = self._load_cache()
        
    def _load_cache(self) -> Dict[str, Any]:
        """Load cache from file"""
        try:
            if self.cache_file.exists():
                with open(self.cache_file, 'rb') as f:
                    return pickle.load(f)
        except Exception as e:
            logging.warning(f"Failed to load cache: {e}")
        return {}
    
    def _save_cache(self):
        """Save cache to file"""
        try:
            with open(self.cache_file, 'wb') as f:
                pickle.dump(self.cache, f)
        except Exception as e:
            logging.warning(f"Failed to save cache: {e}")
    
    def get(self, key: str) -> Optional[Any]:
        """Get cached result"""
        if not self.config.enable_caching:
            return None
        return self.cache.get(key)
    
    def set(self, key: str, value: Any):
        """Cache a result"""
        if not self.config.enable_caching:
            return
        self.cache[key] = value
        self._save_cache()
    
    def get_key(self, operation: str, *args) -> str:
        """Generate cache key"""
        content = f"{operation}_{'_'.join(str(arg) for arg in args)}"
        return hashlib.md5(content.encode()).hexdigest()

# ============================================================================
# CORE ANALYSIS FUNCTIONS (ALL 4 CRITERIA) - ENHANCED
# ============================================================================

class ProteinAnalyzer:
    """Enhanced protein analysis class with improved error handling and caching"""

    _signalp_unavailable_logged = False  # class-level flag

    def __init__(self, config: PipelineConfig, logger: PipelineLogger, cache: CacheManager):
        self.config = config
        self.logger = logger
        self.cache = cache
        self._lock = threading.Lock()
    
    # ------------------------------------------------------------------
    # Remote BLASTp via NCBI (no local blastp needed)
    # ------------------------------------------------------------------
    def _ncbi_blast_put(self, fasta_seq: str) -> Optional[str]:
        """Submit BLASTp job to NCBI, return RID or None."""
        params = {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "nr",   # change to "swissprot" etc. if you prefer
            "QUERY": fasta_seq,
        }
        data = urllib.parse.urlencode(params).encode("utf-8")
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

        try:
            with urllib.request.urlopen(url, data=data, timeout=60) as resp:
                text = resp.read().decode("utf-8")
        except Exception as e:
            self.logger.warning(f"NCBI BLAST Put failed: {e}")
            return None

        # RID is in the response as "RID = <id>"
        match = re.search(r"RID = (\S+)", text)
        if not match:
            self.logger.warning("NCBI BLAST Put: RID not found in response")
            return None
        return match.group(1)

    def _ncbi_blast_check_ready(self, rid: str) -> bool:
        """Poll NCBI BLAST to see if result is ready."""
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo",
        }
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?" + urllib.parse.urlencode(params)

        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                text = resp.read().decode("utf-8")
        except Exception as e:
            self.logger.warning(f"NCBI BLAST Get(SearchInfo) failed: {e}")
            return False

        if "Status=WAITING" in text:
            return False
        if "Status=FAILED" in text or "Status=UNKNOWN" in text:
            self.logger.warning(f"NCBI BLAST status indicates failure/unknown for RID {rid}")
            return False
        if "Status=READY" in text:
            return True
        return False

    def _ncbi_blast_fetch_xml(self, rid: str) -> Optional[str]:
        """Fetch BLAST XML result from NCBI."""
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML",
        }
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?" + urllib.parse.urlencode(params)
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            self.logger.warning(f"NCBI BLAST Get(XML) failed: {e}")
            return None

    def remote_blastp(self, fasta_path: Path) -> Tuple[bool, int, Optional[str]]:
        """
        High-level remote BLASTp:
        - Reads FASTA
        - Submits to NCBI
        - Polls until ready (with a limit)
        - Parses XML to get:
          blast_hit (bool), blast_count (int), top_uniprot (str)
        """
        if not fasta_path or not fasta_path.exists():
            return False, 0, ""

        seq = fasta_path.read_text()
        rid = self._ncbi_blast_put(seq)
        if not rid:
            return False, 0, ""

        # Poll with limited attempts to respect NCBI usage
        max_attempts = 30
        wait_seconds = 20
        for attempt in range(max_attempts):
            if self._ncbi_blast_check_ready(rid):
                break
            time.sleep(wait_seconds)
        else:
            self.logger.warning(f"NCBI BLAST: RID {rid} not ready after {max_attempts} attempts")
            return False, 0, ""

        xml = self._ncbi_blast_fetch_xml(rid)
        if not xml:
            return False, 0, ""

        # Very simple parsing: count <Hit> elements, and extract first accession
        hits = re.findall(r"<Hit>.*?</Hit>", xml, flags=re.DOTALL)
        blast_count = len(hits)
        if blast_count == 0:
            return False, 0, ""

        top_block = hits[0]
        acc_match = re.search(r"<Hit_accession>([^<]+)</Hit_accession>", top_block)
        top_uniprot = acc_match.group(1) if acc_match else ""
        return True, blast_count, top_uniprot
    
    def _validate_file(self, file_path: Path) -> bool:
        """Validate file exists and is readable"""
        try:
            if not file_path.exists():
                self.logger.warning(f"File not found: {file_path}")
                return False
            if file_path.stat().st_size == 0:
                self.logger.warning(f"Empty file: {file_path}")
                return False
            return True
        except Exception as e:
            self.logger.error(f"File validation error: {e}")
            return False
    
    def _retry_with_backoff(self, func, *args, **kwargs):
        """Retry function with exponential backoff"""
        for attempt in range(self.config.max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if attempt == self.config.max_retries - 1:
                    raise e
                wait_time = min((2 ** attempt) * self.config.request_delay + random.uniform(0, 10), 300)
                self.logger.warning(f"Attempt {attempt + 1} failed: {e}. Retrying in {wait_time:.1f}s...")
                time.sleep(wait_time)
    
    def blastp_homology(self, fasta_file: Path) -> Tuple[bool, int, Optional[str]]:
        """
        BLASTP vs UniProt with caching.
        Supports both local BLAST+ and remote NCBI BLAST.
        """
        cache_key = self.cache.get_key("blastp", fasta_file)
        cached_result = self.cache.get(cache_key)
        if cached_result:
            self.logger.debug(f"BLASTP cache hit for {fasta_file}")
            return cached_result

        if not self._validate_file(fasta_file):
            return False, 0, None

        # Check configuration for BLAST mode
        if not self.config.enable_blast:
            self.logger.info(f"BLAST disabled in config; skipping homology search for {fasta_file.name}.")
            result_data = (False, 0, None)
            self.cache.set(cache_key, result_data)
            return result_data

        if self.config.blast_mode == "local":
            # Local BLAST+ mode
            self.logger.info(f"Using local BLAST+ for {fasta_file.name}")
            # Check if blastp is available
            blastp_check = subprocess.run("which blastp", shell=True, capture_output=True)
            if blastp_check.returncode != 0:
                self.logger.warning(f"blastp not found; falling back to remote BLAST for {fasta_file.name}")
                return self._remote_blastp_fallback(fasta_file)
            
            # TODO: Implement local BLAST+ logic here
            self.logger.info(f"Local BLAST+ not implemented yet; using remote BLAST for {fasta_file.name}")
            return self._remote_blastp_fallback(fasta_file)
        
        elif self.config.blast_mode == "remote":
            # Remote NCBI BLAST mode
            self.logger.info(f"Using remote NCBI BLAST for {fasta_file.name}")
            return self.remote_blastp(fasta_file)
        
        else:
            self.logger.warning(f"Unknown BLAST mode '{self.config.blast_mode}'; using remote BLAST for {fasta_file.name}")
            return self.remote_blastp(fasta_file)

    def _remote_blastp_fallback(self, fasta_file: Path) -> Tuple[bool, int, Optional[str]]:
        """Fallback to remote BLAST when local BLAST is not available"""
        self.logger.info(f"Using remote NCBI BLAST as fallback for {fasta_file.name}")
        return self.remote_blastp(fasta_file)

    def predict_secreted_plant(self, fasta_file: Path) -> bool:
        """Enhanced plant secreted prediction with multiple methods"""
        cache_key = self.cache.get_key("secreted", fasta_file)
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            self.logger.debug(f"Secreted prediction cache hit for {fasta_file}")
            return cached_result

        if not self._validate_file(fasta_file):
            return False

        gene_id = fasta_file.stem

        # Read sequence once (used by both SignalP6 and heuristics)
        seq = "".join(
            line.strip()
            for line in fasta_file.read_text().splitlines()
            if not line.startswith(">")
        )
        if not seq:
            self.logger.warning(f"No sequence found in {fasta_file}")
            return False

        # ------------------------------------------------------------------
        # Method 1: Local SignalP 6.0 CLI (replaces broken web API)
        # ------------------------------------------------------------------
        signalp_available = False
        try:
            # Check if signalp6 binary is available
            signalp_check = subprocess.run("which signalp6", shell=True, capture_output=True)
            if signalp_check.returncode == 0:
                signalp_available = True
        except Exception:
            pass

        if signalp_available:
            try:
                if len(seq) < 10:
                    raise ValueError("Sequence too short")

                import tempfile

                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                    tmp.write(f">{gene_id}\n{seq}\n")
                    tmp_path = tmp.name

                try:
                    result = subprocess.run(
                        ['signalp6', '-f', 'summary', tmp_path, '--org', 'eukarya', '--mode', 'fast'],
                        capture_output=True, text=True, timeout=60
                    )

                    # Clean up temp file
                    os.unlink(tmp_path)

                    if result.returncode == 0 and 'SP' in result.stdout and 'YES' in result.stdout.upper():
                        self.cache.set(cache_key, True)
                        self.logger.info(f"   💧 SignalP6: {gene_id} predicted as secreted")
                        return True
                    elif result.returncode != 0:
                        self.logger.debug(f"SignalP6 non-zero exit for {gene_id}: {result.stderr}")

                finally:
                    # Ensure cleanup if something went wrong before os.unlink above
                    try:
                        if os.path.exists(tmp_path):
                            os.unlink(tmp_path)
                    except Exception:
                        pass

            except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
                self.logger.warning(f"Local SignalP6 error for {gene_id}: {e}")
        else:
            if not ProteinAnalyzer._signalp_unavailable_logged:
                self.logger.warning("SignalP6 binary not found in PATH; using heuristics only")
                ProteinAnalyzer._signalp_unavailable_logged = True

        # ------------------------------------------------------------------
        # Method 2: Enhanced heuristics (backup)
        # ------------------------------------------------------------------
        try:
            # seq already read above; reuse it
            if not seq:
                return False

            # Enhanced plant signal peptide patterns
            n_term = seq[:30].upper()

            # Check for common plant signal peptide motifs
            signal_motifs = [
                'MKT', 'MKS', 'MKL', 'MAV', 'MGG', 'MVA', 'MSG', 'MAA',
                'MKC', 'MKN', 'MKQ', 'MKR', 'MKH', 'MKD', 'MKE'
            ]

            has_signal_motif = any(n_term.startswith(motif) for motif in signal_motifs)

            # Basic residue requirements
            basic_residues = sum(1 for aa in n_term if aa in ['R', 'K'])
            has_basic = basic_residues >= 2

            # Hydrophobic stretch analysis
            hydrophobic_window = n_term[:20]
            hydrophobic_count = sum(1 for aa in hydrophobic_window if aa in 'AVILMFYW')
            has_hydrophobic = hydrophobic_count >= 6

            # Transmembrane region check (negative indicator)
            has_transmembrane = any(pattern in n_term[:15] for pattern in ['LL', 'II', 'VV', 'FF'])

            # Length requirement
            min_length = len(seq) > 50

            is_secreted = (has_signal_motif and has_basic and has_hydrophobic and 
                          not has_transmembrane and min_length)

            if is_secreted:
                self.logger.info(f"   💧 Heuristics: {gene_id} predicted as secreted")

            self.cache.set(cache_key, is_secreted)
            return is_secreted

        except Exception as e:
            self.logger.error(f"Heuristic secreted prediction error for {gene_id}: {e}")
            return False

    def get_go_terms(self, uniprot_id: str) -> List[str]:
        """Enhanced GO terms fetching with caching"""
        if not uniprot_id:
            return []
        
        cache_key = self.cache.get_key("go_terms", uniprot_id)
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            return cached_result
        
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
            response = self._retry_with_backoff(
                requests.get, url, timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                xrefs = data.get('uniProtKBCrossReferences', [])
                go_list = [x.get('id') for x in xrefs if x.get('database') == 'GO']
                result = go_list[:5]
                self.cache.set(cache_key, result)
                return result
        except Exception as e:
            self.logger.warning(f"GO terms fetch error for {uniprot_id}: {e}")
        
        self.cache.set(cache_key, [])
        return []

    def alphafold_confidence(self, pdb_file: Path) -> bool:
        """Enhanced AlphaFold confidence analysis"""
        cache_key = self.cache.get_key("alphafold_confidence", pdb_file)
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            return cached_result
        
        if not self._validate_file(pdb_file):
            return False
        
        try:
            content = open(pdb_file).read()
            plddt_matches = re.findall(r'ATOM.{54}([\d.]+)', content)
            if plddt_matches:
                vals = [float(x) for x in plddt_matches if float(x) > 0]
                if vals:
                    avg = float(np.mean(vals))
                    result = avg > self.config.min_plddt
                    self.cache.set(cache_key, result)
                    return result
        except Exception as e:
            self.logger.warning(f"AlphaFold confidence analysis error for {pdb_file}: {e}")
        
        self.cache.set(cache_key, False)
        return False

    def detect_domains(self, foldseek_analysis: Dict[str, Any]) -> bool:
        """Enhanced domain detection.

        Currently uses the number of plant hits as a simple proxy.
        Handles both list- and int-based representations of plant_hits.
        """
        try:
            if not isinstance(foldseek_analysis, dict):
                return False

            plant_hits = foldseek_analysis.get("plant_hits", [])

            # If it's an int, treat as a count directly
            if isinstance(plant_hits, int):
                return plant_hits > 1

            # If it's a list/tuple, use its length
            if isinstance(plant_hits, (list, tuple)):
                return len(plant_hits) > 1

            # Anything else – log at debug and treat as 0
            self.logger.debug(
                "Domain detection: unexpected plant_hits type %s",
                type(plant_hits).__name__,
            )
            return False
        except Exception as e:
            self.logger.warning(f"Domain detection error: {e}")
            return False

    def fpocket_druggable(self, pdb_file: Path) -> Tuple[bool, Optional[float], Optional[float]]:
        """Enhanced fpocket druggability with fallback"""
        cache_key = self.cache.get_key("fpocket", pdb_file)
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            return cached_result
        
        if not self._validate_file(pdb_file):
            return False, None, None
        
        # Check if fpocket is installed
        fpocket_check = subprocess.run("which fpocket", shell=True, capture_output=True)
        if fpocket_check.returncode != 0:
            self.logger.warning(f"fpocket not installed - using structural heuristics for {pdb_file}")
            # Fallback: Use enhanced structural heuristics
            try:
                with open(pdb_file) as f:
                    content = f.read()
                if 'ATOM' in content:
                    plddt_matches = re.findall(r'^ATOM.{60}([\d.]+)', content, re.MULTILINE)
                    if plddt_matches:
                        plddt_values = [float(x) for x in plddt_matches if float(x) > 0]
                        if plddt_values:
                            avg_plddt = sum(plddt_values) / len(plddt_values)
                            is_druggable = avg_plddt > 45
                            result = (is_druggable, avg_plddt, None)
                            self.cache.set(cache_key, result)
                            self.logger.info(f"   💊 Fallback: Avg pLDDT: {avg_plddt:.1f} -> {'druggable' if is_druggable else 'uncertain'}")
                            return result
            except Exception as e:
                self.logger.warning(f"Structural analysis error for {pdb_file}: {e}")
            result = (False, None, None)
            self.cache.set(cache_key, result)
            return result
        
        try:
            gene_id = pdb_file.stem
            result = subprocess.run(
                f"fpocket -f {pdb_file}", 
                shell=True, timeout=self.config.fpocket_timeout, 
                capture_output=True
            )
            
            if result.returncode != 0:
                self.logger.warning(f"fpocket failed for {pdb_file}")
                result = (False, None, None)
                self.cache.set(cache_key, result)
                return result

            grep = subprocess.run(
                "grep -R \"Druggability score\" -n . | head -1",
                shell=True, capture_output=True, text=True
            )
            
            if grep.returncode == 0 and grep.stdout.strip():
                line = grep.stdout.strip().split(":", 1)[1]
                m = re.search(r"([0-9]*\.?[0-9]+)", line)
                if m:
                    score = float(m.group(1))
                    subprocess.run("rm -rf pocket* _fpocket*", shell=True)
                    is_druggable = score > self.config.min_druggability_score
                    result = (is_druggable, None, score)
                    self.cache.set(cache_key, result)
                    self.logger.info(f"   💊 fpocket: Score {score:.3f} -> {'druggable' if is_druggable else 'uncertain'}")
                    return result
            
            subprocess.run("rm -rf pocket* _fpocket*", shell=True)
            
        except Exception as e:
            self.logger.warning(f"fpocket error for {pdb_file}: {e}")
        
        result = (False, None, None)
        self.cache.set(cache_key, result)
        return result

    def get_headers(self) -> Dict[str, str]:
        return {"User-Agent": "Plant-Drug-Target-Pipeline/4.0"}

    def ligand_efficiency(self, pocket_volume: Optional[float], binding_score: Optional[float]) -> Optional[float]:
        """Enhanced ligand efficiency calculation"""
        try:
            if not pocket_volume or not binding_score:
                return None
            vol_norm = pocket_volume / 500.0
            score_norm = float(binding_score)
            le = vol_norm * score_norm
            return float(le)
        except Exception as e:
            self.logger.warning(f"Ligand efficiency calculation error: {e}")
            return None

    def kegg_pathway(self, uniprot_id: str) -> List[str]:
        """Enhanced KEGG pathway mapping"""
        if not uniprot_id:
            return []
        
        cache_key = self.cache.get_key("kegg_pathway", uniprot_id)
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            return cached_result
        
        try:
            url = f"https://rest.kegg.jp/find/genes/{uniprot_id}"
            response = self._retry_with_backoff(
                requests.get, url, timeout=10
            )
            
            if response.status_code != 200:
                return []
                
            lines = response.text.strip().split('\n')
            paths = []
            for line in lines:
                if not line:
                    continue
                parts = line.split('\t')
                if parts:
                    paths.append(parts[0])
            
            self.cache.set(cache_key, paths)
            return paths
            
        except Exception as e:
            self.logger.warning(f"KEGG pathway mapping error for {uniprot_id}: {e}")
            self.cache.set(cache_key, [])
            return []

    def sa_score(self, smiles: str) -> Optional[float]:
        """Enhanced synthetic accessibility score"""
        try:
            from rdkit import Chem
            from rdkit.Chem import QED
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            q = QED.qed(mol)
            sa = (1.0 - q) * 9.0 + 1.0
            return float(sa)
        except ImportError:
            self.logger.warning("RDKit not available for SA score calculation")
            return None
        except Exception as e:
            self.logger.warning(f"SA score calculation error: {e}")
            return None

    def submit_search_job(self, pdb_file_path: str) -> Optional[str]:
        """Enhanced Foldseek submission with better error handling"""
        url = "https://search.foldseek.com/api/ticket"
        
        for attempt in range(self.config.max_retries):
            try:
                with open(pdb_file_path, 'rb') as f:
                    files = {'q': (Path(pdb_file_path).name, f, 'application/octet-stream')}
                    data = [("mode", "3diaa")]
                    for db in self.config.plant_dbs:
                        data.append(("database[]", db))
                    
                    response = self._retry_with_backoff(
                        requests.post, url, files=files, data=data, 
                        headers=self.get_headers(), timeout=self.config.timeout
                    )
                    
                    if response.status_code == 200:
                        result = response.json()
                        job_id = result.get('id') or result.get('jobId')
                        if job_id:
                            self.logger.info(f"Foldseek job submitted: {job_id}")
                        return job_id
                    elif response.status_code == 429:
                        wait_time = min((2 ** attempt) * self.config.request_delay + random.uniform(0, 10), 300)
                        self.logger.warning(f"Rate limited, waiting {wait_time:.0f}s...")
                        time.sleep(wait_time)
                    else:
                        self.logger.warning(f"Foldseek submission failed: {response.status_code}")
                        
            except Exception as e:
                self.logger.warning(f"Foldseek submission error (attempt {attempt + 1}): {e}")
                if attempt == self.config.max_retries - 1:
                    raise e
                time.sleep(5)
        return None

    def check_job_status(self, job_id: str) -> str:
        """Enhanced Foldseek status checking"""
        url = f"https://search.foldseek.com/api/ticket/{job_id}"
        try:
            response = self._retry_with_backoff(
                requests.get, url, headers=self.get_headers(), timeout=self.config.timeout
            )
            if response.status_code == 200:
                data = response.json()
                status = data.get('status', 'UNKNOWN').upper()
                return 'COMPLETE' if status in ('COMPLETE', 'FINISHED') else status
        except Exception as e:
            self.logger.warning(f"Foldseek status check error: {e}")
        return 'ERROR'

    def download_results(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Enhanced Foldseek results download"""
        url = f"https://search.foldseek.com/api/result/{job_id}/0"
        try:
            response = self._retry_with_backoff(
                requests.get, url, headers=self.get_headers(), timeout=self.config.timeout
            )
            if response.status_code == 200:
                return response.json()
        except Exception as e:
            self.logger.warning(f"Foldseek results download error: {e}")
        return None

    def analyze_foldseek_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Enhanced plant-specific Foldseek template analysis"""
        self.logger.info("   🔍 Analyzing Foldseek results...")
        
        if not results or not isinstance(results, dict):
            self.logger.error("   ❌ No Foldseek results")
            return {
                'total_hits': 0,
                'plant_hits': 0,
                'top_tm_score': 0,
                'has_good_template': False,
                'best_plant_hit': None
            }
        
        self.logger.debug(f"   📊 Results keys: {list(results.keys())}")
        
        all_hits = []
        # Enhanced plant codes including POPTR (poplar)
        plant_codes = ['ARATH', 'ORYSJ', 'MAIZE', 'POTTC', 'SOLLC', 'SOLTU', 'POPTR']
        plant_names = ['arabidopsis', 'oryza sativa', 'rice', 'maize', 'potato', 'tomato', 'solanum', 'poplar', 'populus']
        
        db_results = results.get('results', [])
        
        # Check if results is in 'queries' format instead
        if not db_results and 'queries' in results:
            queries = results.get('queries', [])
            for query in queries:
                query_results = query.get('results', [])
                for qr in query_results:
                    alignments = qr.get('alignments', [])
                    if alignments:
                        for hit in alignments:
                            all_hits.append(hit)
        
        # Extract all hits from all databases (original format)
        for db_result in db_results:
            alignments = db_result.get('alignments', [])
            if alignments and alignments[0]:
                for hit in alignments[0][:10]:  # Top 10 hits
                    all_hits.append(hit)
        
        self.logger.info(f"   📈 Found {len(all_hits)} total hits")
        
        top_hits = all_hits[:10] if all_hits else []
        
        # Check for plant hits - use taxName field
        plant_hits = []
        for h in top_hits:
            target = h.get('target', '')
            tax_name = h.get('taxName', '').lower()
            if any(code in str(target).upper() for code in plant_codes) or \
               any(name in tax_name for name in plant_names):
                plant_hits.append(h)
        
        top_hit = top_hits[0] if top_hits else None
        
        # Enhanced score extraction - check multiple fields
        tm_score = 0
        if top_hit:
            for field in ['alntmscore', 'tmscore', 'prob', 'alnprob', 'seqId']:
                if field in top_hit:
                    score = top_hit[field]
                    if field.endswith('Id'):
                        tm_score = score / 100.0
                    else:
                        tm_score = score
                    break
        
        # Also check plant hit for score if main top hit doesn't have good score
        if tm_score < self.config.min_tm_score and plant_hits:
            ph = plant_hits[0]
            for field in ['alntmscore', 'tmscore', 'prob', 'alnprob', 'seqId']:
                if field in ph:
                    score = ph[field]
                    if field.endswith('Id'):
                        tm_score = score / 100.0
                    else:
                        tm_score = score
                    break
        
        # Relaxed threshold for hypothetical proteins
        has_good_template = tm_score > self.config.min_tm_score
        
        self.logger.info(f"   🧬 Top TM-score: {tm_score:.3f} | Plant hits: {len(plant_hits)}")
        
        return {
            'total_hits': len(all_hits),
            'plant_hits': len(plant_hits),
            'top_tm_score': tm_score,
            'has_good_template': has_good_template,
            'best_plant_hit': plant_hits[0] if plant_hits else None
        }

    def full_drug_target_analysis(self, fasta_file: Optional[Path], pdb_file: Path, foldseek_results: Dict[str, Any]) -> Dict[str, Any]:
        """Enhanced complete 4-criteria analysis"""
        gene_id = pdb_file.stem
        
        # 1. BLAST homology
        blast_hit = False
        blast_count = 0
        top_uniprot = ""
        
        if self.config.enable_blast and fasta_file is not None:
            if self.config.blast_mode == "remote":
                blast_hit, blast_count, top_uniprot = self.remote_blastp(fasta_file)
            else:
                blast_hit, blast_count, top_uniprot = self.blastp_homology(fasta_file)  # local blastp version
        
        # 2. Foldseek template
        foldseek_analysis = self.analyze_foldseek_results(foldseek_results)
        foldseek_template = foldseek_analysis.get('has_good_template', False)
        
        # 3. Secreted prediction
        secreted = self.predict_secreted_plant(fasta_file) if fasta_file else False
        
        # 4. Druggable pocket
        druggable_bool, pocket_volume, fpocket_score = self.fpocket_druggable(pdb_file)

        # Additional enhancements
        go_terms = self.get_go_terms(top_uniprot) if top_uniprot else []
        high_confidence = self.alphafold_confidence(pdb_file)
        multi_domain = self.detect_domains(foldseek_analysis)
        le = self.ligand_efficiency(pocket_volume, fpocket_score)
        kegg_paths = self.kegg_pathway(top_uniprot) if top_uniprot else []
        
        # Enhanced DRUG TARGET PRIORITY SCORE (0-5)
        priority_score = sum([
            1 if blast_hit else 0,
            1 if foldseek_template else 0,
            1 if secreted else 0,
            1 if druggable_bool else 0,
            1 if high_confidence else 0
        ])
        
        # Enhanced score with plant_hits bonus
        enhanced_score = priority_score + (foldseek_analysis.get('plant_hits', 0) * 0.5)
        
        return {
            'gene_id': gene_id,
            'blast_hit': blast_hit,
            'blast_count': blast_count,
            'top_uniprot': top_uniprot,
            'foldseek_template': foldseek_template,
            'foldseek_tm': foldseek_analysis.get('top_tm_score', 0),
            'secreted': secreted,
            'druggable': druggable_bool,
            'high_confidence': high_confidence,
            'priority_score': priority_score,
            'enhanced_score': enhanced_score,
            'foldseek_analysis': foldseek_analysis,
            'go_terms': go_terms,
            'multi_domain': multi_domain,
            'pocket_volume': pocket_volume,
            'fpocket_score': fpocket_score,
            'ligand_efficiency': le,
            'kegg_pathways': kegg_paths,
            'analysis_timestamp': datetime.now().isoformat()
        }

    def process_single_protein(self, fasta_file: Optional[Path], pdb_file: Path) -> Dict[str, Any]:
        """Enhanced complete pipeline for one protein with progress tracking"""
        gene_id = pdb_file.stem
        
        self.logger.info(f"\n🔬 [{gene_id}] Running FULL analysis...")
        
        # 1. Submit Foldseek job
        job_id = self.submit_search_job(str(pdb_file))
        if not job_id:
            return {'gene_id': gene_id, 'status': 'FOLDSEEK_FAILED'}
        
        self.logger.info(f"   Job {job_id} submitted. Waiting...")
        
        # 2. Poll Foldseek with progress tracking
        for attempt in range(self.config.max_poll_attempts):
            status = self.check_job_status(job_id)
            if status == 'COMPLETE':
                break
            elif status == 'ERROR':
                return {'gene_id': gene_id, 'status': 'FOLDSEEK_ERROR'}
            self.logger.info(f"   ⏳ Foldseek running... ({attempt*self.config.poll_interval//60}m)")
            time.sleep(self.config.poll_interval)
        
        # 3. Download + analyze ALL 4 criteria
        foldseek_results = self.download_results(job_id)
        if not foldseek_results:
            return {'gene_id': gene_id, 'status': 'NO_FOLDSEEK_RESULTS'}
        
        full_analysis = self.full_drug_target_analysis(fasta_file, pdb_file, foldseek_results)
        full_analysis['foldseek_results'] = foldseek_results
        full_analysis['status'] = 'COMPLETE'
        
        return full_analysis

# ============================================================================
# MAIN PIPELINE - ENHANCED
# ============================================================================

class PipelineManager:
    """Enhanced pipeline manager with parallel processing and checkpointing"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        self.logger = PipelineLogger(config)
        self.cache = CacheManager(config)
        self.analyzer = ProteinAnalyzer(config, self.logger, self.cache)
        
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.results_file = self.output_dir / "all_results.json"
        
    def validate_inputs(self) -> bool:
        """Enhanced input validation"""
        alphafold_dir = Path(self.config.alphafold_output_dir)
        if not alphafold_dir.exists():
            self.logger.error(f"❌ AlphaFold dir missing: {alphafold_dir}")
            return False
        
        pdb_files = list(alphafold_dir.glob("*.pdb"))
        if not pdb_files:
            self.logger.error("❌ No PDB files found!")
            return False
        
        fasta_dir = Path(self.config.fasta_dir)
        if fasta_dir.exists():
            fasta_files = list(fasta_dir.glob("*.fasta"))
            self.logger.info(f"🔍 Found {len(pdb_files)} PDBs + {len(fasta_files)} FASTAs")
        else:
            self.logger.info(f"🔍 Found {len(pdb_files)} PDBs (no FASTA directory)")
        
        return True
    
    def load_checkpoint(self) -> Tuple[List[Dict[str, Any]], List[Path], List[Optional[Path]]]:
        """Load checkpoint data"""
        if not self.checkpoint_file.exists():
            return [], [], []
        
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
            
            completed_results = checkpoint.get('completed_results', [])
            remaining_pdb_files = [Path(p) for p in checkpoint.get('remaining_pdb_files', [])]
            remaining_fasta_files = [Path(p) if p else None for p in checkpoint.get('remaining_fasta_files', [])]
            
            self.logger.info(f"Loaded checkpoint: {len(completed_results)} completed, {len(remaining_pdb_files)} remaining")
            return completed_results, remaining_pdb_files, remaining_fasta_files
            
        except Exception as e:
            self.logger.warning(f"Failed to load checkpoint: {e}")
            return [], [], []
    
    def save_checkpoint(self, completed_results: List[Dict[str, Any]], 
                       remaining_pdb_files: List[Path], remaining_fasta_files: List[Optional[Path]]):
        """Save checkpoint data"""
        try:
            checkpoint = {
                'completed_results': completed_results,
                'remaining_pdb_files': [str(p) for p in remaining_pdb_files],
                'remaining_fasta_files': [str(p) if p else None for p in remaining_fasta_files],
                'timestamp': datetime.now().isoformat()
            }
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint, f, indent=2)
                
        except Exception as e:
            self.logger.warning(f"Failed to save checkpoint: {e}")
    
    def process_protein_batch(self, batch_data: List[Tuple[Optional[Path], Path]]) -> List[Dict[str, Any]]:
        """Process a batch of proteins with progress tracking"""
        batch_results = []
        
        for fasta_file, pdb_file in batch_data:
            try:
                result = self.analyzer.process_single_protein(fasta_file, pdb_file)
                batch_results.append(result)
                
                # Save individual results
                if result.get('status') == 'COMPLETE':
                    result_file = self.output_dir / f"{result['gene_id']}_full_analysis.json"
                    with open(result_file, 'w') as f:
                        json.dump(result, f, indent=2)
                    self._print_status_icons(result)
                
            except Exception as e:
                self.logger.error(f"Error processing {pdb_file}: {e}")
                batch_results.append({
                    'gene_id': pdb_file.stem,
                    'status': 'ERROR',
                    'error_message': str(e)
                })
        
        return batch_results
    
    def _print_status_icons(self, result: Dict[str, Any]):
        """Enhanced visual status display"""
        icons = {
            'blast_hit': '🔍',
            'foldseek_template': '🧬', 
            'secreted': '💧',
            'druggable': '💊',
            'high_confidence': '⭐'
        }
        status = ['❌' if not result.get(k) else v for k, v in icons.items()]
        score = result.get('priority_score', 0)
        enhanced = result.get('enhanced_score', 0)
        
        # Enhanced verdict messages
        verdict = ""
        if score >= 4:
            verdict = " 🟢 WET LAB NOW!"
        elif score >= 3:
            verdict = " 🟡 HIGH VALUE"
        
        extra = []
        if result.get('ligand_efficiency') is not None:
            extra.append(f"LE:{result.get('ligand_efficiency'):.2f}")
        if result.get('fpocket_score') is not None:
            extra.append(f"FP:{result.get('fpocket_score'):.2f}")
        extras = (" | " + ",".join(extra)) if extra else ""
        
        # Show enhanced score indicator if > 4
        bonus_star = "⭐" if enhanced > 4 else ""
        
        self.logger.info(f"   {''.join(status)} [{score}/5]{extras} {verdict} {bonus_star}")

    def _check_external_dependencies(self):
        """Log availability of key external tools."""
        tools: List[str] = []

        # Only require local blastp if BLAST is enabled and mode is local
        if self.config.enable_blast and self.config.blast_mode == "local":
            tools.append("blastp")

        # fpocket remains optional
        tools.append("fpocket")

        for tool in tools:
            path = shutil.which(tool)
            if path:
                self.logger.info(f"✅ Found {tool} at {path}")
            else:
                self.logger.warning(f"⚠️ {tool} not found in PATH; related steps will be skipped or use fallbacks.")

    def run_parallel(self) -> List[Dict[str, Any]]:
        """Run pipeline with parallel processing"""
        alphafold_dir = Path(self.config.alphafold_output_dir)
        pdb_files = list(alphafold_dir.glob("*.pdb"))
        
        fasta_dir = Path(self.config.fasta_dir)
        fasta_files = list(fasta_dir.glob("*.fasta")) if fasta_dir.exists() else []
        
        # Create file mapping
        file_pairs = []
        for pdb_file in pdb_files:
            gene_id = pdb_file.stem
            fasta_file = fasta_dir / f"{gene_id}.fasta" if fasta_dir.exists() else None
            file_pairs.append((fasta_file, pdb_file))
        
        # Load checkpoint
        completed_results, remaining_pdb_files, remaining_fasta_files = self.load_checkpoint()
        
        # If no checkpoint or nothing recorded, start from scratch
        if not completed_results and not remaining_pdb_files:
            # First run - process all files
            remaining_pdb_files = [pair[1] for pair in file_pairs]
            remaining_fasta_files = [pair[0] for pair in file_pairs]
        
        # Process in batches
        all_results = completed_results.copy()
        
        with tqdm(total=len(remaining_pdb_files), desc="Processing proteins") as pbar:
            for i in range(0, len(remaining_pdb_files), self.config.max_workers):
                batch_pdb = remaining_pdb_files[i:i + self.config.max_workers]
                batch_fasta = remaining_fasta_files[i:i + self.config.max_workers]
                batch_data = list(zip(batch_fasta, batch_pdb))
                
                # Process batch
                batch_results = self.process_protein_batch(batch_data)
                all_results.extend(batch_results)
                
                # Update progress
                pbar.update(len(batch_pdb))
                
                # Save checkpoint periodically
                if len(all_results) % self.config.checkpoint_interval == 0:
                    remaining_pdb = remaining_pdb_files[i + self.config.max_workers:]
                    remaining_fasta = remaining_fasta_files[i + self.config.max_workers:]
                    self.save_checkpoint(all_results, remaining_pdb, remaining_fasta)
                
                # Rate limiting: scale sleep with batch size
                if i + self.config.max_workers < len(remaining_pdb_files):
                    delay = self.config.request_delay * (len(batch_pdb) / max(1, self.config.max_workers))
                    self.logger.info(f"💤 Rate limiting for {delay:.1f}s...")
                    time.sleep(delay)
        
        # Save final results
        with open(self.results_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        return all_results

    def run_sequential(self) -> List[Dict[str, Any]]:
        """Run pipeline sequentially"""
        alphafold_dir = Path(self.config.alphafold_output_dir)
        pdb_files = list(alphafold_dir.glob("*.pdb"))
        
        fasta_dir = Path(self.config.fasta_dir)
        fasta_files = list(fasta_dir.glob("*.fasta")) if fasta_dir.exists() else []
        
        # Load checkpoint
        completed_results, remaining_pdb_files, remaining_fasta_files = self.load_checkpoint()
        
        # If no checkpoint or nothing recorded, start from scratch
        if not completed_results and not remaining_pdb_files:
            # First run - process all files
            remaining_pdb_files = pdb_files
            remaining_fasta_files = [fasta_dir / f"{pdb.stem}.fasta" if fasta_dir.exists() else None for pdb in pdb_files]
        
        all_results = completed_results.copy()
        
        for idx, (fasta_file, pdb_file) in enumerate(zip(remaining_fasta_files, remaining_pdb_files)):
            gene_id = pdb_file.stem
            
            self.logger.info(f"\n[{idx + len(completed_results) + 1}/{len(pdb_files)}] {gene_id}")
            
            if idx > 0:
                self.logger.info("💤 Rate limiting...")
                time.sleep(self.config.request_delay)
            
            # FULL PIPELINE
            result = self.analyzer.process_single_protein(fasta_file, pdb_file)
            all_results.append(result)
            
            # Save individual results
            if result.get('status') == 'COMPLETE':
                result_file = self.output_dir / f"{gene_id}_full_analysis.json"
                with open(result_file, 'w') as f:
                    json.dump(result, f, indent=2)
                self._print_status_icons(result)
            
            # Save checkpoint periodically
            if len(all_results) % self.config.checkpoint_interval == 0:
                remaining_pdb = remaining_pdb_files[idx + 1:]
                remaining_fasta = remaining_fasta_files[idx + 1:]
                self.save_checkpoint(all_results, remaining_pdb, remaining_fasta)
        
        # Save final results
        with open(self.results_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        return all_results

    def run(self) -> List[Dict[str, Any]]:
        """Main pipeline execution"""
        self.logger.info("=" * 90)
        self.logger.info("🏆 PLANT PROTEIN DRUG DISCOVERY PIPELINE v4.0")
        self.logger.info("BLAST + Foldseek + Secreted + Druggable → PRIORITIZED TARGETS")
        self.logger.info("=" * 90)

        self._check_external_dependencies()

        if not self.validate_inputs():
            sys.exit(1)

        # Choose execution mode
        if self.config.enable_parallel and self.config.max_workers > 1:
            self.logger.info("🚀 Running in parallel mode")
            all_results = self.run_parallel()
        else:
            self.logger.info("🚀 Running in sequential mode")
            all_results = self.run_sequential()

        # FINAL DRUG TARGET RANKING
        self.generate_drug_target_report(all_results)

        self.logger.info(f"\n🎉 PIPELINE COMPLETE! {self.output_dir.absolute()}")

        # Generate comprehensive visualizations
        self.logger.info("\n📊 Generating visualizations...")
        self.generate_visualizations(all_results)

        return all_results

    def generate_drug_target_report(self, all_results: List[Dict[str, Any]]):
        """Enhanced drug target report generation"""
        
        # Filter perfect candidates (score=4)
        perfect_targets = [r for r in all_results if r.get('priority_score') == 4]
        high_priority = [r for r in all_results if r.get('priority_score') >= 3]
        
        self.logger.info("\n" + "="*90)
        self.logger.info("🏆 FINAL DRUG TARGET RANKING")
        self.logger.info("="*90)
        
        self.logger.info(f"\n🟢 PERFECT TARGETS (4/4 criteria): {len(perfect_targets)}")
        for target in perfect_targets:
            self.logger.info(f"  💎 {target['gene_id']:20s} | "
                           f"TM:{target['foldseek_tm']:.3f} | "
                           f"Blast:{target['blast_count']} | "
                           f"✅ ALL CRITERIA!")
        
        self.logger.info(f"\n🟡 HIGH PRIORITY (3/4): {len(high_priority)}")
        
        # Save JSON report
        report = {
            'perfect_targets': perfect_targets,
            'high_priority': high_priority,
            'all_results': all_results,
            'config': asdict(self.config),
            'timestamp': datetime.now().isoformat()
        }
        report_file = self.output_dir / "DRUG_TARGETS_PERFECT.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Create Excel report with enhanced format
        excel_data = []
        for result in all_results:
            foldseek_analysis = result.get('foldseek_analysis', {})
            best_plant = foldseek_analysis.get('best_plant_hit', {})
            
            excel_data.append({
                'Gene_ID': result.get('gene_id', ''),
                'Priority_Score': result.get('priority_score', 0),
                'Enhanced_Score': result.get('enhanced_score', 0),
                'BLAST_Hit': result.get('blast_hit', False),
                'BLAST_Count': result.get('blast_count', 0),
                'Top_UniProt': result.get('top_uniprot', ''),
                'Foldseek_Template': result.get('foldseek_template', False),
                'Foldseek_TM_Score': result.get('foldseek_tm', 0),
                'Total_Hits': foldseek_analysis.get('total_hits', 0),
                'Plant_Hits': foldseek_analysis.get('plant_hits', 0),
                'Best_Plant_Template': best_plant.get('target', '') if best_plant else '',
                'Best_Plant_Taxonomy': best_plant.get('taxName', '') if best_plant else '',
                'Secreted': result.get('secreted', False),
                'Druggable': result.get('druggable', False),
                'High_Confidence_AF': result.get('high_confidence', False),
                'Multi_Domain': result.get('multi_domain', False),
                'pLDDT_Avg': result.get('pocket_volume', ''),
                'FPocket_Score': result.get('fpocket_score', ''),
                'Ligand_Efficiency': result.get('ligand_efficiency', ''),
                'GO_Terms': ', '.join(result.get('go_terms', [])),
                'KEGG_Pathways': ', '.join(result.get('kegg_pathways', [])),
                'Status': result.get('status', ''),
                'Analysis_Timestamp': result.get('analysis_timestamp', '')
            })
        
        # Create DataFrame and save to Excel
        df = pd.DataFrame(excel_data)
        excel_file = self.output_dir / "DRUG_TARGETS_REPORT.xlsx"
        df.to_excel(excel_file, index=False, engine='openpyxl')
        
        self.logger.info(f"\n📊 JSON report: {report_file}")
        self.logger.info(f"📊 Excel report: {excel_file}")
        self.logger.info("="*90)

    def generate_visualizations(self, all_results: List[Dict[str, Any]]):
        """Enhanced visualization suite"""
        self.logger.info("   📈 Creating visualization suite...")
        
        # Filter complete results
        complete_results = [r for r in all_results if r.get('status') == 'COMPLETE']
        
        if not complete_results:
            self.logger.warning("   ❌ No complete results to visualize")
            return
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Priority Distribution
        self.viz_01_priority_distribution(complete_results)
        
        # 2. Confidence Analysis
        self.viz_02_confidence_all(complete_results)
        
        # 3. Top 10 Candidates
        self.viz_03_top10_candidates(complete_results)
        
        # 4. Score vs Probability
        self.viz_04_score_vs_prob(complete_results)
        
        # 5. Heatmap
        self.viz_05_heatmap(complete_results)
        
        # 6. Statistics
        self.viz_07_statistics(complete_results)
        
        # 7. Comprehensive Report
        self.viz_08_comprehensive_report(complete_results)
        
        self.logger.info("   ✅ All visualizations generated!")

    def viz_01_priority_distribution(self, results: List[Dict[str, Any]]):
        """Priority score distribution"""
        scores = [r.get('priority_score', 0) for r in results]
        
        plt.figure(figsize=self.config.figsize_medium)
        sns.histplot(scores, bins=range(0, 6), kde=False, alpha=0.7)
        plt.title('Drug Target Priority Score Distribution', fontsize=14, fontweight='bold')
        plt.xlabel('Priority Score (0-5)', fontsize=12)
        plt.ylabel('Number of Proteins', fontsize=12)
        plt.grid(axis='y', alpha=0.3)
        
        # Add count labels
        for i, count in enumerate([scores.count(i) for i in range(6)]):
            if count > 0:
                plt.text(i, count + 0.5, str(count), ha='center', fontsize=10)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'viz_01_priority_distribution.png', dpi=self.config.dpi, bbox_inches='tight')
        plt.close()

    def viz_02_confidence_all(self, results: List[Dict[str, Any]]):
        """AlphaFold confidence analysis"""
        confidences = [r.get('high_confidence', False) for r in results]
        confidence_counts = {'High Confidence': sum(confidences), 'Low Confidence': len(confidences) - sum(confidences)}
        
        plt.figure(figsize=self.config.figsize_small)
        colors = ['#2ecc71', '#e74c3c']
        wedges, texts, autotexts = plt.pie(confidence_counts.values(), labels=confidence_counts.keys(), 
                                          autopct='%1.1f%%', colors=colors, startangle=90)
        
        plt.setp(autotexts, size=12, weight='bold')
        plt.setp(texts, size=12)
        plt.title('AlphaFold Structure Confidence Distribution', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'viz_02_confidence_all.png', dpi=self.config.dpi, bbox_inches='tight')
        plt.close()

    def viz_03_top10_candidates(self, results: List[Dict[str, Any]]):
        """Top 10 drug target candidates"""
        # Sort by enhanced score
        sorted_results = sorted(results, key=lambda x: x.get('enhanced_score', 0), reverse=True)[:10]
        
        genes = [r.get('gene_id', 'Unknown')[:15] for r in sorted_results]
        scores = [r.get('enhanced_score', 0) for r in sorted_results]
        priorities = [r.get('priority_score', 0) for r in sorted_results]
        
        plt.figure(figsize=self.config.figsize_large)
        bars = plt.bar(range(len(genes)), scores, alpha=0.7, color='skyblue', edgecolor='navy')
        
        # Add priority score labels
        for i, (bar, priority) in enumerate(zip(bars, priorities)):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1, 
                    f'P:{priority}', ha='center', va='bottom', fontweight='bold')
        
        plt.title('Top 10 Drug Target Candidates', fontsize=14, fontweight='bold')
        plt.xlabel('Gene ID', fontsize=12)
        plt.ylabel('Enhanced Score', fontsize=12)
        plt.xticks(range(len(genes)), genes, rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'viz_03_top10_candidates.png', dpi=self.config.dpi, bbox_inches='tight')
        plt.close()

    def viz_04_score_vs_prob(self, results: List[Dict[str, Any]]):
        """Score vs probability scatter plot"""
        tm_scores = [r.get('foldseek_tm', 0) for r in results]
        priorities = [r.get('priority_score', 0) for r in results]
        
        plt.figure(figsize=self.config.figsize_medium)
        scatter = plt.scatter(tm_scores, priorities, alpha=0.6, s=60, c=priorities, cmap='viridis')
        
        plt.title('TM-Score vs Priority Score', fontsize=14, fontweight='bold')
        plt.xlabel('Foldseek TM-Score', fontsize=12)
        plt.ylabel('Priority Score', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Priority Score', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'viz_04_score_vs_prob.png', dpi=self.config.dpi, bbox_inches='tight')
        plt.close()

    def viz_05_heatmap(self, results: List[Dict[str, Any]]):
        """Criteria satisfaction heatmap"""
        # Create criteria matrix
        criteria = ['blast_hit', 'foldseek_template', 'secreted', 'druggable', 'high_confidence']
        matrix = []
        
        for r in results:
            row = [1 if r.get(c, False) else 0 for c in criteria]
            matrix.append(row)
        
        matrix = np.array(matrix)
        
        plt.figure(figsize=self.config.figsize_large)
        sns.heatmap(matrix.T, 
                    xticklabels=[r.get('gene_id', 'Unknown')[:10] for r in results],
                    yticklabels=['BLAST', 'Template', 'Secreted', 'Druggable', 'High Conf'],
                    cmap='RdYlGn', cbar_kws={'label': 'Criterion Met'})
        
        plt.title('Drug Target Criteria Satisfaction Heatmap', fontsize=14, fontweight='bold')
        plt.xlabel('Proteins', fontsize=12)
        plt.ylabel('Criteria', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'viz_05_heatmap.png', dpi=self.config.dpi, bbox_inches='tight')
        plt.close()

    def viz_07_statistics(self, results: List[Dict[str, Any]]):
        """Statistical summary"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Criteria satisfaction rates
        criteria = ['BLAST', 'Template', 'Secreted', 'Druggable', 'High Conf']
        rates = [
            sum(1 for r in results if r.get('blast_hit', False)) / len(results) * 100,
            sum(1 for r in results if r.get('foldseek_template', False)) / len(results) * 100,
            sum(1 for r in results if r.get('secreted', False)) / len(results) * 100,
            sum(1 for r in results if r.get('druggable', False)) / len(results) * 100,
            sum(1 for r in results if r.get('high_confidence', False)) / len(results) * 100
        ]
        
        ax1.bar(criteria, rates, color=['#3498db', '#2ecc71', '#f39c12', '#e74c3c', '#9b59b6'])
        ax1.set_title('Criteria Satisfaction Rates', fontweight='bold')
        ax1.set_ylabel('Percentage (%)')
        ax1.tick_params(axis='x', rotation=45)
        
        # 2. Priority score distribution
        scores = [r.get('priority_score', 0) for r in results]
        ax2.hist(scores, bins=range(0, 6), alpha=0.7, color='skyblue', edgecolor='black')
        ax2.set_title('Priority Score Distribution', fontweight='bold')
        ax2.set_xlabel('Priority Score')
        ax2.set_ylabel('Count')
        
        # 3. TM-score distribution
        tm_scores = [r.get('foldseek_tm', 0) for r in results if r.get('foldseek_tm', 0) > 0]
        ax3.hist(tm_scores, bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
        ax3.set_title('TM-Score Distribution', fontweight='bold')
        ax3.set_xlabel('TM-Score')
        ax3.set_ylabel('Count')
        
        # 4. Summary statistics
        stats_text = f"""
        Total Proteins: {len(results)}
        Perfect Targets (5/5): {sum(1 for r in results if r.get('priority_score', 0) == 5)}
        High Priority (4+): {sum(1 for r in results if r.get('priority_score', 0) >= 4)}
        Average TM-Score: {np.mean([r.get('foldseek_tm', 0) for r in results]):.3f}
        High Confidence: {sum(1 for r in results if r.get('high_confidence', False))}
        """
        
        ax4.text(0.1, 0.5, stats_text, fontsize=12, verticalalignment='center',
                 bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        ax4.set_title('Summary Statistics', fontweight='bold')
        ax4.axis('off')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'viz_07_statistics.png', dpi=self.config.dpi, bbox_inches='tight')
        plt.close()

    def viz_08_comprehensive_report(self, results: List[Dict[str, Any]]):
        """Generate enhanced HTML comprehensive report"""
        perfect_targets = [r for r in results if r.get('priority_score', 0) >= 4]
        high_priority = [r for r in results if r.get('priority_score', 0) >= 3]
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Plant Drug Target Analysis Report - Enhanced v4.0</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #2c3e50; color: white; padding: 20px; text-align: center; }}
                .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
                .target {{ background-color: #ecf0f1; padding: 10px; margin: 5px 0; border-radius: 3px; }}
                .perfect {{ background-color: #d5f4e6; }}
                .high {{ background-color: #fff3cd; }}
                table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                .metric {{ font-weight: bold; color: #2c3e50; }}
                .footer {{ text-align: center; color: #7f8c8d; font-size: 12px; margin-top: 20px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>🌱 Plant Drug Target Analysis Report - Enhanced v4.0</h1>
                <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>Configuration:</strong> Parallel: {self.config.enable_parallel}, Workers: {self.config.max_workers}, Caching: {self.config.enable_caching}</p>
            </div>
            
            <div class="section">
                <h2>📊 Summary</h2>
                <p><strong>Total Proteins Analyzed:</strong> <span class="metric">{len(results)}</span></p>
                <p><strong>Perfect Targets (4+/5):</strong> <span class="metric">{len(perfect_targets)}</span></p>
                <p><strong>High Priority (3+/5):</strong> <span class="metric">{len(high_priority)}</span></p>
                <p><strong>Average Priority Score:</strong> <span class="metric">{np.mean([r.get('priority_score', 0) for r in results]):.2f}</span></p>
                <p><strong>Average TM-Score:</strong> <span class="metric">{np.mean([r.get('foldseek_tm', 0) for r in results]):.3f}</span></p>
            </div>
            
            <div class="section">
                <h2>🟢 Perfect Targets</h2>
        """
        
        for target in perfect_targets:
            html_content += f"""
                <div class="target perfect">
                    <strong>{target.get('gene_id', 'Unknown')}</strong> | 
                    Score: {target.get('priority_score', 0)}/5 | 
                    Enhanced: {target.get('enhanced_score', 0):.1f} | 
                    TM-Score: {target.get('foldseek_tm', 0):.3f} | 
                    UniProt: {target.get('top_uniprot', 'Novel')} |
                    Status: {target.get('status', 'Unknown')}
                </div>
            """
        
        html_content += """
            </div>
            
            <div class="section">
                <h2>📈 Available Visualizations</h2>
                <ul>
                    <li>viz_01_priority_distribution.png - Priority score distribution</li>
                    <li>viz_02_confidence_all.png - AlphaFold confidence analysis</li>
                    <li>viz_03_top10_candidates.png - Top 10 candidates</li>
                    <li>viz_04_score_vs_prob.png - Score vs probability</li>
                    <li>viz_05_heatmap.png - Criteria satisfaction heatmap</li>
                    <li>viz_07_statistics.png - Statistical summary</li>
                </ul>
            </div>
            
            <div class="section">
                <h2>⚙️ Pipeline Configuration</h2>
                <table>
                    <tr><th>Parameter</th><th>Value</th></tr>
                    <tr><td>Plant DBs</td><td>""" + ", ".join(self.config.plant_dbs) + """</td></tr>
                    <tr><td>Request Delay</td><td>{self.config.request_delay}s</td></tr>
                    <tr><td>Max Retries</td><td>{self.config.max_retries}</td></tr>
                    <tr><td>Min TM Score</td><td>{self.config.min_tm_score}</td></tr>
                    <tr><td>Min PLDDT</td><td>{self.config.min_plddt}</td></tr>
                    <tr><td>Parallel Workers</td><td>{self.config.max_workers}</td></tr>
                    <tr><td>Caching Enabled</td><td>{self.config.enable_caching}</td></tr>
                </table>
            </div>
            
            <div class="footer">
                <p>Enhanced Plant Drug Discovery Pipeline v4.0 - Production Ready</p>
            </div>
        </body>
        </html>
        """
        
        with open(self.output_dir / 'viz_08_comprehensive_report.html', 'w') as f:
            f.write(html_content)

def load_config(config_file: Optional[str] = None) -> PipelineConfig:
    """Load configuration from file or use defaults"""
    if config_file and Path(config_file).exists():
        try:
            with open(config_file, 'r') as f:
                config_dict = json.load(f)
            return PipelineConfig(**config_dict)
        except Exception as e:
            print(f"Warning: Failed to load config file: {e}. Using defaults.")
    
    return PipelineConfig()

def main():
    """Enhanced main function with argument parsing"""
    parser = argparse.ArgumentParser(description='Plant Drug Discovery Pipeline v4.0')
    parser.add_argument('--config', type=str, help='Configuration file path')
    parser.add_argument('--parallel', action='store_true', help='Enable parallel processing')
    parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers')
    parser.add_argument('--no-cache', action='store_true', help='Disable caching')
    parser.add_argument('--output', type=str, help='Output directory')
    parser.add_argument('--alphafold-dir', type=str, help='AlphaFold output directory')
    parser.add_argument('--fasta-dir', type=str, help='FASTA directory')
    parser.add_argument('--blast-mode', type=str, choices=['local', 'remote'],
                        help='BLAST mode: local (blastp binary) or remote (NCBI web API)')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Override config with command line arguments
    if args.parallel:
        config.enable_parallel = True
    if args.workers:
        config.max_workers = args.workers
    if args.no_cache:
        config.enable_caching = False
    if args.output:
        config.output_dir = args.output
    if args.alphafold_dir:
        config.alphafold_output_dir = args.alphafold_dir
    if args.fasta_dir:
        config.fasta_dir = args.fasta_dir
    if getattr(args, "blast_mode", None):
        config.blast_mode = args.blast_mode
    
    # Run pipeline
    manager = PipelineManager(config)
    results = manager.run()
    
    return results

if __name__ == "__main__":
    main()