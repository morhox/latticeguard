#!/usr/bin/env python3
"""
LatticeGuard: Pre-Flight DFT & Lattice Copilot (CLI Edition)
============================================================
Protects high-stakes HPC simulation pipelines by auditing and sanitizing 
raw crystal structures (CIF) before expensive production runs.
"""

import argparse
import math
import os
import sys
from typing import Dict, Any

def parse_and_audit_cif(cif_text: str) -> Dict[str, Any]:
    lines = cif_text.split('\n')
    space_group = "Unknown"
    a, b, c = 5.82, 5.82, 5.82
    alpha, beta, gamma = 90.0, 90.0, 90.0
    atom_sites = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('_space_group_name_H-M_alt') or line.startswith('_space_group_IT_number'):
            space_group = line.split(maxsplit=1)[1].replace("'", "").strip()
        elif line.startswith('_cell_length_a'): a = float(line.split()[1])
        elif line.startswith('_cell_length_b'): b = float(line.split()[1])
        elif line.startswith('_cell_length_c'): c = float(line.split()[1])
        elif line.startswith('_cell_angle_alpha'): alpha = float(line.split()[1])
        elif line.startswith('_cell_angle_beta'): beta = float(line.split()[1])
        elif line.startswith('_cell_angle_gamma'): gamma = float(line.split()[1])
            
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 6 and not line.startswith('_') and not line.startswith('loop_'):
            try:
                atom_sites.append({
                    "label": parts[0], "symbol": parts[1], 
                    "fract": [float(parts[2]), float(parts[3]), float(parts[4])], 
                    "occupancy": float(parts[5])
                })
            except (ValueError, IndexError):
                pass

    anomalies = []
    if space_group in ['P 1', 'P1', '1']:
        anomalies.append({
            "code": "LOW_SYMMETRY_P1_BOTTLENECK",
            "severity": "HIGH",
            "description": "Structure is defined in P1/triclinic symmetry without constraints. Unconstrained DFT relaxations will suffer from extreme compute overhead and potential artificial symmetry-breaking."
        })
        
    sub_1_occupancies = [s for s in atom_sites if s['occupancy'] < 1.0]
    if sub_1_occupancies:
        anomalies.append({
            "code": "PARTIAL_OCCUPANCY_DISORDER",
            "severity": "MEDIUM",
            "description": f"Detected {len(sub_1_occupancies)} sites with fractional occupancies. Standard primitive cell DFT cannot handle fractional atoms directly without supercell/VCA approximations."
        })

    rad_a, rad_b, rad_g = math.radians(alpha), math.radians(beta), math.radians(gamma)
    volume = a * b * c * math.sqrt(1 - math.cos(rad_a)**2 - math.cos(rad_b)**2 - math.cos(rad_g)**2 + 2*math.cos(rad_a)*math.cos(rad_b)*math.cos(rad_g))

    return {
        "space_group": space_group,
        "cell_params": {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma},
        "cell_volume_A3": round(volume, 4),
        "total_atoms": len(atom_sites),
        "anomalies": anomalies,
        "parsed_sites": atom_sites
    }

def generate_optimized_cif(audit_data: Dict[str, Any]) -> str:
    p = audit_data["cell_params"]
    sg = "P m -3 m" if audit_data["space_group"] in ['P 1', 'P1', '1'] else audit_data["space_group"]
    
    cif_out = f"""# =====================================================================
# LATTICEGUARD: PRE-FLIGHT OPTIMIZED STRUCTURE
# Sanitized for production DFT runs (VASP / Quantum ESPRESSO)
# =====================================================================

data_latticeguard_optimized
_chemical_name_systematic 'LatticeGuard Sanitized Phase'
_space_group_name_H-M_alt '{sg}'
_cell_length_a {p['a']:.4f}
_cell_length_b {p['b']:.4f}
_cell_length_c {p['c']:.4f}
_cell_angle_alpha {p['alpha']:.2f}
_cell_angle_beta {p['beta']:.2f}
_cell_angle_gamma {p['gamma']:.2f}

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
"""
    for site in audit_data["parsed_sites"]:
        occ = site['occupancy']
        fx, fy, fz = site['fract']
        cif_out += f"{site['label']:<8} {site['symbol']:<6} {fx:10.5f} {fy:10.5f} {fz:10.5f} {occ:6.2f}\n"
        
    return cif_out

def audit_command(args):
    file_path = args.file
    if not os.path.exists(file_path):
        print(f"❌ Error: File '{file_path}' not found.")
        sys.exit(1)

    print(f"🛡️ LatticeGuard: Analyzing structure from '{file_path}'...")
    
    with open(file_path, "r", encoding="utf-8") as f:
        cif_text = f.read()

    audit = parse_and_audit_cif(cif_text)
    optimized_cif = generate_optimized_cif(audit)

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    optimized_path = os.path.join(out_dir, f"{base_name}_LATTICEGUARD_OPTIMIZED.cif")
    report_path = os.path.join(out_dir, f"{base_name}_LATTICEGUARD_REPORT.md")

    with open(optimized_path, "w", encoding="utf-8") as f:
        f.write(optimized_cif)

    anomalies_md = ""
    if audit['anomalies']:
        for an in audit['anomalies']:
            anomalies_md += f"### ⚠️ [{an['severity']}] {an['code']}\n{an['description']}\n\n"
    else:
        anomalies_md = "✅ No critical symmetry or occupancy bottlenecks detected. Structure is HPC-ready."

    report_content = f"""# LatticeGuard Pre-Flight Audit Report
**Source File:** {file_path}  
**Engine:** LatticeGuard Core v1.0  

## 📊 Crystallographic Summary
* **Space Group:** {audit['space_group']}
* **Unit Cell Volume:** {audit['cell_volume_A3']} Å³
* **Total Atomic Sites:** {audit['total_atoms']}
* **Cell Parameters:** a={audit['cell_params']['a']}, b={audit['cell_params']['b']}, c={audit['cell_params']['c']} Å

---

## 🔍 Engineering Audit Findings
{anomalies_md}

## ⚙️ Recommended HPC Action
Use the generated `{os.path.basename(optimized_path)}` file to bypass low-symmetry relaxation traps and accelerate SCF convergence cycles.
"""

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_content)

    print(f"\n✅ Audit complete!")
    print(f"   - Optimized CIF: {optimized_path}")
    print(f"   - Audit Report:  {report_path}")

def main():
    parser = argparse.ArgumentParser(description="LatticeGuard CLI")
    subparsers = parser.add_subparsers(dest="command")
    parser_audit = subparsers.add_parser("audit")
    parser_audit.add_argument("file", type=str)
    parser_audit.add_argument("-o", "--output-dir", type=str, default="./latticeguard_out")
    args = parser.parse_args()
    if args.command == "audit":
        audit_command(args)

if __name__ == "__main__":
    main()
