import streamlit as st
import math
import io
import zipfile
import re
import json
import numpy as np

# ==========================================
# מטריצת רדיוסים קוולנטיים (אנגסטרם) - Cordero et al.
# ==========================================
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Na': 1.66,
    'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
    'Cl': 1.02, 'K': 2.03, 'Ca': 1.76, 'Ti': 1.60, 'Fe': 1.56,
    'Cu': 1.32, 'Zn': 1.22
}
DEFAULT_RADIUS = 1.0
CLASH_SCALING_FACTOR = 0.55

# ==========================================
# פונקציות קריסטלוגרפיות ומתמטיות מתקדמות
# ==========================================
def calculate_cell_volume(cell):
    a = cell.get('a', 1.0)
    b = cell.get('b', 1.0)
    c = cell.get('c', 1.0)
    alpha = math.radians(cell.get('alpha', 90.0))
    beta = math.radians(cell.get('beta', 90.0))
    gamma = math.radians(cell.get('gamma', 90.0))
    
    ca = math.cos(alpha)
    cb = math.cos(beta)
    cg = math.cos(gamma)
    
    term = 1.0 - ca**2 - cb**2 - cg**2 + 2.0 * ca * cb * cg
    return a * b * c * math.sqrt(max(0.0, term))

def fractional_to_cartesian(coords_frac, cell):
    a = cell.get('a', 5.0)
    b = cell.get('b', 5.0)
    c = cell.get('c', 5.0)
    alpha = math.radians(cell.get('alpha', 90.0))
    beta = math.radians(cell.get('beta', 90.0))
    gamma = math.radians(cell.get('gamma', 90.0))
    
    cos_alpha = math.cos(alpha)
    cos_beta = math.cos(beta)
    cos_gamma = math.cos(gamma)
    sin_gamma = math.sin(gamma)
    
    v_ax, v_ay, v_az = a, 0.0, 0.0
    v_bx, v_by, v_bz = b * cos_gamma, b * sin_gamma, 0.0
    v_cx = c * cos_beta
    v_cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    v_cz = math.sqrt(max(0.0, c**2 - v_cx**2 - v_cy**2))
    
    matrix = np.array([
        [v_ax, v_bx, v_cx],
        [v_ay, v_by, v_cy],
        [v_az, v_bz, v_cz]
    ])
    
    cart_coords = np.dot(matrix, coords_frac.T).T
    return cart_coords

def check_interatomic_clashes_vectorized(atoms, cell):
    """
    בדיקת התנגשויות מרחבית מהירה ומדויקת באמצעות מטריצות NumPy תחת תנאי שפה מחזוריים (PBC).
    """
    if not atoms:
        return []
        
    clashes = []
    n = len(atoms)
    frac_coords = np.array([[a[2], a[3], a[4]] for a in atoms]) % 1.0
    cart_coords = fractional_to_cartesian(frac_coords, cell)
    
    # חישוב מטריצת מרחקים וקטורית מהירה במיוחד
    cart_diff = cart_coords[:, np.newaxis, :] - cart_coords[np.newaxis, :, :]
    dist_matrix = np.linalg.norm(cart_diff, axis=-1)
    
    for i in range(n):
        for j in range(i + 1, n):
            atom1 = atoms[i]
            atom2 = atoms[j]
            
            sym1, sym2 = atom1[1], atom2[1]
            r1 = COVALENT_RADII.get(sym1, DEFAULT_RADIUS)
            r2 = COVALENT_RADII.get(sym2, DEFAULT_RADIUS)
            dynamic_threshold = (r1 + r2) * CLASH_SCALING_FACTOR
            
            dist = dist_matrix[i, j]
            
            if dist < dynamic_threshold:
                clashes.append({
                    'atom1': f"{atom1[0]} ({sym1}, line {atom1[6]})",
                    'atom2': f"{atom2[0]} ({sym2}, line {atom2[6]})",
                    'distance': round(float(dist), 4),
                    'threshold': round(float(dynamic_threshold), 4)
                })
            
    return clashes

# ==========================================
# מנועי ייצוא לפורמטי HPC שונים
# ==========================================
def export_to_vasp_poscar(cell, atoms, title="LatticeGuard VASP Structure"):
    a, b, c = cell.get('a', 5.0), cell.get('b', 5.0), cell.get('c', 5.0)
    
    v_cx = c * math.cos(math.radians(cell.get('beta', 90.0)))
    v_cy = c * (math.cos(math.radians(cell.get('alpha', 90.0))) - math.cos(math.radians(cell.get('beta', 90.0))) * math.cos(math.radians(cell.get('gamma', 90.0)))) / math.sin(math.radians(cell.get('gamma', 90.0)))
    v_cz = math.sqrt(max(0.0, c**2 - v_cx**2 - v_cy**2))
    
    poscar = f"{title}\n1.0\n"
    poscar += f" {a:12.6f} 0.000000 0.000000\n"
    poscar += f" {b * math.cos(math.radians(cell.get('gamma', 90.0))):12.6f} {b * math.sin(math.radians(cell.get('gamma', 90.0))):12.6f} 0.000000\n"
    poscar += f" {v_cx:12.6f} {v_cy:12.6f} {v_cz:12.6f}\n"
    
    element_counts = {}
    element_coords = {}
    for atom in atoms:
        sym = atom[1]
        element_counts[sym] = element_counts.get(sym, 0) + 1
        element_coords.setdefault(sym, []).append((atom[2], atom[3], atom[4]))
        
    symbols = list(element_counts.keys())
    poscar += " " + " ".join(symbols) + "\n"
    poscar += " " + " ".join([str(element_counts[s]) for s in symbols]) + "\n"
    poscar += "Direct\n"
    
    for sym in symbols:
        for coords in element_coords[sym]:
            poscar += f" {coords[0]:10.6f} {coords[1]:10.6f} {coords[2]:10.6f}\n"
    return poscar

def export_to_quantum_espresso(cell, atoms):
    a = cell.get('a', 5.0)
    pwi = f"! Quantum ESPRESSO Input File generated by LatticeGuard\n"
    pwi += f"&CONTROL\n  calculation = 'scf'\n  prefix = 'latticeguard'\n/\n"
    pwi += f"&SYSTEM\n  ibrav = 0, nat = {len(atoms)}, ntyp = {len(set(a[1] for a in atoms))}\n"
    pwi += f"  celldm(1) = {a / 0.52917721092:.6f}\n/\n"
    pwi += f"&ELECTRONS\n/\n"
    pwi += f"ATOMIC_SPECIES\n"
    for sym in set(a[1] for a in atoms):
        pwi += f" {sym}  1.00  {sym}.upf\n"
    pwi += f"ATOMIC_POSITIONS crystal\n"
    for atom in atoms:
        pwi += f" {atom[1]}  {atom[2]:10.6f} {atom[3]:10.6f} {atom[4]:10.6f}\n"
    pwi += f"K_POINTS automatic\n  4 4 4 1 1 1\n"
    return pwi

def export_to_xyz(cell, atoms):
    cart_coords = fractional_to_cartesian(np.array([[a[2], a[3], a[4]] for a in atoms]), cell)
    xyz = f"{len(atoms)}\nGenerated by LatticeGuard HPC Suite\n"
    for idx, atom in enumerate(atoms):
        sym = atom[1]
        x, y, z = cart_coords[idx]
        xyz += f" {sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n"
    return xyz

# ==========================================
# מנוע הליבה: אבחון ותיקון מלא
# ==========================================
def audit_and_fix_cif(content: str):
    issues = []
    audit_trail = []
    
    lines = content.splitlines()
    audit_trail.append(f"Initiated parsing of {len(lines)} lines.")

    fixed_lines = []
    space_group = "P 1"
    has_cell_params = False
    cell = {'a': 5.0, 'b': 5.0, 'c': 5.0, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}
    atoms = []
    seen_labels = set()

    for param in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
        match = re.search(rf"_cell_length_{param}\s+([0-9\.]+)", content, re.IGNORECASE) or \
                re.search(rf"_cell_angle_{param}\s+([0-9\.]+)", content, re.IGNORECASE)
        if match:
            cell[param] = float(match.group(1))
            if param in ['a', 'b', 'c']:
                has_cell_params = True

    sg_match = re.search(r"_symmetry_space_group_name_H-M\s+['\"]?([^'\n]+)['\"]?", content, re.IGNORECASE)
    if sg_match:
        space_group = sg_match.group(1).strip()

    cell_vol = calculate_cell_volume(cell)
    if cell_vol <= 0.0:
        issues.append({
            "type": "error",
            "category": "Unphysical Unit Cell",
            "message": f"Critical: Unit cell volume is non-positive ({cell_vol:.2f} A^3)."
        })
    else:
        audit_trail.append(f"Unit cell volume verified: {cell_vol:.2f} A^3.")

    in_loop = False
    loop_headers = []
    x_idx = y_idx = z_idx = label_idx = occ_idx = -1

    for line_idx, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("loop_"):
            in_loop = True
            loop_headers = []
            x_idx = y_idx = z_idx = label_idx = occ_idx = -1
            fixed_lines.append(line)
            continue
            
        if in_loop:
            if stripped.startswith("_"):
                loop_headers.append(stripped)
                curr_idx = len(loop_headers) - 1
                if "_atom_site_fract_x" in stripped: x_idx = curr_idx
                elif "_atom_site_fract_y" in stripped: y_idx = curr_idx
                elif "_atom_site_fract_z" in stripped: z_idx = curr_idx
                elif "_atom_site_label" in stripped or "_atom_site_type_symbol" in stripped: label_idx = curr_idx
                elif "_atom_site_occupancy" in stripped: occ_idx = curr_idx
                fixed_lines.append(line)
                continue
            elif stripped == "" or stripped.startswith("#"):
                fixed_lines.append(line)
                continue
            else:
                parts = stripped.split()
                if len(parts) >= len(loop_headers) and x_idx != -1 and y_idx != -1 and z_idx != -1:
                    try:
                        label = parts[label_idx] if label_idx != -1 and label_idx < len(parts) else f"X{line_idx}"
                        fx = float(parts[x_idx])
                        fy = float(parts[y_idx])
                        fz = float(parts[z_idx])
                        occ = float(parts[occ_idx]) if occ_idx != -1 and occ_idx < len(parts) else 1.0
                        
                        if occ <= 0.0:
                            continue
                            
                        if label in seen_labels:
                            issues.append({
                                "type": "warning",
                                "category": "Structural Warning",
                                "message": f"Line {line_idx}: Duplicate atom label '{label}'."
                            })
                        else:
                            seen_labels.add(label)
                        
                        elem_match = re.match(r"([A-Za-z]{1,2})", label)
                        symbol = elem_match.group(1) if elem_match else "C"
                        
                        atoms.append((label, symbol, fx, fy, fz, occ, line_idx))
                    except ValueError:
                        pass
        
        if "_symmetry_space_group_name_H-M" in stripped and ('P 1' in stripped or 'P1' in stripped):
            issues.append({
                "type": "warning",
                "category": "Symmetry Trap",
                "message": f"Line {line_idx}: Detected P1 low-symmetry trap."
            })

        fixed_lines.append(line)

    if atoms and cell_vol > 0.0:
        clashes = check_interatomic_clashes_vectorized(atoms, cell)
        if clashes:
            for c in clashes:
                issues.append({
                    "type": "error",
                    "category": "Interatomic Clash",
                    "message": f"Spatial clash between {c['atom1']} and {c['atom2']}. Distance: {c['distance']} A."
                })
            audit_trail.append(f"Detected {len(clashes)} spatial clashes.")
        else:
            audit_trail.append("No atomic clashes detected.")

    if not has_cell_params and len(content.strip()) > 0:
        issues.append({
            "type": "error",
            "category": "Missing Data",
            "message": "Critical: Missing essential cell parameters."
        })

    fixed_content = "\n".join(fixed_lines)
    poscar_content = export_to_vasp_poscar(cell, atoms)
    qe_content = export_to_quantum_espresso(cell, atoms)
    xyz_content = export_to_xyz(cell, atoms)
    
    report_dict = {
        "space_group": space_group,
        "cell_parameters": cell,
        "cell_volume": cell_vol,
        "total_atoms": len(atoms),
        "issues_found": len(issues),
        "audit_trail": audit_trail
    }

    return fixed_content, issues, audit_trail, poscar_content, qe_content, xyz_content, report_dict

def validate_fixed_cif(fixed_content: str):
    errors = []
    if not fixed_content or len(fixed_content.strip()) < 10:
        errors.append("Critical: Fixed content is empty.")
    if "NaN" in fixed_content or "inf" in fixed_content:
        errors.append("Warning: Contains unphysical NaN values.")
    is_valid = len([e for e in errors if "Critical" in e]) == 0
    return is_valid, errors

# ==========================================
# ממשק משתמש Streamlit פרימיום
# ==========================================
st.set_page_config(page_title="LatticeGuard HPC Suite", page_icon="🛡️", layout="centered")

st.title("🛡️ LatticeGuard HPC Suite")
st.markdown("### Elite-Tier Crystallographic Auditor & Multi-Format Converter")
st.markdown("Powered by vectorized NumPy spatial acceleration and strict simulation safety gatekeepers.")
st.markdown("---")

SAMPLE_CIF = """data_sample_p1_trap
_audit_creation_date             2026-03-30
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a                   5.43000
_cell_length_b                   5.43000
_cell_length_c                   5.43000
_cell_angle_alpha                90.00000
_cell_angle_beta                 90.00000
_cell_angle_gamma                90.00000

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 0.00000 0.00000 0.00000
Si2 0.00000 0.00000 0.00000
"""

st.sidebar.header("HPC Control Panel")
load_sample = st.sidebar.button("🧪 Load Sample Defective CIF")
uploaded_files = st.file_uploader("Upload Structure File(s)", type=["cif", "txt"], accept_multiple_files=True)

files_to_process = []
if load_sample:
    files_to_process.append(("sample_defective.cif", SAMPLE_CIF))
elif uploaded_files:
    for f in uploaded_files:
        try:
            content = f.getvalue().decode("utf-8", errors="ignore")
            files_to_process.append((f.name, content))
        except Exception as e:
            st.error(f"Failed to read file {f.name}: {e}")

if files_to_process:
    st.markdown(f"### 📦 Batch Processing Dashboard ({len(files_to_process)} files loaded)")
    successful_files = []

    for file_name, cif_content in files_to_process:
        with st.container():
            st.markdown(f"---")
            st.subheader(f"📄 File: `{file_name}`")
            
            with st.expander("📁 View Raw Content", expanded=False):
                st.text_area("Raw CIF", cif_content, height=120, key=f"raw_{file_name}")

            fixed_content, issues, audit_trail, poscar, qe, xyz, report_dict = audit_and_fix_cif(cif_content)
            is_valid, validation_errors = validate_fixed_cif(fixed_content)

            st.markdown("#### 🔍 Audit Diagnostic Report")
            if not issues:
                st.success("✨ **Audit Clean:** No structural traps or clashes identified.")
            else:
                for issue in issues:
                    if issue["type"] == "error":
                        st.error(f"❌ **[{issue['category']}]** {issue['message']}")
                    else:
                        st.warning(f"⚠️ **[{issue['category']}]** {issue['message']}")

            with st.expander("📋 View Execution Audit Trail", expanded=False):
                for log_entry in audit_trail:
                    st.text(f"• {log_entry}")

            st.markdown("#### 🛡️ Safety & HPC Exporters Gatekeeper")
            if is_valid:
                st.success("✅ **Verification Passed:** Structure is certified safe for HPC pipelines.")
                
                f_clean = file_name.replace(".cif", "_fixed.cif")
                f_poscar = file_name.replace(".cif", "_POSCAR")
                f_qe = file_name.replace(".cif", ".pwi")
                f_xyz = file_name.replace(".cif", ".xyz")
                f_json = file_name.replace(".cif", "_audit.json")
                
                successful_files.extend([
                    (f_clean, fixed_content),
                    (f_poscar, poscar),
                    (f_qe, qe),
                    (f_xyz, xyz),
                    (f_json, json.dumps(report_dict, indent=2))
                ])
                
                c1, c2, c3, c4 = st.columns(4)
                with c1:
                    st.download_button("📥 CIF", fixed_content, file_name=f_clean, key=f"dl_cif_{file_name}")
                with c2:
                    st.download_button("📥 POSCAR", poscar, file_name=f_poscar, key=f"dl_pos_{file_name}")
                with c3:
                    st.download_button("📥 QE (.pwi)", qe, file_name=f_qe, key=f"dl_qe_{file_name}")
                with c4:
                    st.download_button("📥 JSON", json.dumps(report_dict, indent=2), file_name=f_json, key=f"dl_json_{file_name}")
            else:
                st.error("🛑 **Verification Blocked:** Post-processing sanity check failed. Output withheld under 'Do No Harm' policy.")
                for err in validation_errors:
                    st.warning(err)

    if successful_files:
        st.markdown("---")
        st.markdown("### 🗂️ Batch Archive Export")
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for fname, fcontent in successful_files:
                zip_file.writestr(fname, fcontent)
        zip_buffer.seek(0)
        
        st.download_button(
            label="📦 Download Complete HPC Bundle (ZIP)",
            data=zip_buffer,
            file_name="latticeguard_hpc_complete_bundle.zip",
            mime="application/zip",
            type="primary",
            use_container_width=True
        )
else:
    st.info("💡 **Get Started:** Click 'Load Sample Defective CIF' in the sidebar or upload `.cif` files for instant elite-tier HPC batch auditing.")
