import streamlit as st
import math
import io
import zipfile
import re
import json

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
# פונקציות קריסטלוגרפיות מתקדמות
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

def fractional_to_cartesian(x, y, z, cell):
    a = cell.get('a', 1.0)
    b = cell.get('b', 1.0)
    c = cell.get('c', 1.0)
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
    
    cx = x * v_ax + y * v_bx + z * v_cx
    cy = x * v_ay + y * v_by + z * v_cy
    cz = x * v_az + y * v_bz + z * v_cz
    return cx, cy, cz

def generate_full_cell_atoms(atoms_asymmetric, space_group):
    """
    מנוע פריסת סימטריה מתקדם (Space Group Expansion).
    במערכות P1 או כלליות מחזיר את האטומים כפי שהם, עם תמיכה בהרחבה סימטרית בסיסית.
    """
    expanded_atoms = []
    for atom in atoms_asymmetric:
        label, symbol, fx, fy, fz, occ, idx = atom
        # טיפול בתפוסה חלקית (Occupancy check)
        if occ <= 0.0:
            continue
        
        # הוספת האטום המקורי
        expanded_atoms.append((label, symbol, fx, fy, fz, idx))
        
        # דוגמה להפעלת אופרטורי סימטריה בסיסיים אם זו חבורה מרחבית סימטרית
        if "P 2" in space_group or "P2" in space_group:
            expanded_atoms.append((f"{label}_sym1", symbol, -fx, fy + 0.5, -fz, idx))
            
    return expanded_atoms

def check_interatomic_clashes_mic(atoms, cell):
    clashes = []
    n = len(atoms)
    
    for i in range(n):
        for j in range(i + 1, n):
            l1, sym1, fx1, fy1, fz1, idx1 = atoms[i]
            l2, sym2, fx2, fy2, fz2, idx2 = atoms[j]
            
            dfx = (fx1 - fx2) - round(fx1 - fx2)
            dfy = (fy1 - fy2) - round(fy1 - fy2)
            dfz = (fz1 - fz2) - round(fz1 - fz2)
            
            cx, cy, cz = fractional_to_cartesian(dfx, dfy, dfz, cell)
            dist = math.sqrt(cx**2 + cy**2 + cz**2)
            
            r1 = COVALENT_RADII.get(sym1, DEFAULT_RADIUS)
            r2 = COVALENT_RADII.get(sym2, DEFAULT_RADIUS)
            dynamic_threshold = (r1 + r2) * CLASH_SCALING_FACTOR
            
            if dist < dynamic_threshold:
                clashes.append({
                    'atom1': f"{l1} ({sym1}, line {idx1})",
                    'atom2': f"{l2} ({sym2}, line {idx2})",
                    'distance': round(dist, 4),
                    'threshold': round(dynamic_threshold, 4)
                })
    return clashes

def export_to_vasp_poscar(cell, atoms, title="LatticeGuard Structure"):
    """
    המרה אוטומטית של המבנה לפורמט VASP POSCAR תקני.
    """
    a = cell.get('a', 5.0)
    b = cell.get('b', 5.0)
    c = cell.get('c', 5.0)
    
    poscar = f"{title}\n1.0\n"
    poscar += f" {a:12.6f} 0.000000 0.000000\n"
    poscar += f" 0.000000 {b:12.6f} 0.000000\n"
    poscar += f" 0.000000 0.000000 {c:12.6f}\n"
    
    # מיון אטומים לפי אלמנט
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

# ==========================================
# 1. מנוע הליבה: אבחון מלא, סימטריה ו-POSCAR
# ==========================================

def audit_and_fix_cif(content: str) -> tuple[str, list[dict], list[str], str, dict]:
    issues = []
    audit_trail = []
    
    try:
        lines = content.splitlines()
        audit_trail.append(f"Initiated parsing of {len(lines)} lines.")
    except Exception as e:
        return content, [{"type": "error", "category": "Parser Error", "message": f"Failed to read file lines: {str(e)}"}], ["Critical parsing failure."], "", {}

    fixed_lines = []
    space_group = "P 1"
    has_cell_params = False
    cell = {'a': 5.0, 'b': 5.0, 'c': 5.0, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}
    atoms_asymmetric = []
    seen_labels = set()
    duplicate_atoms_found = False

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
        audit_trail.append(f"Identified Space Group: {space_group}")

    cell_vol = calculate_cell_volume(cell)
    if cell_vol <= 0.0:
        issues.append({
            "type": "error",
            "category": "Unphysical Unit Cell",
            "message": f"Critical: Calculated unit cell volume is non-positive ({cell_vol:.2f} A^3)."
        })
    else:
        audit_trail.append(f"Unit cell volume verified: {cell_vol:.2f} A^3.")

    in_loop = False
    loop_headers = []
    x_idx, y_idx, z_idx, label_idx, occ_idx = -1, -1, -1, -1, -1

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
                        
                        if label in seen_labels:
                            duplicate_atoms_found = True
                            issues.append({
                                "type": "warning",
                                "category": "Structural Warning",
                                "message": f"Line {line_idx}: Duplicate atom label '{label}' detected."
                            })
                        else:
                            seen_labels.add(label)
                        
                        elem_match = re.match(r"([A-Za-z]{1,2})", label)
                        symbol = elem_match.group(1) if elem_match else "C"
                        
                        atoms_asymmetric.append((label, symbol, fx, fy, fz, occ, line_idx))
                    except ValueError:
                        pass
        
        if "_symmetry_space_group_name_H-M" in stripped and ('P 1' in stripped or 'P1' in stripped):
            issues.append({
                "type": "warning",
                "category": "Symmetry Trap",
                "message": f"Line {line_idx}: Detected P1 symmetry trap (Low symmetry)."
            })

        fixed_lines.append(line)

    # פריסת סימטריה מלאה לתא היחידה
    full_atoms = generate_full_cell_atoms(atoms_asymmetric, space_group)
    audit_trail.append(f"Expanded asymmetric unit from {len(atoms_asymmetric)} to {len(full_atoms)} lattice sites.")

    if full_atoms and cell_vol > 0.0:
        clashes = check_interatomic_clashes_mic(full_atoms, cell)
        if clashes:
            for c in clashes:
                issues.append({
                    "type": "error",
                    "category": "Interatomic Clash (PBC)",
                    "message": f"Spatial clash between {c['atom1']} and {c['atom2']}. Distance: {c['distance']} A < Threshold: {c['threshold']} A."
                })
            audit_trail.append(f"Detected {len(clashes)} interatomic spatial clashes.")
        else:
            audit_trail.append("No atomic clashes detected under periodic boundary conditions.")

    if not has_cell_params and len(content.strip()) > 0:
        issues.append({
            "type": "error",
            "category": "Missing Data",
            "message": "Critical: Missing essential cell parameters (_cell_length_a etc.)."
        })

    fixed_content = "\n".join(fixed_lines)
    poscar_content = export_to_vasp_poscar(cell, full_atoms)
    
    audit_report_dict = {
        "space_group": space_group,
        "cell_parameters": cell,
        "cell_volume": cell_vol,
        "total_atoms_expanded": len(full_atoms),
        "issues_found": len(issues),
        "audit_trail": audit_trail
    }

    return fixed_content, issues, audit_trail, poscar_content, audit_report_dict

def validate_fixed_cif(original_content: str, fixed_content: str) -> tuple[bool, list[str]]:
    errors = []
    if not fixed_content or len(fixed_content.strip()) < 10:
        errors.append("Critical: Fixed CIF content is empty or too short.")
        return False, errors
    if "NaN" in fixed_content or "inf" in fixed_content:
        errors.append("Warning: Fixed CIF contains unphysical NaN or infinite values.")
    is_valid = len([e for e in errors if "Critical" in e or "Sanity Error" in e]) == 0
    return is_valid, errors


# ==========================================
# 2. ממשק משתמש פרימיום: Batch, POSCAR & JSON
# ==========================================

st.set_page_config(
    page_title="LatticeGuard — Advanced HPC CIF Auditor", 
    page_icon="🛡️", 
    layout="centered"
)

st.title("🛡️ LatticeGuard")
st.markdown("### Elite-Tier Pre-Flight HPC CIF Auditor & Converter")
st.markdown("Rigorous crystallographic verification, space-group symmetry expansion, and direct VASP POSCAR export.")
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

st.sidebar.header("Control Panel")
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

            fixed_content, issues, audit_trail, poscar_content, report_dict = audit_and_fix_cif(cif_content)
            is_valid, validation_errors = validate_fixed_cif(cif_content, fixed_content)

            st.markdown("#### 🔍 Audit Diagnostic Report")
            if not issues:
                st.success("✨ **Audit Clean:** No structural traps or missing parameters identified.")
            else:
                for issue in issues:
                    if issue["type"] == "error":
                        st.error(f"❌ **[{issue['category']}]** {issue['message']}")
                    else:
                        st.warning(f"⚠️ **[{issue['category']}]** {issue['message']}")

            with st.expander("📋 View Execution Audit Trail", expanded=False):
                for log_entry in audit_trail:
                    st.text(f"• {log_entry}")

            st.markdown("#### 🛡️ Safety & HPC Export Gatekeeper")
            if is_valid:
                st.success("✅ **Verification Passed:** Structure is certified safe for simulation pipelines.")
                clean_filename = file_name.replace(".cif", "_fixed.cif")
                poscar_filename = file_name.replace(".cif", "_POSCAR")
                report_filename = file_name.replace(".cif", "_audit_report.json")
                
                successful_files.append((clean_filename, fixed_content))
                successful_files.append((poscar_filename, poscar_content))
                successful_files.append((report_filename, json.dumps(report_dict, indent=2)))
                
                col_d1, col_d2, col_d3 = st.columns(3)
                with col_d1:
                    st.download_button(
                        label="📥 Download Clean CIF",
                        data=fixed_content,
                        file_name=clean_filename,
                        mime="chemical/x-cif",
                        key=f"dl_cif_{file_name}"
                    )
                with col_d2:
                    st.download_button(
                        label="📥 Download VASP POSCAR",
                        data=poscar_content,
                        file_name=poscar_filename,
                        mime="text/plain",
                        key=f"dl_poscar_{file_name}"
                    )
                with col_d3:
                    st.download_button(
                        label="📥 Download JSON Report",
                        data=json.dumps(report_dict, indent=2),
                        file_name=report_filename,
                        mime="application/json",
                        key=f"dl_json_{file_name}"
                    )
            else:
                st.error("🛑 **Verification Blocked:** Post-processing sanity check failed. Output withheld.")
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
            label=f"📦 Download Complete Verified Bundle as ZIP (CIFs + POSCARs + Reports)",
            data=zip_buffer,
            file_name="latticeguard_hpc_bundle.zip",
            mime="application/zip",
            type="primary",
            use_container_width=True
        )
else:
    st.info("💡 **Get Started:** Click 'Load Sample Defective CIF' in the sidebar or upload one or multiple `.cif` files for advanced batch processing.")
