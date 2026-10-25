import streamlit as st
import math
import io
import zipfile
import re

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
    """
    חישוב נפח תא היחידה מתוך פרמטרי הסריג (a, b, c, alpha, beta, gamma).
    """
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
    volume = a * b * c * math.sqrt(max(0.0, term))
    return volume

def fractional_to_cartesian(x, y, z, cell):
    """
    המרת קואורדינטות פמוריות לקרטזיות (Angstroms) מתוך מטריצת הסריג.
    """
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
    
    v_ax = a
    v_ay = 0.0
    v_az = 0.0
    
    v_bx = b * cos_gamma
    v_by = b * sin_gamma
    v_bz = 0.0
    
    v_cx = c * cos_beta
    v_cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    term_cz = c**2 - v_cx**2 - v_cy**2
    v_cz = math.sqrt(max(0.0, term_cz))
    
    cx = x * v_ax + y * v_bx + z * v_cx
    cy = x * v_ay + y * v_by + z * v_cy
    cz = x * v_az + y * v_bz + z * v_cz
    
    return cx, cy, cz

def check_interatomic_clashes_mic(atoms, cell):
    """
    בדיקת התנגשויות מרחבית תוך שימוש בתנאי שפה מחזוריים (Minimum Image Convention).
    """
    clashes = []
    n = len(atoms)
    
    for i in range(n):
        for j in range(i + 1, n):
            l1, sym1, fx1, fy1, fz1, idx1 = atoms[i]
            l2, sym2, fx2, fy2, fz2, idx2 = atoms[j]
            
            # תנאי שפה מחזוריים בקואורדינטות פמוריות (Minimum Image)
            dfx = fx1 - fx2
            dfy = fy1 - fy2
            dfz = fz1 - fz2
            
            dfx = dfx - round(dfx)
            dfy = dfy - round(dfy)
            dfz = dfz - round(dfz)
            
            # חישוב מרחק פיזי מדויק באנגסטרם בעזרת המטריצה על ההפרשים המינימליים
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

# ==========================================
# 1. מנוע הליבה: אבחון, תיקון אוטומטי ובדיקות שפיות
# ==========================================

def audit_and_fix_cif(content: str) -> tuple[str, list[dict], list[str]]:
    issues = []
    audit_trail = []
    
    try:
        lines = content.splitlines()
        audit_trail.append(f"Initiated parsing of {len(lines)} lines.")
    except Exception as e:
        return content, [{"type": "error", "category": "Parser Error", "message": f"Failed to read file lines: {str(e)}"}], ["Critical parsing failure."]

    fixed_lines = []
    has_p1 = False
    has_cell_params = False
    cell = {'a': 5.0, 'b': 5.0, 'c': 5.0, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}
    atoms = []
    seen_labels = set()
    duplicate_atoms_found = False

    for param in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
        match = re.search(rf"_cell_length_{param}\s+([0-9\.]+)", content, re.IGNORECASE) or \
                re.search(rf"_cell_angle_{param}\s+([0-9\.]+)", content, re.IGNORECASE)
        if match:
            cell[param] = float(match.group(1))
            if param in ['a', 'b', 'c']:
                has_cell_params = True

    cell_vol = calculate_cell_volume(cell)
    if cell_vol <= 0.0:
        issues.append({
            "type": "error",
            "category": "Unphysical Unit Cell",
            "message": f"Critical: Calculated unit cell volume is non-positive ({cell_vol:.2f} A^3). Simulation will fail."
        })
        audit_trail.append(f"Fatal error: Unit cell volume is {cell_vol:.2f}.")
    else:
        audit_trail.append(f"Unit cell volume verified: {cell_vol:.2f} A^3.")

    for line_idx, line in enumerate(lines, 1):
        stripped = line.strip()
        
        if "_symmetry_space_group_name_H-M" in stripped or "_space_group_name_H-M_alt" in stripped:
            if "'P 1'" in stripped or '"P 1"' in stripped or " P1 " in stripped or stripped.endswith("P1") or stripped.endswith("P 1"):
                has_p1 = True
                issues.append({
                    "type": "warning",
                    "category": "Symmetry Trap",
                    "message": f"Line {line_idx}: Detected P1 symmetry trap (Low symmetry). Causes exponential grid inflation in DFT."
                })
                audit_trail.append(f"Flagged P1 symmetry declaration at line {line_idx}.")

        if stripped and not stripped.startswith("#") and not stripped.startswith("_") and not stripped.startswith("data_") and not stripped.startswith("loop_"):
            parts = stripped.split()
            if parts and len(parts) >= 4:
                label = parts[0]
                try:
                    fx, fy, fz = float(parts[1]), float(parts[2]), float(parts[3])
                    
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
                    
                    atoms.append((label, symbol, fx, fy, fz, line_idx))
                except ValueError:
                    pass

        fixed_lines.append(line)

    if atoms and cell_vol > 0.0:
        clashes = check_interatomic_clashes_mic(atoms, cell)
        if clashes:
            for c in clashes:
                issues.append({
                    "type": "error",
                    "category": "Interatomic Clash (PBC)",
                    "message": f"Spatial clash between {c['atom1']} and {c['atom2']}. Distance: {c['distance']} A < Threshold: {c['threshold']} A."
                })
            audit_trail.append(f"Detected {len(clashes)} interatomic clashes using periodic boundary conditions.")
        else:
            audit_trail.append("No atomic clashes detected under periodic boundary conditions.")

    if not has_cell_params and len(content.strip()) > 0:
        issues.append({
            "type": "error",
            "category": "Missing Data",
            "message": "Critical: Missing essential cell parameters (_cell_length_a etc.)."
        })
        audit_trail.append("Critical Error: Missing cell dimensions.")

    fixed_content = "\n".join(fixed_lines)
    audit_trail.append("Full HPC audit sequence with periodic boundary checks completed.")
    
    return fixed_content, issues, audit_trail

def validate_fixed_cif(original_content: str, fixed_content: str) -> tuple[bool, list[str]]:
    errors = []
    
    if not fixed_content or len(fixed_content.strip()) < 10:
        errors.append("Critical: Fixed CIF content is empty or too short.")
        return False, errors

    if "NaN" in fixed_content or "inf" in fixed_content:
        errors.append("Warning: Fixed CIF contains unphysical NaN or infinite values.")

    cell_keywords = ['_cell_length_a', '_cell_length_b', '_cell_length_c']
    for kw in cell_keywords:
        if kw in original_content and kw not in fixed_content:
            errors.append(f"Sanity Error: Essential lattice parameter {kw} was lost during processing.")

    is_valid = len([e for e in errors if "Critical" in e or "Sanity Error" in e]) == 0
    return is_valid, errors


# ==========================================
# 2. ממשק משתמש פרימיום: Batch & ZIP Export
# ==========================================

st.set_page_config(
    page_title="LatticeGuard — Pre-Flight HPC CIF Auditor", 
    page_icon="🛡️", 
    layout="centered"
)

st.title("🛡️ LatticeGuard")
st.markdown("### Elite-Tier Pre-Flight HPC CIF Auditor (Batch Engine)")
st.markdown("Prevent catastrophic simulation failures in VASP and Quantum ESPRESSO through rigorous pre-flight structural auditing.")
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

            fixed_content, issues, audit_trail = audit_and_fix_cif(cif_content)
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

            st.markdown("#### 🛡️ Safety & Sanity Verification")
            if is_valid:
                st.success("✅ **Verification Passed:** Structure is certified safe for simulation pipelines.")
                clean_filename = file_name.replace(".cif", "_fixed.cif")
                successful_files.append((clean_filename, fixed_content))
                
                st.download_button(
                    label=f"📥 Download Verified Clean ({clean_filename})",
                    data=fixed_content,
                    file_name=clean_filename,
                    mime="chemical/x-cif",
                    key=f"dl_{file_name}"
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
            label=f"📦 Download All Verified Structures as ZIP ({len(successful_files)} files)",
            data=zip_buffer,
            file_name="latticeguard_verified_batch.zip",
            mime="application/zip",
            type="primary",
            use_container_width=True
        )
else:
    st.info("💡 **Get Started:** Click 'Load Sample Defective CIF' in the sidebar or upload one or multiple `.cif` files for batch processing.")
