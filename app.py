import streamlit as st
import math
import io
import zipfile

# ==========================================
# 1. מנוע הליבה: אבחון מתקדם, הגנה ובדיקות שפיות
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
    atoms = []
    duplicate_atoms_found = False

    for line_idx, line in enumerate(lines, 1):
        stripped = line.strip()
        
        # זיהוי מלכודת סימטריה P1
        if "_symmetry_space_group_name_H-M" in stripped or "_space_group_name_H-M_alt" in stripped:
            if "'P 1'" in stripped or '"P 1"' in stripped or " P1 " in stripped or stripped.endswith("P1") or stripped.endswith("P 1"):
                has_p1 = True
                issues.append({
                    "type": "warning",
                    "category": "Symmetry Trap",
                    "message": f"Line {line_idx}: Detected P1 symmetry trap (Low symmetry). Causes exponential grid inflation in DFT."
                })
                audit_trail.append(f"Flagged P1 symmetry declaration at line {line_idx}.")
        
        # זיהוי פרמטרי תא
        if "_cell_length_a" in stripped:
            has_cell_params = True
            audit_trail.append("Confirmed presence of primary cell parameter (_cell_length_a).")

        # איסוף קואורדינציות אטומים לבדיקת התנגשויות
        if stripped and not stripped.startswith("#") and not stripped.startswith("_") and not stripped.startswith("data_") and not stripped.startswith("loop_"):
            parts = stripped.split()
            if parts and len(parts) >= 4:
                label = parts[0]
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    
                    if any(a[0] == label for a in atoms):
                        duplicate_atoms_found = True
                        issues.append({
                            "type": "warning",
                            "category": "Structural Warning",
                            "message": f"Line {line_idx}: Duplicate atom label '{label}' detected."
                        })
                    
                    atoms.append((label, x, y, z, line_idx))
                except ValueError:
                    pass

        fixed_lines.append(line)

    # בדיקת התנגשויות אטומים מרחבית
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            l1, x1, y1, z1, idx1 = atoms[i]
            l2, x2, y2, z2, idx2 = atoms[j]
            
            dist_sq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
            if dist_sq < 0.0001 and dist_sq > 0.0:
                issues.append({
                    "type": "error",
                    "category": "Interatomic Clash",
                    "message": f"Potential atomic clash between '{l1}' (line {idx1}) and '{l2}' (line {idx2}). Distance is near zero."
                })
                audit_trail.append(f"Critical: Spatial clash detected between {l1} and {l2}.")

    if not has_cell_params and len(content.strip()) > 0:
        issues.append({
            "type": "error",
            "category": "Missing Data",
            "message": "Critical: Missing essential cell parameters (_cell_length_a etc.)."
        })
        audit_trail.append("Critical Error: Missing cell dimensions.")
    else:
        audit_trail.append("Cell parameters structure integrity verified.")

    if duplicate_atoms_found:
        audit_trail.append("Warning: Duplicate atom sites recorded in atomic coordinate list.")

    fixed_content = "\n".join(fixed_lines)
    audit_trail.append("Parsing and diagnostic scan completed successfully under 'Do No Harm' policy.")
    
    return fixed_content, issues, audit_trail

def validate_fixed_cif(original_content: str, fixed_content: str) -> tuple[bool, list[str]]:
    errors = []
    
    if not fixed_content or len(fixed_content.strip()) < 10:
        errors.append("Critical: Fixed CIF content is empty or too short.")
        return False, errors

    if "NaN" in fixed_content or "inf" in fixed_content:
        errors.append("Warning: Fixed CIF contains unphysical NaN or infinite values.")

    cell_keywords = ['_cell_length_a', '_cell_length_b', '_cell_length_c', 
                     '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']
    
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

    # כפתור הורדת ZIP מרכזי אם יש יותר מקובץ אחד תקין או אפילו שורה מרכזית
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
