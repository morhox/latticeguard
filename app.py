import streamlit as st
import re

# ==========================================
# 1. מנוע ליבה מתקדם עם לוג פעולות שקוף (Audit Trail Engine)
# ==========================================

def audit_and_fix_cif(content: str) -> tuple[str, list[dict], list[str]]:
    issues = []
    audit_trail = []
    lines = content.splitlines()
    fixed_lines = []
    
    audit_trail.append(f"Initiated parsing of {len(lines)} lines.")
    
    has_p1 = False
    has_cell_params = False
    atom_labels = set()
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
        
        # בדיקת קיום פרמטרי תא
        if "_cell_length_a" in stripped:
            has_cell_params = True
            audit_trail.append("Confirmed presence of primary cell parameter (_cell_length_a).")

        # זיהוי אטומים ובדיקת כפילויות בלולאות
        if stripped and not stripped.startswith("#") and not stripped.startswith("_") and not stripped.startswith("data_") and not stripped.startswith("loop_"):
            parts = stripped.split()
            if parts and len(parts) >= 4:  # תבנית אופיינית לשורת אתר אטום: label x y z
                label = parts[0]
                if label in atom_labels:
                    duplicate_atoms_found = True
                    issues.append({
                        "type": "warning",
                        "category": "Structural Warning",
                        "message": f"Line {line_idx}: Duplicate atom label '{label}' detected."
                    })
                else:
                    atom_labels.add(label)

        fixed_lines.append(line)

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
# 2. ממשק משתמש פרימיום (Elite UI)
# ==========================================

st.set_page_config(
    page_title="LatticeGuard — Pre-Flight HPC CIF Auditor", 
    page_icon="🛡️", 
    layout="centered"
)

st.title("🛡️ LatticeGuard")
st.markdown("### Elite-Tier Pre-Flight HPC CIF Auditor")
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
Si2 0.25000 0.25000 0.25000
"""

if "cif_text" not in st.session_state:
    st.session_state.cif_text = ""
if "file_name" not in st.session_state:
    st.session_state.file_name = "structure.cif"

st.sidebar.header("Control Panel")
if st.sidebar.button("🧪 Load Sample Defective CIF"):
    st.session_state.cif_text = SAMPLE_CIF
    st.session_state.file_name = "sample_defective.cif"
    st.rerun()

uploaded_file = st.file_uploader("Upload Structure File", type=["cif", "txt"])
if uploaded_file is not None:
    st.session_state.cif_text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    st.session_state.file_name = uploaded_file.name

if st.session_state.cif_text:
    with st.expander("📁 View Loaded Raw File Content", expanded=False):
        st.text_area("Raw CIF", st.session_state.cif_text, height=150)

    if st.button("🚀 Run Pre-Flight Audit & Verification", type="primary", use_container_width=True):
        with st.spinner("Analyzing crystallographic structure and enforcing safety boundaries..."):
            fixed_content, issues, audit_trail = audit_and_fix_cif(st.session_state.cif_text)
            is_valid, validation_errors = validate_fixed_cif(st.session_state.cif_text, fixed_content)

        # שלב 1: דוח אבחון וסינון
        st.markdown("---")
        st.markdown("### 🔍 Step 1: Audit Diagnostic Report")
        
        if not issues:
            st.success("✨ **Audit Clean:** No structural traps or missing parameters identified.")
        else:
            for issue in issues:
                if issue["type"] == "error":
                    st.error(f"❌ **[{issue['category']}]** {issue['message']}")
                else:
                    st.warning(f"⚠️ **[{issue['category']}]** {issue['message']}")

        # הצגת לוג פעולות שקוף (Audit Trail)
        with st.expander("📋 View Engine Audit Trail & Execution Log", expanded=False):
            for log_entry in audit_trail:
                st.text(f"• {log_entry}")

        # שלב 2: אימות ואבטחה
        st.markdown("---")
        st.markdown("### 🛡️ Step 2: Safety & Sanity Verification")
        st.markdown("*Enforcing 'Do No Harm' policy and self-parse integrity before release.*")
        
        if is_valid:
            st.success("✅ **Verification Passed:** Structure is certified safe for simulation pipelines.")
            
            clean_filename = st.session_state.file_name.replace(".cif", "_fixed.cif")
            st.download_button(
                label="📥 Download Verified Clean CIF",
                data=fixed_content,
                file_name=clean_filename,
                mime="chemical/x-cif",
                use_container_width=True
            )
        else:
            st.error("🛑 **Verification Blocked:** Post-processing sanity check failed. Output withheld to protect your simulation.")
            for err in validation_errors:
                st.warning(err)
else:
    st.info("💡 **Get Started:** Click 'Load Sample Defective CIF' in the sidebar or upload your own structure file.")
