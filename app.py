import streamlit as st
from latticeguard import audit_and_fix_cif, validate_fixed_cif

# הגדרת עמוד Streamlit
st.set_page_config(
    page_title="LatticeGuard — Pre-Flight HPC CIF Auditor", 
    page_icon="🛡️", 
    layout="centered"
)

st.title("🛡️ LatticeGuard — Pre-Flight HPC CIF Auditor")
st.markdown("Catch symmetry traps, missing cell parameters, and structural errors **before** burning expensive HPC hours on VASP or Quantum ESPRESSO.")

# קובץ דוגמה לבדיקה מהירה בממשק
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

# ניהול מצב (Session State) לטקסט הקובץ
if "cif_text" not in st.session_state:
    st.session_state.cif_text = ""
if "file_name" not in st.session_state:
    st.session_state.file_name = "structure.cif"

# צד ה-Sidebar: טעינת דוגמה מהירה
st.sidebar.header("Quick Actions")
if st.sidebar.button("🧪 Load Sample Defective CIF"):
    st.session_state.cif_text = SAMPLE_CIF
    st.session_state.file_name = "sample_defective.cif"
    st.rerun()

# העלאת קובץ מצד המשתמש
uploaded_file = st.file_uploader("Upload your `.cif` structure file", type=["cif", "txt"])
if uploaded_file is not None:
    st.session_state.cif_text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    st.session_state.file_name = uploaded_file.name

# תצוגת תוכן והרצת אודיט
if st.session_state.cif_text:
    st.subheader("📁 Loaded File Preview")
    st.text_area("Raw Content", st.session_state.cif_text, height=150)

    if st.button("🚀 Run Pre-Flight Audit & Fix", type="primary"):
        with st.spinner("Analyzing crystallographic structure and verifying safety boundaries..."):
            # שלב 1: ניתוח ותיקון
            fixed_content, issues = audit_and_fix_cif(st.session_state.cif_text)
            
            # שלב 2: בדיקת שפיות ואמינות (Sanity Check)
            is_valid, validation_errors = validate_fixed_cif(st.session_state.cif_text, fixed_content)

        st.markdown("---")
        st.subheader("🔍 Audit Diagnostic Report")
        
        if not issues:
            st.success("✨ No structural traps or errors detected. Structure looks clean!")
        else:
            for issue in issues:
                if issue["type"] == "error":
                    st.error(f"❌ {issue['message']}")
                else:
                    st.warning(f"⚠️ {issue['message']}")

        st.markdown("---")
        st.subheader("🛡️ Safety & Sanity Verification")
        
        # שלב 3: הצגת כפתור הורדה רק אם הבדיקה עברה בהצלחה מלאה
        if is_valid:
            st.success("✅ Pre-flight verification passed! The processed file is certified safe for simulation pipelines.")
            
            clean_filename = st.session_state.file_name.replace(".cif", "_fixed.cif")
            st.download_button(
                label="📥 Download Verified Clean CIF",
                data=fixed_content,
                file_name=clean_filename,
                mime="chemical/x-cif"
            )
        else:
            st.error("🛑 Safety Warning: Post-processing verification failed. Download blocked to protect your simulation.")
            for err in validation_errors:
                st.warning(err)
else:
    st.info("💡 Tip: Click 'Load Sample Defective CIF' in the sidebar or upload your own file to get started.")
