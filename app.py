import streamlit as st
import os
import tempfile

# ייבוא פונקציות הליבה מהקובץ שיצרנו קודם
from latticeguard import parse_and_audit_cif, generate_optimized_cif

st.set_page_config(
    page_title="LatticeGuard: Pre-Flight Copilot",
    page_icon="🛡️",
    layout="wide"
)

st.title("🛡️ LatticeGuard: Pre-Flight DFT & Lattice Copilot")
st.markdown("""
**Gatekeeper for High-Stakes HPC Simulations.** 
Upload your raw crystal structure (`.cif`) to catch low-symmetry traps, occupancy disorders, and convergence bottlenecks *before* burning compute budget on cluster runs.
""")

uploaded_file = st.file_uploader("Upload Crystal Structure (.cif)", type=["cif"])

if uploaded_file is not None:
    # שמירת הקובץ שהועלה זמנית
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        tmp_path = tmp_file.name

    with open(tmp_path, "r", encoding="utf-8") as f:
        cif_text = f.read()

    # הרצת האודיט
    audit = parse_and_audit_cif(cif_text)
    optimized_cif = generate_optimized_cif(audit)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📊 Crystallographic Summary")
        st.metric("Space Group", audit['space_group'])
        st.metric("Unit Cell Volume", f"{audit['cell_volume_A3']} Å³")
        st.metric("Total Atomic Sites", audit['total_atoms'])
        
        p = audit['cell_params']
        st.text(f"Cell Parameters:\na = {p['a']} Å\nb = {p['b']} Å\nc = {p['c']} Å\nα={p['alpha']}°, β={p['beta']}°, γ={p['gamma']}°")

    with col2:
        st.subheader("🔍 Engineering Audit Findings")
        if audit['anomalies']:
            for an in audit['anomalies']:
                if an['severity'] == 'HIGH':
                    st.error(f"**[{an['severity']}] {an['code']}**\n\n{an['description']}")
                else:
                    st.warning(f"**[{an['severity']}] {an['code']}**\n\n{an['description']}")
        else:
            st.success("✅ No critical symmetry or occupancy bottlenecks detected. Structure is HPC-ready.")

    st.markdown("---")
    st.subheader("⚙️ Optimized HPC Asset")
    st.text_area("Sanitized CIF File Preview", optimized_cif, height=200)

    # כפתור הורדה לקובץ המנוקה
    st.download_button(
        label="📥 Download Sanitized CIF for HPC",
        data=optimized_cif,
        file_name=f"{uploaded_file.name.split('.')[0]}_LATTICEGUARD_OPTIMIZED.cif",
        mime="text/plain"
    )
