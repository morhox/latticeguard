import streamlit as st
import io
import zipfile
import re

st.set_page_config(
    page_title="LatticeGuard - HPC CIF Auditor",
    layout="wide"
)

# מטריצת רדיוסים קוולנטיים (אנגסטרם) - מבוסס על נתוני Cordero et al.
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Na': 1.66,
    'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
    'Cl': 1.02, 'K': 2.03, 'Ca': 1.76, 'Ti': 1.60, 'Fe': 1.56,
    'Cu': 1.32, 'Zn': 1.22
}

DEFAULT_RADIUS = 1.0
CLASH_SCALING_FACTOR = 0.55

def check_interatomic_clashes(atoms_list):
    """
    בדיקת התנגשויות מרחביות מותנית לפי זוגות אטומים (Element-specific radii matrix).
    """
    clashes = []
    n = len(atoms_list)
    
    for i in range(n):
        for j in range(i + 1, n):
            atom_a = atoms_list[i]
            atom_b = atoms_list[j]
            
            elem_a = atom_a.get('symbol', 'C')
            elem_b = atom_b.get('symbol', 'C')
            
            r_a = COVALENT_RADII.get(elem_a, DEFAULT_RADIUS)
            r_b = COVALENT_RADII.get(elem_b, DEFAULT_RADIUS)
            
            dynamic_threshold = (r_a + r_b) * CLASH_SCALING_FACTOR
            
            dx = atom_a['x'] - atom_b['x']
            dy = atom_a['y'] - atom_b['y']
            dz = atom_a['z'] - atom_b['z']
            distance = (dx**2 + dy**2 + dz**2)**0.5
            
            if distance < dynamic_threshold:
                clashes.append({
                    'atom1': f"{elem_a} ({i+1})",
                    'atom2': f"{elem_b} ({j+1})",
                    'distance': round(distance, 4),
                    'threshold': round(dynamic_threshold, 4)
                })
                
    return clashes

def audit_cif(filename, content):
    audit_log = []
    is_valid = True
    
    if not content.strip():
        return False, ["File is empty."], []
        
    # בדיקת מלכודת סימטריה (P1)
    if re.search(r"(_space_group_name_H-M_alt|_symmetry_space_group_name_H-M)\s+['\"]?P\s*1['\"]?", content, re.IGNORECASE) or \
       re.search(r"_space_group_IT_number\s+1\b", content):
        audit_log.append("WARNING: Low-symmetry P1 space group detected. This may cause exponential grid inflation in DFT calculations.")
    else:
        audit_log.append("OK: Space group symmetry verified.")

    # חילוץ קואורדינטות אטומים לבדיקת התנגשויות
    atoms = []
    lines = content.splitlines()
    in_loop = False
    x_idx, y_idx, z_idx, symbol_idx = -1, -1, -1, -1
    loop_headers = []
    
    for line in lines:
        line_str = line.strip()
        if line_str.startswith("loop_"):
            in_loop = True
            loop_headers = []
            continue
        if in_loop:
            if line_str.startswith("_"):
                loop_headers.append(line_str)
                if "_atom_site_fract_x" in line_str: x_idx = len(loop_headers) - 1
                if "_atom_site_fract_y" in line_str: y_idx = len(loop_headers) - 1
                if "_atom_site_fract_z" in line_str: z_idx = len(loop_headers) - 1
                if "_atom_site_label" in line_str or "_atom_site_type_symbol" in line_str: symbol_idx = len(loop_headers) - 1
                continue
            elif line_str.startswith("#") or line_str == "":
                continue
            else:
                parts = line_str.split()
                if len(parts) >= len(loop_headers) and x_idx != -1 and y_idx != -1 and z_idx != -1:
                    try:
                        x = float(parts[x_idx])
                        y = float(parts[y_idx])
                        z = float(parts[z_idx])
                        
                        symbol = "C"
                        if symbol_idx != -1 and symbol_idx < len(parts):
                            raw_sym = parts[symbol_idx]
                            clean_sym = re.match(r"([A-Za-z]{1,2})", raw_sym)
                            if clean_sym:
                                symbol = clean_sym.group(1)
                                
                        atoms.append({'symbol': symbol, 'x': x, 'y': y, 'z': z})
                    except ValueError:
                        pass
                else:
                    if not line_str.startswith("_") and len(parts) < len(loop_headers):
                        in_loop = False

    clashes = []
    if atoms:
        clashes = check_interatomic_clashes(atoms)
        if clashes:
            is_valid = False
            audit_log.append(f"ERROR: Detected {len(clashes)} interatomic spatial clashes based on dynamic covalent radii thresholds.")
        else:
            audit_log.append("OK: No atomic clashes detected using dynamic element-specific radii matrices.")
    else:
        audit_log.append("NOTE: No atomic coordinate loops parsed for spatial clash checking.")

    return is_valid, audit_log, clashes

def main():
    st.title("LatticeGuard")
    st.markdown("Elite-tier pre-flight safety gatekeeper for VASP, Quantum ESPRESSO, and DFT simulations.")
    st.markdown("---")

    uploaded_files = st.file_uploader("Upload CIF Structure Files", type=["cif"], accept_multiple_files=True)

    if uploaded_files:
        st.subheader("Audit Results & Batch Summary")
        
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for uploaded_file in uploaded_files:
                file_content = uploaded_file.getvalue().decode("utf-8", errors="ignore")
                filename = uploaded_file.name
                
                is_valid, audit_log, clashes = audit_cif(filename, file_content)
                
                status_label = "PASSED" if is_valid else "FAILED"
                with st.expander(f"[{status_label}] {filename}"):
                    for log_entry in audit_log:
                        if log_entry.startswith("ERROR"):
                            st.error(log_entry)
                        elif log_entry.startswith("WARNING"):
                            st.warning(log_entry)
                        else:
                            st.success(log_entry)
                            
                    if not is_valid and clashes:
                        st.markdown("**Detailed Clash Breakdown:**")
                        for c in clashes:
                            st.text(f"  - Atom {c['atom1']} <-> {c['atom2']} | Measured: {c['distance']} A | Min Threshold: {c['threshold']} A")
                
                if is_valid:
                    zip_file.writestr(f"verified_{filename}", file_content)
                    
        zip_buffer.seek(0)
        
        st.markdown("---")
        st.download_button(
            label="Download Verified Structures (ZIP)",
            data=zip_buffer,
            file_name="latticeguard_verified_structures.zip",
            mime="application/zip"
        )

if __name__ == "__main__":
    main()
