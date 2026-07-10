import streamlit as st
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, AllChem, DataStructs
import pubchempy as pcp
import ssl
import io

# SSL 우회 (NCBI, PubChem 접속용)
ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "your_email@example.com"

st.set_page_config(page_title="신약개발 툴킷", page_icon="💊", layout="wide")

# ===== 사이드바 메뉴 =====
st.sidebar.title("💊 신약개발 툴킷")
menu = st.sidebar.radio(
    "기능을 선택하세요:",
    [
        "🏠 홈",
        "🧪 분자량 계산기",
        "🧬 DNA 서열 분석",
        "🧬🔬 DNA 서열 비교",
        "🧬📁 DNA 파일 분석",
        "🧬🌐 NCBI 검색",
        "💊 화합물 분석",
        "💊🔬 화합물 유사도 비교",
        "💊🌐 PubChem 검색",
    ]
)

# ===== 공통 함수들 =====

atomic_weights = {
    "H": 1.008, "C": 12.011, "O": 15.999, "N": 14.007,
    "Na": 22.990, "Cl": 35.453, "S": 32.06,
}

def parse_formula(formula):
    result = {}
    i = 0
    while i < len(formula):
        element = formula[i]
        i += 1
        if i < len(formula) and formula[i].islower():
            element += formula[i]
            i += 1
        num = ""
        while i < len(formula) and formula[i].isdigit():
            num += formula[i]
            i += 1
        count = int(num) if num else 1
        result[element] = result.get(element, 0) + count
    return result


def analyze_dna(dna_string):
    dna_string = dna_string.upper().strip()
    valid_bases = set("ATGC")
    if not set(dna_string).issubset(valid_bases):
        return None
    seq = Seq(dna_string)
    count_a, count_t = dna_string.count("A"), dna_string.count("T")
    count_g, count_c = dna_string.count("G"), dna_string.count("C")
    gc_content = (count_g + count_c) / len(seq) * 100
    result = {
        "length": len(seq), "count_a": count_a, "count_t": count_t,
        "count_g": count_g, "count_c": count_c, "gc_content": gc_content,
        "complement": seq.complement(), "reverse_complement": seq.reverse_complement(),
        "rna": seq.transcribe(),
    }
    try:
        result["protein"] = seq.translate()
    except Exception:
        result["protein"] = "번역 불가"
    return result


# ===== 페이지별 화면 =====

if menu == "🏠 홈":
    st.title("💊 신약개발 툴킷")
    st.write("DNA 분석부터 화합물 검색까지, 신약개발 기초 도구 모음입니다.")
    st.info("왼쪽 사이드바에서 원하는 기능을 선택하세요!")
    
    st.write("### 포함된 기능")
    st.markdown("""
    - 🧪 **분자량 계산기**: 화학식으로 분자량 계산
    - 🧬 **DNA 서열 분석**: 염기 구성, GC함량, 단백질 번역
    - 🧬🔬 **DNA 서열 비교**: 두 서열의 유사도 및 돌연변이 위치 찾기
    - 🧬📁 **DNA 파일 분석**: FASTA 파일 업로드 분석
    - 🧬🌐 **NCBI 검색**: 실제 유전자 데이터베이스 검색
    - 💊 **화합물 분석**: SMILES 기반 Lipinski's Rule 평가
    - 💊🔬 **화합물 유사도 비교**: Tanimoto 유사도 계산
    - 💊🌐 **PubChem 검색**: 약물 이름으로 자동 검색
    """)


elif menu == "🧪 분자량 계산기":
    st.title("🧪 분자량 계산기")
    formula = st.text_input("화학식을 입력하세요 (예: H2O, C6H12O6, NaCl):")
    
    if st.button("계산하기"):
        if formula:
            elements = parse_formula(formula)
            total_weight = 0
            details = []
            unknown = None
            
            for element, count in elements.items():
                if element not in atomic_weights:
                    unknown = element
                    break
                weight = atomic_weights[element] * count
                total_weight += weight
                details.append((element, count, atomic_weights[element], weight))
            
            if unknown:
                st.error(f"⚠️ '{unknown}' 원소 정보가 없습니다.")
            else:
                st.success(f"✅ 총 분자량: **{total_weight:.3f} g/mol**")
                for element, count, weight_per, total_w in details:
                    st.write(f"- {element}: {count}개 × {weight_per} = {total_w:.3f}")


elif menu == "🧬 DNA 서열 분석":
    st.title("🧬 DNA 서열 분석기")
    dna_input = st.text_input("DNA 서열을 입력하세요 (예: ATGGCCATTGTAATG):")
    
    if st.button("분석하기"):
        if dna_input:
            result = analyze_dna(dna_input)
            if result is None:
                st.error("⚠️ A, T, G, C만 입력해주세요.")
            else:
                st.success(f"서열 길이: {result['length']} 염기")
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("A", result["count_a"])
                col2.metric("T", result["count_t"])
                col3.metric("G", result["count_g"])
                col4.metric("C", result["count_c"])
                
                st.write("### GC 함량")
                st.progress(result["gc_content"] / 100)
                st.write(f"**{result['gc_content']:.2f}%**")
                
                st.write("### 결과")
                st.code(f"상보서열: {result['complement']}\n"
                        f"역상보서열: {result['reverse_complement']}\n"
                        f"RNA: {result['rna']}\n"
                        f"단백질: {result['protein']}")


elif menu == "🧬🔬 DNA 서열 비교":
    st.title("🧬🔬 DNA 서열 비교기")
    col1, col2 = st.columns(2)
    with col1:
        seq1_input = st.text_input("첫 번째 서열:")
    with col2:
        seq2_input = st.text_input("두 번째 서열:")
    
    if st.button("비교하기"):
        if seq1_input and seq2_input:
            s1, s2 = seq1_input.upper().strip(), seq2_input.upper().strip()
            valid = set("ATGC")
            if not set(s1).issubset(valid) or not set(s2).issubset(valid):
                st.error("⚠️ A, T, G, C만 입력해주세요.")
            else:
                min_len = min(len(s1), len(s2))
                matches = sum(1 for i in range(min_len) if s1[i] == s2[i])
                diffs = [i+1 for i in range(min_len) if s1[i] != s2[i]]
                similarity = matches / min_len * 100
                
                st.metric("서열 유사도", f"{similarity:.2f}%")
                st.progress(similarity / 100)
                st.write(f"일치: {matches}개 / 차이: {len(diffs)}개")
                if diffs:
                    st.code(", ".join(str(d) for d in diffs[:20]))


elif menu == "🧬📁 DNA 파일 분석":
    st.title("🧬📁 DNA 파일 분석기")
    uploaded_file = st.file_uploader("FASTA 파일 업로드", type=["fasta", "fa", "txt"])
    
    if uploaded_file:
        stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        records = list(SeqIO.parse(stringio, "fasta"))
        
        if records:
            names = [r.id for r in records]
            selected = st.selectbox("서열 선택:", names)
            record = next(r for r in records if r.id == selected)
            dna_seq = str(record.seq)
            
            if st.button("분석하기"):
                result = analyze_dna(dna_seq)
                if result:
                    st.success(f"길이: {result['length']} 염기")
                    st.write(f"GC 함량: {result['gc_content']:.2f}%")
                    st.code(f"단백질: {result['protein']}")


elif menu == "🧬🌐 NCBI 검색":
    st.title("🧬🌐 NCBI 유전자 검색기")
    search_term = st.text_input("검색어 (예: human insulin):")
    
    if st.button("검색하기"):
        if search_term:
            with st.spinner("검색 중..."):
                handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5)
                record = Entrez.read(handle)
                handle.close()
                ids = record["IdList"]
            
            if ids:
                st.session_state["ncbi_ids"] = ids
                st.success(f"{len(ids)}개 결과 발견!")
    
    if "ncbi_ids" in st.session_state:
        selected_id = st.selectbox("서열 선택:", st.session_state["ncbi_ids"])
        if st.button("가져와서 분석"):
            with st.spinner("가져오는 중..."):
                handle = Entrez.efetch(db="nucleotide", id=selected_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
            
            st.write(f"### {record.description}")
            seq_str = str(record.seq)[:3000]
            seq_str = seq_str[:len(seq_str) - len(seq_str) % 3]
            result = analyze_dna(seq_str)
            if result:
                st.write(f"GC 함량: {result['gc_content']:.2f}%")
                st.code(f"단백질(일부): {str(result['protein'])[:200]}")


elif menu == "💊 화합물 분석":
    st.title("💊 화합물 분석기 (SMILES)")
    smiles = st.text_input("SMILES 입력 (예: CC(=O)OC1=CC=CC=C1C(=O)O):")
    
    if st.button("분석하기"):
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                st.error("⚠️ 올바른 SMILES가 아닙니다.")
            else:
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hd = Descriptors.NumHDonors(mol)
                ha = Descriptors.NumHAcceptors(mol)
                
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("분자량", f"{mw:.1f}")
                col2.metric("LogP", f"{logp:.2f}")
                col3.metric("H공여체", hd)
                col4.metric("H수용체", ha)
                
                passed = sum([mw<=500, hd<=5, ha<=10, logp<=5])
                st.write(f"### Lipinski's Rule: {passed}/4 통과")
                
                img = Draw.MolToImage(mol, size=(400, 400))
                st.image(img, caption="분자 구조")


elif menu == "💊🔬 화합물 유사도 비교":
    st.title("💊🔬 화합물 유사도 비교기")
    col1, col2 = st.columns(2)
    with col1:
        smiles1 = st.text_input("화합물 1 SMILES:")
    with col2:
        smiles2 = st.text_input("화합물 2 SMILES:")
    
    if st.button("비교하기"):
        if smiles1 and smiles2:
            mol1, mol2 = Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)
            if mol1 is None or mol2 is None:
                st.error("⚠️ 올바른 SMILES가 아닙니다.")
            else:
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                
                st.metric("구조 유사도", f"{similarity*100:.1f}%")
                st.progress(similarity)
                
                if similarity >= 0.7:
                    st.success("🟢 매우 유사한 구조")
                elif similarity >= 0.4:
                    st.warning("🟡 어느 정도 유사")
                else:
                    st.error("🔴 구조적으로 많이 다름")


elif menu == "💊🌐 PubChem 검색":
    st.title("💊🌐 PubChem 화합물 검색기")
    drug_name = st.text_input("약물 이름 (예: aspirin):")
    
    if st.button("검색하기"):
        if drug_name:
            with st.spinner("검색 중..."):
                compounds = pcp.get_compounds(drug_name, 'name')
            
            if compounds:
                compound = compounds[0]
                st.write(f"### {compound.iupac_name}")
                st.write(f"분자식: {compound.molecular_formula}")
                st.code(f"SMILES: {compound.canonical_smiles}")
                
                mol = Chem.MolFromSmiles(compound.canonical_smiles)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    st.metric("분자량", f"{mw:.2f}")
                    img = Draw.MolToImage(mol, size=(400, 400))
                    st.image(img, caption=drug_name)
            else:
                st.error("⚠️ 검색 결과가 없습니다.")