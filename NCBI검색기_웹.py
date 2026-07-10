import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import io
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

# NCBI 접속시 필요한 이메일 (본인 이메일로 바꿔도 되고, 그냥 둬도 됨)
Entrez.email = "your_email@example.com"

def search_ncbi(search_term, database="nucleotide", max_results=5):
    """NCBI에서 유전자 검색하는 함수"""
    try:
        handle = Entrez.esearch(db=database, term=search_term, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        st.error(f"검색 중 오류 발생: {e}")
        return []


def fetch_sequence(seq_id, database="nucleotide"):
    """NCBI에서 실제 서열 데이터 가져오는 함수"""
    try:
        handle = Entrez.efetch(db=database, id=seq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        st.error(f"서열을 가져오는 중 오류 발생: {e}")
        return None


def analyze_dna(dna_string):
    """DNA 서열을 분석하는 함수"""
    dna_string = dna_string.upper().strip()
    valid_bases = set("ATGC")
    
    if not set(dna_string).issubset(valid_bases):
        return None
    
    seq = Seq(dna_string)
    count_a = dna_string.count("A")
    count_t = dna_string.count("T")
    count_g = dna_string.count("G")
    count_c = dna_string.count("C")
    gc_content = (count_g + count_c) / len(seq) * 100
    
    result = {
        "length": len(seq),
        "count_a": count_a, "count_t": count_t,
        "count_g": count_g, "count_c": count_c,
        "gc_content": gc_content,
        "complement": seq.complement(),
        "reverse_complement": seq.reverse_complement(),
        "rna": seq.transcribe(),
    }
    
    try:
        result["protein"] = seq.translate()
    except Exception:
        result["protein"] = "번역 불가 (길이가 3의 배수 아님)"
    
    return result


# ===== 웹 화면 구성 =====

st.title("🧬🌐 NCBI 유전자 검색기")
st.write("실제 NCBI 데이터베이스에서 유전자 서열을 검색하고 분석해요.")

st.info("예시 검색어: human insulin, SARS-CoV-2 spike protein, BRCA1")

search_term = st.text_input("검색어를 입력하세요:")

if st.button("검색하기"):
    if search_term:
        with st.spinner("NCBI에서 검색 중..."):
            ids = search_ncbi(search_term)
        
        if ids:
            st.success(f"{len(ids)}개의 결과를 찾았습니다!")
            st.session_state["search_ids"] = ids
        else:
            st.warning("검색 결과가 없습니다.")

# 검색 결과가 있으면 선택할 수 있게
if "search_ids" in st.session_state and st.session_state["search_ids"]:
    selected_id = st.selectbox("분석할 서열 ID를 선택하세요:", st.session_state["search_ids"])
    
    if st.button("이 서열 가져와서 분석하기"):
        with st.spinner("서열 데이터를 가져오는 중..."):
            record = fetch_sequence(selected_id)
        
        if record:
            st.write(f"### {record.description}")
            
            full_seq = str(record.seq)
            
            # 서열이 너무 길면 앞부분만 분석 (RNA 번역 정확도를 위해 3의 배수로 자름)
            max_len = 3000
            if len(full_seq) > max_len:
                st.warning(f"서열이 너무 깁니다 (총 {len(full_seq)}염기). 앞부분 {max_len}염기만 분석합니다.")
                full_seq = full_seq[:max_len]
            
            trimmed_len = len(full_seq) - (len(full_seq) % 3)
            full_seq = full_seq[:trimmed_len]
            
            result = analyze_dna(full_seq)
            
            if result:
                st.success(f"✅ 분석한 서열 길이: {result['length']} 염기")
                
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("A", result["count_a"])
                col2.metric("T", result["count_t"])
                col3.metric("G", result["count_g"])
                col4.metric("C", result["count_c"])
                
                st.write("### GC 함량")
                st.progress(result["gc_content"] / 100)
                st.write(f"**{result['gc_content']:.2f}%**")
                
                st.write("### 단백질 번역 (일부)")
                st.code(str(result["protein"])[:200] + "...")